#!/usr/bin/env Rscript

print("Run scMEGA on paired RNA/ATAC object to infer GRN using only genes in GMs + TFs")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(Seurat)
library(Signac)
library(scMEGA)
library(SummarizedExperiment)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(VennDiagram)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("", "--verbose"), action = "store", type = "logical", help = "Verbose", default = FALSE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    addArchRThreads(threads = 1)
    
    # full data
    data_path = "./output/NF-downstream_analysis/Processing/FullData/scMEGA/MEGA_cell_pairing_and_chromvar/rds_files/" # paired object with chromvar run
    
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/"
    rds_path = "./rds_files/"
    csv_path = "./csv_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  dir.create(csv_path, recursive = T)
}

set.seed(42)

############################## EDITED FUNCTIONS #######################################

TrajectoryHeatmap_updated <- function (trajectory, varCutOff = 0.9, maxFeatures = 25000, scaleRows = TRUE, 
                                       rowOrder = NULL, limits = c(-1.5, 1.5), labelRows = FALSE, 
                                       pal = NULL, labelMarkers = NULL, labelTop = 50, name = "Heatmap", 
                                       returnMatrix = FALSE) 
{
  mat <- assay(trajectory)
  rSNA <- rowSums(is.na(mat))
  if (sum(rSNA > 0) > 0) {
    message("Removing rows with NA values...")
    mat <- mat[rSNA == 0, ]
  }
  varQ <- ArchR:::.getQuantiles(matrixStats::rowVars(mat))
  orderedVar <- FALSE
  if (is.null(rowOrder)) {
    mat <- mat[order(varQ, decreasing = TRUE), ]
    orderedVar <- TRUE
    if (is.null(varCutOff) & is.null(maxFeatures)) {
      n <- nrow(mat)
    }
    else if (is.null(varCutOff)) {
      n <- maxFeatures
    }
    else if (is.null(maxFeatures)) {
      n <- (1 - varCutOff) * nrow(mat)
    }
    else {
      n <- min((1 - varCutOff) * nrow(mat), maxFeatures)
    }
    n <- min(n, nrow(mat))
    mat <- mat[head(seq_len(nrow(mat)), n), ]
  }
  if (!is.null(labelTop) & labelTop > 0) {
    if (orderedVar) {
      idxLabel <- rownames(mat)[seq_len(labelTop)]
    }
    else {
      idxLabel <- rownames(mat)[order(varQ, decreasing = TRUE)][seq_len(labelTop)]
    }
  }
  else {
    idxLabel <- NULL
  }
  if (scaleRows) {
    # Subset the matrix to include only rows with non-zero standard deviations
    row_sds <- matrixStats::rowSds(mat)
    mat <- mat[row_sds != 0, ]
    # then can scale and shouldnt create any NAs
    mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), 
                 `/`)
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
  }
  if (nrow(mat) == 0) {
    stop("No Features Remaining!")
  }
  if (is.null(pal)) {
    pal <- ArchR::paletteContinuous(set = "blueYellow", n = 100)
  }
  if (!is.null(rowOrder)) {
    idx <- rowOrder
  }
  else {
    idx <- order(apply(mat, 1, which.max))
  }
  if (!is.null(idxLabel)) {
    customRowLabel <- match(idxLabel, rownames(mat[idx, ]))
  }
  else {
    customRowLabel <- NULL
  }
  ht <- ArchR:::.ArchRHeatmap(mat = mat[idx, ], scale = FALSE,
                              limits = c(min(mat), max(mat)), color = pal, clusterCols = FALSE,
                              clusterRows = FALSE, labelRows = labelRows, labelCols = FALSE,
                              customRowLabel = customRowLabel, showColDendrogram = TRUE,
                              name = name, draw = FALSE)
  if (returnMatrix) {
    return(mat[idx, ])
  }
  else {
    return(ht)
  }
  return(mat)
}

SelectTFs_updated <- function (object, tf.assay = "chromvar", rna.assay = "RNA", atac.assay = "ATAC", 
                               trajectory.name = "Trajectory", groupEvery = 1, return.heatmap = TRUE,
                               p.cutoff = 0.01, cor.cutoff = 0.3) 
  {
  trajMM <- GetTrajectory_updated(object, assay = tf.assay, trajectory.name = trajectory.name, 
                                  groupEvery = groupEvery, slot = "data", smoothWindow = 7, 
                                  log2Norm = FALSE)
  rownames(trajMM) <- object@assays[[atac.assay]]@motifs@motif.names
  trajRNA <- GetTrajectory_updated(object, assay = rna.assay, trajectory.name = trajectory.name, 
                                   groupEvery = groupEvery, slot = "data", smoothWindow = 7, 
                                   log2Norm = TRUE)
  
  df.cor <- GetCorrelation(trajMM, trajRNA) # correlation between TF expression RNAseq and binding activity ATACseq using cor.test
  if (is.null(p.cutoff) & is.null(cor.cutoff)){
    print("No filtering being applied to TFs!")
  } else{
    print("TFs being filtered on correlation of expression and chromvar score!")
    df.cor <- df.cor[df.cor$adj_p < p.cutoff & df.cor$correlation > 
                       cor.cutoff, ]
  }
  
  matMM <- suppressMessages(TrajectoryHeatmap_updated(trajMM, varCutOff = 0, 
                                              pal = paletteContinuous(set = "solarExtra"), limits = c(-2, 
                                                                                                      2), name = "TF activity", returnMatrix = TRUE))
  df_tf_time_point <- data.frame(tfs = rownames(matMM), time_point = seq(1, 
                                                                         100, length.out = nrow(matMM)))
  rownames(df_tf_time_point) <- df_tf_time_point$tfs
  df_tf_time_point <- df_tf_time_point[df.cor$tfs, ]
  df.cor$time_point <- df_tf_time_point$time_point
  df.cor <- df.cor[order(df.cor$time_point), ]
  trajMM <- trajMM[df.cor$tfs, ]
  trajRNA <- trajRNA[df.cor$tfs, ]
  if (return.heatmap) {
    ht <- suppressMessages(CorrelationHeatmap(trajectory1 = trajMM, 
                                              trajectory2 = trajRNA, name1 = "TF activity", name2 = "Gene expression"))
    res <- list(tfs = df.cor, heatmap = ht)
  }
  else {
    res <- list(tfs = df.cor)
  }
  return(res)
}

SelectGenes_updated <- function (object, atac.assay = "ATAC", rna.assay = "RNA", var.cutoff.gene = 0.9, 
                                 trajectory.name = "Trajectory", distance.cutoff = 2000, groupEvery = 1, 
                                 cor.cutoff = 0, fdr.cutoff = 1e-04, return.heatmap = TRUE, 
                                 labelTop1 = 10, labelTop2 = 10, genome = "hg38") 
{
  trajRNA <- GetTrajectory_updated(object, assay = rna.assay, trajectory.name = trajectory.name, 
                                   groupEvery = groupEvery, slot = "data", smoothWindow = 7, 
                                   log2Norm = TRUE)
  trajATAC <- GetTrajectory_updated(object, assay = atac.assay, groupEvery = groupEvery, 
                                    trajectory.name = trajectory.name, slot = "data", smoothWindow = 7, 
                                    log2Norm = TRUE)
  groupMatRNA <- suppressMessages(TrajectoryHeatmap_updated(trajRNA, 
                                                    varCutOff = var.cutoff.gene, pal = paletteContinuous(set = "horizonExtra"), 
                                                    limits = c(-2, 2), returnMatrix = TRUE))
  groupMatATAC <- suppressMessages(TrajectoryHeatmap_updated(trajATAC, 
                                                     varCutOff = 0, maxFeatures = nrow(trajATAC), pal = paletteContinuous(set = "solarExtra"), 
                                                     limits = c(-2, 2), name = "Chromatin accessibility", 
                                                     returnMatrix = TRUE))
  message("Linking cis-regulatory elements to genes...")
  df.p2g <- PeakToGene(peak.mat = groupMatATAC, gene.mat = groupMatRNA, 
                       genome = genome)
  df.p2g <- df.p2g %>% subset(distance > distance.cutoff) %>% 
    subset(Correlation > cor.cutoff & FDR < fdr.cutoff)
  trajATAC <- trajATAC[df.p2g$peak, ]
  trajRNA <- trajRNA[df.p2g$gene, ]
  if (return.heatmap) {
    ht <- suppressMessages(CorrelationHeatmap(trajectory1 = trajATAC, 
                                              trajectory2 = trajRNA, name1 = "Chromatin accessibility", 
                                              name2 = "Gene expression", labelTop1 = labelTop1, 
                                              labelTop2 = labelTop2, labelRows1 = FALSE, labelRows2 = FALSE))
    res <- list(p2g = df.p2g, heatmap = ht)
  }
  else {
    res <- list(p2g = df.p2g)
  }
  return(res)
}

GetTFGeneCorrelation_updated <- function (object, tf.use = NULL, gene.use = NULL, tf.assay = "chromvar", 
                                          gene.assay = "RNA", atac.assay = "ATAC", trajectory.name = "Trajectory", 
                                          groupEvery = 1) 
{
  trajMM <- GetTrajectory_updated(object, assay = tf.assay, slot = "data", 
                                  trajectory.name = trajectory.name, groupEvery = groupEvery, 
                                  smoothWindow = 7, log2Norm = FALSE)
  trajRNA <- GetTrajectory_updated(object, assay = gene.assay, slot = "data", 
                                   trajectory.name = trajectory.name, groupEvery = groupEvery, 
                                   smoothWindow = 7, log2Norm = TRUE)
  rownames(trajMM) <- object@assays[[atac.assay]]@motifs@motif.names
  tf_activity <- suppressMessages(TrajectoryHeatmap_updated(trajMM, 
                                                    varCutOff = 0, pal = paletteContinuous(set = "solarExtra"), 
                                                    limits = c(-2, 2), name = "TF activity", returnMatrix = TRUE))
  gene_expression <- suppressMessages(TrajectoryHeatmap_updated(trajRNA, 
                                                        varCutOff = 0, pal = paletteContinuous(set = "solarExtra"), 
                                                        limits = c(-2, 2), name = "Gene expression", returnMatrix = TRUE))
  if (!is.null(tf.use)) {
    tf_activity <- tf_activity[tf.use, ]
  }
  if (!is.null(gene.use)) {
    sel_genes <- intersect(rownames(gene_expression), gene.use)
    gene_expression <- gene_expression[sel_genes, ]
  }
  df.cor <- t(cor(t(tf_activity), t(gene_expression))) %>% 
    as.data.frame()
  if (!is.null(tf.use)) {
    df.cor <- df.cor[, tf.use]
  }
  df.cor$gene <- rownames(df.cor)
  df.cor <- df.cor %>% tidyr::pivot_longer(!gene, names_to = "tf", 
                                           values_to = "correlation") %>% dplyr::select(c(tf, gene, correlation))
  df.cor$t_stat <- (df.cor$correlation/sqrt((pmax(1 - df.cor$correlation^2, 
                                                  1e-17, na.rm = TRUE))/(ncol(tf_activity) - 2)))
  df.cor$p_value <- 2 * pt(-abs(df.cor$t_stat), ncol(tf_activity) - 
                             2)
  df.cor$fdr <- p.adjust(df.cor$p_value, method = "fdr")
  return(df.cor)
}


GetTrajectory_updated <- function (object = NULL, trajectory.name = "Trajectory", assay = NULL, 
                                   slot = "counts", groupEvery = 1, log2Norm = TRUE, scaleTo = 10000, 
                                   smoothWindow = 11) 
{
  if (is.null(assay) | !assay %in% Assays(object)) {
    stop("Please provide an available assay!")
  }
  if (!(trajectory.name %in% colnames(object@meta.data))) {
    stop(glue::glue("Cannot find trajecotry {trajectory.name}!"))
  }
  trajectory <- object@meta.data[trajectory.name]
  trajectory <- trajectory[!is.na(trajectory[, 1]), , drop = FALSE]
  breaks <- seq(0, 100, groupEvery)
  if (!all(is.numeric(trajectory[, 1]))) {
    stop("Trajectory must be a numeric. Did you add the trajectory with addTrajectory?")
  }
  if (!all(trajectory[, 1] >= 0 & trajectory[, 1] <= 100)) {
    stop("Trajectory values must be between 0 and 100. Did you add the trajectory with addTrajectory?")
  }
  groupList <- lapply(seq_along(breaks), function(x) {
    if (x == 1) {
      NULL
    }
    else {
      rownames(trajectory)[which(trajectory[, 1] > breaks[x - 
                                                            1] & trajectory[, 1] <= breaks[x])]
    }
  })[-1]
  names(groupList) <- paste0("T.", breaks[-length(breaks)], 
                             "_", breaks[-1])
  message("Creating Trajectory Group Matrix..")
  
  ### LayerData fails due to seurat package version issue so have hacked this section myself to extract and deal with
  # both ATAC/RNA count data but also motif data (NB motif data must be in a slot called 'chromvar')
  #data.use <- GetAssayData(object, assay = assay, slot = slot)
  #data.use <- LayerData(object, assay = assay)
  if (assay %in% c("chromvar", "target")){
    data.use <- object[[assay]]@data
  } else {
    data.use <- object[[assay]]$counts
  }
  
  # Adapted due to errors when there is only one cell in a group
  groupMat <- lapply(1:length(groupList), function(x) {
    cell_names <- groupList[[x]]
    if (length(cell_names) == 1) {
      mat <- data.use[, cell_names]
    } else {
      mat <- Matrix::rowMeans(data.use[, cell_names]) # rowMeans doesnt work if there is only one cell in group
    }
  }) %>% Reduce(cbind, .)
  
  # Added this to avoid the NAs that arise when there are no cells in a group
  groupMat[is.na(groupMat)] <- 0
  
  colnames(groupMat) <- names(groupList)
  if (!is.null(scaleTo)) {
    if (any(groupMat < 0)) {
      message("Some values are below 0, this could be the Motif activity matrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!")
    }
    else {
      groupMat <- t(t(groupMat)/colSums(groupMat)) * scaleTo
    }
  }
  
  # Added this because NaNs arrive now
  groupMat[is.nan(groupMat)] <- 0
  
  if (log2Norm) {
    if (any(groupMat < 0)) {
      message("Some values are below 0, this could be a Motif activity matrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!")
    }
    else {
      groupMat <- log2(groupMat + 1)
    }
  }
  if (!is.null(smoothWindow)) {
    message("Smoothing...")
    smoothGroupMat <- as.matrix(t(apply(groupMat, 1, function(x) centerRollMean(x, 
                                                                                k = smoothWindow))))
    colnames(smoothGroupMat) <- paste0(colnames(groupMat))
    colnames(groupMat) <- paste0(colnames(groupMat))
    seTrajectory <- SummarizedExperiment(assays = SimpleList(smoothMat = as.matrix(smoothGroupMat), 
                                                             mat = as.matrix(groupMat)))
  }
  else {
    colnames(groupMat) <- paste0(colnames(groupMat))
    seTrajectory <- SummarizedExperiment(assays = SimpleList(mat = as.matrix(groupMat)))
  }
  return(seTrajectory)
}


centerRollMean <- function(v = NULL, k = NULL){
  o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
  if(k%%2==0){
    o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else if(k%%2==1){
    o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else{
    stop("Error!")
  }
  o2
}

CorrelationHeatmap_updated <- function (trajectory1, trajectory2, name1 = NULL, name2 = NULL, 
          labelRows1 = TRUE, labelRows2 = TRUE, labelTop1 = 50, labelTop2 = 50, 
          limits1 = c(-2, 2), limits2 = c(-2, 2)) 
{
  trajCombined <- trajectory1
  assay(trajCombined, withDimnames = FALSE) <- t(apply(assay(trajectory2), 
                                                       1, scale)) + t(apply(assay(trajectory1), 1, scale))
  combinedMat <- TrajectoryHeatmap_updated(trajCombined, returnMatrix = TRUE, 
                                   varCutOff = 0)
  rowOrder <- match(rownames(combinedMat), rownames(trajectory1))
  ht1 <- TrajectoryHeatmap_updated(trajectory1, pal = paletteContinuous(set = "solarExtra"), 
                           varCutOff = 0, maxFeatures = nrow(trajectory1), rowOrder = rowOrder, 
                           limits = limits1, labelRows = labelRows1, labelTop = labelTop1, 
                           name = name1)
  ht2 <- TrajectoryHeatmap_updated(trajectory2, pal = paletteContinuous(set = "horizonExtra"), 
                           varCutOff = 0, maxFeatures = nrow(trajectory2), rowOrder = rowOrder, 
                           limits = limits2, labelRows = labelRows2, labelTop = labelTop2, 
                           name = name2)
  ht <- ht1 + ht2
  return(ht)
}


############################## Read in seurat object #######################################

print("reading in data...")

## read in paired seurat object
obj.pair <- readRDS(paste0(data_path, "rds_files/paired_object_chromvar.RDS"))
obj.pair

## read in P2G linkage df
P2G <- read_csv(paste0(data_path, "Peak_to_gene_linkage_df_250000_distance.csv"))

############################## Create trajectories from lineage probabilities #######################################

head(obj.pair@meta.data)

# then need them to be between 0 and 100 - check for placodal
summary(obj.pair@meta.data$lineage_placodal_probability)
obj.pair@meta.data$lineage_placodal_probability <- obj.pair@meta.data$lineage_placodal_probability * 100
summary(obj.pair@meta.data$lineage_placodal_probability)
# hist(obj.pair@meta.data$lineage_placodal_probability, breaks = 100)

# then run for neural and NC
obj.pair@meta.data$lineage_NC_probability <- obj.pair@meta.data$lineage_NC_probability * 100
obj.pair@meta.data$lineage_neural_probability <- obj.pair@meta.data$lineage_neural_probability * 100

############################## Save seurat object #######################################

print("saving whole seurat object:")

saveRDS(obj.pair, paste0(rds_path, "paired_object_chromvar.RDS"), compress = FALSE)

############################## Extract all known TF names #######################################

print("extracting TF names...")

TF_names <- obj.pair@assays$ATAC@motifs@motif.names
names(TF_names) <- NULL
TF_names <- unlist(TF_names)
length(TF_names) # 746 known TFs
fileConn <- file(paste0(csv_path, "Known_TF_names.txt"))
writeLines(TF_names, fileConn)
close(fileConn)

############################## GENE MODULES #######################################

# set up list of genes
GM12 <- c("AKR1D1", "ATP1A1", "ATP1B1", "ATP2B1", "B3GNT7", "CCDC25", "CGNL1", "DAG1", "EMILIN2", "ENSGALG00000004814", "EPCAM", "FAM184B", "FAM89A", "GATA2", "GATA3", "IRAK2", "IRF6", "MAP7", "METRNL", "MPZL3", "NET1", "PLEKHA5", "POGLUT2", "PPFIBP1", "RGN", "SLC16A10", "SLC25A4", "TSPAN13", "TTC39A", "UNC5B", "Z-TJP2")
GM14 <- c("BRINP1", "CITED4", "DLX5", "DLX6", "ENSGALG00000023936", "ENSGALG00000040010", "FN1", "HESX1", "HIF1A", "KCNAB1", "LAMB1", "LRP11", "NFKB1", "PAX6", "PITX1", "PITX2", "SEM1", "SFRP1", "SH3D19", "SHISA2", "SIX3", "SPON1", "SST", "WDR1", "Z-HAPLN1")
GM13 <- c("AKAP12", "ASS1", "BASP1", "CD99", "ENSGALG00000011296", "ENSGALG00000041054", "ENSGALG00000042443", "EYA2", "FERMT2", "LGMN", "METTL24", "NR2F2", "NUCKS1", "SIX1")
# GM23 <- c("ASS1", "BMP4", "BMP6", "CLDN3", "CSRP2", "DLX5", "DLX6", "EMILIN2", "ENSGALG00000001885", "ENSGALG00000023936", "ENSGALG00000040010", "ENSGALG00000042443", "ENSGALG00000050334", "ENSGALG00000052786", "EYA2", "FABP3", "FAM184B", "FAM89A", "FN1", "GATA2", "GATA3", "HAS2", "KCNAB1", "KRT18", "KRT19", "KRT7", "LAMB1", "NEDD9", "NET1", "PITX1", "SIX1", "SPON1", "TFAP2A", "TUBAL3", "UNC5B", "Z-HAPLN1")
gm_genes <- unique(c(GM12, GM14, GM13))
length(GM12)
length(GM14)
length(GM13)

# plot venn of nodes
myCol <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(GM12, GM14, GM13),
  category.names = c("GM12", "GM14", "GM13"),
  filename = paste0(plot_path, 'Gene_modules.png'),
  output=TRUE, disable.logging = TRUE,
  # Output features
  imagetype="png",
  height = 600, 
  width = 600, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.7,
  cat.fontface = "bold"
)


######################################################################################
##############################    PLACODAL     #######################################
######################################################################################

print("Starting GRNi for placodal lineage...")

temp_csv_path = "./csv_files/placodal_lineage/"
dir.create(temp_csv_path, recursive = T)
# temp_lineage_plot_path = "./plots/placodal_lineage/lineage_dynamics_plots/"
# dir.create(temp_lineage_plot_path, recursive = T)

# set trajectory
trajectory <- "lineage_placodal_probability"

# remove cells that have a placodal probability of less than 25%
# obj.pair[["RNA"]] <- as(obj.pair[["RNA"]], "Assay5") # debugging seurat's subset function
obj.traj <- subset(obj.pair, subset = lineage_placodal_probability > 25)
obj.traj

# save subsetted object
saveRDS(obj.traj, file = paste0(rds_path, "Placodal_traj_obj.RDS"))

############################## Select nodes #######################################

temp_plot_path = "./plots/placodal_lineage/node_selection/"
dir.create(temp_plot_path, recursive = T)

############  SOURCE NODES

print("Selecting source nodes...")

# extract the source nodes - not filtering based on correlation of gene expression and binding
res <- SelectTFs_updated(object = obj.traj, trajectory.name = trajectory, return.heatmap = TRUE,
                         groupEvery = 2,
                         p.cutoff = NULL, cor.cutoff = NULL)

# save the selected source nodes
df.tfs <- res$tfs
write.csv(df.tfs, file = paste0(temp_csv_path, "TF_correlations_all.csv"), row.names = FALSE)

# how many source nodes are there
nrow(df.tfs) # 311
unique(df.tfs$tfs)

# plot TF activity dynamics across trajectory
ht <- res$heatmap
png(paste0(temp_plot_path, 'TFs_all_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
draw(ht)
graphics.off()

# sense check that all source nodes are also a known TF
if (sum(df.tfs$tfs %in% TF_names) != length(df.tfs$tfs)){
  stop("Problem! Not all source nodes are known TFs!")
}

############  TARGET NODES

print("Selecting target nodes...")

## the scMEGA function 'SelectGenes' uses the human or mouse genome to connect genes to enhancers
## instead just use the previously calculated Peak to gene linkage 250000 distance from ArchR

# rename colnames so match output of SelectGenes
names(P2G)[names(P2G) == 'PeakID'] <- 'peak'
names(P2G)[names(P2G) == 'gene_name'] <- 'gene'
head(P2G)
dim(P2G)

# filter peak-enhancer interactions on GMs and also threshold correlation
df.p2g <- P2G %>% 
  dplyr::filter(FDR < 0.01) %>%
  dplyr::filter(gene %in% unique(c(gm_genes, df.tfs$tfs)))
dim(df.p2g)

# save final target node df
write.csv(df.p2g, file = paste0(temp_csv_path, "Target_genes_with_matched_enhancers.csv"), row.names = FALSE)

# plot heatmaps of target nodes and their connected peaks
trajRNA <- GetTrajectory_updated(obj.traj, assay = "RNA", trajectory.name = trajectory, 
                                 groupEvery = 2, slot = "data", smoothWindow = 7, log2Norm = TRUE)
trajATAC <- GetTrajectory_updated(obj.traj, assay = "ATAC", groupEvery = 2, 
                                  trajectory.name = trajectory, slot = "data", smoothWindow = 7, 
                                  log2Norm = TRUE)
trajATAC <- trajATAC[df.p2g$peak, ]
trajRNA <- trajRNA[df.p2g$gene, ]
ht <- CorrelationHeatmap_updated(trajectory1 = trajATAC, trajectory2 = trajRNA, 
                                 name1 = "Chromatin accessibility", name2 = "Gene expression", 
                                 labelTop1 = 50, labelTop2 = 50, labelRows1 = FALSE, labelRows2 = FALSE)

png(paste0(temp_plot_path, 'Genes_peaks_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
draw(ht)
graphics.off()

## how many peaks is each target node associated with?
npeaks <- as.data.frame(table(df.p2g$gene))
png(paste0(temp_plot_path, 'Target_node_nPeaks.png'), height = 10, width = 15, units = 'cm', res = 400)
print(hist(npeaks$Freq, breaks = 100))
graphics.off()
summary(npeaks$Freq)

############  TOTAL NODES NUMBERS

print("Total node numbers...")

## plot node numbers
df <- data.frame(
  SourceNodes = nrow(df.tfs),
  TargetNodes = length(unique(df.p2g$gene)),
  BothNodes = sum(unique(df.p2g$gene) %in% df.tfs$tfs)
)
png(paste0(temp_plot_path, 'Node_numbers.png'), height = 8, width = 30, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

## node numbers from each GM
df <- data.frame(
  GM12_nodes = sum(unique(df.p2g$gene) %in% GM12),
  GM14_nodes = sum(unique(df.p2g$gene) %in% GM14),
  GM13_nodes = sum(unique(df.p2g$gene) %in% GM13)
)
png(paste0(temp_plot_path, 'GM_node_numbers.png'), height = 8, width = 30, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

############################## Extract TF by peak matrix #######################################

temp_plot_path = "./plots/TF_by_peaks/"
dir.create(temp_plot_path, recursive = T)

print("Motif matrix...")

# peak by TF matrix to indicate presence of binding sites
motif.matching <- obj.traj@assays$ATAC@motifs@data
colnames(motif.matching) <- obj.traj@assays$ATAC@motifs@motif.names
motif.matching <- motif.matching[unique(df.p2g$peak), unique(df.tfs$tfs)]
dim(motif.matching)

# distribution of motifs by TF
n_hits_per_TF <- colSums(motif.matching)
png(paste0(temp_plot_path, 'Motif_hits_per_TF.png'), height = 8, width = 10, units = 'cm', res = 400)
hist(n_hits_per_TF, breaks = 100)
graphics.off()
summary(n_hits_per_TF)
print(n_hits_per_TF[order(n_hits_per_TF)])

# distribution of peaks by TF
n_hits_per_peak <- rowSums(motif.matching)
png(paste0(temp_plot_path, 'Motif_hits_per_peak.png'), height = 8, width = 10, units = 'cm', res = 400)
hist(n_hits_per_peak, breaks = 100)
graphics.off()
summary(n_hits_per_peak)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
print(paste0("mode: ", getmode(n_hits_per_peak)))

############################## Building quantiative GRN #######################################
# by using correlation between accessibility of selected TF targets and expression of selected genes over trajectory

temp_plot_path = "./plots/placodal_lineage/full_network/"
dir.create(temp_plot_path, recursive = T)

print("Building quantiative GRN...")

# GRN inference
df.tf.gene <- GetTFGeneCorrelation_updated(object = obj.traj,
                                    tf.use = df.tfs$tfs,
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar",
                                    gene.assay = "RNA",
                                    groupEvery = 2,
                                    trajectory.name = trajectory)

# check that all the selected TFs and genes are in the correlation matrix
if ( sum(unique(df.tfs$tf) %in% unique(df.tf.gene$tf)) != length(unique(df.tfs$tf)) ){
  stop("Not all selected TFs are in TF-gene correlation matrix!")
}
if ( sum(unique(df.p2g$gene) %in% unique(df.tf.gene$gene)) != length(unique(df.p2g$gene)) ){
  stop("Not all selected genes are in TF-gene correlation matrix!")
}

# save the correlation matrix
dim(df.tf.gene) # 189837 (ie 2403 genes x 79 TFs)
write.csv(df.tf.gene, file = paste0(temp_csv_path, "TF_to_gene_correlations.csv"), row.names = FALSE)

# plot TF-gene correlation heatmap
ht <- GRNHeatmap(df.tf.gene, tf.timepoint = df.tfs$time_point, km = 1)
png(paste0(temp_plot_path, 'TF_gene_corr_heatmap.png'), height = 30, width = 115, units = 'cm', res = 400)
ht
graphics.off()

############################## Building enhancer-based GRN #######################################
# To associate genes to TFs, we will use the peak-to-gene links and TF binding sites information
# if a gene is regulated by a peak AND this peak is bound by a TF, THEN we say this gene is a target of this TF

print("Building enhancer GRN...")

# take the TF-gene correlation, peak-TF binding prediction and peaks-to-genes linkage to build network
df.grn <- GetGRN(motif.matching = motif.matching,
                 df.cor = df.tf.gene,
                 df.p2g = df.p2g)
nrow(df.grn) # 148,113 interactions found
head(df.grn)

# check which tfs and genes in final network
if ( sum(unique(df.tfs$tf) %in% unique(df.grn$tf)) != length(unique(df.tfs$tf)) ){
  stop("Not all selected TFs are in TF-gene correlation matrix!")
}
print("How many selected genes didn't make it to final GRN: ")
print(table(unique(df.p2g$gene) %in% unique(df.grn$gene)))

############################## Final full GRN #######################################

print("Final full GRN...")

# full network numbers
df <- data.frame(
  nSource = length(unique(df.grn$tf)),
  nTarget = length(unique(df.grn$gene)),
  nBoth = sum(unique(df.grn$tf) %in% unique(df.grn$gene)),
  nInteractions = nrow(df.grn),
  nPositiveInteractions = length(which(df.grn$correlation > 0)),
  nNegativeInteractions = length(which(df.grn$correlation < 0))
           )
png(paste0(temp_plot_path, 'Network_all_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

## node numbers from each GM
df <- data.frame(
  GM12_nodes = sum(unique(df.p2g$gene) %in% GM12),
  GM14_nodes = sum(unique(df.p2g$gene) %in% GM14),
  GM13_nodes = sum(unique(df.p2g$gene) %in% GM13)
)
png(paste0(temp_plot_path, 'Unfiltered_GRN_GM_node_numbers.png'), height = 8, width = 30, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# add column for positive/negative correlations and abs correlation
df.grn <- df.grn %>%
  dplyr::mutate(direction = ifelse(correlation > 0, "P", "N")) %>%
  mutate(abscorr = abs(correlation))

# save full network (but shouldn't really use this)
write_tsv(df.grn, file = paste0(temp_csv_path, "GRN_initial.txt"))

############################## Filter full GRN #######################################

temp_plot_path = "./plots/placodal_lineage/filtered_network/"
dir.create(temp_plot_path, recursive = T)

print("Filtering final GRN..")

# filter network
df.grn <- df.grn %>%
  dplyr::filter(fdr < 0.01) # only keep significant correlations between TF binding and target gene expression (FDR < 0.01)

# filtered network numbers
df <- data.frame(
  nSource = length(unique(df.grn$tf)),
  nTarget = length(unique(df.grn$gene)),
  nBoth = sum(unique(df.grn$tf) %in% unique(df.grn$gene)),
  nInteractions = nrow(df.grn),
  nPositiveInteractions = length(which(df.grn$correlation > 0)),
  nNegativeInteractions = length(which(df.grn$correlation < 0))
)
png(paste0(temp_plot_path, 'Network_filtered_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

## node numbers from each GM
df <- data.frame(
  GM12_nodes = sum(unique(df.p2g$gene) %in% GM12),
  GM14_nodes = sum(unique(df.p2g$gene) %in% GM14),
  GM13_nodes = sum(unique(df.p2g$gene) %in% GM13)
)
png(paste0(temp_plot_path, 'Filtered_GRN_GM_node_numbers.png'), height = 8, width = 30, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# plot venn of nodes
myCol <- brewer.pal(2, "Pastel2")
venn.diagram(
  x = list(unique(df.grn$tf), unique(df.grn$gene)[!is.na(unique(df.grn$gene))]),
  category.names = c("Source node", "Target node"),
  filename = paste0(temp_plot_path, 'Nodes.png'),
  output=TRUE, disable.logging = TRUE,
  # Output features
  imagetype="png",
  height = 600, 
  width = 600, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0,
  cat.fontface = "bold"
)

# plot piechart of interactions
pie <- c(length(which(df.grn$correlation > 0)), length(which(df.grn$correlation < 0)))
png(paste0(temp_plot_path, 'Interactions_pie.png'), height = 8, width = 8, units = 'cm', res = 400)
pie(pie, labels = c("Positive interactions", "Negative interactions"))
graphics.off()

# plot TF-gene correlation heatmap
df.tf.gene.subset <- df.tf.gene %>%
  dplyr::filter(tf %in% unique(df.grn$tf))
df.tfs.subset <- df.tfs %>%
  dplyr::filter(tfs %in% unique(df.grn$tf))
ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = 1)

png(paste0(temp_plot_path, 'TF_gene_corr_heatmap.png'), height = 30, width = 115, units = 'cm', res = 400)
ht
graphics.off()

## how many peaks supports each interaction
summary(df.grn$n_peaks)
png(paste0(temp_plot_path, 'Interactions_nPeaks.png'), height = 10, width = 15, units = 'cm', res = 400)
print(hist(df.grn$n_peaks, breaks = 100))
graphics.off()

png(paste0(temp_plot_path, 'Interactions_nPeaks_log10.png'), height = 10, width = 15, units = 'cm', res = 400)
print(hist(log10(df.grn$n_peaks), breaks = 100))
graphics.off()

subset <- df.grn %>% dplyr::filter(n_peaks < 40)
png(paste0(temp_plot_path, 'Interactions_nPeaks_under_40.png'), height = 10, width = 15, units = 'cm', res = 400)
print(hist(subset$n_peaks, breaks = 100))
graphics.off()

# save filtered network
write_tsv(df.grn, file = paste0(temp_csv_path, "GRN_filtered.txt"))

############################## Subset GRN to only include positive TFs and all target nodes they interact with #######################################
# only keep source nodes with overall positive correlate with target nodes

temp_plot_path = "./plots/placodal_lineage/filtered_network_pos_corr/"
dir.create(temp_plot_path, recursive = T)

# filter correlation matrix
summarized_tf_gene_corr <- df.tf.gene %>%
  group_by(tf) %>%
  summarize(mean_correlation = mean(correlation)) %>%
  arrange(desc(mean_correlation)) %>%
  filter(mean_correlation > 0) # filter

# how many TFs remaining
length(summarized_tf_gene_corr$tf) # 152

# filter grn to only include these positive corr TFs
sum(summarized_tf_gene_corr$tf %in% df.grn$tf) # 2 have already been removed
df.grn.pos <- df.grn[df.grn$tf %in% summarized_tf_gene_corr$tf, ]
nrow(df.grn.pos) # 33529

# TF network numbers
df <- data.frame(
  nSource = length(unique(df.grn.pos$tf)),
  nTarget = length(unique(df.grn.pos$gene)),
  nBoth = sum(unique(df.grn.pos$tf) %in% unique(df.grn.pos$gene)),
  nInteractions = nrow(df.grn.pos),
  nPositiveInteractions = length(which(df.grn.pos$correlation > 0)),
  nNegativeInteractions = length(which(df.grn.pos$correlation < 0))
)
png(paste0(temp_plot_path, 'Network_filtered_positive_source_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

## node numbers from each GM
df <- data.frame(
  GM12_nodes = sum(unique(df.p2g$gene) %in% GM12),
  GM14_nodes = sum(unique(df.p2g$gene) %in% GM14),
  GM13_nodes = sum(unique(df.p2g$gene) %in% GM13)
)
png(paste0(temp_plot_path, 'Pos_corr_GRN_GM_node_numbers.png'), height = 8, width = 30, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# plot venn of nodes
myCol <- brewer.pal(2, "Pastel2")
venn.diagram(
  x = list(unique(df.grn.pos$tf), unique(df.grn.pos$gene)[!is.na(unique(df.grn.pos$gene))]),
  category.names = c("Source node", "Target node"),
  filename = paste0(temp_plot_path, 'Nodes_pos_corr.png'),
  output=TRUE, disable.logging = TRUE,
  # Output features
  imagetype="png",
  height = 600,
  width = 600,
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[1:2],
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0,
  cat.fontface = "bold"
)

# plot piechart of interactions
pie <- c(length(which(df.grn.pos$correlation > 0)), length(which(df.grn.pos$correlation < 0)))
png(paste0(temp_plot_path, 'Interactions_pie_pos_corr.png'), height = 8, width = 8, units = 'cm', res = 400)
pie(pie, labels = c("Positive interactions", "Negative interactions"))
graphics.off()

# plot TF-gene correlation heatmap
df.tf.gene.subset <- df.tf.gene %>%
  dplyr::filter(tf %in% unique(df.grn.pos$tf))
df.tfs.subset <- df.tfs %>%
  dplyr::filter(tfs %in% unique(df.grn.pos$tf))
ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = 1)

png(paste0(temp_plot_path, 'TF_gene_corr_heatmap.png'), height = 30, width = 115, units = 'cm', res = 400)
ht
graphics.off()

# save network
write_tsv(df.grn.pos, file = paste0(temp_csv_path, "GRN_filtered_pos_corr.txt"))

############################## Subset filtered GRN to only include TFs #######################################
# only keep nodes which are known TFs (they might not regulate any other nodes in this network)

temp_plot_path = "./plots/placodal_lineage/filtered_network_TFs/"
dir.create(temp_plot_path, recursive = T)

# filter grn to only include known TFs
df.grn.TFs <- df.grn[df.grn$gene %in% TF_names, ]
nrow(df.grn.TFs) # 3,587 interactions

# TF network numbers
df <- data.frame(
  nSource = length(unique(df.grn.TFs$tf)),
  nTarget = length(unique(df.grn.TFs$gene)),
  nBoth = sum(unique(df.grn.TFs$tf) %in% unique(df.grn.TFs$gene)),
  nInteractions = nrow(df.grn.TFs),
  nPositiveInteractions = length(which(df.grn.TFs$correlation > 0)),
  nNegativeInteractions = length(which(df.grn.TFs$correlation < 0))
)
png(paste0(temp_plot_path, 'Network_filtered_TFs_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# save network
write_tsv(df.grn.TFs, file = paste0(temp_csv_path, "GRN_filtered_TFs.txt"))

############################## Subset positive corr GRN to only include TFs #######################################
# only keep TF nodes which overall positive correlate with target nodes

temp_plot_path = "./plots/placodal_lineage/filtered_network_pos_corr_TFs/"
dir.create(temp_plot_path, recursive = T)

# filter grn to only include thes positive TFs
df.grn.pos.TFs <- df.grn.pos[df.grn.pos$gene %in% TF_names, ]
nrow(df.grn.pos.TFs) # 6,707 interactions

# TF network numbers
df <- data.frame(
  nSource = length(unique(df.grn.pos.TFs$tf)),
  nTarget = length(unique(df.grn.pos.TFs$gene)),
  nBoth = sum(unique(df.grn.pos.TFs$tf) %in% unique(df.grn.pos.TFs$gene)),
  nInteractions = nrow(df.grn.pos.TFs),
  nPositiveInteractions = length(which(df.grn.pos.TFs$correlation > 0)),
  nNegativeInteractions = length(which(df.grn.pos.TFs$correlation < 0))
)
png(paste0(temp_plot_path, 'Network_filtered_positive_source_TFs_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# save network
write_tsv(df.grn.pos.TFs, file = paste0(temp_csv_path, "GRN_filtered_pos_corr_TFs.txt"))

############################## Make node metadata table to import into Cytoscape #######################################

print("Making node metdata...")

# make df with all nodes
node_metadata <- data.frame(
  node = unique(c(df.grn$tf, df.grn$gene))
)

# add whether node is a TF or not
node_metadata <- node_metadata %>%
  dplyr::mutate(TF = ifelse(node %in% TF_names, "TF", "G"))

# add ordering of source nodes
df.tfs <- df.tfs[order(df.tfs$time_point), ]
df.timepoint <- df.tfs %>% dplyr::select(c(tfs, time_point))
node_metadata <- merge(node_metadata, df.timepoint, by.x = "node", by.y = "tfs", all = TRUE)

# add which GM node is part of (GMs calculated on ss8)
df.gm <- data.frame(
  node <- c(GM12, GM14, GM13),
  ss8_GM <- c( rep("GM12", length(GM12)), rep("GM14", length(GM14)), rep("GM13", length(GM13)) )
)
colnames(df.gm) <- c("node", "ss8_GM")
node_metadata <- merge(node_metadata, df.gm, by = "node", all = TRUE)

# check final
head(node_metadata)
nrow(node_metadata) # 1,327

# save
write_tsv(node_metadata, file = paste0(temp_csv_path, "Node_metadata.txt"))