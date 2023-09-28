#!/usr/bin/env Rscript

print("Run scMEGA on paired RNA/ATAC object to infer GRN")

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
    data_path = "./input/rds_files/"
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

SelectTFs_updated <- function (object, tf.assay = "chromvar", rna.assay = "RNA", atac.assay = "ATAC", 
                               trajectory.name = "Trajectory", groupEvery = 1, p.cutoff = 0.01, 
                               cor.cutoff = 0.3, return.heatmap = TRUE) 
{
  trajMM <- GetTrajectory_updated(object, assay = tf.assay, trajectory.name = trajectory.name, 
                                  groupEvery = groupEvery, slot = "data", smoothWindow = 7, 
                                  log2Norm = FALSE)
  rownames(trajMM) <- object@assays[[atac.assay]]@motifs@motif.names
  trajRNA <- GetTrajectory_updated(object, assay = rna.assay, trajectory.name = trajectory.name, 
                                   groupEvery = groupEvery, slot = "data", smoothWindow = 7, 
                                   log2Norm = TRUE)
  df.cor <- GetCorrelation(trajMM, trajRNA) # correlation between TF expression RNAseq and binding activity ATACseq
  df.cor <- df.cor[df.cor$adj_p < p.cutoff & df.cor$correlation > 
                     cor.cutoff, ]
  matMM <- suppressMessages(TrajectoryHeatmap(trajMM, varCutOff = 0, 
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
  groupMatRNA <- suppressMessages(TrajectoryHeatmap(trajRNA, 
                                                    varCutOff = var.cutoff.gene, pal = paletteContinuous(set = "horizonExtra"), 
                                                    limits = c(-2, 2), returnMatrix = TRUE))
  groupMatATAC <- suppressMessages(TrajectoryHeatmap(trajATAC, 
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
  tf_activity <- suppressMessages(TrajectoryHeatmap(trajMM, 
                                                    varCutOff = 0, pal = paletteContinuous(set = "solarExtra"), 
                                                    limits = c(-2, 2), name = "TF activity", returnMatrix = TRUE))
  gene_expression <- suppressMessages(TrajectoryHeatmap(trajRNA, 
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

### New function to plot network using igraph
# colours nodes by timepoint
# colours edges by positive/negative and thickness is by correlation strength
# edit so can optionally colour nodes

PlotTFNetwork <- function(df.grn, tfs.timepoint){
  # extract data to plot
  links <- df.grn %>% 
    dplyr::select(c("tf", "gene", "correlation")) %>%
    dplyr::mutate(abs_corr = abs(correlation))
  
  # timepoint expression information of TFs
  nodes <- rownames_to_column(as.data.frame(tfs.timepoint), var = "tf")
  length(unique(nodes$tf))
  nodes <- nodes %>% 
    dplyr::filter(tf %in% c(links$tf, links$gene)) %>% # remove TFs/genes with no links
    mutate(tfs.timepoint = round(tfs.timepoint, digits = 0))
  length(unique(nodes$tf))
  
  # colours for latent time
  cols <- data.frame(col = paletteContinuous(set = "beach", n = 100),
                     val = seq(1:100))
  cols_matched <- merge(cols, nodes, by.x = "val", by.y = "tfs.timepoint")
  
  # Turn it into igraph object
  network <- graph_from_data_frame(d = links, vertices = nodes, directed = T)
  
  # colours for positive or negative interactions
  E(network)$color <- ifelse(E(network)$correlation > 0,'red','blue')
  
  # make plot
  p <- plot(network, edge.arrow.size=E(network)$abs_corr*0.75, edge.width = E(network)$abs_corr*5,
            vertex.color = cols_matched$col,
            vertex.label.color = "black", vertex.label.family = "Helvetica", vertex.label.cex = 0.8,
            layout = layout_nicely(network))
  
  return(p)
}



############################## Read in seurat object #######################################

print("reading in data...")

## read in paired seurat object
obj.pair <- readRDS(paste0(data_path, "paired_object_chromvar.RDS"))
obj.pair

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

saveRDS(obj.pair, paste0(rds_path, "paired_object_chromvar.RDS"), compress = FALSE)


############################## Extract all known TF names #######################################

TF_names <- obj.pair@assays$ATAC@motifs@motif.names
names(TF_names) <- NULL
TF_names <- unlist(TF_names)
length(TF_names) # 746 known TFs
fileConn <- file(paste0(rds_path, "Known_TF_names.txt"))
writeLines(TF_names, fileConn)
close(fileConn)

######################################################################################
##############################    PLACODAL     #######################################
######################################################################################

temp_plot_path = "./plots/placodal_lineage/"
dir.create(temp_plot_path, recursive = T)
temp_csv_path = "./csv_files/placodal_lineage/"
dir.create(temp_csv_path, recursive = T)
temp_lineage_plot_path = "./plots/placodal_lineage/lineage_dynamics_plots/"
dir.create(temp_lineage_plot_path, recursive = T)

# set trajectory
trajectory <- "lineage_placodal_probability"

# remove cells that have a placodal probability of less than 25%
obj.traj <- subset(obj.pair, subset = lineage_placodal_probability > 25)
obj.traj

############################## Select nodes #######################################

############  SOURCE NODES

print("Selecting source nodes...")

# select the source nodes - dont worry about correlating activity and gex
res <- SelectTFs_updated(object = obj.traj, trajectory.name = trajectory, return.heatmap = TRUE,
                         groupEvery = 2,
                         p.cutoff = 1, cor.cutoff = 0)

# save the selected source nodes
df.tfs <- res$tfs
write.csv(df.tfs, file = paste0(temp_csv_path, "TF_correlations.csv"), row.names = FALSE)

# how many source nodes are there
nrow(df.tfs) # 155
unique(df.tfs$tfs)
# [1] "SOX2"    "ARNT2"   "POU4F3"  "POU4F2"  "HIF1A"   "PAX7"    "GABPA"  
# [8] "POU4F1"  "OLIG2"   "SMAD5"   "ZIC5"    "MAFK"    "FOS"     "TBX3"   
# [15] "TBX20"   "TEF"     "OVOL2"   "TCF7"    "BHLHE40" "MAFG"    "HES5"   
# [22] "ZIC3"    "ASCL1"   "MSANTD3" "TBXT"    "HLF"     "PROX1"   "ZIC1"   
# [29] "RUNX3"   "ZBTB26"  "ETV3"    "JUN"     "RBPJ"    "MEF2A"   "CREB3L1"
# [36] "NFKB1"   "TGIF1"   "ELF1"    "TGIF2"   "DMRTA2"  "LIN54"   "ESR2"   
# [43] "NKX2-3"  "HMBOX1"  "STAT1"   "TP63"    "SOX21"   "NR5A1"   "NR2C1"  
# [50] "ONECUT1" "JUND"    "HES6"    "HOXA9"   "IRF6"    "NFE2L1"  "NR3C2"  
# [57] "PITX1"   "HESX1"   "BACH2"   "EBF1"    "EBF3"    "ZEB1"    "SNAI1"  
# [64] "TCF3"    "TBR1"    "BACH1"   "EOMES"   "TBX2"    "HOXD4"   "NFE2"   
# [71] "ETV4"    "TFAP2A"  "TFAP2B"  "GATA4"   "TFAP2C"  "GATA2"   "SNAI2"  
# [78] "GATA3"   "GATA5"   "GATA6"   "TFAP2E"  "HOXA2"   "SP1"     "MEOX1"  
# [85] "DRGX"    "HOXA1"   "GBX2"    "EMX1"    "HOXA5"   "DLX6"    "EMX2"   
# [92] "EVX1"    "GBX1"    "KLF6"    "DLX5"    "EN2"     "MIXL1"   "LMX1A"  
# [99] "KLF10"   "KLF11"   "RAX2"    "LHX9"    "FOXP2"   "MSX2"    "BARX1"  
# [106] "PRRX2"   "LMX1B"   "MGA"     "KLF5"    "KLF3"    "ZBTB6"   "REL"    
# [113] "JDP2"    "PPARG"   "SREBF1"  "NFKB2"   "KLF15"   "SREBF2"  "CUX1"   
# [120] "TEAD1"   "TEAD4"   "TEAD3"   "SIX1"    "SP2"     "TFCP2"   "FOXK2"  
# [127] "HOXA7"   "HOXD8"   "FOXP3"   "FEV"     "ATF3"    "ELK4"    "FOXK1"  
# [134] "FOXG1"   "CREM"    "ATF2"    "ZBTB7A"  "FOXN3"   "BATF"    "THRB"   
# [141] "BATF3"   "PITX2"   "E2F1"    "PLAG1"   "CEBPG"   "GMEB2"   "EGR1"   
# [148] "EHF"     "RELA"    "RARA"    "PRDM4"   "ZNF652"  "HEY1"    "ETV1"   
# [155] "NKX2-5" 

# how do these source nodes overlap with previous predictions
Alex_TFs <- c("GATA2", "SIX1", "DLX5", "TFAP2A", "IRF6", "DLX6", 
              "GATA3", "TFAP2C", "EPAS1", "PITX1", "PITX2", "RARB", "HRKB1", "IRX1" )
Chromvar_TFs <- c("TEAD3", "TEAD4", "TEAD2", "TEAD1", "TFAP2A", "TFAP2C", 
                  "Rbpjl", "Ptf1a", "SNAI2", "SNAI1", "SNAI3", "TCF12")   

print("Overlap with Alex's TFs:")
print(Alex_TFs[Alex_TFs %in% unique(df.tfs$tfs)])
# [1] "GATA2"  "SIX1"   "DLX5"   "TFAP2A" "IRF6"   "DLX6"   "GATA3"  "TFAP2C"
# [9] "PITX1"  "PITX2" 

print("Overlap with ChromVar TFs:")
print(Chromvar_TFs[Chromvar_TFs %in% unique(df.tfs$tfs)])
# "TEAD3"  "TEAD4"  "TEAD1"  "TFAP2A" "TFAP2C" "SNAI2"  "SNAI1" 

# plot TF activity dynamics across trajectory
ht <- res$heatmap
png(paste0(temp_plot_path, 'TFs_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
draw(ht)
graphics.off()

# sense check that all source nodes are also a known TF
if (sum(df.tfs$tfs %in% TF_names) != length(df.tfs$tfs)){
  stop("Problem! Not all source nodes are known TFs!")
}

############  TARGET NODES

print("Selecting target nodes...")

# select target nodes that vary with the trajectory AND correlate with peaks
res <- SelectGenes_updated(obj.traj, trajectory.name = trajectory, groupEvery = 2,
                           var.cutoff.gene = 0.7, # how much gene expression has to vary across trajectory
                           cor.cutoff = 0.7, fdr.cutoff = 1e-04) # how much peaks and genes need to correlate

# save target nodes
df.p2g.var <- res$p2g
write.csv(df.p2g.var, file = paste0(temp_csv_path, "Variable_genes_with_matched_enhancers.csv"), row.names = FALSE)

# plot the dynamics of these genes across trajectory
ht <- res$heatmap
png(paste0(temp_plot_path, 'Genes_peaks_variable_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
draw(ht)
graphics.off()

# select genes correlate with peaks (don't have to be variable across trajectory)
res <- SelectGenes_updated(obj.traj, trajectory.name = trajectory, groupEvery = 2,
                           var.cutoff.gene = 0.01, # how much gene expression has to vary across trajectory
                           cor.cutoff = 0.7, fdr.cutoff = 1e-04) # how much peaks and genes need to correlate

# save target nodes
df.p2g <- res$p2g
write.csv(df.p2g, file = paste0(temp_csv_path, "Genes_with_matched_enhancers.csv"), row.names = FALSE)

# plot the dynamics of these genes across trajectory
ht <- res$heatmap
png(paste0(temp_plot_path, 'Genes_peaks_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
draw(ht)
graphics.off()

# for final target node df, combine the variable ones with all TFs from non variable one
df.p2g.tfs <- df.p2g %>% 
  dplyr::filter(gene %in% TF_names) %>% # filter p2g to only include TFs
  dplyr::filter(!gene %in% df.p2g.var$gene)
nrow(df.p2g.tfs) # 1,307 new interactions added
length(unique(df.p2g.tfs$gene)) # 67 new nodes added

df.p2g.final <- rbind(df.p2g.var, df.p2g.tfs)
nrow(df.p2g.final) # 74,125 final interactions

# save final target node df
write.csv(df.p2g.final, file = paste0(temp_csv_path, "Final_target_genes_with_matched_enhancers.csv"), row.names = FALSE)

# plot heatmaps of target nodes and their connected peaks
trajRNA <- GetTrajectory_updated(obj.traj, assay = "RNA", trajectory.name = trajectory, 
                                 groupEvery = 2, slot = "data", smoothWindow = 7, 
                                 log2Norm = TRUE)
trajATAC <- GetTrajectory_updated(obj.traj, assay = "ATAC", groupEvery = 2, 
                                  trajectory.name = trajectory, slot = "data", smoothWindow = 7, 
                                  log2Norm = TRUE)
trajATAC <- trajATAC[df.p2g.final$peak, ]
trajRNA <- trajRNA[df.p2g.final$gene, ]
ht <- suppressMessages(CorrelationHeatmap(trajectory1 = trajATAC, trajectory2 = trajRNA, 
                                          name1 = "Chromatin accessibility", name2 = "Gene expression", 
                                          labelTop1 = 50, labelTop2 = 50, labelRows1 = FALSE, labelRows2 = FALSE))

png(paste0(temp_plot_path, 'Genes_peaks_heatmap_final.png'), height = 30, width = 45, units = 'cm', res = 400)
draw(ht)
graphics.off()

## plot node numbers
df <- data.frame(
  SourceNodes = nrow(df.tfs),
  P2G_genes = length(unique(df.p2g$gene)),
  P2G_genes_var = length(unique(df.p2g.var$gene)),
  P2G_genes_tfs = sum(unique(df.p2g$gene) %in% TF_names),
  P2G_genes_var_tfs = sum(unique(df.p2g.var$gene) %in% TF_names),
  TargetNodes = length(unique(df.p2g.final$gene)),
  BothNodes = sum(unique(df.p2g.final$gene) %in% df.tfs$tfs)
)
png(paste0(temp_plot_path, 'Node_numbers.png'), height = 8, width = 30, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

############################## Building quantiative GRN #######################################
# by using correlation between accessibility of selected TF targets and expression of selected genes over trajectory

print("Building quantiative GRN...")

# GRN inference
df.tf.gene <- GetTFGeneCorrelation_updated(object = obj.traj,
                                    tf.use = df.tfs$tfs,
                                    gene.use = unique(df.p2g.final$gene),
                                    tf.assay = "chromvar",
                                    gene.assay = "RNA",
                                    groupEvery = 2,
                                    trajectory.name = trajectory)

# check that all the selected TFs and genes are in the correlation matrix
if ( sum(unique(df.tfs$tf) %in% unique(df.tf.gene$tf)) != length(unique(df.tfs$tf)) ){
  stop("Not all selected TFs are in TF-gene correlation matrix!")
}
if ( sum(unique(df.p2g.final$gene) %in% unique(df.tf.gene$gene)) != length(unique(df.p2g.final$gene)) ){
  stop("Not all selected genes are in TF-gene correlation matrix!")
}

# save the correlation matrix
dim(df.tf.gene) # 189837 (ie 2403 genes x 79 TFs)
write.csv(df.tf.gene, file = paste0(temp_csv_path, "TF_to_gene_correlations.csv"), row.names = FALSE)

# plot TF-gene correlation heatmap
ht <- GRNHeatmap(df.tf.gene, tf.timepoint = df.tfs$time_point, km = 1)
png(paste0(temp_plot_path, 'TF_gene_corr_heatmap.png'), height = 30, width = 60, units = 'cm', res = 400)
ht
graphics.off()

############################## Building enhancer-based GRN #######################################
# To associate genes to TFs, we will use the peak-to-gene links and TF binding sites information
# if a gene is regulated by a peak AND this peak is bound by a TF, THEN we say this gene is a target of this TF

print("Building enhancer GRN...")

# peak by TF matrix to indicate precence of binding sites
motif.matching <- obj.pair@assays$ATAC@motifs@data
colnames(motif.matching) <- obj.pair@assays$ATAC@motifs@motif.names
motif.matching <- motif.matching[unique(df.p2g$peak), unique(df.tfs$tf)]

# take the TF-gene correlation, peak-TF binding prediction and peaks-to-genes linkage to build network
df.grn <- GetGRN(motif.matching = motif.matching,
                 df.cor = df.tf.gene,
                 df.p2g = df.p2g.final)
nrow(df.grn) # 148,113 interactions found

# check which tfs and genes in final network
if ( sum(unique(df.tfs$tf) %in% unique(df.grn$tf)) != length(unique(df.tfs$tf)) ){
  stop("Not all selected TFs are in TF-gene correlation matrix!")
}
print("How many selected genes didn't make it to final GRN: ")
print(table(unique(df.p2g.final$gene) %in% unique(df.grn$gene)))

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

# add column for positive/negative correlations and abs correlation
df.grn <- df.grn %>%
  dplyr::mutate(direction = ifelse(correlation > 0, "P", "N")) %>%
  mutate(abscorr = abs(correlation))

# save full network (but shouldn't really use this)
write.csv(df.grn, file = paste0(temp_csv_path, "GRN_unfiltered.csv"), row.names = FALSE)

############################## Filter full GRN #######################################

print("Filtering final GRN..")

# distribution of nPeaks
summary(df.grn$n_peaks)
png(paste0(temp_plot_path, 'nPeaks_initial_distribution.png'), height = 8, width = 12, units = 'cm', res = 400)
hist(df.grn$n_peaks, breaks = 80)
graphics.off()

# filter network
df.grn <- df.grn %>%
  dplyr::filter(fdr < 0.01) %>% # only keep significant correlations between TF binding and target gene expression (FDR < 0.01)
  dplyr::filter(n_peaks > 3) # only keep interactions where at least 3 connected peaks have the TF binding site

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

# save filtered network
write.csv(df.grn, file = paste0(temp_csv_path, "GRN_filtered.csv"), row.names = FALSE)
write_tsv(df.grn, file = paste0(temp_csv_path, "GRN_filtered.txt"))

############################## Subset GRN to only include source nodes #######################################
# only keep nodes that regulate other nodes

print("Making source GRN...")

# filter grn to only include source TFs
df.grn.source <- df.grn[df.grn$gene %in% df.grn$tf, ]
nrow(df.grn.source) # 1872 interactions

# source network numbers
df <- data.frame(
  nSource = length(unique(df.grn.source$tf)),
  nTarget = length(unique(df.grn.source$gene)),
  nBoth = sum(unique(df.grn.source$tf) %in% unique(df.grn.source$gene)),
  nInteractions = nrow(df.grn.source),
  nPositiveInteractions = length(which(df.grn.source$correlation > 0)),
  nNegativeInteractions = length(which(df.grn.source$correlation < 0))
)
png(paste0(temp_plot_path, 'Network_filtered_source_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# save network
write.csv(df.grn.source, file = paste0(temp_csv_path, "GRN_filtered_source_nodes.csv"), row.names = FALSE)
write_tsv(df.grn.source, file = paste0(temp_csv_path, "GRN_filtered_source_nodes.txt"))

############################## Subset GRN to only include TFs #######################################
# only keep nodes which are known TFs (they might not regulate any other nodes in this network)

print("Making TF GRN...")

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
write_tsv(df.grn.TFs, file = paste0(temp_csv_path, "GRN_filtered_TF_nodes.txt"))

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

# check final
head(node_metadata)
nrow(node_metadata) # 1,327

# save
write.csv(node_metadata, file = paste0(temp_csv_path, "Node_metadata.csv"), row.names = FALSE)
write_tsv(node_metadata, file = paste0(temp_csv_path, "Node_metadata.txt"))

############################## Plot GRN maps using scMEGA #######################################

print("Plotting GRNs with scMEGA...")

# set order for cols
tfs.timepoint <- df.tfs$time_point
names(tfs.timepoint) <- df.tfs$tfs

# plot whole GRN
png(paste0(temp_plot_path, 'Network_all_importance_MEGA.png'), height = 10, width = 35, units = 'cm', res = 400)
p <- GRNPlot(df.grn,
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 42,
             plot.importance = TRUE,
             min.importance = 2,
             remove.isolated = FALSE)
graphics.off()

png(paste0(temp_plot_path, 'Network_all_MEGA.png'), height = 30, width = 45, units = 'cm', res = 400)
print(p)
graphics.off()

# only plot TFs
p <- GRNPlot(df.grn,
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             plot.importance = TRUE,
             genes.use = df.tfs$tfs,
             remove.isolated = TRUE)

png(paste0(temp_plot_path, 'Network_TFs_MEGA.png'), height = 30, width = 45, units = 'cm', res = 400)
print(p)
graphics.off()`

############################## Plot whole GRN map using igraphs #######################################

print("Plotting GRNs with igraph...")

# Plot source network - its different everytime, cant seem to fix randomness
png(paste0(temp_plot_path, 'Network_TFs_igraph.png'), height = 30, width = 45, units = 'cm', res = 400)
PlotTFNetwork(df.grn.source, tfs.timepoint)
graphics.off()


# ############################## Plot lineage dynamics #######################################

# print("Plotting lineage dynamics...")

# obj.temp <- AddTargetAssay_updated(obj.traj, df.grn = df.grn)
# 
# # plot dynamics of each TF in network
# for (TF in df.tfs$tfs){
#   print(TF)
# 
#   png(paste0(temp_lineage_plot_path, TF, '_dynamics_plot.png'), height = 15, width = 20, units = 'cm', res = 400)
#   print(PseudotimePlot_updated(obj.temp, trajectory.name = trajectory,
#                          tf.use = TF))
#   graphics.off()
# }
