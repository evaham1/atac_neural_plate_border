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
  if (assay == "chromvar"){
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
                                          groupEvery = 1, var.cutoff.gene = NULL) 
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
                                                        varCutOff = var.cutoff.gene, pal = paletteContinuous(set = "solarExtra"), 
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

############################## Read in seurat object #######################################

print("reading in data...")

## read in paired seurat object
obj.pair <- readRDS(paste0(data_path, "paired_object_chromvar.RDS"))
obj.pair

############################## Create trajectories from lineage probabilities #######################################

head(obj.pair@meta.data)

# then need them to be between 0 and 100
summary(obj.pair@meta.data$rna_lineage_placodal_probability)
obj.pair@meta.data$rna_lineage_placodal_probability <- obj.pair@meta.data$rna_lineage_placodal_probability * 100
summary(obj.pair@meta.data$rna_lineage_placodal_probability)
# hist(obj.pair@meta.data$rna_lineage_placodal_probability, breaks = 100)

obj.pair@meta.data$rna_lineage_NC_probability <- obj.pair@meta.data$rna_lineage_NC_probability * 100
obj.pair@meta.data$rna_lineage_neural_probability <- obj.pair@meta.data$rna_lineage_neural_probability * 100

############################## Save seurat object #######################################

saveRDS(obj.pair, paste0(rds_path, "paired_object_chromvar.RDS"), compress = FALSE)

######################################################################################
##############################    PLACODAL     #######################################
######################################################################################

############################## Filter genes and TFs #######################################

temp_plot_path = "./plots/placodal_lineage/"
dir.create(temp_plot_path, recursive = T)
temp_csv_path = "./csv_files/placodal_lineage/"
dir.create(temp_csv_path, recursive = T)

# select genes that vary with the trajectory and correlate with peaks
res <- SelectGenes_updated(obj.pair, trajectory.name = "rna_lineage_placodal_probability", groupEvery = 2,
                           var.cutoff.gene = 0.3, # how much gene expression has to vary across trajectory
                           cor.cutoff = 0.3, fdr.cutoff = 1e-04) # how much peaks and genes need to correlate

# save the most variable genes
df.p2g <- res$p2g
dim(df.p2g)
length(unique(df.p2g$gene))
write.csv(df.p2g, file = paste0(temp_csv_path, "Variable_genes_and_matched_enhancers.csv"), row.names = FALSE)

# plot the dynamics of these genes across trajectory
ht <- res$heatmap
png(paste0(temp_plot_path, 'Genes_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
draw(ht)
graphics.off()

# select the TFs that correlate in 'activity' as defined by chromvar and expression
res <- SelectTFs_updated(object = obj.pair, trajectory.name = "rna_lineage_placodal_probability", return.heatmap = TRUE,
                 groupEvery = 2,
                 p.cutoff = 0.01, cor.cutoff = 0.3)

# save the selected TFs
df.tfs <- res$tfs
dim(df.tfs)
write.csv(df.tfs, file = paste0(temp_csv_path, "TF_correlations.csv"), row.names = FALSE)
unique(df.tfs$tfs)

# plot TF activity dynamics across trajectory
ht <- res$heatmap
png(paste0(temp_plot_path, 'TFs_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
draw(ht)
graphics.off()

############################## Building quantiative GRN #######################################
# by using correlation between accessibility of selected TF targets and expression of selected genes over trajectory

# GRN inference
df.tf.gene <- GetTFGeneCorrelation_updated(object = obj.pair, 
                                    tf.use = df.tfs$tfs, 
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar", 
                                    gene.assay = "RNA",
                                    var.cutoff.gene = 0.3,
                                    groupEvery = 2,
                                    trajectory.name = "rna_lineage_placodal_probability")

# check that all the selected TFs and genes are in the correlation matrix
if ( sum(unique(df.tfs$tf) %in% unique(df.tf.gene$tf)) != length(unique(df.tfs$tf)) ){
  stop("Not all selected TFs are in TF-gene correlation matrix!")
}
if ( sum(unique(df.p2g$gene) %in% unique(df.tf.gene$gene)) != length(unique(df.p2g$gene)) ){
  stop("Not all selected genes are in TF-gene correlation matrix!")
}

# save the correlation matrix
dim(df.tf.gene)
write.csv(df.tf.gene, file = paste0(temp_csv_path, "TF_to_gene_correlations.csv"), row.names = FALSE)

# plot TF-gene correlation heatmap
ht <- GRNHeatmap(df.tf.gene, tf.timepoint = df.tfs$time_point, km = 1)
png(paste0(temp_plot_path, 'TF_gene_corr_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
ht
graphics.off()

############################## Building enhancer-based GRN #######################################
# To associate genes to TFs, we will use the peak-to-gene links and TF binding sites information 
# if a gene is regulated by a peak AND this peak is bound by a TF, THEN we say this gene is a target of this TF

# peak by TF matrix to indicate prescence of binding sites
motif.matching <- obj.pair@assays$ATAC@motifs@data
colnames(motif.matching) <- obj.pair@assays$ATAC@motifs@motif.names
motif.matching <- motif.matching[unique(df.p2g$peak), unique(df.tfs$tf)]

# take the TF-gene correlation, peak-TF binding prediction and peaks-to-genes linkage to build network
df.grn <- GetGRN(motif.matching = motif.matching, 
                 df.cor = df.tf.gene, 
                 df.p2g = df.p2g)
dim(df.grn)

# check which tfs and genes in final network
if ( sum(unique(df.tfs$tf) %in% unique(df.grn$tf)) != length(unique(df.tfs$tf)) ){
  stop("Not all selected TFs are in TF-gene correlation matrix!")
}
print("How many selected genes didn't make it to final GRN: ")
print(table(unique(df.p2g$gene) %in% unique(df.grn$gene)))

# save network
write.csv(df.grn, file = paste0(temp_csv_path, "GRN_data.csv"), row.names = FALSE)


# 
# ######################################################################################
# ##############################    NC     #######################################
# ######################################################################################
# 
# ############################## Filter genes and TFs #######################################
# 
# temp_plot_path = "./plots/NC_lineage/"
# dir.create(temp_plot_path, recursive = T)
# temp_csv_path = "./csv_files/NC_lineage/"
# dir.create(temp_csv_path, recursive = T)
# 
# # select genes that vary with the trajectory and correlate with peaks
# res <- SelectGenes_updated(obj.pair, trajectory.name = "rna_lineage_NC_probability", groupEvery = 2,
#                            var.cutoff.gene = 0.3, # how much gene expression has to vary across trajectory
#                            cor.cutoff = 0.3, fdr.cutoff = 1e-04) # how much peaks and genes need to correlate
# 
# # save the most variable genes
# df.p2g <- res$p2g
# dim(df.p2g)
# write.csv(df.p2g, file = paste0(temp_csv_path, "Variable_genes_and_matched_enhancers.csv"), row.names = FALSE)
# 
# # plot the dynamics of these genes across trajectory
# ht <- res$heatmap
# png(paste0(temp_plot_path, 'Genes_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
# draw(ht)
# graphics.off()
# 
# # select the TFs that correlate in 'activity' and expression
# res <- SelectTFs_updated(object = obj.pair, trajectory.name = "rna_lineage_NC_probability", return.heatmap = TRUE,
#                          groupEvery = 2,
#                          p.cutoff = 0.01, cor.cutoff = 0.3)
# 
# # save the (top 100?) TFs that correlate
# df.cor <- res$tfs
# dim(df.cor)
# write.csv(df.cor, file = paste0(temp_csv_path, "TF_correlations.csv"), row.names = FALSE)
# 
# # plot TF activity dynamics across trajectory
# ht <- res$heatmap
# png(paste0(temp_plot_path, 'TFs_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
# draw(ht)
# graphics.off()
# 
# ############################## GRN inference #######################################
# 
# # GRN inference
# tf.gene.cor <- GetTFGeneCorrelation_updated(object = obj.pair, 
#                                             tf.use = df.cor$tfs, 
#                                             gene.use = unique(df.p2g$gene),
#                                             tf.assay = "chromvar", 
#                                             gene.assay = "RNA",
#                                             trajectory.name = "rna_lineage_NC_probability")
# 
# ## plot TF-gene correlation heatmap
# ht <- GRNHeatmap(tf.gene.cor, tf.timepoint = df.cor$time_point)
# png(paste0(temp_plot_path, 'TF_gene_corr_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
# ht
# graphics.off()
# 
# # To associate genes to TFs, we will use the peak-to-gene links and TF binding sites information. Specifically, 
# # if a gene is regulated by a peak and this peak is bound by a TF, then we consider this gene to be a target of this TF.
# motif.matching <- obj.pair@assays$ATAC@motifs@data
# colnames(motif.matching) <- obj.pair@assays$ATAC@motifs@motif.names
# motif.matching <- motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]
# 
# # take the TF-gene correlation, peak-TF binding prediction and peaks-to-genes linkage to build network
# df.grn <- GetGRN(motif.matching = motif.matching, 
#                  df.cor = tf.gene.cor, 
#                  df.p2g = df.p2g)
# dim(df.grn)
# 
# # save network
# write.csv(df.grn, file = paste0(temp_csv_path, "GRN_data.csv"), row.names = FALSE)
# 
# ############################## GRN visualisation #######################################
# 
# # define colors for nodes representing TFs (i.e., regulators)
# df.cor <- df.cor[order(df.cor$time_point), ]
# tfs.timepoint <- df.cor$time_point
# names(tfs.timepoint) <- df.cor$tfs
# 
# # plot whole GRN
# p <- GRNPlot(df.grn,
#              tfs.timepoint = tfs.timepoint,
#              show.tf.labels = TRUE,
#              seed = 42, 
#              plot.importance = TRUE,
#              min.importance = 2,
#              remove.isolated = FALSE)
# 
# png(paste0(temp_plot_path, 'Network_all.png'), height = 30, width = 45, units = 'cm', res = 400)
# print(p)
# graphics.off()
# 
# # filter grn
# df.grn2 <- df.grn %>%
#   subset(correlation > 0.7) %>%
#   dplyr::select(c(tf, gene, correlation)) %>%
#   rename(weights = correlation)
# 
# # which TFs are left in the GRN
# unique(df.grn2$tf)
# 
# # plot filtered GRN
# p <- GRNPlot(df.grn,
#              tfs.timepoint = tfs.timepoint,
#              show.tf.labels = TRUE,
#              plot.importance = TRUE,
#              min.importance = 2,
#              remove.isolated = TRUE)
# 
# png(paste0(temp_plot_path, 'Network_filtered.png'), height = 30, width = 45, units = 'cm', res = 400)
# print(p)
# graphics.off()
# 
# ######################################################################################
# ##############################    NEURAL     #######################################
# ######################################################################################
# 
# ############################## Filter genes and TFs #######################################
# 
# temp_plot_path = "./plots/neural_lineage/"
# dir.create(temp_plot_path, recursive = T)
# temp_csv_path = "./csv_files/neural_lineage/"
# dir.create(temp_csv_path, recursive = T)
# 
# # select genes that vary with the trajectory and correlate with peaks
# res <- SelectGenes_updated(obj.pair, trajectory.name = "rna_lineage_neural_probability", groupEvery = 2,
#                            var.cutoff.gene = 0.3, # how much gene expression has to vary across trajectory
#                            cor.cutoff = 0.3, fdr.cutoff = 1e-04) # how much peaks and genes need to correlate
# 
# # save the most variable genes
# df.p2g <- res$p2g
# dim(df.p2g)
# write.csv(df.p2g, file = paste0(temp_csv_path, "Variable_genes_and_matched_enhancers.csv"), row.names = FALSE)
# 
# # plot the dynamics of these genes across trajectory
# ht <- res$heatmap
# png(paste0(temp_plot_path, 'Genes_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
# draw(ht)
# graphics.off()
# 
# # select the TFs that correlate in 'activity' and expression
# res <- SelectTFs_updated(object = obj.pair, trajectory.name = "rna_lineage_neural_probability", return.heatmap = TRUE,
#                          groupEvery = 2,
#                          p.cutoff = 0.01, cor.cutoff = 0.3)
# 
# # save the (top 100?) TFs that correlate
# df.cor <- res$tfs
# dim(df.cor)
# write.csv(df.cor, file = paste0(temp_csv_path, "TF_correlations.csv"), row.names = FALSE)
# 
# # plot TF activity dynamics across trajectory
# ht <- res$heatmap
# png(paste0(temp_plot_path, 'TFs_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
# draw(ht)
# graphics.off()
# 
# ############################## GRN inference #######################################
# 
# # GRN inference
# tf.gene.cor <- GetTFGeneCorrelation_updated(object = obj.pair, 
#                                             tf.use = df.cor$tfs, 
#                                             gene.use = unique(df.p2g$gene),
#                                             tf.assay = "chromvar", 
#                                             gene.assay = "RNA",
#                                             trajectory.name = "rna_lineage_neural_probability")
# 
# ## plot TF-gene correlation heatmap
# ht <- GRNHeatmap(tf.gene.cor, tf.timepoint = df.cor$time_point)
# png(paste0(temp_plot_path, 'TF_gene_corr_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
# ht
# graphics.off()
# 
# # To associate genes to TFs, we will use the peak-to-gene links and TF binding sites information. Specifically, 
# # if a gene is regulated by a peak and this peak is bound by a TF, then we consider this gene to be a target of this TF.
# motif.matching <- obj.pair@assays$ATAC@motifs@data
# colnames(motif.matching) <- obj.pair@assays$ATAC@motifs@motif.names
# motif.matching <- motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]
# 
# # take the TF-gene correlation, peak-TF binding prediction and peaks-to-genes linkage to build network
# df.grn <- GetGRN(motif.matching = motif.matching, 
#                  df.cor = tf.gene.cor, 
#                  df.p2g = df.p2g)
# dim(df.grn)
# 
# # save network
# write.csv(df.grn, file = paste0(temp_csv_path, "GRN_data.csv"), row.names = FALSE)
# 
# ############################## GRN visualisation #######################################
# 
# # define colors for nodes representing TFs (i.e., regulators)
# df.cor <- df.cor[order(df.cor$time_point), ]
# tfs.timepoint <- df.cor$time_point
# names(tfs.timepoint) <- df.cor$tfs
# 
# # plot whole GRN
# p <- GRNPlot(df.grn,
#              tfs.timepoint = tfs.timepoint,
#              show.tf.labels = TRUE,
#              seed = 42, 
#              plot.importance = TRUE,
#              min.importance = 2,
#              remove.isolated = FALSE)
# 
# png(paste0(temp_plot_path, 'Network_all.png'), height = 30, width = 45, units = 'cm', res = 400)
# print(p)
# graphics.off()
# 
# # filter grn
# df.grn2 <- df.grn %>%
#   subset(correlation > 0.7) %>%
#   dplyr::select(c(tf, gene, correlation)) %>%
#   rename(weights = correlation)
# 
# # which TFs are left in the GRN
# unique(df.grn2$tf)
# 
# # plot filtered GRN
# p <- GRNPlot(df.grn,
#              tfs.timepoint = tfs.timepoint,
#              show.tf.labels = TRUE,
#              plot.importance = TRUE,
#              min.importance = 2,
#              remove.isolated = TRUE)
# 
# png(paste0(temp_plot_path, 'Network_filtered.png'), height = 30, width = 45, units = 'cm', res = 400)
# print(p)
# graphics.off()
