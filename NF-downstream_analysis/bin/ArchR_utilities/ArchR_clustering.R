#!/usr/bin/env Rscript

print("clustering ArchR")

#### need to fix plotting for transfer labels object: stage_clusters, stage_scHelper_cell_type, stage_cluster_labels

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)
library(clustree)
library(presto)
library(Seurat)
library(gtools)
library(ComplexHeatmap)
library(scHelper)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
    make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
    make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
    make_option(c("", "--stage_clust_res"), action = "store", type = "double", help = "clustering resolution for stage data", default = 1),
    make_option(c("", "--full_clust_res"), action = "store", type = "double", help = "clustering resolution for full data", default = 2),
    make_option(c("", "--clustree_stage"), action = "store", type = "logical", help = "whether to run clustree plot on stage data", default = FALSE),
    make_option(c("", "--clustree_full"), action = "store", type = "logical", help = "whether to run clustree plot on full data", default = FALSE),
    make_option(c("", "--stage_clustree_by"), action = "store", type = "double", help = "clustering res intervals for clustree for stages", default = 0.1),
    make_option(c("", "--full_clustree_by"), action = "store", type = "double", help = "clustering res intervals for clustree for full data", default = 0.2),
    make_option(c("", "--GeneScore_heatmaps_stage"), action = "store", type = "logical", help = "whether to run gene score heatmaps on stage data", default = FALSE),
    make_option(c("", "--GeneScore_heatmaps_full"), action = "store", type = "logical", help = "whether to run gene score heatmaps on full data", default = FALSE),
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
    
    # already clustered
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Clustering/rds_files/"
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

set.seed(42)

############################## Read in ArchR project #######################################

# If files are not in rds_files subdirectory look in input dir 
label <- unique(sub('_.*', '', list.files(data_path)))
print(label) 

if (length(label) == 0){
  data_path = "./input/"
  label <- sub('_.*', '', list.files(data_path))
  print(label)
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
} else {
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
}

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getPeakSet(ArchR)
getAvailableMatrices(ArchR)

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

#################################################################################
############################## PROCESSING #######################################
# Dimensionality reduction
ArchR <- addIterativeLSI(ArchR, force = TRUE, seed = 1)
print("iterative LSI ran")

# Run UMAP
ArchR <- addUMAP(ArchR, force = TRUE, seed = 1)
print("UMAP added")

# Cluster
if (length(unique(ArchR$stage)) == 1){
  ArchR <- addClusters(ArchR, name = "clusters", resolution = opt$stage_clust_res, force = TRUE, seed = 1)
} else {
  ArchR <- addClusters(ArchR, name = "clusters", resolution = opt$full_clust_res, force = TRUE, seed = 1)
}
print("clustering ran")

# Plot Clustree (optional)
if (length(unique(ArchR$stage)) == 1){
  if (isTRUE(opt$clustree_stage)) {
    print("running clustree plot...")
    png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
      plot(ArchR_ClustRes(ArchR, by = opt$stage_clustree_by))
      graphics.off()
    print("clustree plot ran") }

} else {
  if (isTRUE(opt$clustree_full)) {
    print("running clustree plot...")
    png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
      plot(ArchR_ClustRes(ArchR, by = opt$full_clustree_by))
      graphics.off()
    print("clustree plot ran") }
}

################## Save clustered ArchR project #################################
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)
print("ArchR object saved")

#######################################################################################
############################ CELL COUNTS PER CLUSTER ##################################
plot_path_temp = "./plots/cell_counts/"
dir.create(plot_path_temp, recursive = T)

# Plot number of cells in each cluster
cluster_cell_counts <- as.data.frame(table(substr(ArchR$clusters, 2, nchar(ArchR$clusters))))
cluster_cell_counts <- cluster_cell_counts %>% 
  dplyr::rename(Cell_count = Freq, Cluster_number = Var1) %>%
  dplyr::mutate(Cluster_number = as.numeric(as.character(Cluster_number))) %>%
  dplyr::arrange(Cluster_number)
print("Cluster cell counts: ")
print(cluster_cell_counts)

cluster_cell_counts_totals <- cluster_cell_counts %>%
  rbind(c("Total", sum(cluster_cell_counts$Cell_count)))
print("Cluster cell counts with totals: ")
print(cluster_cell_counts_totals)
png(paste0(plot_path_temp, 'cluster_cell_counts_table.png'), height = 25, width = 10, units = 'cm', res = 400)
grid.arrange(tableGrob(cluster_cell_counts_totals, rows=NULL, theme = ttheme_minimal()))
graphics.off()

p<-ggplot(data=cluster_cell_counts, aes(x=`Cluster_number`, y=`Cell_count`)) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks = round(seq(min(cluster_cell_counts$Cluster_number), max(cluster_cell_counts$Cluster_number), by = 1),1))

png(paste0(plot_path_temp, 'cluster_cell_counts_barchart.png'), height = 10, width = 20, units = 'cm', res = 400)
print(p)
graphics.off()

# Plot contribution of each stage to each cluster
if (length(unique(ArchR$stage)) > 1){
  png(paste0(plot_path_temp, "cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
  ArchRCellCountsHeatmap(ArchR = ArchR, group1 = "clusters", group2 = "stage")
  graphics.off()
}

print("cell counts calculated")

###############################################################################################
################################### cluster UMAPS #############################################
plot_path_temp = "./plots/UMAPs/"
dir.create(plot_path_temp, recursive = T)

p1 <- plotEmbedding(ArchR, 
                    name = "stage",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    pal = stage_colours, randomize = TRUE)
p2 <- plotEmbedding(ArchR, 
                    name = "clusters",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0,
                    randomize = TRUE)

png(paste0(plot_path_temp, "UMAPs.png"), width=60, height=40, units = 'cm', res = 200)
ggAlignPlots(p1, p2, type = "h")
graphics.off()

png(paste0(plot_path_temp, 'UMAP_clusters.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "clusters", plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1), baseSize = 0, 
              labelSize = 10, legendSize = 0, randomize = TRUE, labelAsFactors = FALSE)
graphics.off()

if ( !(is.null(ArchR$stage_clusters)) ) {
  
  png(paste0(plot_path_temp, 'UMAP_stage_clusters.png'), height = 20, width = 20, units = 'cm', res = 400)
  plotEmbedding(ArchR, name = "stage_clusters", plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1), baseSize = 0, 
              labelSize = 10, legendSize = 0, randomize = TRUE, labelAsFactors = FALSE)
  graphics.off()

}

#################################################################################
############################ GENE SCORE PLOTS ###################################

plot_path = "./plots/Gene_score_plots/"
dir.create(plot_path, recursive = T)

##########    Feature plots

ArchR <- addImputeWeights(ArchR, seed = 1)

# set genes of interest
TFs <- c("SIX1", "EYA2", "IRF6", "DLX5", "DLX6", "GATA2", "GATA3", 
         "TFAP2A", "TFAP2B", "TFAP2C", 
         "PITX1", "PITX2",
         "PAX7", "MSX1", "ETS1", 
         "SOX2", "SOX9", "SOX8", "SOX10", "SOX5", "SOX21", "SOX3",
         "NKX6-2", "CSRNP1", "SNAI2", "LMX1B", "ZEB2",
         "EPAS1", "BMP4", "YEATS4", "HOXB1", "EOMES", "ADMP")

ArchR <- addImputeWeights(ArchR)

# Plot ridge plot of each TF deviation
for (TF in TFs){
  print(TF)
  
  # Plot distribution of GeneScore values for each cluster
  png(paste0(plot_path, TF, '_gene_score_ridge_plot.png'), height = 12, width = 10, units = 'cm', res = 400)
  print(plotGroups(ArchR, groupBy = "clusters", name = TF, colorBy = "GeneScoreMatrix") + 
    theme_ArchR(baseSize = 17, plotMarginCm = 0.5))
  graphics.off()
  
  # Plot GeneScore values on UMAP
  png(paste0(plot_path, TF, '_gene_score_UMAP.png'), height = 12, width = 14, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = TF,
                plotAs = "points", size = 1.8,
                colorBy = "GeneScoreMatrix", continuousSet = "horizon") + 
    theme_ArchR(legendTextSize = 12, baseSize = 16, plotMarginCm = 0.5))
  graphics.off()
  
}

print("Feature plots done")

##########    Heatmaps (optional)

run_heatmaps <- ifelse(length(unique(ArchR$stage)) == 1 & isTRUE(opt$GeneScore_heatmaps_stage) | length(unique(ArchR$stage)) > 1 & isTRUE(opt$GeneScore_heatmaps_full),
       TRUE, FALSE)

if (isTRUE(run_heatmaps)) {

  seMarker <- getMarkerFeatures(
    ArchRProj = ArchR, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "clusters")
  seMarker <- ArchRAddUniqueIdsToSe(seMarker, ArchR, matrix_type = "GeneScoreMatrix")

  # prepare for plotting
  normalised_matrix <- ArchR_ExtractMeansFromSe(seMarker, Log2norm = TRUE, scaleTo = 10^4) # extract means df from se object and log2norm all features in each cell group

  # heatmap palette 
  pal = viridis::magma(100)

  # Heatmap of positive markers which pass cutoff thresholds
  ids <- ArchR_ExtractIds(seMarker, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE) # extract ids
  if (length(ids) < 2){
    print(paste0(length(ids), " features passed cutoff - not enough to make heatmap"))
  } else {
    print(paste0(length(ids), " features passed cutoff - now plotting heatmap"))
    subsetted_matrix <- normalised_matrix[ids, ]
    
    png(paste0(plot_path_temp, 'diff_cutoff_heatmap.png'), height = 40, width = 20, units = 'cm', res = 400)
    print(ArchR_PlotMarkerHeatmap(subsetted_matrix, pal = pal, labelRows = TRUE))
    graphics.off()
  }

  # Heatmap of positive markers top 10 per cell group
  ids <- ArchR_ExtractIds(seMarker, cutOff = "FDR <= 0.05 & Log2FC >= 0", top_n = TRUE, n = 10) # extract ids
  subsetted_matrix <- normalised_matrix[ids, ]
  
  png(paste0(plot_path_temp, 'diff_top10_heatmap.png'), height = 40, width = 20, units = 'cm', res = 400)
  print(ArchR_PlotMarkerHeatmap(subsetted_matrix, labelRows = TRUE, pal = pal, cluster_columns = FALSE, cluster_rows = FALSE))
  graphics.off()
  
}

#################################################################################
############################ CELL LABEL PLOTS ###################################

###### schelper cell type colours
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak', 
                              'PGC', 'BI', 'meso', 'endo')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                                "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                                "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                       'vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo')

if ( !(is.null(ArchR$scHelper_cell_type)) ) {

  atac_scHelper_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type)]

  plot_path_temp <- "./plots/RNA_label_plots/"
  dir.create(plot_path_temp, recursive = T)
  
  png(paste0(plot_path_temp, 'UMAP_scHelper_cell_type.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "scHelper_cell_type", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_cols, labelAsFactors = FALSE))
  graphics.off()

  png(paste0(plot_path_temp, 'UMAP_scHelper_cell_type_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "scHelper_cell_type", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_cols))
  graphics.off()

}

if ( !(is.null(ArchR$stage_scHelper_cell_type)) ) {

  atac_scHelper_cols <- scHelper_cell_type_colours[unique(ArchR$stage_scHelper_cell_type)]

  plot_path_temp <- "./plots/RNA_label_plots/"
  dir.create(plot_path_temp, recursive = T)
  
  png(paste0(plot_path_temp, 'UMAP_stage_scHelper_cell_type.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "stage_scHelper_cell_type", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_cols, labelAsFactors = FALSE))
  graphics.off()

  png(paste0(plot_path_temp, 'UMAP_stage_scHelper_cell_type_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "stage_scHelper_cell_type", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_cols))
  graphics.off()

}

if ( !(is.null(ArchR$cluster_labels)) ) {

  plot_path_temp <- "./plots/RNA_label_plots/"
  dir.create(plot_path_temp, recursive = T)
  
  png(paste0(plot_path_temp, 'UMAP_cluster_labels.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "cluster_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_cols, labelAsFactors = FALSE))
  graphics.off()

  png(paste0(plot_path_temp, 'UMAP_cluster_labels_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "cluster_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_cols))
  graphics.off()

}

if ( !(is.null(ArchR$stage_cluster_labels)) ) {

  plot_path_temp <- "./plots/RNA_label_plots/"
  dir.create(plot_path_temp, recursive = T)
  
  png(paste0(plot_path_temp, 'UMAP_stage_cluster_labels.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "stage_cluster_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_cols, labelAsFactors = FALSE))
  graphics.off()

  png(paste0(plot_path_temp, 'UMAP_stage_cluster_labels_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "stage_cluster_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_cols))
  graphics.off()

}