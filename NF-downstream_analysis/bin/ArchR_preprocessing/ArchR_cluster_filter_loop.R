#!/usr/bin/env Rscript

print("loops until all poor quality clusters are removed")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)
library(clustree)
library(scHelper)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("", "--clust_res"), action = "store", type = "double", help = "clustering resolution for stage data", default = 1),
  make_option(c("", "--clustree_by"), action = "store", type = "double", help = "clustering res intervals for clustree", default = 0.1),
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
    
    data_path = "./output/NF-downstream_analysis/Upstream_processing/FILTERING/ss8/rds_files/" # already clustered
    plot_path = "./output/NF-downstream_analysis/Upstream_processing/FILTERING/ss8/rds_files/"
    rds_path = "./output/NF-downstream_analysis/Upstream_processing/FILTERING/ss8/rds_files/"

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
label <- sub('_.*', '', list.files(data_path))
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

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

#################################################################################
############################## LOOP #############################################

quantiles = c(0.2, 0.8)
iteration = 0
base_plot_path <- plot_path

outliers = "outliers"

while ( length(outliers) > 0 ) {
  
  iteration <- iteration + 1
  print(paste0("iteration ", iteration))
  
  plot_path <- paste0(base_plot_path, "iteration_", iteration, "/")
  print(plot_path)
  dir.create(plot_path, recursive = T)
  
  ##############  PROCESSING  ##################
  # Dimensionality reduction
  ArchR <- addIterativeLSI(ArchR, force = TRUE)
  print("iterative LSI ran")

  # Run UMAP
  ArchR <- addUMAP(ArchR, force = TRUE)
  print("UMAP added")

  # Clustree
  png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
  print(ArchR_ClustRes(ArchR, by = opt$clustree_by, starting_res = -opt$clustree_by))
  graphics.off()

  # Cluster
  ArchR <- addClusters(ArchR, name = "clusters", resolution = opt$clust_res, force = TRUE)
  print("clustering ran")

  ##############  PLOTS  ##################

  # Labelled UMAP
  png(paste0(plot_path, 'UMAP_clusters.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "clusters", plotAs = "points", size = 1.8, baseSize = 0, 
                labelSize = 10, legendSize = 0, randomize = TRUE, labelAsFactors = FALSE))
  graphics.off()
  
  # QC UMAPs
  png(paste0(plot_path, "UMAP_nFrags.png"), width=20, height=20, units = 'cm', res = 200)
  print(plotEmbedding(ArchR, name = "nFrags",
                plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                baseSize = 20, labelSize = 0, legendSize = 20, randomize = TRUE))
  graphics.off()
  png(paste0(plot_path, "UMAP_NucleosomeRatio.png"), width=20, height=20, units = 'cm', res = 200)
  print(plotEmbedding(ArchR, name = "NucleosomeRatio",
                plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                baseSize = 20, labelSize = 0, legendSize = 20, randomize = TRUE))
  graphics.off()
  png(paste0(plot_path, "UMAP_TSSEnrichment.png"), width=20, height=20, units = 'cm', res = 200)
  print(plotEmbedding(ArchR, name = "TSSEnrichment",
                plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                baseSize = 20, labelSize = 0, legendSize = 20, randomize = TRUE))
  graphics.off()

  # Violin Plots
  metrics = "nFrags"
  p <- plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", name = metrics, plotAs = "Violin", baseSize = 25, alpha = 0.4, addBoxPlot = TRUE)
  p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1]), linetype = "dashed", color = "red")
  p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2]), linetype = "dashed", color = "red")
  png(paste0(plot_path, "VlnPlot_thresholds_nFrags.png"), width=50, height=20, units = 'cm', res = 200)
  print(p)
  graphics.off()
  metrics = "TSSEnrichment"
  p <- plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", name = metrics, plotAs = "Violin", baseSize = 25, alpha = 0.4, addBoxPlot = TRUE)
  p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1]), linetype = "dashed", color = "red")
  p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2]), linetype = "dashed", color = "red")
  png(paste0(plot_path, "VlnPlot_thresholds_TSSEnrichment.png"), width=50, height=20, units = 'cm', res = 200)
  print(p)
  graphics.off()
  metrics = "NucleosomeRatio"
  p <- plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", name = metrics, plotAs = "Violin", baseSize = 25, alpha = 0.4, addBoxPlot = TRUE)
  p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1]), linetype = "dashed", color = "red")
  p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2]), linetype = "dashed", color = "red")
  png(paste0(plot_path, "VlnPlot_thresholds_NucleosomeRatio.png"), width=50, height=20, units = 'cm', res = 200)
  print(p)
  graphics.off()
  
  ##############  OUTLIERS  ##################
  outliers_TSS <- ArchR_IdentifyOutliers(ArchR, group_by = 'clusters', metrics = "TSSEnrichment", intersect_metrics = FALSE, quantiles = quantiles)
  outliers_nucleosome <- ArchR_IdentifyOutliers(ArchR, group_by = 'clusters', metrics = "NucleosomeRatio", intersect_metrics = FALSE, quantiles = quantiles)
  outliers <- as.character(unique(unlist(c(outliers_TSS, outliers_nucleosome))))
  print(outliers)
  
  idxSample <- BiocGenerics::which(!(ArchR$clusters %in% outliers))
  cellsSample <- ArchR$cellNames[idxSample]
  ArchR_filtered <- ArchR[cellsSample, ]
  
  # which clusters filtered
  unfiltered <- as.data.frame(table(ArchR$clusters))
  colnames(unfiltered) <- c("Cluster_ID", "Unfiltered_cell_count")
  unfiltered <- unfiltered %>% mutate(Cluster_ID = as.numeric(gsub('^.', '', Cluster_ID))) %>%
    arrange(Cluster_ID)
  filtered <- as.data.frame(table(ArchR_filtered$clusters)) %>% 
    mutate(Var1 = as.numeric(gsub('^.', '', Var1)))
  cell_counts <- merge(unfiltered, filtered, by.x = "Cluster_ID", by.y = "Var1", all = TRUE)
  cell_counts[is.na(cell_counts)] <- 0
  
  png(paste0(plot_path, 'cluster_cell_counts_table.png'), height = 30, width = 15, units = 'cm', res = 400)
  grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
               tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
  graphics.off()
  
  # how many cells filtered total
  unfiltered <- table(ArchR$stage)
  filtered <- table(ArchR_filtered$stage)
  cell_counts <- as_tibble(rbind(unfiltered, filtered))
  cell_counts <- cbind(cell_counts, Total = rowSums(cell_counts))
  
  png(paste0(plot_path, 'stage_cell_counts_table.png'), height = 10, width = 10, units = 'cm', res = 400)
  grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
               tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
  graphics.off()
  
  # set filtered ArchR for next iteration
  ArchR <- ArchR_filtered
  
  if (length(outliers) == 0){print("no more outliers found!")}

}


##############  SAVE  ##################
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)