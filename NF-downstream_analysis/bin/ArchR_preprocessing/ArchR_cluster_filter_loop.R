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
    
    plot_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ss8_Save-ArchR/ArchR_clustering/plots/"
    rds_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ss8_Save-ArchR/ArchR_clustering/rds_files/"
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ArchR_split/rds_files/"
    
    # already clustered
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/QC_MED/ss8/prefiltering/clustering/rds_files/"
    
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

############################### FUNCTIONS - adapted from scHelper #################################################
ArchR_IdentifyOutliers <- function(ArchR, group_by = 'clusters', metrics, intersect_metrics = TRUE, quantiles){
  outlier <- list()
  if(!length(quantiles) == 2){
    stop('quantiles must be an array of length == 2')
  }
  for(metric in metrics){
    min = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1])
    max = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2])
    
    outlier[[metric]] <- as.tibble(getCellColData(ArchR)) %>%
      group_by((!!as.symbol(group_by))) %>%
      summarise(median = median((!!as.symbol(metric)))) %>%
      filter(median > max | median < min) %>%
      pull(!!as.symbol(group_by))
  }
  
  if(intersect_metrics){
    if(length(Reduce(intersect, outlier)) == 0){
      cat('No outliers detected!')
    } else {
      return(Reduce(intersect, outlier))
    }
  } else{
    if(length(as.character(unique(unlist(outlier)))) == 0){
      cat('No outliers detected!')
    } else {
      return(as.character(unique(unlist(outlier))))
    }
  }
}

ArchR_ClustRes <- function(ArchR, by = 0.1, starting_res = 0){
  plots <- list()
  resolutions <- c(seq(starting_res, starting_res+9*by, by=by))
  cluster_df <- data.frame(ArchR$cellNames)
  
  if(length(ArchR@reducedDims) == 0){stop("Carry out dimensionality reduction before clustering")}
  
  for(res in resolutions[2:length(resolutions)]){
    ArchR_clustered <- addClusters(input = ArchR, name = "clusters", force = TRUE, resolution = res)
    plots[[paste(res)]] <- plotEmbedding(ArchR_clustered, name = "clusters") +
      ggtitle(paste("resolution = ", res))
    title <- paste0("clustering_res_", res)
    cluster_df <- cluster_df %>% mutate(!!title := ArchR_clustered@cellColData$clusters)
  }
  
  plots[["clustree"]] <- clustree(cluster_df, prefix = "clustering_res_")
  lay <- rbind(c(1,1,1,2,3,4),
               c(1,1,1,5,6,7),
               c(1,1,1,8,9,10))
  lay <- rbind(c(10,10,10,1,2,3),
               c(10,10,10,4,5,6),
               c(10,10,10,7,8,9))
  plots2 <- gridExtra::arrangeGrob(grobs = plots, layout_matrix = lay)
  return(gridExtra::grid.arrange(plots2))
}

# function to make heatmap showing contribution of cell groups to other cell groups
cell_counts_heatmap <- function(ArchR = ArchR, group1 = "scHelper_cell_type_new", group2 = "clusters") {
  group1_data <- getCellColData(ArchR, select = group1)[,1]
  group2_data <- getCellColData(ArchR, select = group2)[,1]
  cM <- confusionMatrix(paste0(group2_data), paste0(group1_data))
  cM <- cM / Matrix::rowSums(cM)
  
  p <- pheatmap::pheatmap(
    mat = cM,
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
  )
}

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
  print(paste0("iteration ", iteration, "/"))
  
  plot_path <- paste0(base_plot_path, "iteration_", iteration)
  dir.create(plot_path, recursive = T)
  
  ##############  PROCESSING  ##################
  # Dimensionality reduction
  ArchR <- addIterativeLSI(ArchR, force = TRUE)
  print("iterative LSI ran")

  # Run UMAP
  ArchR <- addUMAP(ArchR, force = TRUE)
  print("UMAP added")

  # # Clustree
  # png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
  # print(ArchR_ClustRes(ArchR, by = opt$clustree_by, starting_res = -opt$clustree_by))
  # graphics.off()

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
  p <- plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", name = metrics, plotAs = "Violin", baseSize = 20, alpha = 0.4, addBoxPlot = TRUE)
  p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1]), linetype = "dashed", color = "red")
  p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2]), linetype = "dashed", color = "red")
  png(paste0(plot_path, "VlnPlot_thresholds_nFrags.png"), width=50, height=20, units = 'cm', res = 200)
  print(p)
  graphics.off()
  metrics = "TSSEnrichment"
  p <- plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", name = metrics, plotAs = "Violin", baseSize = 20, alpha = 0.4, addBoxPlot = TRUE)
  p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1]), linetype = "dashed", color = "red")
  p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2]), linetype = "dashed", color = "red")
  png(paste0(plot_path, "VlnPlot_thresholds_TSSEnrichment.png"), width=50, height=20, units = 'cm', res = 200)
  print(p)
  graphics.off()
  metrics = "NucleosomeRatio"
  p <- plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", name = metrics, plotAs = "Violin", baseSize = 20, alpha = 0.4, addBoxPlot = TRUE)
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
  
  # how many cells filtered
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
  
  # set filtered ArchR for next iteration
  ArchR <- ArchR_filtered
  
  if (length(outliers) == 0){print("no more outliers found!")}

}


##############  SAVE  ##################
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)