#!/usr/bin/env Rscript

print("nice heatmaps")
# calculates differential peaks between clusters, scHelper_cell_types_old and stage -> plots on heatmaps

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
library(plyr)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-g", "--group_by"), action = "store", type = "character", help = "How to group cells to call peaks", default = "clusters",),
  make_option(c("-m", "--matrix"), action = "store", type = "character", help = "Matrix to use", default = "PeakMatrix",),
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
    
    #ss8
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/FILTERING/ss8/postfiltering/peak_calling/rds_files/"
    plot_path = "./output/NF-downstream_analysis/ArchR_preprocessing/FILTERING/ss8/postfiltering/heatmaps/plots"
  
    # stage_clusters on full data
    #data_path = "./output/NF-downstream_analysis/ArchR_integration/FullData/7_peak_calling/rds_files/"
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
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

############################### FUNCTIONS ####################################

### Function to add unique ids to se peak object so can subset properly
add_unique_ids_to_se <- function(seMarker, ArchR) {
  tmp_peaks = data.frame(ArchR@peakSet)
  tmp_diff_peaks = data.frame(rowData(seMarker))
  diff_peaks_join_peakset = left_join(tmp_diff_peaks, tmp_peaks, 
                                      by = c("seqnames" = "seqnames", "start" = "start", "end" = "end"))
  diff_peaks_join_peakset$gene_name = paste(diff_peaks_join_peakset$nearestGene, diff_peaks_join_peakset$distToTSS,sep="_")
  diff_peaks_join_peakset$unique_id = paste0(diff_peaks_join_peakset$seqnames, ":", diff_peaks_join_peakset$start, "-", diff_peaks_join_peakset$end)
  
  rowData(seMarker) = diff_peaks_join_peakset
  return(seMarker)
}

### Function to extract means from se object into matrix for plotting
extract_means_from_se <- function(seMarker) {
  mat <- as.data.frame(SummarizedExperiment::assays(seMarker)[["Mean"]])
  rownames(mat) <- rowData(seMarker)$unique_id
  
  return(mat)
}

### Function to Log2normalise matrix (must be run before subsetting!)
Log2norm <- function(mat, scaleTo = 10^4) {
  mat <- log2(t(t(mat)/colSums(mat)) * scaleTo + 1) # normalising means for depth of cluster
  return(mat)
}

### Function to extract IDs to be plotted, either by cut off or cut off + top n features
extract_ids <- function(seMarker, cutOff = "FDR <= 1 & Log2FC >= 0", top_n = TRUE, n = 10, group_name = "clusters") {
  
  markerList <- getMarkers(seMarker, cutOff = cutOff) # extract features that pass threshold
  
  df <- data.frame() # merged all features into a df
  for (i in 1:length(names(markerList))) {
    print(i)
    df_i <- as.data.frame(markerList[i])
    df <- rbind(df, df_i)
  }
  
  if (top_n == FALSE){
    ids <- df$unique_id
  } else {
    df <- df %>%
      group_by(group_name) %>%
      top_n(n, Log2FC) %>%
      dplyr::arrange(Log2FC, .by_group = TRUE)
    ids <- unique(df$unique_id)
  }
  
  return(ids)
}

### Function to subset normalised matrix using IDs
subset_matrix <- function(mat, ids) {
  subsetted_matrix <- mat[ids, ]
  return(subsetted_matrix)
}

### Function to plot marker heatmap
marker_heatmap <- function(mat, pal = NULL, 
                           labelRows = TRUE, clusterRows = TRUE, showRowDendrogram = TRUE, fontSizeRows = 12,
                           labelCols = TRUE, clusterCols = TRUE, showColDendrogram = TRUE, fontSizeCols = 12) {
  
  # scale each feature independently and add min/max limits
  limits <- c(-2, 2) # could make this user-defined
  mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), 
               `/`)
  mat[mat > max(limits)] <- max(limits)
  mat[mat < min(limits)] <- min(limits)
  
  # colours - set default if NULL
  if (is.null(pal) == TRUE) {
    pal <- paletteContinuous(set = "solarExtra", n = 100)
  }
  
  # legend
  legend <- list(at = c(0, 1),
                 labels = c(round(min(limits),2), round(max(limits),2)),
                 color_bar = "continuous",
                 legend_direction = "horizontal",
                 legend_width = unit(3, "cm"),
                 title = "Z-scores"
  )
  
  Heatmap(
    matrix = mat,
    col = pal,
    heatmap_legend_param = legend,
    top_annotation = topAnno, 
    # add raster stuff?
    
    #Column Options
    show_column_names = labelCols,
    cluster_columns = clusterCols,
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    column_names_gp = gpar(fontsize = fontSizeCols),
    column_names_max_height = unit(100, "mm"),
    column_split = colData$stage,
    
    #Row Options
    show_row_names = labelRows,
    cluster_rows = clusterRows,
    show_row_dend = showRowDendrogram,
    clustering_method_rows = "ward.D2",
    row_names_gp = gpar(fontsize = fontSizeRows)
    #row_split = row_split_params
  )
  
  return(Heatmap)
  
}

################################################################################








