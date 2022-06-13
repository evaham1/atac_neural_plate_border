#!/usr/bin/env Rscript

print("differential accessibility between stages")

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
library(ComplexHeatmap)

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
    
    # test input folder made
    data_path = "./output/NF-downstream_analysis/ArchR_integration/transfer_labels/peak_call/rds_files/"
    plot_path = "./output/NF-downstream_analysis/ArchR_integration/transfer_labels_late_peaks/plots/"
    
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
}

############################### FUNCTIONS ####################################

### Function to add unique ids to se peak object so can subset properly
add_unique_ids_to_se <- function(seMarker, ArchR, matrix_type) {
  
  if (matrix_type == "PeakMatrix") {
    tmp_peaks = data.frame(ArchR@peakSet)
    tmp_diff_peaks = data.frame(rowData(seMarker))
    diff_peaks_join_peakset = left_join(tmp_diff_peaks, tmp_peaks, 
                                        by = c("seqnames" = "seqnames", "start" = "start", "end" = "end"))
    diff_peaks_join_peakset$gene_name = paste(diff_peaks_join_peakset$nearestGene, diff_peaks_join_peakset$distToTSS,sep="_")
    diff_peaks_join_peakset$unique_id = paste0(diff_peaks_join_peakset$seqnames, ":", diff_peaks_join_peakset$start, "-", diff_peaks_join_peakset$end)
    rowData(seMarker) = diff_peaks_join_peakset
  }
  
  if (matrix_type == "GeneScoreMatrix") {
    rowData <- as.data.frame(rowData(seMarker))
    
    duplicated_gene_names <- rowData$name[duplicated(rowData$name)]
    duplicated_genes <- rowData[which(rowData$name %in% duplicated_gene_names), ]
    duplicated_genes <- duplicated_genes %>% group_by(name) %>% 
      arrange(name) %>% mutate(unique_id = paste0(name, "-", rowid(name)))
    join_df = left_join(rowData, duplicated_genes,
                        by = c("seqnames" = "seqnames", "start" = "start", "end" = "end", "name" = "name", "idx" = "idx", "strand" = "strand"))
    
    join_df <- join_df %>% mutate(unique_id = coalesce(unique_id, name))
    
    rowData(seMarker) = join_df
  }
  
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
                           labelRows = FALSE, fontSizeRows = 12, clusterRows = TRUE,
                           labelCols = TRUE, fontSizeCols = 12, clusterCols = TRUE) {
  
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
  
  # order rows by eucladian distance
  if (clusterRows == TRUE){
    dist_mat <- dist(mat, method = 'euclidean')
    hclust_avg <- hclust(dist_mat, method = 'average')
    ordered_features <- hclust_avg$labels[c(hclust_avg$order)]
    mat <- mat[match(ordered_features, rownames(mat)), ]
  }
  
  # order columns by eucladian distance
  if (clusterCols == TRUE){
    dist_mat <- dist(t(mat), method = 'euclidean')
    hclust_avg <- hclust(dist_mat, method = 'average')
    ordered_cell_groups <- hclust_avg$labels[c(hclust_avg$order)]
    mat <- mat[ , match(ordered_cell_groups, colnames(mat))]
  }
  
  Heatmap(
    matrix = mat,
    col = pal,
    heatmap_legend_param = list(title = "z-scores"),
    #top_annotation = topAnno, 
    # add raster stuff?
    
    #Column Options
    cluster_columns = FALSE,
    show_column_names = labelCols,
    column_names_gp = gpar(fontsize = fontSizeCols),
    column_names_max_height = unit(100, "mm"),
    #column_split = colData$stage,
    
    #Row Options
    cluster_rows = FALSE,
    show_row_names = labelRows,
    row_names_gp = gpar(fontsize = fontSizeRows)
    #row_split = row_split_params
  )
  
}

###########################################################################################
############################## Read in ArchR projects #####################################

stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

pal = paletteContinuous(set = "solarExtra", n = 100)

# Read in all data
files <- list.files(data_path, full.names = TRUE)
print(files)

full_data <- grep("FullData", files, invert = F, value = TRUE)
FullData <- loadArchRProject(path = full_data, force = TRUE, showLogo = FALSE)

getAvailableMatrices(FullData)
FullData@peakSet


##########################################################################
############################## PLOTS #####################################

plot_path = "./plots/full/"
dir.create(plot_path, recursive = T)

##### Prepare FULL data for plotting 
Full_se <- getMarkerFeatures(
  ArchRProj = FullData, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage")
Full_se <- add_unique_ids_to_se(Full_se, FullData, matrix_type = "PeakMatrix")

matrix <- extract_means_from_se(Full_se) # extract means df from se object
normalised_matrix <- Log2norm(matrix) # log2norm across all features in each cell group

###### Plot all differential peaks
ids <- extract_ids(Full_se, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE)
print(paste0("all peaks: ", length(ids)))

if (length(ids) > 4){
  matrix <- extract_means_from_se(Full_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)

  png(paste0(plot_path, 'full_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()
}

ids <- extract_ids(Full_se, cutOff = "FDR <= 0.05 & Log2FC >= 0", top_n = TRUE, n = 20)
print(paste0("all peaks top 20: ", length(ids)))

if (length(ids) > 4){
  matrix <- extract_means_from_se(Full_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)

  png(paste0(plot_path, 'full_heatmap_top.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
  graphics.off()
}

###### Individual volcano plots

plot_path = "./plots/volcano_plots/"
dir.create(plot_path, recursive = T)

png(paste0(plot_path, 'HH5_volcano.png'), height = 20, width = 20, units = 'cm', res = 400)
print(markerPlot(seMarker = Full_se, name = "HH5", cutOff = "FDR <= 0.05 & Log2FC >= 0", plotAs = "Volcano"))
graphics.off()

png(paste0(plot_path, 'HH6_volcano.png'), height = 20, width = 20, units = 'cm', res = 400)
print(markerPlot(seMarker = Full_se, name = "HH6", cutOff = "FDR <= 0.05 & Log2FC >= 0", plotAs = "Volcano"))
graphics.off()

png(paste0(plot_path, 'HH7_volcano.png'), height = 20, width = 20, units = 'cm', res = 400)
print(markerPlot(seMarker = Full_se, name = "HH7", cutOff = "FDR <= 0.05 & Log2FC >= 0", plotAs = "Volcano"))
graphics.off()

png(paste0(plot_path, 'ss4_volcano.png'), height = 20, width = 20, units = 'cm', res = 400)
print(markerPlot(seMarker = Full_se, name = "ss4", cutOff = "FDR <= 0.05 & Log2FC >= 0", plotAs = "Volcano"))
graphics.off()

png(paste0(plot_path, 'ss8_volcano.png'), height = 20, width = 20, units = 'cm', res = 400)
print(markerPlot(seMarker = Full_se, name = "ss8", cutOff = "FDR <= 0.05 & Log2FC >= 0", plotAs = "Volcano"))
graphics.off()

###### Individual stage heatmaps

plot_path = "./plots/heatmaps/"
dir.create(plot_path, recursive = T)

# create matrix across clusters
clusters_se <- getMarkerFeatures(
  ArchRProj = FullData, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters")
clusters_se <- add_unique_ids_to_se(clusters_se, FullData, matrix_type = "PeakMatrix")

matrix <- extract_means_from_se(Full_se) # extract means df from se object
normalised_matrix <- Log2norm(matrix) # log2norm across all features in each cell group

##

se <- getMarkerFeatures(
  ArchRProj = FullData, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage",
  useGroups = "ss8")
se <- add_unique_ids_to_se(se, FullData, matrix_type = "PeakMatrix")
ids <- extract_ids(se, cutOff = "FDR <= 0.05 & Log2FC >= 0", top_n = FALSE)
print(paste0("ss8: ", length(ids)))

if (length(ids) > 4){
  matrix <- extract_means_from_se(Full_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)
  png(paste0(plot_path, 'ss8_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()

  matrix <- extract_means_from_se(clusters_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)
  png(paste0(plot_path, 'ss8_clusters_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()
}

##

se <- getMarkerFeatures(
  ArchRProj = FullData, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage",
  useGroups = "ss4")
se <- add_unique_ids_to_se(se, FullData, matrix_type = "PeakMatrix")
ids <- extract_ids(se, cutOff = "FDR <= 0.05 & Log2FC >= 0", top_n = FALSE)
print(paste0("ss4: ", length(ids)))

if (length(ids) > 4){
  matrix <- extract_means_from_se(Full_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)
  png(paste0(plot_path, 'ss4_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()

  matrix <- extract_means_from_se(clusters_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)
  png(paste0(plot_path, 'ss4_clusters_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()
}

##

se <- getMarkerFeatures(
  ArchRProj = FullData, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage",
  useGroups = "HH7")
se <- add_unique_ids_to_se(se, FullData, matrix_type = "PeakMatrix")
ids <- extract_ids(se, cutOff = "FDR <= 0.05 & Log2FC >= 0", top_n = FALSE)
print(paste0("HH7: ", length(ids)))

if (length(ids) > 4){
  matrix <- extract_means_from_se(Full_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)
  png(paste0(plot_path, 'HH7_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()

  matrix <- extract_means_from_se(clusters_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)
  png(paste0(plot_path, 'HH7_clusters_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()
}

##

se <- getMarkerFeatures(
  ArchRProj = FullData, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage",
  useGroups = "HH6")
se <- add_unique_ids_to_se(se, FullData, matrix_type = "PeakMatrix")
ids <- extract_ids(se, cutOff = "FDR <= 0.05 & Log2FC >= 0", top_n = FALSE)
print(paste0("HH6: ", length(ids)))

if (length(ids) > 4){
  matrix <- extract_means_from_se(Full_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)
  png(paste0(plot_path, 'HH6_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()

  matrix <- extract_means_from_se(clusters_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)
  png(paste0(plot_path, 'HH6_clusters_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()
}

##

se <- getMarkerFeatures(
  ArchRProj = FullData, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage",
  useGroups = "HH5")
se <- add_unique_ids_to_se(se, FullData, matrix_type = "PeakMatrix")
ids <- extract_ids(se, cutOff = "FDR <= 0.05 & Log2FC >= 0", top_n = FALSE)
print(paste0("HH5: ", length(ids)))

if (length(ids) > 4){
  matrix <- extract_means_from_se(Full_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)
  png(paste0(plot_path, 'HH5_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()

  matrix <- extract_means_from_se(clusters_se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)
  png(paste0(plot_path, 'HH5_clusters_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
  print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()
}