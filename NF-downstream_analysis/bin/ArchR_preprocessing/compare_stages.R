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
library(ComplexHeatmap)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
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
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/test_input/"
    plot_path = "./output/NF-downstream_analysis/ArchR_preprocessing/compare_stages/"
    
    # stage_clusters on full data
    #data_path = "./output/NF-downstream_analysis/ArchR_integration/FullData/7_peak_calling/rds_files/"
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/"
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
                           labelRows = FALSE, fontSizeRows = 12,
                           labelCols = TRUE, fontSizeCols = 12) {
  
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
  dist_mat <- dist(mat, method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'average')
  ordered_features <- hclust_avg$labels[c(hclust_avg$order)]
  mat <- mat[match(ordered_features, rownames(mat)), ]
  
  # order columns by eucladian distance
  dist_mat <- dist(t(mat), method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'average')
  ordered_cell_groups <- hclust_avg$labels[c(hclust_avg$order)]
  mat <- mat[ , match(ordered_cell_groups, colnames(mat))]
  
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
############################## Read in ArchR projects #######################################

stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

# Read in all data
files <- list.files(data_path, full.names = TRUE)
print(files)
stages_data <- grep("FullData", files, invert = T, value = TRUE) # source data from which labels are extracted
print(paste0("Stages data: ", stages_data))

HH5 <- loadArchRProject(path = stages_data[1], force = TRUE, showLogo = FALSE)
print(HH5)
HH6 <- loadArchRProject(path = stages_data[2], force = TRUE, showLogo = FALSE)
print(HH6)
HH7 <- loadArchRProject(path = stages_data[3], force = TRUE, showLogo = FALSE)
print(HH7)
ss4 <- loadArchRProject(path = stages_data[4], force = TRUE, showLogo = FALSE)
print(ss4)
ss8 <- loadArchRProject(path = stages_data[5], force = TRUE, showLogo = FALSE)
print(ss8)

####################### Calculate diff features between clusters in each stage ##########################

HH5_se <- getMarkerFeatures(
  ArchRProj = HH5, 
  useMatrix = opt$matrix, 
  groupBy = "clusters")
HH6_se <- getMarkerFeatures(
  ArchRProj = HH6, 
  useMatrix = opt$matrix, 
  groupBy = "clusters")
HH7_se <- getMarkerFeatures(
  ArchRProj = HH7, 
  useMatrix = opt$matrix, 
  groupBy = "clusters")
ss4_se <- getMarkerFeatures(
  ArchRProj = ss4, 
  useMatrix = opt$matrix, 
  groupBy = "clusters")
ss8_se <- getMarkerFeatures(
  ArchRProj = ss8, 
  useMatrix = opt$matrix, 
  groupBy = "clusters")

HH5_se <- add_unique_ids_to_se(HH5_se, HH5, matrix_type = opt$matrix)
HH6_se <- add_unique_ids_to_se(HH6_se, HH6, matrix_type = opt$matrix)
HH7_se <- add_unique_ids_to_se(HH7_se, HH7, matrix_type = opt$matrix)
ss4_se <- add_unique_ids_to_se(ss4_se, ss4, matrix_type = opt$matrix)
ss8_se <- add_unique_ids_to_se(ss8_se, ss8, matrix_type = opt$matrix)


###################### Plots showing distribution of FDR and Logf2c values #############################

HH5_df <- data.frame(LogFC = c(t(assays(HH5_se)$Log2FC)), FDR = c(t(assays(HH5_se)$FDR)), stage = "HH5", stringsAsFactors=FALSE)
HH6_df <- data.frame(LogFC = c(t(assays(HH6_se)$Log2FC)), FDR = c(t(assays(HH6_se)$FDR)), stage = "HH6", stringsAsFactors=FALSE)
HH7_df <- data.frame(LogFC = c(t(assays(HH7_se)$Log2FC)), FDR = c(t(assays(HH7_se)$FDR)), stage = "HH7", stringsAsFactors=FALSE)
ss4_df <- data.frame(LogFC = c(t(assays(ss4_se)$Log2FC)), FDR = c(t(assays(ss4_se)$FDR)), stage = "ss4", stringsAsFactors=FALSE)
ss8_df <- data.frame(LogFC = c(t(assays(ss8_se)$Log2FC)), FDR = c(t(assays(ss8_se)$FDR)), stage = "ss8", stringsAsFactors=FALSE)

png(paste0(plot_path, 'HH5_FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(HH5_df, aes(x = -LogFDR, y = LogFC)) + 
  geom_point(alpha = 0.5) + 
  geom_segment(aes(x = 2, xend = 2, y = 1, yend = max(LogFC), colour = "black")) +
  geom_segment(aes(x = 2, xend = max(abs(LogFDR)), y = 1, yend = 1, colour = "black"))
graphics.off()

png(paste0(plot_path, 'HH6_FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(HH6_df, aes(x = -LogFDR, y = LogFC)) + 
  geom_point(alpha = 0.5) + 
  geom_segment(aes(x = 2, xend = 2, y = 1, yend = max(LogFC), colour = "black")) +
  geom_segment(aes(x = 2, xend = max(abs(LogFDR)), y = 1, yend = 1, colour = "black"))
graphics.off()

png(paste0(plot_path, 'HH7_FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(HH7_df, aes(x = -LogFDR, y = LogFC)) + 
  geom_point(alpha = 0.5) + 
  geom_segment(aes(x = 2, xend = 2, y = 1, yend = max(LogFC), colour = "black")) +
  geom_segment(aes(x = 2, xend = max(abs(LogFDR)), y = 1, yend = 1, colour = "black"))
graphics.off()

png(paste0(plot_path, 'ss4_FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(ss4_df, aes(x = -LogFDR, y = LogFC)) + 
  geom_point(alpha = 0.5) + 
  geom_segment(aes(x = 2, xend = 2, y = 1, yend = max(LogFC), colour = "black")) +
  geom_segment(aes(x = 2, xend = max(abs(LogFDR)), y = 1, yend = 1, colour = "black"))
graphics.off()

png(paste0(plot_path, 'ss8_FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(ss8_df, aes(x = -LogFDR, y = LogFC)) + 
  geom_point(alpha = 0.5) + 
  geom_segment(aes(x = 2, xend = 2, y = 1, yend = max(LogFC), colour = "black")) +
  geom_segment(aes(x = 2, xend = max(abs(LogFDR)), y = 1, yend = 1, colour = "black"))
graphics.off()

## all stages together

df <- do.call("rbind", list(HH5_df, HH6_df, HH7_df, ss4_df, ss8_df))

unique(df$stage)

## Boxplots
png(paste0(plot_path, 'Log2FC_boxplot.png'), height = 20, width = 20, units = 'cm', res = 400)
ggplot(df, aes(x=stage, y=LogFC)) + 
  geom_boxplot()
graphics.off()

png(paste0(plot_path, 'FDR_boxplot.png'), height = 20, width = 20, units = 'cm', res = 400)
ggplot(df, aes(x=stage, y=FDR)) + 
  geom_boxplot()
graphics.off()

df_cut <- df %>% group_by(stage) %>% filter(FDR < 0.05)

png(paste0(plot_path, 'FDR_0.05_cutoff_boxplot.png'), height = 20, width = 20, units = 'cm', res = 400)
ggplot(df_cut, aes(x=stage, y=FDR)) + 
  geom_boxplot()
graphics.off()

df <- df %>% mutate(LogFDR = as.numeric(log10(FDR))) %>%
  mutate(Passed = as.factor(ifelse(FDR < 0.01 & LogFC > 1,"passed", "failed")))

# randomly shuffle rows so plots point in a random order
set.seed(42)
rows <- sample(nrow(df))
plot_data <- df[rows, ]

plot_data <- plot_data[complete.cases(plot_data), ]

# scatter plot of -10log9FDR) VS Log2FC coloured by stage and significant features highlighted
png(paste0(plot_path, 'FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(plot_data, aes(x = -LogFDR, y = LogFC, color = stage, shape = Passed)) + 
  geom_point(alpha = 0.5) + 
  scale_color_manual(values=stage_colours) +
  scale_shape_manual(values=c(16, 17)) +
  geom_segment(aes(x = 2, xend = 2, y = 1, yend = max(LogFC), colour = "black")) +
  geom_segment(aes(x = 2, xend = max(abs(LogFDR)), y = 1, yend = 1, colour = "black"))
graphics.off()

###################### Plots showing how many features pass different thresholds #############################

sequence <- seq(from = 1, to = 0.01, by = -0.01)

############. LogFC >= 0
number_of_features <- c()
for (i in sequence){
  cutOff <- paste0("FDR <= ", i, " & Log2FC >= 0")
  ids <- extract_ids(ss8_se, cutOff = cutOff, top_n = FALSE) # extract ids
  print(length(ids))
  number_of_features <- c(number_of_features, length(ids))
}
ss8_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "ss8")

number_of_features <- c()
for (i in sequence){
  cutOff <- paste0("FDR <= ", i, " & Log2FC >= 0")
  ids <- extract_ids(ss4_se, cutOff = cutOff, top_n = FALSE) # extract ids
  print(length(ids))
  number_of_features <- c(number_of_features, length(ids))
}
ss4_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "ss4")

number_of_features <- c()
for (i in sequence){
  cutOff <- paste0("FDR <= ", i, " & Log2FC >= 0")
  ids <- extract_ids(HH7_se, cutOff = cutOff, top_n = FALSE) # extract ids
  print(length(ids))
  number_of_features <- c(number_of_features, length(ids))
}
HH7_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH7")

number_of_features <- c()
for (i in sequence){
  cutOff <- paste0("FDR <= ", i, " & Log2FC >= 0")
  ids <- extract_ids(HH6_se, cutOff = cutOff, top_n = FALSE) # extract ids
  print(length(ids))
  number_of_features <- c(number_of_features, length(ids))
}
HH6_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH6")

number_of_features <- c()
for (i in sequence){
  cutOff <- paste0("FDR <= ", i, " & Log2FC >= 0")
  ids <- extract_ids(HH5_se, cutOff = cutOff, top_n = FALSE) # extract ids
  print(length(ids))
  number_of_features <- c(number_of_features, length(ids))
}
HH5_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH5")

df <- do.call("rbind", list(HH5_df, HH6_df, HH7_df, ss4_df, ss8_df))

png(paste0(plot_path, 'changing_FDR_cutoffs_logFC_0_linegraph.png'), height = 20, width = 20, units = 'cm', res = 400)
ggplot(df, aes(x=FDR_cutoff, y=number_of_features, group=stage, colour=stage)) +
  geom_line() + scale_color_manual(values=stage_colours) + theme_minimal()
graphics.off()

df_cut <- df %>% group_by(stage) %>% filter(FDR_cutoff < 0.1)
png(paste0(plot_path, 'changing_FDR_cutoffs_zoom_logFC_0_linegraph.png'), height = 20, width = 20, units = 'cm', res = 400)
ggplot(df_cut, aes(x=FDR_cutoff, y=number_of_features, group=stage, colour=stage)) +
  geom_line() + scale_color_manual(values=stage_colours) + theme_minimal()
graphics.off()


############. LogFC >= 1
number_of_features <- c()
for (i in sequence){
  cutOff <- paste0("FDR <= ", i, " & Log2FC >= 1")
  ids <- extract_ids(ss8_se, cutOff = cutOff, top_n = FALSE) # extract ids
  print(length(ids))
  number_of_features <- c(number_of_features, length(ids))
}
ss8_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "ss8")

number_of_features <- c()
for (i in sequence){
  cutOff <- paste0("FDR <= ", i, " & Log2FC >= 1")
  ids <- extract_ids(ss4_se, cutOff = cutOff, top_n = FALSE) # extract ids
  print(length(ids))
  number_of_features <- c(number_of_features, length(ids))
}

ss4_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "ss4")

number_of_features <- c()
for (i in sequence){
  cutOff <- paste0("FDR <= ", i, " & Log2FC >= 1")
  ids <- extract_ids(HH7_se, cutOff = cutOff, top_n = FALSE) # extract ids
  print(length(ids))
  number_of_features <- c(number_of_features, length(ids))
}
HH7_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH7")

number_of_features <- c()
for (i in sequence){
  cutOff <- paste0("FDR <= ", i, " & Log2FC >= 1")
  ids <- extract_ids(HH6_se, cutOff = cutOff, top_n = FALSE) # extract ids
  print(length(ids))
  number_of_features <- c(number_of_features, length(ids))
}
HH6_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH6")

number_of_features <- c()
for (i in sequence){
  cutOff <- paste0("FDR <= ", i, " & Log2FC >= 1")
  ids <- extract_ids(HH5_se, cutOff = cutOff, top_n = FALSE) # extract ids
  print(length(ids))
  number_of_features <- c(number_of_features, length(ids))
}
HH5_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH5")

df <- do.call("rbind", list(HH5_df, HH6_df, HH7_df, ss4_df, ss8_df))

png(paste0(plot_path, 'changing_FDR_cutoffs_logFC_1_linegraph.png'), height = 20, width = 20, units = 'cm', res = 400)
ggplot(df, aes(x=FDR_cutoff, y=number_of_features, group=stage, colour=stage)) +
  geom_line() + scale_color_manual(values=stage_colours) + theme_minimal()
graphics.off()

df_cut <- df %>% group_by(stage) %>% filter(FDR_cutoff < 0.1)
png(paste0(plot_path, 'changing_FDR_cutoffs_zoom_logFC_1_linegraph.png'), height = 20, width = 20, units = 'cm', res = 400)
ggplot(df_cut, aes(x=FDR_cutoff, y=number_of_features, group=stage, colour=stage)) +
  geom_line() + scale_color_manual(values=stage_colours) + theme_minimal()
graphics.off()