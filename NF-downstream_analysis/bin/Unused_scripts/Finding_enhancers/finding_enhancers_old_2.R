#!/usr/bin/env Rscript

print("finding enhancers")
# look for peaks which are differentially upregulated in the NC or PPR clusters at ss8/ss4 and are also open earlier

##### NB: THIS SCRIPT RUNS DIFFERENTLY INTERACTIVELY AND THROUGH PIPELINE!!!!
##### DIFFERENT PEAKS AT END, DIFFERENT NUMBERS AND IDENTITIES

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
    data_path = "./output/NF-downstream_analysis/ArchR_peak_exploration/transfer_labels/peak_call/rds_files/"
    plot_path = "./output/NF-downstream_analysis/ArchR_peak_exploration/transfer_labels_late_peaks/plots/"
    se_path = "./NF-downstream_analysis/work/b8/6306d2907d1376bc51d5f92d5e48d7/rds_files/Full_se.RDS"
    
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

# takes matrix of average cut sites per cluster, sums the average across clusters within each stage
# then filters peaks that have an aggregated count score above the threshold in every stage
open_across_stages_test <- function(m, threshold_type = "min", threshold_HH5 = 0, threshold_HH6 = 0,
                                    threshold_HH7 = 0, threshold_ss4 = 0, threshold_ss8 = 0){
  df <- as.data.frame(m)
  sums <- data.frame(HH5_sum = rowSums(dplyr::select(df, starts_with("HH5"))),
                     HH6_sum = rowSums(dplyr::select(df, starts_with("HH6"))),
                     HH7_sum = rowSums(dplyr::select(df, starts_with("HH7"))),
                     ss4_sum = rowSums(dplyr::select(df, starts_with("ss4"))),
                     ss8_sum = rowSums(dplyr::select(df, starts_with("ss8"))))
  rownames(sums) <- rownames(m)
  
  if (!(threshold_type %in% c("min","max"))){
    print("threshold_type must either be 'min' or 'max'")
  }
  if (threshold_type == "min"){
    filtered <- filter(sums, HH5_sum > threshold_HH5, HH6_sum > threshold_HH6, 
                       HH7_sum > threshold_HH7, ss4_sum > threshold_ss4, ss8_sum > threshold_ss8)
  }
  if (threshold_type == "max"){
    filtered <- filter(sums, HH5_sum < threshold_HH5, HH6_sum < threshold_HH6, 
                       HH7_sum < threshold_HH7, ss4_sum < threshold_ss4, ss8_sum < threshold_ss8)
  }
  
  return(rownames(filtered))
}

######### Function to plot tracks of ids of interest - only works with one cluster done at a time!
plot_browser_tracks <- function(ArchR, se, cutOff = "FDR <= 0.01 & Log2FC >= 1", extend = 50000, 
                                groupBy = "clusters", ids = ids, plot_path = plot_path, prefix = "_enhancer_") {
  
  # extract granges objects and add unique ids
  gr <- getMarkers(se, cutOff = cutOff, returnGR = TRUE)
  gr <- gr[[1]]
  values(gr) <- DataFrame(unique_id = paste0(gr@seqnames, ":", gr@ranges))
  
  # extend granges object so can see around the peak
  gr_extended <- extendGR(gr = gr, upstream = extend, downstream = extend)
  
  # get coordinates of granges object that want to plot based on ids
  gr_subset <- gr[(elementMetadata(gr_extended)[, "unique_id"] %in% ids)]
  
  # plot granges objects
  for (row in 1:length(gr_subset$unique_id)){
    print(row)
    name <- str_replace(gr_subset$unique_id[row], ":", "-")
    print(name)
    p <- plotBrowserTrack(ArchR, region = gr_subset[row], groupBy = groupBy)
    png(paste0(plot_path, prefix, name, '.png'), height = 50, width = 50, units = 'cm', res = 400)
    grid::grid.draw(p)
    graphics.off()
  }
}

## function to create a granges object from a peak id and optionally extend it
make_gr_object <- function(id, extend = TRUE, extend_by = 50000){
  seqnames <- sub("\\:.*", "", id)
  range <- sub(".*\\:", "", id)
  start <- sub("\\-.*", "", range)
  end <- sub(".*\\-", "", range)
  df <- data.frame(seqnames=seqnames, start=start, end=end,
                   strand="*")
  
  gr <- makeGRangesFromDataFrame(df)
  values(gr) <- DataFrame(unique_id = id)
  
  if(extend == TRUE){
    gr <- extendGR(gr = gr, upstream = extend_by, downstream = extend_by)
  }
  
  return(gr)
}

#################################################################################
############################## Read in data #####################################

stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

pal = paletteContinuous(set = "solarExtra", n = 100)

# Read in ArchR project data
files <- list.files(data_path, full.names = TRUE)
print(files)

full_data <- grep("FullData", files, invert = F, value = TRUE)
FullData <- loadArchRProject(path = full_data, force = TRUE, showLogo = FALSE)

getAvailableMatrices(FullData)
FullData@peakSet

# Read in calculated summarised exp object across all data
se_data <- grep("Full_se", files, invert = F, value = TRUE)
print(se_data)
Full_se <- readRDS(se_data)

# when running interactively:
#Full_se <- readRDS(se_path)

#############################################################################
###############################   PPR   #####################################

##################### ss8 PPR vs ss8 everything else #######################
# diff accessible in PPR ss8 clusters (C7+C8) vs other cells at ss8
# then filter on annotation 
# then filter so open at all stages (early peaks) 
# OR only open from HH7 onwards (late peaks)

plot_path <- "./plots/PPR/diff_within_ss8"
dir.create(plot_path, recursive = T)

### Step 1: differentially accessible in PPR
se <- getMarkerFeatures(
  ArchRProj = FullData, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("ss8_C7", "ss8_C8"),
  bgdGroups = c("ss8_C1", "ss8_C2", "ss8_C3", "ss8_C4", "ss8_C5", "ss8_C6", "ss8_C9", "ss8_C10"))
se <- add_unique_ids_to_se(se, FullData, matrix_type = "PeakMatrix")
unique_ids <- unique(extract_ids(se, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE))

matrix <- extract_means_from_se(Full_se)
normalised_matrix <- Log2norm(matrix)
subsetted_matrix <- subset_matrix(normalised_matrix, unique_ids)

png(paste0(plot_path, 'diff_accessible.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
graphics.off()

print(paste0("length of unique_ids: ", length(unique_ids))) #10481
id_data <- as.data.frame(rowData(se)[which(rowData(se)$unique_id %in% unique_ids), ])
print(dim(id_data)) #10481 x 21

### Step 2: filter out peaks in genes
filtered_id_data <- id_data[which(id_data$peakType %in% c("Distal", "Intronic")), ]
filtered_ids <- filtered_id_data$unique_id
print(paste0("length of filtered_ids: ", length(filtered_ids))) # 6841

subsetted_matrix <- subset_matrix(normalised_matrix, filtered_ids)

png(paste0(plot_path, 'diff_accessible_annot_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
graphics.off()

id_data <- id_data %>% mutate(annotation_filter = ifelse(unique_id %in% filtered_ids == TRUE, "T", "F"))

### Step 3: early peaks (open at every stage from HH5)
subsetted_raw_matrix <- subset_matrix(matrix, filtered_ids)
early_peaks <- open_across_stages_test(subsetted_raw_matrix, threshold_type = "min", threshold_HH5 = 1.5,
                                      threshold_HH6 = 1.5, threshold_HH7 = 1.5, threshold_ss4 = 1.5, threshold_ss8 = 1.5)
print(paste0("length of early_peaks: ", length(early_peaks))) # 20

subsetted_matrix <- subset_matrix(normalised_matrix, early_peaks)

png(paste0(plot_path, 'diff_accessible_annot_filtered_open_from_HH5.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(early_filter = ifelse(unique_id %in% early_peaks == TRUE, "T", "F"))

### Step 3: late peaks (open only from HH7)
late_peaks <- open_across_stages_test(subsetted_raw_matrix, threshold_type = "max", threshold_HH5 = 0.015,
                                      threshold_HH6 = 0.015, threshold_HH7 = 10, threshold_ss4 = 10, threshold_ss8 = 10)
print(paste0("length of late_peaks: ", length(late_peaks))) # 23

subsetted_matrix <- subset_matrix(normalised_matrix, late_peaks)

png(paste0(plot_path, 'diff_accessible_annot_filtered_open_from_HH7.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(late_filter = ifelse(unique_id %in% late_peaks == TRUE, "T", "F"))

### Step 4: export peaks and visualise them
write.csv(id_data, file = paste0(plot_path, "ss8_PPR_putative_enhancers_table.csv"))

# make genome browser plots for open peaks
plot_path <- paste0(plot_path, "browser_tracks/early_peaks/")
dir.create(plot_path, recursive = T)
for (id in early_peaks){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(FullData, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}

plot_path <- paste0(plot_path, "browser_tracks/late_peaks/")
dir.create(plot_path, recursive = T)
for (id in late_peaks){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(FullData, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}


################ PPR: HH7,ss4,ss8 PPR vs everything else ####################
# diff accessible in PPR clusters at HH7, ss4 and ss8 vs all other cells
# then filter on annotation, 
# then filter so open at all stages (early peaks) 
# OR only open from HH7 onwards (late peaks)

plot_path <- "./plots/PPR/diff_HH7_ss4_ss8"
dir.create(plot_path, recursive = T)

### Step 1: differentially accessible in PPR
se <- getMarkerFeatures(
  ArchRProj = FullData, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("HH7_C4", "ss4_C2", "ss4_C3", "ss8_C7", "ss8_C8"))
se <- add_unique_ids_to_se(se, FullData, matrix_type = "PeakMatrix")
unique_ids <- unique(extract_ids(se, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE))

matrix <- extract_means_from_se(Full_se)
normalised_matrix <- Log2norm(matrix)
subsetted_matrix <- subset_matrix(normalised_matrix, unique_ids)

png(paste0(plot_path, 'diff_accessible.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
graphics.off()

print(paste0("length of unique_ids: ", length(unique_ids))) #7883
id_data <- as.data.frame(rowData(se)[which(rowData(se)$unique_id %in% unique_ids), ])
print(dim(id_data))

### Step 2: filter out peaks in genes
filtered_id_data <- id_data[which(id_data$peakType %in% c("Distal", "Intronic")), ]
filtered_ids <- filtered_id_data$unique_id
print(paste0("length of filtered_ids: ", length(filtered_ids))) # 7139

subsetted_matrix <- subset_matrix(normalised_matrix, filtered_ids)

png(paste0(plot_path, 'diff_accessible_annot_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
graphics.off()

id_data <- id_data %>% mutate(annotation_filter = ifelse(unique_id %in% filtered_ids == TRUE, "T", "F"))

### Step 3: early peaks (open at every stage from HH5)
subsetted_raw_matrix <- subset_matrix(matrix, filtered_ids)
early_peaks <- open_across_stages_test(subsetted_raw_matrix, threshold_type = "min", threshold_HH5 = 1,
                                      threshold_HH6 = 1, threshold_HH7 = 1, threshold_ss4 = 1, threshold_ss8 = 1)
print(paste0("length of early_peaks: ", length(early_peaks))) # 39

subsetted_matrix <- subset_matrix(normalised_matrix, early_peaks)

png(paste0(plot_path, 'diff_accessible_annot_filtered_open_from_HH5.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(open_early_filter = ifelse(unique_id %in% early_peaks == TRUE, "T", "F"))

### Step 3: late peaks (only open from HH7)
late_peaks <- open_across_stages_test(subsetted_raw_matrix, threshold_type = "max", threshold_HH5 = 0.01,
                                        threshold_HH6 = 0.01, threshold_HH7 = 10, threshold_ss4 = 10, threshold_ss8 = 10)
print(paste0("length of late_peaks: ", length(late_peaks))) # 7

subsetted_matrix <- subset_matrix(normalised_matrix, late_peaks)

png(paste0(plot_path, 'diff_accessible_annot_filtered_open_from_HH7.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(closed_early_filter = ifelse(unique_id %in% late_peaks == TRUE, "T", "F"))

### Step 4: export peaks and visualise them
write.csv(id_data, file = paste0(plot_path, "HH7_ss4_ss8_PPR_putative_enhancers_table.csv"))

# make genome browser plots for open peaks
plot_path <- paste0(plot_path, "browser_tracks/early_peaks/")
dir.create(plot_path, recursive = T)
for (id in early_peaks){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(FullData, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}

plot_path <- paste0(plot_path, "browser_tracks/late_peaks/")
dir.create(plot_path, recursive = T)
for (id in late_peaks){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(FullData, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}

#############################################################################
##############################    NC    #####################################

# diff accessible in NC ss8 clusters (C1) vs other clusters at ss8 (although we have NC mix...)
# filter on annotation, 
# then filter so open at all stages OR only open from HH7 onwards

plot_path <- "./plots/NC/diff_within_ss8"
dir.create(plot_path, recursive = T)

### Step 1: differentially accessible in NC
se <- getMarkerFeatures(
  ArchRProj = FullData, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("ss8_C1"),
  bgdGroups = c("ss8_C2", "ss8_C3", "ss8_C4", "ss8_C5", "ss8_C6", "ss8_C7", "ss8_C8", "ss8_C9", "ss8_C10"))
se <- add_unique_ids_to_se(se, FullData, matrix_type = "PeakMatrix")
unique_ids <- unique(extract_ids(se, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE))

matrix <- extract_means_from_se(Full_se)
normalised_matrix <- Log2norm(matrix)
subsetted_matrix <- subset_matrix(normalised_matrix, unique_ids)

png(paste0(plot_path, 'diff_accessible.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
graphics.off()

print(paste0("length of unique_ids: ", length(unique_ids))) #4312
id_data <- as.data.frame(rowData(se)[which(rowData(se)$unique_id %in% unique_ids), ])
print(dim(id_data)) #4312 x 21

### Step 2: filter out peaks in genes - FILTER_1
filtered_id_data <- id_data[which(id_data$peakType %in% c("Distal", "Intronic")), ]
filtered_ids <- filtered_id_data$unique_id
print(paste0("length of filtered_ids: ", length(filtered_ids))) # 3938

subsetted_matrix <- subset_matrix(normalised_matrix, filtered_ids)

png(paste0(plot_path, 'diff_accessible_annot_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
graphics.off()

id_data <- id_data %>% mutate(filter_1 = ifelse(unique_id %in% filtered_ids == TRUE, "T", "F"))

### Step 3: earlyish peaks (open from HH7)
subsetted_raw_matrix <- subset_matrix(matrix, filtered_ids)
earlyish_peaks <- open_across_stages_test(subsetted_raw_matrix, threshold_type = "min", threshold_HH5 = 0.00001,
                                      threshold_HH6 = 0.00001, threshold_HH7 = 1, threshold_ss4 = 1, threshold_ss8 = 1)
print(paste0("length of earlyish_peaks: ", length(earlyish_peaks))) # 35

subsetted_matrix <- subset_matrix(normalised_matrix, earlyish_peaks)

png(paste0(plot_path, 'diff_accessible_annot_filtered_open_from_HH7.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(filter_2 = ifelse(unique_id %in% earlyish_peaks == TRUE, "T", "F"))

### Step 3: late peaks (open from ss4)
late_peaks <- open_across_stages_test(subsetted_raw_matrix, threshold_type = "max", threshold_HH5 = 0.01,
                                        threshold_HH6 = 0.01, threshold_HH7 = 0.01, threshold_ss4 = 10, threshold_ss8 = 10)
print(paste0("length of late_peaks: ", length(late_peaks))) # 12

subsetted_matrix <- subset_matrix(normalised_matrix, late_peaks)

png(paste0(plot_path, 'diff_accessible_annot_filtered_open_from_ss4.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(closed_early_filter = ifelse(unique_id %in% late_peaks == TRUE, "T", "F"))

### Step 4: export peaks and visualise them
write.csv(id_data, file = paste0(plot_path, "ss8_NC_putative_enhancers_table.csv"))

# make genome browser plots
plot_path <- paste0(plot_path, "browser_tracks/earlyish_peaks/")
dir.create(plot_path, recursive = T)
for (id in earlyish_peaks){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(FullData, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}

plot_path <- paste0(plot_path, "browser_tracks/late_peaks/")
dir.create(plot_path, recursive = T)
for (id in late_peaks){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(FullData, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}


# diff accessible in NC HH7,ss4,ss8 clusters, 
# filter on annotation, 
# then filter for earlyish (open from HH7)
# or late (from ss4)

plot_path <- "./plots/NC/diff_HH7_ss4_ss8"
dir.create(plot_path, recursive = T)

### Step 1: differentially accessible in NC
se <- getMarkerFeatures(
  ArchRProj = FullData, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("HH7_C5", "HH7_C6", "ss4_C6", "ss8_C1"))
se <- add_unique_ids_to_se(se, FullData, matrix_type = "PeakMatrix")
unique_ids <- unique(extract_ids(se, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE))

matrix <- extract_means_from_se(Full_se)
normalised_matrix <- Log2norm(matrix)
subsetted_matrix <- subset_matrix(normalised_matrix, unique_ids)

png(paste0(plot_path, 'diff_accessible.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
graphics.off()

print(paste0("length of unique_ids: ", length(unique_ids))) #5825
id_data <- as.data.frame(rowData(se)[which(rowData(se)$unique_id %in% unique_ids), ])
print(dim(id_data)) #7531 x 21

### Step 2: filter out peaks in genes - FILTER_1
filtered_id_data <- id_data[which(id_data$peakType %in% c("Distal", "Intronic")), ]
filtered_ids <- filtered_id_data$unique_id
print(paste0("length of filtered_ids: ", length(filtered_ids))) # 5320

subsetted_matrix <- subset_matrix(normalised_matrix, filtered_ids)

png(paste0(plot_path, 'diff_accessible_annot_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
graphics.off()

id_data <- id_data %>% mutate(filter_1 = ifelse(unique_id %in% filtered_ids == TRUE, "T", "F"))

### Step 3: earlyish peaks (open from HH7)
subsetted_raw_matrix <- subset_matrix(matrix, filtered_ids)
earlyish_peaks <- open_across_stages_test(subsetted_raw_matrix, threshold_type = "min", threshold_HH5 = 0.00001,
                                          threshold_HH6 = 0.00001, threshold_HH7 = 1, threshold_ss4 = 1, threshold_ss8 = 1)
print(paste0("length of earlyish_peaks: ", length(earlyish_peaks))) # 35

subsetted_matrix <- subset_matrix(normalised_matrix, earlyish_peaks)

png(paste0(plot_path, 'diff_accessible_annot_filtered_open_from_HH7.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(filter_2 = ifelse(unique_id %in% earlyish_peaks == TRUE, "T", "F"))

### Step 3: late peaks (open from ss4)
late_peaks <- open_across_stages_test(subsetted_raw_matrix, threshold_type = "max", threshold_HH5 = 0.01,
                                      threshold_HH6 = 0.01, threshold_HH7 = 0.01, threshold_ss4 = 10, threshold_ss8 = 10)
print(paste0("length of late_peaks: ", length(late_peaks))) # 12

subsetted_matrix <- subset_matrix(normalised_matrix, late_peaks)

png(paste0(plot_path, 'diff_accessible_annot_filtered_open_from_ss4.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(closed_early_filter = ifelse(unique_id %in% late_peaks == TRUE, "T", "F"))

### Step 4: export peaks and visualise them
write.csv(id_data, file = paste0(plot_path, "ss8_NC_putative_enhancers_table.csv"))

# make genome browser plots
plot_path <- paste0(plot_path, "browser_tracks/earlyish_peaks/")
dir.create(plot_path, recursive = T)
for (id in earlyish_peaks){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(FullData, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}

plot_path <- paste0(plot_path, "browser_tracks/late_peaks/")
dir.create(plot_path, recursive = T)
for (id in late_peaks){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(FullData, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}