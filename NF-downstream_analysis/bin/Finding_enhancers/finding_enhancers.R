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

############################## Set colours #######################################

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

###### schelper cell type colours ~ 29 cell states
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

pal = paletteContinuous(set = "solarExtra", n = 100)

###########################################################################################
############################## Read in data #####################################

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

getAvailableMatrices(ArchR)
ArchR@peakSet

# Read in calculated summarised exp object across all data
se_data <- readRDS(paste0(data_path, label, "_SE.RDS"))
print(se_data)

# when running interactively:
#Full_se <- readRDS(se_path)

# Extract matrix for heatmap plotting
matrix <- extract_means_from_se(se_data)
normalised_matrix <- Log2norm(matrix)

#############################################################################
###############################   PPR   #####################################
# diff accessible in PPR clusters at different stages
# then filter on annotation 
# then filter to remove ap differences

## Diff accessibility test to find ap differences in PPR clusters at ss4 and ss8 so can filter these out
ss8_ap_filter_se <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("ss8_C7"),
  bgdGroups = c("ss8_C8"))
ss8_ap_filter_se <- add_unique_ids_to_se(ss8_ap_filter_se, ArchR, matrix_type = "PeakMatrix")
ss8_ap_filter_ids <- unique(extract_ids(ss8_ap_filter_se, cutOff = "FDR <= 0.1", top_n = FALSE))
print(paste0("length of ss8_ap_filter_ids: ", length(ss8_ap_filter_ids))) # 5675

ss4_ap_filter_se <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("ss4_C2"),
  bgdGroups = c("ss4_C3"))
ss4_ap_filter_se <- add_unique_ids_to_se(ss4_ap_filter_se, ArchR, matrix_type = "PeakMatrix")
ss4_ap_filter_ids <- unique(extract_ids(ss4_ap_filter_se, cutOff = "FDR <= 0.1", top_n = FALSE))
print(paste0("length of ss4_ap_filter_ids: ", length(ss4_ap_filter_ids))) # 5675

##################### ss8 #######################
plot_path <- "./plots/PPR/diff_ss8/"
dir.create(plot_path, recursive = T)

### Step 1: differentially accessible in ss8 PPR clusters vs everything else
se <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("ss8_C7", "ss8_C8"))
se <- add_unique_ids_to_se(se, ArchR, matrix_type = "PeakMatrix")
ids_1 <- unique(extract_ids(se, cutOff = "FDR <= 0.01 & Log2FC >= 3", top_n = FALSE))
print(paste0("length of ids_1: ", length(ids_1))) 

subsetted_matrix <- subset_matrix(normalised_matrix, ids_1)
png(paste0(plot_path, '1_diff_accessible.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- as.data.frame(rowData(se)[which(rowData(se)$unique_id %in% ids_1), ])
print(dim(id_data)) #55 x 21

### Step 2: filter out peaks in genes
annot_id_data <- id_data[which(id_data$peakType %in% c("Distal", "Intronic")), ]
annot_keep_ids <- annot_id_data$unique_id
print(paste0("length of annot_keep_ids: ", length(annot_keep_ids))) # 51

ids_2 <- intersect(ids_1, annot_keep_ids)
print(paste0("length of ids_2: ", length(ids_2))) # 51

subsetted_matrix <- subset_matrix(normalised_matrix, ids_2)
png(paste0(plot_path, '2_diff_accessible_annot_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- id_data %>% mutate(annotation_filter = ifelse(unique_id %in% ids_2 == TRUE, "T", "F"))

### Step 3: filter out peaks that are differentially accessible between aPPR and pPPR
ap_filter_ids <- ids_2[ids_2 %in% ss8_ap_filter_ids]
print(paste0("length of ss8_ap_filter_ids that are in ids_2: ", length(ap_filter_ids))) # 22
ids_3 <- setdiff(ids_2, ap_filter_ids)
print(paste0("length of ids_3: ", length(ids_3))) # 29

subsetted_matrix <- subset_matrix(normalised_matrix, ids_3)
png(paste0(plot_path, '3_diff_accessible_annot_filtered_ap_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(ap_filter = ifelse(unique_id %in% ids_3 == TRUE, "T", "F"))

### Step 3: export peaks and visualise them
write.csv(id_data, file = paste0(plot_path, "putative_enhancers_table.csv"))

# make genome browser plots for open peaks
plot_path <- paste0(plot_path, "browser_tracks/")
dir.create(plot_path, recursive = T)
for (id in ids_3){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(ArchR, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}

##################### ss4 + ss8 #######################
plot_path <- "./plots/PPR/diff_ss4_ss8/"
dir.create(plot_path, recursive = T)

### Step 1: differentially accessible in ss8 PPR clusters vs everything else
se <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("ss4_C2", "ss4_C3", "ss8_C7", "ss8_C8"))
se <- add_unique_ids_to_se(se, ArchR, matrix_type = "PeakMatrix")
ids_1 <- unique(extract_ids(se, cutOff = "FDR <= 0.01 & Log2FC >= 5", top_n = FALSE))
print(paste0("length of ids_1: ", length(ids_1))) # 52

subsetted_matrix <- subset_matrix(normalised_matrix, ids_1)
png(paste0(plot_path, '1_diff_accessible.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- as.data.frame(rowData(se)[which(rowData(se)$unique_id %in% ids_1), ])
print(dim(id_data)) #52 x 21

### Step 2: filter out peaks in genes
annot_id_data <- id_data[which(id_data$peakType %in% c("Distal", "Intronic")), ]
annot_keep_ids <- annot_id_data$unique_id
print(paste0("length of annot_keep_ids: ", length(annot_keep_ids))) # 51

ids_2 <- intersect(ids_1, annot_keep_ids)
print(paste0("length of ids_2: ", length(ids_2))) # 51

subsetted_matrix <- subset_matrix(normalised_matrix, ids_2)
png(paste0(plot_path, '2_diff_accessible_annot_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- id_data %>% mutate(annotation_filter = ifelse(unique_id %in% ids_2 == TRUE, "T", "F"))

### Step 3: filter out peaks that are differentially accessible between aPPR and pPPR
ap_filter_ids <- ids_2[ids_2 %in% ss8_ap_filter_ids]
print(paste0("length of ss8_ap_filter_ids that are in ids_2: ", length(ap_filter_ids))) # 22
ids_3 <- setdiff(ids_2, ap_filter_ids)
print(paste0("length of ids_3 after filtering for ap differences at ss8: ", length(ids_3))) # 29

ap_filter_ids <- ids_2[ids_2 %in% ss4_ap_filter_ids]
print(paste0("length of ss4_ap_filter_ids that are in ids_2: ", length(ap_filter_ids))) # 22
ids_3 <- setdiff(ids_3, ap_filter_ids)
print(paste0("length of ids_3 after filtering for ap differences at ss4: ", length(ids_3))) # 29

subsetted_matrix <- subset_matrix(normalised_matrix, ids_3)
png(paste0(plot_path, '3_diff_accessible_annot_filtered_ap_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(ap_filter = ifelse(unique_id %in% ids_3 == TRUE, "T", "F"))

### Step 3: export peaks and visualise them
write.csv(id_data, file = paste0(plot_path, "putative_enhancers_table.csv"))

# make genome browser plots for open peaks
plot_path <- paste0(plot_path, "browser_tracks/")
dir.create(plot_path, recursive = T)
for (id in ids_3){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(ArchR, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}

##################### HH7 + ss4 + ss8 #######################
plot_path <- "./plots/PPR/diff_HH7_ss4_ss8/"
dir.create(plot_path, recursive = T)

### Step 1: differentially accessible in ss8 PPR clusters vs everything else
se <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("HH7_C4", "ss4_C2", "ss4_C3", "ss8_C7", "ss8_C8"))
se <- add_unique_ids_to_se(se, ArchR, matrix_type = "PeakMatrix")
ids_1 <- unique(extract_ids(se, cutOff = "FDR <= 0.01 & Log2FC >= 5", top_n = FALSE))
print(paste0("length of ids_1: ", length(ids_1))) # 52

subsetted_matrix <- subset_matrix(normalised_matrix, ids_1)
png(paste0(plot_path, '1_diff_accessible.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- as.data.frame(rowData(se)[which(rowData(se)$unique_id %in% ids_1), ])
print(dim(id_data)) #52 x 21

### Step 2: filter out peaks in genes
annot_id_data <- id_data[which(id_data$peakType %in% c("Distal", "Intronic")), ]
annot_keep_ids <- annot_id_data$unique_id
print(paste0("length of annot_keep_ids: ", length(annot_keep_ids))) # 51

ids_2 <- intersect(ids_1, annot_keep_ids)
print(paste0("length of ids_2: ", length(ids_2))) # 51

subsetted_matrix <- subset_matrix(normalised_matrix, ids_2)
png(paste0(plot_path, '2_diff_accessible_annot_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- id_data %>% mutate(annotation_filter = ifelse(unique_id %in% ids_2 == TRUE, "T", "F"))

### Step 3: filter out peaks that are differentially accessible between aPPR and pPPR
ap_filter_ids <- ids_2[ids_2 %in% ss8_ap_filter_ids]
print(paste0("length of ss8_ap_filter_ids that are in ids_2: ", length(ap_filter_ids))) # 22
ids_3 <- setdiff(ids_2, ap_filter_ids)
print(paste0("length of ids_3 after filtering for ap differences at ss8: ", length(ids_3))) # 29

ap_filter_ids <- ids_2[ids_2 %in% ss4_ap_filter_ids]
print(paste0("length of ss4_ap_filter_ids that are in ids_2: ", length(ap_filter_ids))) # 22
ids_3 <- setdiff(ids_3, ap_filter_ids)
print(paste0("length of ids_3 after filtering for ap differences at ss4: ", length(ids_3))) # 29

subsetted_matrix <- subset_matrix(normalised_matrix, ids_3)
png(paste0(plot_path, '3_diff_accessible_annot_filtered_ap_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(ap_filter = ifelse(unique_id %in% ids_3 == TRUE, "T", "F"))

### Step 3: export peaks and visualise them
write.csv(id_data, file = paste0(plot_path, "putative_enhancers_table.csv"))

# make genome browser plots for open peaks
plot_path <- paste0(plot_path, "browser_tracks/")
dir.create(plot_path, recursive = T)
for (id in ids_3){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(ArchR, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}

##################### HH7 + ss4 + ss8 - open from HH5 #######################
plot_path <- "./plots/PPR/diff_HH7_ss4_ss8_open_from_HH5/"
dir.create(plot_path, recursive = T)

### Step 1: differentially accessible in ss8 PPR clusters vs everything else
se <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("HH7_C4", "ss4_C2", "ss4_C3", "ss8_C7", "ss8_C8"))
se <- add_unique_ids_to_se(se, ArchR, matrix_type = "PeakMatrix")
ids_1 <- unique(extract_ids(se, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE))
print(paste0("length of ids_1: ", length(ids_1)))

subsetted_matrix <- subset_matrix(normalised_matrix, ids_1)
png(paste0(plot_path, '1_diff_accessible.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- as.data.frame(rowData(se)[which(rowData(se)$unique_id %in% ids_1), ])
print(dim(id_data))

### Step 2: filter out peaks in genes
annot_id_data <- id_data[which(id_data$peakType %in% c("Distal", "Intronic")), ]
annot_keep_ids <- annot_id_data$unique_id
print(paste0("length of annot_keep_ids: ", length(annot_keep_ids))) # 51

ids_2 <- intersect(ids_1, annot_keep_ids)
print(paste0("length of ids_2: ", length(ids_2))) # 51

subsetted_matrix <- subset_matrix(normalised_matrix, ids_2)
png(paste0(plot_path, '2_diff_accessible_annot_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- id_data %>% mutate(annotation_filter = ifelse(unique_id %in% ids_2 == TRUE, "T", "F"))

### Step 3: filter out peaks that are differentially accessible between aPPR and pPPR
ap_filter_ids <- ids_2[ids_2 %in% ss8_ap_filter_ids]
print(paste0("length of ss8_ap_filter_ids that are in ids_2: ", length(ap_filter_ids)))
ids_3 <- setdiff(ids_2, ap_filter_ids)
print(paste0("length of ids_3 after filtering for ap differences at ss8: ", length(ids_3)))

ap_filter_ids <- ids_2[ids_2 %in% ss4_ap_filter_ids]
print(paste0("length of ss4_ap_filter_ids that are in ids_2: ", length(ap_filter_ids)))
ids_3 <- setdiff(ids_3, ap_filter_ids)
print(paste0("length of ids_3 after filtering for ap differences at ss4: ", length(ids_3)))

subsetted_matrix <- subset_matrix(normalised_matrix, ids_3)
png(paste0(plot_path, '3_diff_accessible_annot_filtered_ap_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- id_data %>% mutate(ap_filter = ifelse(unique_id %in% ids_3 == TRUE, "T", "F"))

### Step 3: filter to only include those that open early
subsetted_raw_matrix <- subset_matrix(matrix, ids_3)
ids_4 <- open_across_stages_test(subsetted_raw_matrix, threshold_type = "min", threshold_HH5 = 1,
                                 threshold_HH6 = 1, threshold_HH7 = 1, threshold_ss4 = 1, threshold_ss8 = 1)
print(paste0("length of ids_4: ", length(ids_4)))

subsetted_matrix <- subset_matrix(normalised_matrix, ids_4)

png(paste0(plot_path, '4_diff_accessible_annot_filtered_ap_filtered_open_from_HH5.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(early_filter = ifelse(unique_id %in% ids_4 == TRUE, "T", "F"))

### Step 4: export peaks and visualise them
write.csv(id_data, file = paste0(plot_path, "putative_enhancers_table.csv"))

# make genome browser plots for open peaks
plot_path <- paste0(plot_path, "browser_tracks/")
dir.create(plot_path, recursive = T)
for (id in ids_4){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(ArchR, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}

#############################################################################
###############################   NC   #####################################
# diff accessible in NC clusters at different stages
# then filter on annotation 

##################### ss8 #######################
plot_path <- "./plots/NC/diff_ss8/"
dir.create(plot_path, recursive = T)

### Step 1: differentially accessible in ss8 NC clusters vs everything else
se <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("ss8_C1"))
se <- add_unique_ids_to_se(se, ArchR, matrix_type = "PeakMatrix")
ids_1 <- unique(extract_ids(se, cutOff = "FDR <= 0.01 & Log2FC >= 5.5", top_n = FALSE))
print(paste0("length of ids_1: ", length(ids_1))) # 33

subsetted_matrix <- subset_matrix(normalised_matrix, ids_1)
png(paste0(plot_path, '1_diff_accessible.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- as.data.frame(rowData(se)[which(rowData(se)$unique_id %in% ids_1), ])
print(dim(id_data)) #33 x 21

### Step 2: filter out peaks in genes
annot_id_data <- id_data[which(id_data$peakType %in% c("Distal", "Intronic")), ]
annot_keep_ids <- annot_id_data$unique_id
print(paste0("length of annot_keep_ids: ", length(annot_keep_ids))) # 51

ids_2 <- intersect(ids_1, annot_keep_ids)
print(paste0("length of ids_2: ", length(ids_2))) # 51

subsetted_matrix <- subset_matrix(normalised_matrix, ids_2)
png(paste0(plot_path, '2_diff_accessible_annot_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(annotation_filter = ifelse(unique_id %in% ids_2 == TRUE, "T", "F"))

### Step 3: export peaks and visualise them
write.csv(id_data, file = paste0(plot_path, "putative_enhancers_table.csv"))

# make genome browser plots for open peaks
plot_path <- paste0(plot_path, "browser_tracks/")
dir.create(plot_path, recursive = T)
for (id in ids_2){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(ArchR, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}

##################### ss4 + ss8 #######################
plot_path <- "./plots/NC/diff_ss4_ss8/"
dir.create(plot_path, recursive = T)

### Step 1: differentially accessible in ss8 NC clusters vs everything else
se <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("ss4_C6", "ss8_C1"))
se <- add_unique_ids_to_se(se, ArchR, matrix_type = "PeakMatrix")
ids_1 <- unique(extract_ids(se, cutOff = "FDR <= 0.01 & Log2FC >= 2", top_n = FALSE))
print(paste0("length of ids_1: ", length(ids_1))) # 37

subsetted_matrix <- subset_matrix(normalised_matrix, ids_1)
png(paste0(plot_path, '1_diff_accessible.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- as.data.frame(rowData(se)[which(rowData(se)$unique_id %in% ids_1), ])
print(dim(id_data)) #37 x 21

### Step 2: filter out peaks in genes
annot_id_data <- id_data[which(id_data$peakType %in% c("Distal", "Intronic")), ]
annot_keep_ids <- annot_id_data$unique_id
print(paste0("length of annot_keep_ids: ", length(annot_keep_ids))) # 35

ids_2 <- intersect(ids_1, annot_keep_ids)
print(paste0("length of ids_2: ", length(ids_2))) # 35

subsetted_matrix <- subset_matrix(normalised_matrix, ids_2)
png(paste0(plot_path, '2_diff_accessible_annot_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- id_data %>% mutate(annotation_filter = ifelse(unique_id %in% ids_2 == TRUE, "T", "F"))

### Step 3: filter to only include those that open at ss4
subsetted_raw_matrix <- subset_matrix(matrix, ids_2)
ids_3 <- open_across_stages_test(subsetted_raw_matrix, threshold_type = "min", threshold_HH5 = 0.001,
                                 threshold_HH6 = 0.001, threshold_HH7 = 0.001, threshold_ss4 = 1, threshold_ss8 = 1)
print(paste0("length of ids_3: ", length(ids_3))) #31

subsetted_matrix <- subset_matrix(normalised_matrix, ids_3)
png(paste0(plot_path, '3_diff_accessible_annot_filtered_open_from_ss4.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(early_filter = ifelse(unique_id %in% ids_3 == TRUE, "T", "F"))

### Step 3: export peaks and visualise them
write.csv(id_data, file = paste0(plot_path, "putative_enhancers_table.csv"))

# make genome browser plots for open peaks
plot_path <- paste0(plot_path, "browser_tracks/")
dir.create(plot_path, recursive = T)
for (id in ids_3){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(ArchR, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}

##################### HH7 + ss4 + ss8 #######################
plot_path <- "./plots/NC/diff_HH7_ss4_ss8/"
dir.create(plot_path, recursive = T)

### Step 1: differentially accessible in ss8 NC clusters vs everything else
se <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "stage_clusters",
  useGroups = c("HH7_C5", "ss4_C6", "ss8_C1"))
se <- add_unique_ids_to_se(se, ArchR, matrix_type = "PeakMatrix")
ids_1 <- unique(extract_ids(se, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE))
print(paste0("length of ids_1: ", length(ids_1))) # 37

subsetted_matrix <- subset_matrix(normalised_matrix, ids_1)
png(paste0(plot_path, '1_diff_accessible.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- as.data.frame(rowData(se)[which(rowData(se)$unique_id %in% ids_1), ])
print(dim(id_data)) #37 x 21

### Step 2: filter out peaks in genes
annot_id_data <- id_data[which(id_data$peakType %in% c("Distal", "Intronic")), ]
annot_keep_ids <- annot_id_data$unique_id
print(paste0("length of annot_keep_ids: ", length(annot_keep_ids))) # 35

ids_2 <- intersect(ids_1, annot_keep_ids)
print(paste0("length of ids_2: ", length(ids_2))) # 35

subsetted_matrix <- subset_matrix(normalised_matrix, ids_2)
png(paste0(plot_path, '2_diff_accessible_annot_filtered.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = FALSE))
graphics.off()

id_data <- id_data %>% mutate(annotation_filter = ifelse(unique_id %in% ids_2 == TRUE, "T", "F"))

### Step 3: filter to only include those that open at HH7
subsetted_raw_matrix <- subset_matrix(matrix, ids_2)
ids_3 <- open_across_stages_test(subsetted_raw_matrix, threshold_type = "min", threshold_HH5 = 0.001,
                                 threshold_HH6 = 0.001, threshold_HH7 = 1, threshold_ss4 = 1, threshold_ss8 = 1)
print(paste0("length of ids_3: ", length(ids_3))) #36

subsetted_matrix <- subset_matrix(normalised_matrix, ids_3)
png(paste0(plot_path, '3_diff_accessible_annot_filtered_open_from_HH7.png'), height = 20, width = 30, units = 'cm', res = 400)
print(marker_heatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

id_data <- id_data %>% mutate(early_filter = ifelse(unique_id %in% ids_3 == TRUE, "T", "F"))

### Step 3: export peaks and visualise them
write.csv(id_data, file = paste0(plot_path, "putative_enhancers_table.csv"))

# make genome browser plots for open peaks
plot_path <- paste0(plot_path, "browser_tracks/")
dir.create(plot_path, recursive = T)
for (id in ids_3){
  print(id)
  gr <- make_gr_object(id = id, extend = TRUE, extend_by = 10000)
  p <- plotBrowserTrack(ArchR, region = gr, groupBy = "stage_clusters", baseSize = 20, facetbaseSize = 20,
                        plotSummary = c("bulkTrack", "featureTrack", "geneTrack"), sizes = c(10, 1.5, 2),
                        title = paste0("Peak ID:", id))
  
  name <- str_replace(id, ":", "-")
  png(paste0(plot_path, name, '_extended_by_10000.png'), height = 50, width = 50, units = 'cm', res = 400)
  grid::grid.draw(p)
  graphics.off()
}