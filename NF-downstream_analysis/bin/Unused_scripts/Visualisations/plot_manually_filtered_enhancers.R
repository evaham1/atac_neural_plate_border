#!/usr/bin/env Rscript

# make visualisations of enhancers that were manually filtered (need to hard coded in)

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
    data_path = "./NF-downstream_analysis/work/17/d51f87f947475ee9e5e230be2c5426/rds_files/" #need to update
    plot_path = "./output/NF-downstream_analysis/Downstream_processing/transfer_labels/plot_manually_filtered_enhancers/plots/"
    se_path = "./NF-downstream_analysis/work/17/d51f87f947475ee9e5e230be2c5426/rds_files/" 
    
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

# Extract matrix for heatmap plotting
matrix <- extract_means_from_se(Full_se)
normalised_matrix <- Log2norm(matrix)


##################################################################################################
##########################    Read in manually selected peaks   ##################################

NC_peak_ids <- c("chr1:4741943-4742443", "chr1:36224641-36225141",
                 "chr18:8888066-8888566", "chr2:41875394-41875894",
                 "chr2:42661874-42662374", "chr2:140273671-140274171",
                 "chr4:70459983-70460483", "chr8:7216956-7217456",
                 "chr3:60034724-60035224", "chr4:51680702-51681202",
                 "chr9:5834645-5835145", "chr1:57370425-57370925",
                 "chr1:61674271-61674771", "chr4:39407032-39407532")

PPR_peak_ids <- c("chr2:10216477-10216977", "chr2:11936117-11936617",
                  "chr2:45931424-45931924", "chr21:5992393-5992893",
                  "chr9:16689582-16690082", "chr1:44512244-44512744",
                  "chr12:1723295-1723795", "chr4:69910791-69911291",
                  "chr8:25914832-25915332", "chr1:88390937-88391437",
                  "chr2:149406064-149406564", "chr4:5183754-5184254",
                  "chr4:11263569-11264069")


###########################################################################
##########################    Heatmaps   ##################################

## NC
NC_matrix <- subset_matrix(normalised_matrix, NC_peak_ids)

png(paste0(plot_path, 'NC/NC_heatmap_clustered.png'), height = 20, width = 40, units = 'cm', res = 400)
print(marker_heatmap(NC_matrix, pal = pal, clusterCols = TRUE, labelRows = TRUE))
graphics.off()

png(paste0(plot_path, 'NC/NC_heatmap_unclustered.png'), height = 20, width = 40, units = 'cm', res = 400)
print(marker_heatmap(NC_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

## PPR
PPR_matrix <- subset_matrix(normalised_matrix, PPR_peak_ids)

png(paste0(plot_path, 'PPR/PPR_heatmap_clustered.png'), height = 20, width = 40, units = 'cm', res = 400)
print(marker_heatmap(PPR_matrix, pal = pal, clusterCols = TRUE, labelRows = TRUE))
graphics.off()

png(paste0(plot_path, 'PPR/PPR_heatmap_unclustered.png'), height = 20, width = 40, units = 'cm', res = 400)
print(marker_heatmap(PPR_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
graphics.off()

####################################################################################
##########################    Genome Browser Plots   ###############################

## create new function for this and add here

######################################################################################
##########################    Annotation Piechart   ##################################
Full_peaks_data <- getPeakSet(FullData)

## NC
NC_peaks_data <- Full_peaks_data[which(Full_peaks_data$name %in% gsub(":", "-", NC_peak_ids)), ]

counts <- as.data.frame(table(NC_peaks_data$peakType))
colnames(counts) <- c("Peak Annotation", "Number of peaks")
counts <- counts %>% mutate('Percentage' = round(100*(`Number of peaks`/sum(counts$`Number of peaks`))))

png(paste0(plot_path, 'NC/peak_counts_per_type.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Peak Counts per type", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path, 'NC/peak_counts_per_type_piechart.png'), height = 10, width = 20, units = 'cm', res = 400)
ggplot(counts, aes(x="", y=`Number of peaks`, fill=`Peak Annotation`)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  geom_text(aes(label = paste0(Percentage, "%")),
            position = position_stack(vjust = 0.5))
graphics.off()

## PPR
PPR_peaks_data <- Full_peaks_data[which(Full_peaks_data$name %in% gsub(":", "-", PPR_peak_ids)), ]

counts <- as.data.frame(table(PPR_peaks_data$peakType))
colnames(counts) <- c("Peak Annotation", "Number of peaks")
counts <- counts %>% mutate('Percentage' = round(100*(`Number of peaks`/sum(counts$`Number of peaks`))))

png(paste0(plot_path, 'PPR/peak_counts_per_type.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Peak Counts per type", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path, 'PPR/peak_counts_per_type_piechart.png'), height = 10, width = 20, units = 'cm', res = 400)
ggplot(counts, aes(x="", y=`Number of peaks`, fill=`Peak Annotation`)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  geom_text(aes(label = paste0(Percentage, "%")),
            position = position_stack(vjust = 0.5))
graphics.off()

##########################################################################################
##########################    Distance to nearest TSSs   ##################################

## NC
NC_peaks_data <- as.data.frame(NC_peaks_data, row.names = NULL)

png(paste0(plot_path, 'NC/dist_to_TSS.png'), height = 30, width = 20, units = 'cm', res = 400)
ggplot(NC_peaks_data, aes(x = peakType, y = distToTSS)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  geom_hline(aes(yintercept = 50000), color = "black", linetype = "dashed", size = 1)
graphics.off()

NC_peaks_data_minus_outlier <- NC_peaks_data[which(NC_peaks_data$distToTSS < 50000), ]

png(paste0(plot_path, 'NC/dist_to_TSS_minus_outlier.png'), height = 30, width = 20, units = 'cm', res = 400)
ggplot(NC_peaks_data_minus_outlier, aes(x = peakType, y = distToTSS)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  geom_hline(aes(yintercept = 500), color = "black", linetype = "dashed", size = 1)
graphics.off()

## PPR
PPR_peaks_data <- as.data.frame(PPR_peaks_data, row.names = NULL)

png(paste0(plot_path, 'PPR/dist_to_TSS.png'), height = 30, width = 20, units = 'cm', res = 400)
ggplot(PPR_peaks_data, aes(x = peakType, y = distToTSS)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  geom_hline(aes(yintercept = 50000), color = "black", linetype = "dashed", size = 1)
graphics.off()
