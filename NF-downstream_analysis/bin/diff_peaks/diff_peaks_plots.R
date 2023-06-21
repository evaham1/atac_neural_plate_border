#!/usr/bin/env Rscript

print("late differences across all timepoints")
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
library(scHelper)

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

### Function to Log2normalise matrix (must be run before subsetting!)
Log2norm <- function(mat, scaleTo = 10^4) {
  mat <- log2(t(t(mat)/colSums(mat)) * scaleTo + 1) # normalising means for depth of cluster
  return(mat)
}

### Function to subset normalised matrix using IDs
subset_matrix <- function(mat, ids) {
  subsetted_matrix <- mat[ids, ]
  return(subsetted_matrix)
}

###########################################################################################
############################## Read in data #####################################

# Retrieve object label
label <- unique(sub('_.*', '', list.files(data_path)))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
print('data read in')
print(ArchR)

# check that gene score matrix and gene integration matrix are available
getAvailableMatrices(ArchRProj = ArchR)

# read in se object
se <- readRDS(paste0(data_path, label, "_SE.RDS"))

print("data read in!")

############################################################################################
############################## COLOURS #######################################

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

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
# set colour palettes for UMAPs
atac_scHelper_old_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_old)]

# set pal for heatmaps
pal = paletteContinuous(set = "solarExtra", n = 100)

##########################################################################
############################## PLOTS #####################################

##### Prepare data for plotting 
matrix <- scHelper::ArchR_ExtractMeansFromSe(se) # extract means df from se object
normalised_matrix <- Log2norm(matrix) # log2norm across all features in each cell group


##### Plot all differential peaks
ids <- scHelper::ArchR_ExtractIds(se, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE)
print(paste0("all peaks: ", length(ids)))

if (length(ids) > 4){
  matrix <- ArchR_ExtractMeansFromSe(se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)

  png(paste0(plot_path, 'full_heatmap.png'), height = 70, width = 60, units = 'cm', res = 400)
  print(scHelper::ArchR_PlotMarkerHeatmap(subsetted_matrix, pal = pal, clusterCols = FALSE))
  graphics.off()
}

ids <- scHelper::ArchR_ExtractIds(se, cutOff = "FDR <= 0.05 & Log2FC >= 0", top_n = TRUE)
print(paste0("all peaks top 10: ", length(ids)))

if (length(ids) > 4){
  matrix <- scHelper::ArchR_ExtractMeansFromSe(se)
  normalised_matrix <- Log2norm(matrix)
  subsetted_matrix <- subset_matrix(normalised_matrix, ids)

  png(paste0(plot_path, 'full_heatmap_top10.png'), height = 40, width = 30, units = 'cm', res = 400)
  print(scHelper::ArchR_PlotMarkerHeatmap(subsetted_matrix, pal = pal, clusterCols = FALSE, labelRows = TRUE))
  graphics.off()
}

############################# Add peak information to markersPeaks object #######################################

tmp_peaks = data.frame(ArchR@peakSet)
tmp_diff_peaks = data.frame(rowData(se))

diff_peaks_join_peakset = left_join(tmp_diff_peaks, tmp_peaks, by = c("seqnames" = "seqnames", "start" = "start", "end" = "end"))
diff_peaks_join_peakset$name = paste(diff_peaks_join_peakset$nearestGene, diff_peaks_join_peakset$distToTSS,sep="_")
diff_peaks_join_peakset$unique_id = paste(diff_peaks_join_peakset$seqnames, diff_peaks_join_peakset$start, diff_peaks_join_peakset$end, sep=":")
rowData(se) = diff_peaks_join_peakset

############################# Promoter Peaks #######################################

markersPeaks_promoter <- subset(se, rowData(se)$peakType == "Promoter")
marker_tables_promoter = markersPeaks_promoter %>% getMarkers(cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
marker_tables_promoter_tmp = marker_tables_promoter %>%  as.data.frame()
write.csv(marker_tables_promoter_tmp, paste0(plot_path, "all_differentially_expressed_peaks_promoter.csv"), row.names = FALSE)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks_promoter,
  cutOff = "FDR <= 0.3 & Log2FC >= 0.5",
  nLabel = 3)
png(paste0(plot_path, 'diff_promoter_peak_cutoff_heatmap.png'), height = 50, width = 40, units = 'cm', res = 400)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
graphics.off()

############################# Distal Peaks #######################################

markersPeaks_distal <- subset(se, rowData(se)$peakType == "Distal")
marker_tables_distal = markersPeaks_distal %>% getMarkers(cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
marker_tables_distal_tmp = marker_tables_distal %>%  as.data.frame()
write.csv(marker_tables_distal_tmp, paste0(plot_path, "all_differentially_expressed_peaks_distal.csv"), row.names = FALSE)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks_distal,
  cutOff = "FDR <= 0.3 & Log2FC >= 0.5",
  nLabel = 3)
png(paste0(plot_path, 'diff_distal_peak_cutoff_heatmap.png'), height = 50, width = 40, units = 'cm', res = 400)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
graphics.off()

###################### Plots showing distribution of FDR and Logf2c values #############################

png(paste0(plot_path, 'Log2FC_boxplot.png'), height = 20, width = 20, units = 'cm', res = 400)
boxplot(assays(se)$Log2FC)
graphics.off()

png(paste0(plot_path, 'FDR_boxplot.png'), height = 20, width = 20, units = 'cm', res = 400)
boxplot(assays(se)$FDR)
graphics.off()

df <- data.frame(LogFC = c(t(assays(se)$Log2FC)), FDR = c(t(assays(se)$FDR)), 
                 LogFDR = log10(c(t(assays(se)$FDR))), stringsAsFactors=FALSE)

df <- df %>% mutate(Passed = as.factor(ifelse(FDR < 0.01 & LogFC > 1,"passed", "failed")))

set.seed(42)
rows <- sample(nrow(df))
df <- df[rows, ]

png(paste0(plot_path, 'FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(df, aes(x = -LogFDR, y = LogFC, color = Passed, shape = Passed)) + 
  geom_point() + 
  scale_color_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(16, 17))
graphics.off()