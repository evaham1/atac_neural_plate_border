#!/usr/bin/env Rscript

print("differential peaks ArchR - just calculate no visualisations")
# calculates differential peaks between clusters, scHelper_cell_types_old and stage + few basic heatmaps

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
    
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/QC_HIGH/ss8/prefiltering/peak_calling/rds_files/"
    
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

### function to turn summarised experiment object into a df that can be printed to see diff features
extract_features_table <- function(markersPeaks) {
  markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 1") # keep all peaks
  df <- data.frame()
  for (i in 1:length(names(markerList))) {
    print(i)
    df_i <- as.data.frame(markerList[i])
    df <- rbind(df, df_i)
  }
  return(df)
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

getPeakSet(ArchR)
getAvailableMatrices(ArchR)

############################# Calculate Diff Peaks between cell populations #######################################

markersPeaks <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = opt$group_by)
markersPeaks

peaks_table <- extract_features_table(markersPeaks)
write.csv(peaks_table, paste0(plot_path, "all_differentially_expressed_peaks.csv"), row.names = FALSE)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0",
  nLabel = 3)
png(paste0(plot_path, 'diff_peak_cutoff_heatmap_4.png'), height = 50, width = 40, units = 'cm', res = 400)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
graphics.off()

############################# Add peak information to markersPeaks object #######################################

tmp_peaks = data.frame(ArchR@peakSet)
tmp_diff_peaks = data.frame(rowData(markersPeaks))

diff_peaks_join_peakset = left_join(tmp_diff_peaks, tmp_peaks, by = c("seqnames" = "seqnames", "start" = "start", "end" = "end"))
diff_peaks_join_peakset$name = paste(diff_peaks_join_peakset$nearestGene, diff_peaks_join_peakset$distToTSS,sep="_")
diff_peaks_join_peakset$unique_id = paste(diff_peaks_join_peakset$seqnames, diff_peaks_join_peakset$start, diff_peaks_join_peakset$end, sep=":")
rowData(markersPeaks) = diff_peaks_join_peakset

############################# Promoter Peaks #######################################

markersPeaks_promoter <- subset(markersPeaks, rowData(markersPeaks)$peakType == "Promoter")
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

markersPeaks_distal <- subset(markersPeaks, rowData(markersPeaks)$peakType == "Distal")
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