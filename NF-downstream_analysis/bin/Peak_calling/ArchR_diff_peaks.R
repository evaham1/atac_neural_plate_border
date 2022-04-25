#!/usr/bin/env Rscript

print("differential peaks ArchR")

############################## Load libraries #######################################
library(getopt)
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
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    
    data_path = "./output/NF-downstream_analysis/ArchR_peak_calling/ss8/rds_files/"
    plot_path = "./output/NF-downstream_analysis/ArchR_peak_calling/ss8/2_differential_peaks/plots/"
    
    
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

# ############################## Subset to speed up interactive work #######################################
# ArchR_big <- ArchR
# 
# idxSample <- BiocGenerics::which(ArchR_big$scHelper_cell_type_old %in% c("NC", "dNC", "aPPR", "pPPR"))
# cellsSample <- ArchR_big$cellNames[idxSample]
# 
# ArchR <- subsetArchRProject(ArchR_big, cellsSample, outputDirectory = "ArchRSubset",
#   dropCells = TRUE, logFile = NULL, force = TRUE)
# 
# ArchR <- loadArchRProject(path = "./ArchRSubset", force = TRUE, showLogo = TRUE)


############################# CLUSTER Marker Peaks #######################################

markersPeaks <- getMarkerFeatures(
    ArchRProj = ArchR, 
    useMatrix = "PeakMatrix", 
    groupBy = "clusters")
markersPeaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.3 & Log2FC >= 0.5",
  transpose = TRUE)

png(paste0(plot_path, 'clusters_diff_peak_FDR<0.3_Log2FC>0.5_heatmap.png'), height = 30, width = 40, units = 'cm', res = 400)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
graphics.off()

############################# SC_HELPER_CELL_TYPE Marker Peaks #######################################

ArchR$scHelper_cell_type_old <- as.character(ArchR$scHelper_cell_type_old)

markersPeaks <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "scHelper_cell_type_old")
markersPeaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.3 & Log2FC >= 0.5",
  transpose = TRUE)

  png(paste0(plot_path, 'stage_diff_peak_FDR<0.3_Log2FC>0.5_heatmap.png'), height = 30, width = 40, units = 'cm', res = 400)
  print(draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot"))
  graphics.off()


############################# STAGE Marker Peaks #######################################

if (length(unique(ArchR$stage)) > 1) {
  markersPeaks <- getMarkerFeatures(
    ArchRProj = ArchR, 
    useMatrix = "PeakMatrix", 
    groupBy = "stage")
  markersPeaks

  markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
  markerList

  heatmapPeaks <- plotMarkerHeatmap(
    seMarker = markersPeaks, 
    cutOff = "FDR <= 0.3 & Log2FC >= 0.5",
    transpose = TRUE)

  png(paste0(plot_path, 'schelper_diff_peak_FDR<0.3_Log2FC>0.5_heatmap.png'), height = 30, width = 40, units = 'cm', res = 400)
  print(draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot"))
  graphics.off()

}

# ############################## Individual group plots #######################################

# ### make for loop so make MA/volcano plot for each group in group_by
# for (i in )
# pma <- plotMarkers(seMarker = markersPeaks, name = "aPPR", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
# pma

# pv <- plotMarkers(seMarker = markersPeaks, name = "C6", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
# pv


# ############################## Browser tracks #######################################

# p <- plotBrowserTrack(
#   ArchRProj = ArchR, 
#   groupBy = "clusters", 
#   geneSymbol = c("SIX1"),
#   features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["C6"],
#   upstream = 50000,
#   downstream = 50000
# )

# grid::grid.newpage()
# grid::grid.draw(p$SIX1)