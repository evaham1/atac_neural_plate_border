#!/usr/bin/env Rscript

print("calculates SE object by running 'getMarkerFeatures'")
# calculates SE object by running 'getMarkerFeatures', this is useful to save time downstream!
# outputs the se object and the original ArchR object

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
  make_option(c("-g", "--group_by"), action = "store", type = "character", help = "How cells are grouped when peaks called", default = "clusters",),
  make_option(c("-v", "--verbose"), action = "store", type = "logical", help = "Verbose", default = FALSE)
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
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/rds_files/"
    rds_path = "./rds_files/"
    ncores = opt$cores
    
    # addArchRThreads(threads = ncores)
    addArchRThreads(threads = 1) 
    
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

###########################################################################################
############################## Read in ArchR project #####################################

# If files are not in rds_files subdirectory look in input dir
label <- unique(sub('_.*', '', list.files(data_path)))
print(label)

if (length(label) == 0){
  data_path = "./input/"
  label <- unique(sub('_.*', '', list.files(data_path)))
  print(label)
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
} else {
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
}

getAvailableMatrices(ArchR)
ArchR@peakSet

###########################################################################################
############################## Write out ArchR project #####################################

paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)
print("ArchR object saved")

###########################################################################################
########################## Calculate se across all clusters ###############################

print("Calculating se across all cell groups...")
print(paste0("Cells grouped by: ", opt$group_by))

se <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = opt$group_by)
se <- scHelper::ArchRAddUniqueIdsToSe(se, ArchR, matrix_type = "PeakMatrix")

saveRDS(se, file = paste0(rds_path, label, "_SE.RDS"))

print("se RDS saved")

###########################################################################################
########################## Extract table of differentially accessible peaks ###############################

peaks_table <- extract_features_table(se)
write.csv(peaks_table, paste0(plot_path, "all_differentially_accessible_peaks.csv"), row.names = FALSE)

print("diff peaks table saved")