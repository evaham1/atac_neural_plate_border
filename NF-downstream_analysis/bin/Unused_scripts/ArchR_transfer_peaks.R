#!/usr/bin/env Rscript

print("transfer peak set from one ArchR object onto another ArchR object")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)
library(clustree)
library(ComplexHeatmap)
library(BSgenome.Ggallus.UCSC.galGal6)
library(scHelper)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-g", "--group_by"), action = "store", type = "character", help = "How to group cells to call peaks", default = "clusters",),
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
    
    # peaks not called, to save time make smaller data
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Clustering/rds_files/"
    
    # rds_path
    rds_path = "./output/NF-downstream_analysis/Processing/ss8/Peak_call/test/rds_files/"
    
    # once read in ArchR run this to speed it up
    ArchR_backup <- ArchR
    idxSample <- BiocGenerics::which(ArchR$clusters %in% c("C10"))
    cellsSample <- ArchR$cellNames[idxSample]
    ArchR <- ArchR[cellsSample, ]
    
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

############################## Read in ArchR project #######################################

# Extract stage name
label <- unique(sub('_.*', '', list.files(data_path)))
print(paste0("Label: ", label))

# Read in target ArchR object (in rds_files because it came from previous process)
print("reading in target ArchR object in rds_files folder")
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getPeakSet(ArchR)
getAvailableMatrices(ArchR)

# Read in donor ArchR object (in input folder)
print("reading in donor ArchR object in input folder")
ArchR_donor <- loadArchRProject(path = paste0("./input/FullData_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# see what is in the ArchR donor object already
print("ArchR donor object info: ")
print(ArchR_donor)
getPeakSet(ArchR_donor)
getAvailableMatrices(ArchR_donor)

# how cells are grouped for pseudobulk replicates + peak calling
print(paste0("Cells grouped by: ", opt$group_by))

##########################################################################################
############################## Transfer peak set onto ArchR ###############################

# extract peakset granges from donor ArchR object
peakset_granges <- getPeakSet(ArchR_donor)

# add the peak set to the target ArchR object and calculate the matrix
ArchR_peaks <- addPeakSet(ArchR, peakset_granges, force = TRUE)
ArchR_peaks <- addPeakMatrix(ArchR_peaks, force = TRUE)

#################################################################################
############################## Save ArchR project ###############################

output_directory <- paste0(rds_path, label, "_Save-ArchR")
print(paste0("Output directory of object with peak set transferred onto it: ", output_directory))
saveArchRProject(ArchRProj = ArchR_peaks, outputDirectory = output_directory, load = FALSE)

########################################################################
##################### How many cut sites per peak ######################

plot_path_temp <- paste0(plot_path, "cut_sites_per_peak/")
dir.create(plot_path_temp, recursive = T)

peak_data <- getMatrixFromProject(ArchR_peaks, useMatrix = "PeakMatrix")
peak_matrix <- t(assays(peak_data)[[1]])
colnames(peak_matrix) <- rowData(peak_data)$name

png(paste0(plot_path_temp, 'cutsites_per_peak_histogram.png'), height = 20, width = 40, units = 'cm', res = 400)
hist(colSums(peak_matrix), breaks = 1000, main = "Histogram of cut sites per peak", 
     xlab = "Number of cut sites", ylab = "Frequency")
graphics.off()

png(paste0(plot_path_temp, 'cutsites_per_peak_boxplot.png'), height = 30, width = 15, units = 'cm', res = 400)
boxplot(colSums(peak_matrix), breaks = 1000, main = "Boxplot of cut sites per peak", 
        xlab = "Number of cut sites", ylab = "Frequency")
graphics.off()

cutsites <- summary(colSums(peak_matrix))
table <- data.frame(Stats = names(cutsites), Value = as.vector(cutsites))
png(paste0(plot_path_temp, 'cutsites_per_peak_table.png'), height = 30, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Cut sites per peak", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(table, rows=NULL, theme = ttheme_minimal()))
graphics.off()

print("cut sites per peak calculated")