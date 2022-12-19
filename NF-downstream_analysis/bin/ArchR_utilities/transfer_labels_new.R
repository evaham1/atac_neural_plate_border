#!/usr/bin/env Rscript

print("ArchR transfer cluster IDs from stage data to full data")
# transfer labels new

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(gridExtra)
library(grid)
library(parallel)

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
    addArchRThreads(threads = 1) 
    
    data_path = "./output/NF-downstream_analysis/test_inputs/test_input_transfer_labels/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/"
    rds_path = "./rds_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Read in ArchR projects #######################################
# Read in all data
files <- list.files(data_path, full.names = TRUE)
print(paste0("All input files: ", files))

stages_data <- grep("FullData", files, invert = T, value = TRUE) # source data from which labels are extracted
print(paste0("Stages data: ", stages_data))

full_data <- grep("FullData", files, invert = F, value = TRUE) # destination data where labels will be transfered onto
print(paste0("Destination data: ", full_data))
ArchR_full <- loadArchRProject(path = full_data, force = FALSE, showLogo = TRUE)

################## Transfer info from stages onto full data #######################################
# extract cell IDs from each of the stages datasets
all_cell_ids <- c()
all_cluster_ids <- c()
for (i in 1:5) {
  stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
  print(paste0("Iteration: ", i, ", Stage: ", stage))
  
  ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
  cell_ids <- ArchR$cellNames
  print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
  all_cell_ids <- c(all_cell_ids, cell_ids)
  
  stage <- unique(ArchR$stage)
  cluster_ids <- paste0(stage, "_", ArchR$clusters)
  print(length(cluster_ids))
  all_cluster_ids <- c(all_cluster_ids, cluster_ids)
}

print(paste0("Length of all cell ids: ", length(all_cell_ids)))
print(paste0("Length of all cluster ids: ", length(all_cluster_ids)))

data <- as.vector(all_cluster_ids)
cells <- as.vector(all_cell_ids)

ArchR <- addCellColData(ArchRProj = ArchR_full, 
                        data = data,
                        cells = cells, 
                        name = "stage_clusters",
                        force = TRUE)
ArchR <- ArchR_full
ArchR@cellColData[name] <- NA
ArchR@cellColData[cells, name] <- data
ArchR@cellColData[cells, "stage"]

ArchR@cellColData$stage
sum(is.na(ArchR$cellNames))
length(ArchR$cellNames)

print(table(ArchR$stage_clusters))

# save filtered data
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, "FullData_Save-ArchR"), load = FALSE)


#####################################################################################
############################## Visualisations #######################################
###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

p1 <- plotEmbedding(ArchR, 
                    name = "stage",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    pal = stage_colours, randomize = TRUE)
p2 <- plotEmbedding(ArchR, 
                    name = "stage_clusters",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0,
                    randomize = TRUE)

png(paste0(plot_path, "UMAPs.png"), width=60, height=40, units = 'cm', res = 200)
ggAlignPlots(p1, p2, type = "h")
graphics.off()