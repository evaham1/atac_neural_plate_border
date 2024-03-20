#!/usr/bin/env Rscript

print("Transfer latent time and lineage probabilities from single cell ATAC object to metacell ATAC metadata by averaging them per metacell")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(Seurat)
library(SummarizedExperiment)

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
    
    # data_path = "./output/NF-downstream_analysis/Processing/FullData/Single_cell_integration/" # ATAC
    # data_path = "./output/NF-downstream_analysis/rna_objects/" # RNA
    # rds_path = "./output/NF-downstream_analysis/Processing/FullData/Transfer_latent_time/"
    
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

set.seed(42)

############################## Read in data #######################################

# read in ATAC single cell data
ArchR <- loadArchRProject(path = paste0(data_path, "rds_files/", "TransferLabels_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# read in ATAC metacell metadata csv
metadata <- read.csv(paste0(data_path, "Combined_SEACell_integrated_metadata.csv"), row.names = 'ATAC')

# Add stage to metadata using SEACell IDs
substrRight <- function(x, n){
  sapply(x, function(xx)
    substr(xx, (nchar(xx)-n+1), nchar(xx))
  )
}
metadata <- metadata %>% mutate(stage = substrRight(rownames(metadata), 3))
metadata <- metadata[,-1]

# Change cell names to match matrix
rownames(metadata) <- gsub('-', '_', rownames(metadata))

# Check metadata
print(head(metadata))

print("Data read in!")