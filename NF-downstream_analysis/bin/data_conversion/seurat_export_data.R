#!/usr/bin/env Rscript

print("Export data from seurat object so can downgrade it for MEGA compatability")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(Seurat)
library(Signac)

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
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/rds_files/"
    rds_path = "./rds_files/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

set.seed(42)

############################## Read in seurat objects #######################################

## read in paired seurat object
obj.pair <- readRDS(paste0(data_path, "paired_object_chromvar.RDS"))
obj.pair

############################## Extract data #######################################

rna_data <- obj.pair$RNA
atac_data <- obj.pair$ATAC
chromvar_data <- obj.pair$chromvar
metadata <- obj.pair[[]]

############################## Save data #######################################
saveRDS(rna_data, paste0(rds_path, "rna_data.RDS"), compress = FALSE)
saveRDS(atac_data, paste0(rds_path, "atac_data.RDS"), compress = FALSE)
saveRDS(chromvar_data, paste0(rds_path, "chromvar_data.RDS"), compress = FALSE)
saveRDS(metadata, paste0(rds_path, "metadata.RDS"), compress = FALSE)