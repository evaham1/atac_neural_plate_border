#!/usr/bin/env Rscript

### script to investigate the loops made by HiCDCPlus
print("script to investigate the loops made by HiCDCPlus")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(parallel)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(HiCDCPlus)

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
    
    plot_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output/plots/"
    rds_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output/rds_files/"
    data_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output/rds_files/" # for HiCDCPlus output
    data_path = "./output/NF-hichip-downstream/bins/genes_intersect" # for filtering by genes
    data_path = "./output/NF-hichip-downstream/bins/peaks_intersect" # for filtering by peaks
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


######################################   Read in data   ####################################################

# read in HiCDCPlus output


# read in genes intersected with bins


# read in peaks intersected with bins


######################################   Investigate loops   ####################################################

