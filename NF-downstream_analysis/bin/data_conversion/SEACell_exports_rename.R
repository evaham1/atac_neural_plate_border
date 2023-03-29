#!/usr/bin/env Rscript

print("take .csv files from SEACell computation (cell_metadata.csv, feature_metadata.csv and summarised_by_metacell_counts.csv)")
print("and rename them based on stage extracted from cell_metadata.csv")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(tidyverse)
library(plyr)
library(dplyr)
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
    
    # local interactive paths
    data_path = "./local_test_data/SEACells_from_camp/SEACELLS_ATAC_WF/2_SEACells_computation/exported_data/"
    rds_path = "./local_test_data/SEACells_outputs_renamed/csv_files/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    rds_path = "./csv_files/"
    data_path = "./input/exported_data/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(rds_path, recursive = T)
}


############################## Read in .csv files #######################################
input_files <- list.files(path = data_path, pattern = "*.csv", full.names = TRUE)
print(input_files)

############################## Extract stage name #######################################
cell_metadata <- read.csv(paste0(data_path, "cell_metadata.csv"))
stage <- substr(cell_metadata$index[1], 8, 10)
print(paste0("Stage detected: ", stage))

############################## Read in and rewrite .csv files #######################################

write.csv(cell_metadata, paste0(rds_path, stage, '_cell_metadata.csv'))

feature_metadata <- read.csv(paste0(data_path, "Feature_metadata.csv"))
write.csv(feature_metadata, paste0(rds_path, stage, '_feature_metadata.csv'))

counts <- read.csv(paste0(data_path, "Summarised_by_metacells_counts.csv"))
write.csv(counts, paste0(rds_path, stage, '_summarised_by_metacells_counts.csv'))



