#!/usr/bin/env Rscript

print("make TxDB file from gtf")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(gridExtra)
library(grid)
library(parallel)
library(Rsamtools)

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

    data_path = "./output/NF-downstream_analysis/Upstream_processing//PREPROCESSING/edit_gtf/"

  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    data_path = "./input/"
    rds_path = "./txdb/"
    
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(rds_path, recursive = T)
}

set.seed(42)


#################### R script to take gtf and turn into txdb #########################

txdb <- makeTxDbFromGFF(paste0(data_path, "galgal6/tag_chroms.gtf"))
print("txdb made")

print(genes(txdb))

saveDb(txdb, paste0(rds_path, "txdb"))