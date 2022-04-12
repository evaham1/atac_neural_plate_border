#!/usr/bin/env Rscript

print("ArchR_subset")

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

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
    make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
    make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
    make_option(c("-m", "--meta_col1"), action = "store", type = "character", help = "Name of first metadata column containing groups to subset", default = NULL),
    #make_option(c("", "--meta_col2"), action = "store", type = "character", help = "Name of second metadata column containing groups to subset", default = NULL),
    #make_option(c("-o", "--output"), action = "store", type = "character", help = "Name of output ArchR project", default = 'ArchR_subset'),
    make_option(c("-g", "--groups1"), action = "store", type = "character", help = "Classifications of cells (within meta_col1) to subset from dataset. \
    If multiple classifications are used to subest, must be provided as a comma separated list i.e. --groups celltype1,celltype2", default = NULL),
    #make_option(c("", "--groups2"), action = "store", type = "character", help = "Classifications of cells (within meta_col2) to subset from dataset.", default = NULL),
    #make_option(c("-i", "--invert1"), action = "store", type = "logical", help = "Boolean for whether to invert groups1 selection", default = FALSE),
    #make_option(c("", "--invert2"), action = "store", type = "logical", help = "Boolean for whether to invert groups2 selection", default = FALSE),
    make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    
    #plot_path = "./output/NF-downstream_analysis/4_ArchR_filter_clusters/plots/"
    #rds_path = "./output/NF-downstream_analysis/4_ArchR_filter_clusters/rds_files/"
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ss8/ArchR_integration/rds_files/"

    addArchRThreads(threads = 1)
    
    opt$invert1 = FALSE
    opt$invert2 = FALSE
    opt$groups1 = NULL
    opt$groups2 = NULL
    opt$meta_col1 = NULL
    opt$meta_col2 = NULL
  
    
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

#### TEMPORARY: just remove contamination
opt$meta_col1 = "scHelper_cell_type_old"
opt$groups1 = "NC"


# Set up options
if(!is.null(opt$groups1)){
  opt$groups1 = strsplit(opt$groups1, ',')[[1]]}
#if(!is.null(opt$groups2)){
#  opt$groups2 = strsplit(opt$groups2, ',')[[1]]}

if(is.na(opt$meta_col1)){
  stop("meta_col1 parameter must be provided. See script usage (--help)")}
if(is.null(opt$groups1)){
  stop("groups1 parameter must be provided. See script usage (--help)")}


############################## Function to subset ArchR project #######################################
subset_ArchR <- function(ArchR_object, meta_col, groups){
  idxPass <- which(ArchR_object@cellColData[,opt$meta_col] %in% opt$groups)
  cellsPass <- ArchR$cellNames[idxPass]
  ArchR_filtered <- ArchR[cellsPass, ]
  return(ArchR_filtered)
}

############################## Read in ArchR project #######################################
# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# see all available metadata cols which could be used to subset
print(colnames(ArchR@cellColData))

############################## Subset ArchR object #######################################
### will need to add extra functionality here 
ArchR_subset <- subset_ArchR(ArchR, meta_col = opt$meta_col1, groups = opt$groups1)

############################## Save new ArchR object #######################################
saveArchRProject(ArchRProj = ArchR_subset, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)