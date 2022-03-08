#!/usr/bin/env Rscript

print("ArchR split stages")

############################## Load libraries #######################################
library(getopt)
library(future)
library(tidyverse)
library(grid)
library(gridExtra)
library(clustree)
library(GenomeInfoDb)
library(ggplot2)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(parallel)
library(ArchR)
library(GenomicFeatures)

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
    
    setwd("~/NF-downstream_analysis")
    ncores = 8
    
    plot_path = "../output/NF-downstream_analysis/5_ArchR_clustering_postfiltering/plots/"
    rds_path = "../output/NF-downstream_analysis/5_ArchR_clustering_postfiltering/rds_files/"
    data_path = "../output/NF-downstream_analysis/6_ArchR_stage_split/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 155* 1024^3)
    addArchRThreads(threads = ncores) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

addArchRThreads(threads = 16) 

############################## Function to split samples ################################

## need to make this generic (not just split by stage) need to overcome passing variable to $ issue
split_ArchR_by_stage <- function(ArchR_project){
  split_names <- unique(ArchR_project$stage)
  split_ArchR <- c()
  for (i in split_names){
    idxSample <- BiocGenerics::which(ArchR_project$stage %in% i)
    cellsSample <- ArchR_project$cellNames[idxSample]
    split_ArchR <- c(split_ArchR, ArchR_project[cellsSample, ])
  }
  return(split_ArchR)
}

############################## Read in ArchR project #####################################
ArchR <- loadArchRProject(path = paste0(data_path, "./rds_files/Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

############################## Split ArchR project #######################################
split_ArchR <- split_ArchR_by_stage(ArchR)
names(split_ArchR) <- unique(ArchR$stage)

# save RDS object for each stage/run
for(split in names(split_ArchR)){
  dir.create(paste0(rds_path, split, "_Save-ArchR"), recursive = T)
  saveArchRProject(ArchRProj = split_ArchR[[split]], outputDirectory = paste0(rds_path, split, "_Save-ArchR"), load = FALSE)
}