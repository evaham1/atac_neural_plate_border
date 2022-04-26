#!/usr/bin/env Rscript

print("ArchR split stages")

############################## Load libraries #######################################
library(getopt)
library(future)
library(parallel)
library(ArchR)

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
    data_path = "../output/NF-downstream_analysis/6_ArchR_stage_split/rds_files/"

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

############################## Read in ArchR project #######################################
# If files are not in rds_files subdirectory look in input dir
label <- sub('_.*', '', list.files(data_path))
print(label)

if (length(label) == 0){
  data_path = "./input/"
  label <- sub('_.*', '', list.files(data_path))
  print(label)
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
} else {
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
}

############################## Split ArchR project #######################################
split_ArchR <- split_ArchR_by_stage(ArchR)
print("ArcR object split")

names(split_ArchR) <- unique(ArchR$stage)
print(names(split_ArchR))

# save RDS object for each stage/run
for(split in names(split_ArchR)){
  dir.create(paste0(rds_path, split, "_Save-ArchR"), recursive = T)
  saveArchRProject(ArchRProj = split_ArchR[[split]], outputDirectory = paste0(rds_path, split, "_Save-ArchR"), load = FALSE)
}