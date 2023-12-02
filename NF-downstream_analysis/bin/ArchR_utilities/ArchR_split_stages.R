#!/usr/bin/env Rscript

print("ArchR split stages")

############################## Load libraries #######################################
library(getopt)
library(future)
library(parallel)
library(ArchR)
library(scHelper)

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
    
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Peak_call/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/FullData/Split_stages/rds_files/" 
    rds_path = "./output/NF-downstream_analysis/Processing/FullData/Split_stages/plots/"
    

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

set.seed(42)

############################## Read in ArchR project #######################################
# If files are not in rds_files subdirectory look in input dir
label <- unique(sub('_.*', '', list.files(data_path)))
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
split_ArchR <- ArchR_SplitbyStage(ArchR)
print("ArcR object split")

names(split_ArchR) <- unique(ArchR$stage)
print(names(split_ArchR))

# save RDS object for each stage/run
for(split in names(split_ArchR)){
  dir.create(paste0(rds_path, split, "_Save-ArchR"), recursive = T)
  saveArchRProject(ArchRProj = split_ArchR[[split]], outputDirectory = paste0(rds_path, split, "_Save-ArchR"), load = FALSE)
}