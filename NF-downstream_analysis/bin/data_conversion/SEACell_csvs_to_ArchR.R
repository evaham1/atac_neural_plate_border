#!/usr/bin/env Rscript

print("inputs: 1) original ArchR object, 2) csv mapping single cell IDs to SEACell IDs, 3) csv mapping SEACell IDs to integrate cell type label")
print("script uses csvs to map SEACell integrated labels on single cell data, adds metadata to ArchR and visualise on UMAP")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
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
    addArchRThreads(threads = 1) 
    
    ### Interactively will have to read from a few different places
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Peak_call/rds_files/" # for ArchR object
    data_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_ATAC_WF/2_SEACells_computation/exported_data/" # for single cell to seacell mapping
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Integrated_SEACells/" # integrated seacell mapping

    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    rds_path = "./rds_files/"
    plot_path = "./plots/"
    data_path = "./input/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(rds_path, recursive = T)
}

############################## Read data #######################################

# List all files in input
input_paths <- list.files(data_path)
print(input_paths)

# 1) Read in ArchR object
ArchR_path <- list.files(path = data_path, pattern = "*_Save-ArchR", full.names = TRUE)
print(paste0("ArchR path: ", ArchR_path))

ArchR <- loadArchRProject(path = ArchR_path, force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# 2) Read in ATAC single cell to SEACell mapping
SEACell_map_path <- list.files(path = paste0(data_path, "csv_files/"), pattern = "*_cell_metadata.csv", full.names = TRUE)
#SEACell_map_path <- list.files(path = paste0(data_path, ""), pattern = "*Cell_metadata.csv", full.names = TRUE) #interactive
print(paste0("SEACell map path: ", SEACell_map_path))

SEACell_map <- read.csv(SEACell_map_path)
print(head(SEACell_map))
print(dim(SEACell_map))

# 3) Read in ATAC SEACell to RNA SEACell mapping
Integration_map_path <- list.files(path = paste0(data_path, "rds_files/"), pattern = "*_mappings_cell_type.csv", full.names = TRUE)
print(paste0("Integration map path: ", Integration_map_path))

Integration_map <- read.csv(Integration_map_path)
print(head(Integration_map))
print(dim(Integration_map))

############################## Merge maps #######################################

length(unique(SEACell_map$SEACell))
length(unique(Integration_map$ATAC))
length(unique(Integration_map$RNA))

sum(unique(SEACell_map$SEACell) %in% unique(Integration_map$ATAC))

df_new <- merge(SEACell_map, Integration_map, by.x = "SEACell", by.y = "ATAC", all.x = TRUE)
head(df_new)

dim(df_new)
df_new$scHelper_cell_type

############ Add SEACell integrated scHelper_cell_type to ArchR object ###################


############################## Visualise on UMAPs #######################################

############################## Save ArchR #######################################


