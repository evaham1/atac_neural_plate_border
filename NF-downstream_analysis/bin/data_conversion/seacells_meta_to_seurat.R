#!/usr/bin/env Rscript

print("importing summarised by metacells counts into seurat object")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(parallel)
library(Seurat)
library(dplyr)
library(tibble)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-i", "--input"), action = "store", type = "character", help = "Name of input file to convert"),
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
    data_path = "./local_test_data/SEACells_from_camp/SEACELLS_RNA_WF/exported_SEACells_data/rds_files/"
    opt$input = "AnnData_summarised_by_metacells_peak_counts.csv"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    label = "AnnData_summarised_by_metacells.h5ad"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

################### Functions ##########################

## function to aggregate matrix from seurat object and summarise by cell groupings
summarise_seurat_data <- function(seurat, data_slot = "counts", category = "SEACell"){
  
  # extract data into dataframe
  df <- GetAssayData(object = seurat, slot = data_slot)
  df <- as.data.frame(t(as.data.frame(df)))
  
  # convert cell ids to category ids
  category_ids <- select(seurat@meta.data, category)[,1]
  df <- df %>% mutate(category = category_ids)
  
  # aggregate df based on category
  df_summarised <- aggregate(. ~ category, df, sum)
  
  # format df so can be added back to seurat object
  df_summarised <- t(column_to_rownames(df_summarised, var = "category"))
  
  return(df_summarised)
}

#################### Read in data and add metadata to seurat object #########################

# Read in metdata
metacell_metadata <- read.csv(paste0(data_path, "AnnData_metacells_assigned_cell_metadata.csv"))
metacell_dictionary <- select(metacell_metadata, c("index", "SEACell"))

# Read in seurat object
seurat <- readRDS("~/local_test_data/seurat_from_camp/ss8_clustered_data.RDS")

# Reorder seacells metadata to match cell order in seurat object
metacell_dictionary <- metacell_dictionary[match(rownames(seurat_secell@meta.data), metacell_dictionary$index),]

# Add seacells metadata to seurat object
seurat <- AddMetaData(seurat, metacell_dictionary$SEACell, col.name = "SEACell")

# Save seurat object
saveRDS(seurat, paste0(rds_path, "seurat.RDS"), compress = FALSE)

#################### Create new seurat object with only metacells #########################

# Filter seurat object to only include SEACells
seacells_seurat <- seurat[,colnames(seurat) %in% unique(metacell_dictionary$SEACell)]

#################### Add up counts across metacells #########################

###### RNA slot: 3 assays: counts (raw), data (normalised), scale.data -> only add up raw 'counts'
DefaultAssay(object = seacells_seurat) <- "RNA"
DefaultAssay(object = seacells_seurat)

summarised_RNA_counts <- summarise_seurat_data(seurat = seurat, data_slot = "counts", category = "SEACell")
seacells_seurat <- SetAssayData(object = seacells_seurat, slot = "counts", new.data = summarised_RNA_counts, assay = "RNA")

###### Integrated slot: 2 assays: data, scale.data -> only add up 'data'
DefaultAssay(object = seacells_seurat) <- "integrated"
DefaultAssay(object = seacells_seurat)

summarised_integrated_data <- summarise_seurat_data(seurat = seurat, data_slot = "data", category = "SEACell")
seacells_seurat <- SetAssayData(object = seacells_seurat, slot = "data", new.data = summarised_integrated_data, assay = "integrated")


#####################################################################################
############################ Re-process added up counts #############################
#####################################################################################




