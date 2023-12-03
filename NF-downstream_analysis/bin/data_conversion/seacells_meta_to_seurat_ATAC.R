#!/usr/bin/env Rscript

print("This script to takes the ArchR exported gene scores matrix and cell metadata seacell assignments")
print("It 1) summarises the gene scores matrix by seacell assignments")
print("and 2) uses this matrix to make a seacells ATAC gene score seurat object")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(parallel)
library(Seurat)
library(dplyr)
library(tibble)
library(scHelper)
library(ggplot2)
library(future)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-m", "--metadata_file_name"), action = "store", type = "character", help = "Name of csv file which assigns cell ids to metacell ids", default = "exported_data/Cell_metadata.csv"),
  make_option(c("-g", "--matrix_file_name"), action = "store", type = "character", help = "Name of csv file with the gene score matrix", default = "exported_data/gene_scores.csv"),
  make_option(c("", "--verbose"), action = "store", type = "logical", help = "Verbose", default = TRUE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    data_path = "./local_test_data/test_inputs/test_input_seacells_meta_to_seurat_ATAC//"
    rds_path = "./local_test_data/convert_seacells_ATAC_to_seurat/rds_files/"
    plot_path = "./local_test_data/convert_seacells_ATAC_to_seurat/plots/"
    
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

set.seed(42)

# ## function to aggregate matrix from seurat object and summarise by cell groupings
# SEACells_SummariseSeuratData <- function(seurat, data_slot = "counts", category = "SEACell"){
  
#   # extract data into dataframe
#   df <- GetAssayData(object = seurat, slot = data_slot)
#   df <- as.data.frame(t(as.data.frame(df)))
  
#   # convert cell ids to category ids
#   category_ids <- select(seurat@meta.data, category)[,1]
#   df <- df %>% mutate(category = category_ids)
  
#   # aggregate df based on category
#   df_summarised <- aggregate(. ~ category, df, sum)
  
#   # format df so can be added back to seurat object
#   df_summarised <- t(column_to_rownames(df_summarised, var = "category"))
  
#   return(df_summarised)
# }


#######################################################################################
######################    Read in data   #########################
#######################################################################################

# Read in metadata (use metadata_file_name to find correct object)
metacell_metadata <- read.csv(paste0(data_path, opt$metadata_file_name))
metacell_dictionary <- select(metacell_metadata, c("index", "SEACell"))
metacell_dictionary$index <- gsub('#', '-', metacell_dictionary$index)

print("Metacell IDs read in")
print(head(metacell_dictionary))
print(paste0("Number of cells in dictionary: ", dim(metacell_dictionary))[1])

# Read in exported gene scores matrix
gene_score_matrix <- read.csv(paste0(data_path, opt$matrix_file_name))
colnames(gene_score_matrix) <- gsub('\\.', '-', colnames(gene_score_matrix))

print("Gene score matrix read in")
print(gene_score_matrix[1:3, 1:3])
print(dim(gene_score_matrix))

if (dim(metacell_dictionary)[1] == dim(gene_score_matrix)[2]-1){
  print("Number of cells in dictionary and matrix match!")
} else {
  print("ERROR: Number of cells in dictionary and matrix DO NOT match!!")
}

#####################################################################################
######################    Create summarised seurat object   #########################
#####################################################################################

#################### Add up gene score counts across metacells #########################

summarised_gene_score_counts <- SEACells_SummariseGeneScoreMatrix(matrix = gene_score_matrix, dictionary = metacell_dictionary)

dim(summarised_gene_score_counts) # 7547  202
sum(is.na(summarised_gene_score_counts)) # 0 NA values

print("gene score counts summarised by metacells")
print(summarised_gene_score_counts[1:3, 1:3])
rownames(summarised_gene_score_counts) <- gene_score_matrix$X

#################### Create metdata: just stage #########################

stage <- substr(metacell_dictionary$index[1], 8, 10)

#################### Create new seurat object #########################

# Create object using summarised RNA counts and newly created cell metadata
seacells_seurat <- CreateSeuratObject(
  project = "seacells_seurat",
  counts = summarised_gene_score_counts,
  assay = "RNA",
  names.field = 2,
  names.delim = "-",
  meta.data = NULL
)

# Add stage metadata
seacells_seurat <- AddMetaData(seacells_seurat, stage, col.name = "stage")
print(seacells_seurat@meta.data$stage)

# Print object
print(seacells_seurat)

print("created summarised seurat object!")

#################### Save new seurat object #########################

## save seacells seurat object
saveRDS(seacells_seurat, paste0(rds_path, "seacells_seurat.RDS"), compress = FALSE)
