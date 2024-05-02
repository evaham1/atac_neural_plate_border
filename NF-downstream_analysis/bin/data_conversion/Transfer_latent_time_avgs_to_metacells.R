#!/usr/bin/env Rscript

print("Transfer latent time and lineage probabilities from single cell ATAC object to metacell ATAC metadata by averaging them per metacell")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(Seurat)
library(SummarizedExperiment)

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
    
    # data_path = "./output/NF-downstream_analysis/Processing/FullData/Single_cell_integration/" # ATAC
    # data_path = "./output/NF-downstream_analysis/rna_objects/" # RNA
    # rds_path = "./output/NF-downstream_analysis/Processing/FullData/Transfer_latent_time/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/"
    rds_path = "./rds_files/"
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

############################## Read in data #######################################

# read in ATAC single cell data
ArchR <- loadArchRProject(path = paste0(data_path, "rds_files/", "TransferLabels_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# read in ATAC metacell metadata csv
metadata <- read.csv(paste0(data_path, "Combined_SEACell_integrated_metadata.csv"), row.names = 'ATAC')

# Add stage to metadata using SEACell IDs
substrRight <- function(x, n){
  sapply(x, function(xx)
    substr(xx, (nchar(xx)-n+1), nchar(xx))
  )
}
metadata <- metadata %>% mutate(stage = substrRight(rownames(metadata), 3))
metadata <- metadata[,-1]

# make metacell IDs a column so can merge with more metadata
metadata$SEACell_ID <- row.names(metadata)
row.names(metadata) <- NULL

# Check metadata
print(head(metadata))

print("Data read in!")


############################## Add latent time and lineage probs to metacells #######################################

# extract data
single_cell_map <- as.data.frame(ArchR@cellColData) %>% dplyr::select(c("rna_latent_time", "rna_lineage_neural_probability", 
                                                         "rna_lineage_NC_probability", "rna_lineage_placodal_probability",
                                                         "in_rna_lineage", "SEACell_ID"))

# summarise by averaging the latent time and lineage probabilities, and calculate the percentage of single cells in metacell that are in a lineage
metacells_latent_time <- single_cell_map %>%
  dplyr::group_by(SEACell_ID) %>%
  dplyr::summarise(across(c("rna_latent_time", "rna_lineage_neural_probability", "rna_lineage_NC_probability", "rna_lineage_placodal_probability"), mean), 
            across(c("in_rna_lineage"), ~mean(. == TRUE) * 100, .names = "percentage_{.col}"))
        
head(metacells_latent_time)

# add that to exisiting metacell metadata
new_metadata <- merge(metadata, metacells_latent_time, by = "SEACell_ID")

# Write out unaltered SEACell metadata
write.csv(new_metadata, paste0(rds_path, "Combined_SEACell_integrated_metadata_latent_time.csv"), col.names = TRUE)