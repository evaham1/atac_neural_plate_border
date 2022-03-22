#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
library(Seurat)
library(future)
library(tidyverse)
library(grid)
library(gridExtra)
library(clustree)
library(scHelper)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    #plot_path = "./output/NF-downstream_analysis_stacas//merged/plots/"
    #rds_path = "./output/NF-downstream_analysis_stacas/seurat/merged/rds_files/"
    #data_path = "./output/NF-scRNAseq_alignment/cellranger/count/filtered_feature_bc_matrix"
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 16gb
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

# Make dataframe with stage and replicate info extracted from path
print(list.files(data_path, recursive = FALSE, full.names = TRUE))

input <- list.files(data_path, recursive = FALSE, full.names = TRUE)
input <- data.frame(sample = sub('.*/', '', input), path = input)
print(input)

# Init list of seurat objects then merge
seurat_list <- apply(input, 1, function(x) readRDS(file = x[["path"]]))
names(seurat_list) <- substr(input$sample, 1, 3)
seurat_all <- merge(x = seurat_list[[1]], y=seurat_list[-1], add.cell.ids = names(seurat_list), project = "chick.10x")

# Save RDS output
saveRDS(seurat_all, paste0(rds_path, "merged_seurat.RDS"), compress = FALSE)