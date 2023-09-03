#!/usr/bin/env Rscript

print("Run scMEGA motif annotation and chromvar")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(Seurat)
library(Signac)
library(scMEGA)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Ggallus.UCSC.galGal6)
library(SummarizedExperiment)
library(igraph)
library(ggraph)
library(BiocParallel)

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
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/rds_files/"
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

############################## Read in seurat objects #######################################

## read in paired seurat object
obj.pair <- readRDS(paste0(data_path, "TransferLabel_paired_object.RDS"))

############################## Add motif information #######################################

# download motif database
motifList <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))

# add motif information to ATAC data
obj.motifs <- AddMotifs(
  object = obj.pair,
  genome = BSgenome.Ggallus.UCSC.galGal6,
  pfm = motifList,
  assay = "ATAC"
)

############################## Run chromvar #######################################

## NB if I try to rename motifList by TF names like I do when running chromvar in ArchR this fails

BiocParallel::register(SerialParam())

# run chromvar
obj_chromvar <- RunChromVAR(
  object = obj.motifs,
  genome = BSgenome.Ggallus.UCSC.galGal6,
  assay = 'ATAC'
)

# save chromvar object
saveRDS(obj_chromvar, paste0(rds_path, "paired_object_chromvar.RDS"), compress = FALSE)