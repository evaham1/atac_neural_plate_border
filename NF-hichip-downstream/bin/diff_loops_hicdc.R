#!/usr/bin/env Rscript

### script to use HiCDCPlus to find differential significant loops from HiChip data
print("script to use HiCDCPlus to find differential significant loops from HiChip data")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(parallel)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(HiCDCPlus)

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
    
    plot_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_diff_interactions/plots/"
    rds_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_diff_interactions/rds_files/"
    data_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output/rds_files/" 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################################################################################################
################################## Find Differential loops #################################################
############################################################################################################

######################################   Read in loops   ####################################################

print("reading in loops...")

## read in loops for all WE samples
# interactions <- data.table::fread(paste0(data_path, "HiCDC_output_filtered.txt"))
# head(interactions)
# dim(interactions)

## read in loops for all NF samples
# interactions <- data.table::fread(paste0(data_path, "HiCDC_output_filtered.txt"))
# head(interactions)
# dim(interactions)

## turn into index file
# #save index file---union of significants at 50kb
# colnames(indexfile)<-c('chr','startI','startJ')
# data.table::fwrite(indexfile,
#             paste0(outdir,'/GSE131651_analysis_indices.txt.gz'),
#             sep='\t',row.names=FALSE,quote=FALSE)

####################################   Run hicdcdiff   ###################################################

# #Differential analysis using modified DESeq2 (see ?hicdcdiff)

# hicdcdiff(input_paths=list(NSD2=c(paste0(outdir,'/GSE131651_NSD2_LOW_arima_example.txt.gz'),
#                  paste0(outdir,'/GSE131651_NSD2_HIGH_arima_example.txt.gz')),
# TKO=c(paste0(outdir,'/GSE131651_TKOCTCF_new_example.txt.gz'),
# paste0(outdir,'/GSE131651_NTKOCTCF_new_example.txt.gz'))),
# filter_file=paste0(outdir,'/GSE131651_analysis_indices.txt.gz'),
# output_path=paste0(outdir,'/diff_analysis_example/'),
# fitType = 'mean',
# binsize=50000,
# diagnostics=TRUE)

#### LOOK FOR CONSISTENCY BETWEEN SAMPLES, PCA PLOT??? - venn diagram of all interactions across 6 samples

#### LOOK IF DIFF INTERACTIONS ARE FOUND IN AT LEAST 2 OF 3 REPLICATES