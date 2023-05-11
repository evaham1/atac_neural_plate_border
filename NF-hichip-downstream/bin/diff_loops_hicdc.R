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

install.packages("DESeq2")

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-m", "--Dmin"), action = "store", type = "integer", help = "Minimum distance (included) to check for significant interactions"),
  make_option(c("-n", "--Dmax"), action = "store", type = "integer", help = "Maximum distance (included) to check for significant interactions"),
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
    
    plot_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_diff_interactions/plots/"
    rds_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_diff_interactions/rds_files/"
    
    ## Created folder with all HiCDC+ outputs to test interactively
    data_path = "./local_test_data/all_HiCDC_outputs/"
    
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

# read in all file paths
print("All samples:")
files <- list.files(data_path, full.names = TRUE)
print(files)

# read in loops for all WE samples
print("WE samples:")
WE_samples <-  grep("WE", files, invert = T, value = TRUE)
print(WE_samples)

# read in loops for all NF samples
print("NF samples:")
NF_samples <-  grep("NF", files, invert = T, value = TRUE)
print(NF_samples)

######################################   Create index file   ####################################################
# union of significant interactions, chr, startI, startJ

indexfile <- data.frame()
for (file in files) {
  # read in hicDC+ output
  output <- data.table::fread(file)
  # add unique interactions to indexfile
  indexfile <- unique(rbind(indexfile, output[,c('chrI','startI','startJ')]))
}

print(head(indexfile))
dim(indexfile)

#save index file
colnames(indexfile) <- c('chr','startI','startJ')
data.table::fwrite(indexfile,
                   paste0(rds_path,'/Indexfile.txt.gz'),
                   sep='\t',row.names=FALSE,quote=FALSE)

####################################   Run hicdcdiff   ###################################################
# Differential analysis using modified DESeq2 (see ?hicdcdiff)

hicdcdiff(input_paths = list(WE = WE_samples, NF = NF_samples),
          filter_file = paste0(rds_path,'/Indexfile.txt.gz'),
          output_path=paste0(rds_path,'/HicDCDiff_output/'),
          # Bins
          bin_type = 'Bins-uniform',
          binsize = 5000,
          # Interaction sizes
          Dmin = opt$Dmin,
          Dmax = opt$Dmax,
          # fitType in DESeq2::estimateDispersions
          fitType = 'mean',
          # What objects are made
          diagnostics=TRUE, # generates diagnostic plots of normalisation factors and MA plots
          DESeq.save = TRUE # saves DESEq objects for each chromosome
          )

#### LOOK FOR CONSISTENCY BETWEEN SAMPLES, PCA PLOT??? - venn diagram of all interactions across 6 samples

#### LOOK IF DIFF INTERACTIONS ARE FOUND IN AT LEAST 2 OF 3 REPLICATES