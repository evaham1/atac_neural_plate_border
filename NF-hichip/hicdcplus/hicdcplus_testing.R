#!/usr/bin/env Rscript

##### Move to docker image - run on archr docker ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("HiCDCPlus")
######

print("script to use HiCDCPlus to call significant loop from HiChip data")

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
    
    plot_path = "./output/NF-hichip/HicDC/plots/"
    rds_path = "./output/NF-hichip/HicDC/rds_files/"
    data_path = "./output/NF-hichip/NF-hic/"
    
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


## Finding Significant Interactions from HiChIP

## Create genomic features
# finds all restriction enzyme cutsites of a given genome and genome version and computes GC content, mappability (if a relevant .bigWig file is provided) and effective fragment length for uniform bin or across specified multiples of restriction enzyme cutsites of given pattern(s).

construct_features(output_path = paste0(rds_path,"test"),
                   gen = "Ggallus", gen_ver = "galGal6", # BSgenome.Ggallus.UCSC.galGal6, same as used for ArchR preprocessing
                   sig = c("GATC"), # this is the cut site of restriction enzyme Mboi which was used to make HiChip data
                   bin_type = "Bins-uniform", binsize = 5000
                   )
bintolen <- data.table::fread(paste0(rds_path,"test_bintolen.txt.gz"))
head(bintolen,20)

# Create gi_list from this bintolen object
gi_list<-generate_bintolen_gi_list(
  bintolen_path=paste0(rds_path,"/test_bintolen.txt.gz"),
  gen = "Ggallus", gen_ver = "galGal6")
head(gi_list)

#add .hic counts
hicfile_path = ""
valid_pair_path = paste0(data_path, "NF_HiChip_r1.allValidPairs")

gi_list <- add_hicpro_allvalidpairs_counts(gi_list, allvalidpairs_path = valid_pair_path)

test <- read.table(valid_pair_path,
  sep="\t", header=FALSE)

#expand features for modeling
gi_list<-expand_1D_features(gi_list)

#run HiC-DC+ on 2 cores
set.seed(1010) #HiC-DC downsamples rows for modeling
gi_list<-HiCDCPlus_parallel(gi_list,ncore=2)
head(gi_list)



