#!/usr/bin/env Rscript

### script to make bins for HiCDCPlus analysis + to intersect with peak and gene beds
print("script to make bins for HiCDCPlus analysis + to intersect with peak and gene beds")

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
    
    plot_path = "./output/NF-hichip-downstream/1_bins/plots/"
    rds_path = "./output/NF-hichip-downstream/1_bins/rds_files/"
    data_path = ""
    
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

############################## Create bins #######################################

# finds all restriction enzyme cutsites of a given genome and genome version and computes GC content, mappability (if a relevant .bigWig file is provided) and effective fragment length for uniform bin or across specified multiples of restriction enzyme cutsites of given pattern(s).

# construct features
construct_features(output_path = paste0(rds_path, "bins"),
                   gen = "Ggallus", gen_ver = "galGal6", # BSgenome.Ggallus.UCSC.galGal6, same as used for ArchR preprocessing
                   sig = c("GATC"), # this is the cut site of restriction enzyme Mboi which was used to make HiChip data
                   bin_type = "Bins-uniform", 
                   binsize = 5000 # resolution = 5kb
                   )

# read in bintolen object to check
bintolen <- data.table::fread(paste0(rds_path,"bins_bintolen.txt.gz"))
head(bintolen,20)

# Write out bins into simplified bed file so they can be intersected with peaks and genes
split <- function(str, col) { 
  return (strsplit(as.character(str), '-')[[1]][col])
}

df <- data.frame(chrom = unlist(lapply(bintolen$bins, FUN = split, col=1)),
                 chromStart = unlist(lapply(bintolen$bins, FUN = split, col=2)),
                 chromEnd = unlist(lapply(bintolen$bins, FUN = split, col=3)),
                 name = bintolen$bins
)
head(df)

write.table(df, file = paste0(rds_path, "bins.bed"), sep="\t", row.names=F, col.names=F, quote = F)