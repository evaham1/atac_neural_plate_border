#!/usr/bin/env Rscript

### script to use HiCDCPlus to call significant loop from HiChip data
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
    data_path = "./output/NF-hichip/HicDC/"
    
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

############################################################################################################
############################################# Call loops ###################################################
############################################################################################################

######################################   Read in bins   ####################################################

# Read in bintolen object with the generated bins
bintolen <- data.table::fread(paste0(rds_path, "rds_files/bins_bintolen.txt.gz"))
head(bintolen,20)

# Create gi_list from this bintolen object
gi_list<-generate_bintolen_gi_list(
  bintolen_path=paste0(rds_path, "rds_files/bins_bintolen.txt.gz"),
  gen = "Ggallus", gen_ver = "galGal6")

# Check gi_list
gi_list_validate(gi_list) # passes without errors
head(gi_list)
# GInteractions object with 1457234 interactions and 1 metadata column:
#   seqnames1           ranges1     seqnames2           ranges2 |         D
# <Rle>         <IRanges>         <Rle>         <IRanges> | <integer>
#   [1]     chr13            0-5000 ---     chr13            0-5000 |         0
# [2]     chr13            0-5000 ---     chr13        5000-10000 |      5000
# [3]     chr13            0-5000 ---     chr13       10000-15000 |     10000
# [4]     chr13            0-5000 ---     chr13       15000-20000 |     15000
# [5]     chr13            0-5000 ---     chr13       20000-25000 |     20000

####################################   Read in ValidPairs   ###################################################

# Add Validpairs
valid_pair_path = paste0(data_path, "edited_ValidPairs.txt")
gi_list_with_valid_pairs <- add_hicpro_allvalidpairs_counts(gi_list, allvalidpairs_path = valid_pair_path)

# Check file now
gi_list_validate(gi_list_with_valid_pairs)
head(gi_list_with_valid_pairs)
# length(gi_list_with_valid_pairs) # 42 - each element in list is a chromosome
# gi_list_with_valid_pairs[[1]]
# mcols(gi_list[[1]])
# An object of class "GInteractions"
# Slot "anchor1":
#   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# [41] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

# Slot "anchor2":
#   [1]   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
# [21]  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40

# Slot "regions":
#   GRanges object with 39522 ranges and 2 metadata columns:
#   seqnames              ranges strand |        gc       len
# <Rle>           <IRanges>  <Rle> | <numeric> <numeric>
#   [1]     chr1              0-5000      * |  0.547500      1771
# [2]     chr1          5000-10000      * |  0.562083      2279
# [3]     chr1         10000-15000      * |  0.430000      1000
# [4]     chr1         15000-20000      * |  0.424536      3604
# [5]     chr1         20000-25000      * |  0.583636      4801
# ...      ...                 ...    ... .       ...       ...
# [39518]     chr1 197585000-197590000      * |  0.520000      3961
# [39519]     chr1 197590000-197595000      * |  0.544423      4006
# [39520]     chr1 197595000-197600000      * |  0.556905      3978
# [39521]     chr1 197600000-197605000      * |  0.525588      3376
# [39522]     chr1 197605000-197608386      * |  0.563571      3135
# -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths
# 
# Slot "NAMES":
#   NULL
# 
# Slot "elementMetadata":
#   DataFrame with 15768122 rows and 2 columns
# D    counts
# <integer> <numeric>
#   1                0         0
# 2             5000         0
# 3            10000         0
# 4            15000         0
# 5            20000         0
# ...            ...       ...
# 15768118      5000         0
# 15768119      9193         0
# 15768120         0         5
# 15768121      4193         6
# 15768122         0         1
# 
# Slot "metadata":
#   list()

####################################   Run HiCDC+   ##################################################

# Expand features for modeling - adds 2D features in metadata handle? what does that mean?
expanded_gi_list_with_valid_pairs <- expand_1D_features(gi_list_with_valid_pairs)

# Run HiC-DC+
set.seed(1010) #HiC-DC downsamples rows for modeling
# finds significant interactions in HiC-DC readable matrix and expresses statistical significance of counts with p-val, q-val, FDR corrected pval (mu)
# ncore defaults to parallel::detectCores()-1
expanded_gi_list_with_valid_pairs_HiCDC <- HiCDCPlus_parallel(expanded_gi_list_with_valid_pairs,
                                                              covariates = NULL,
                                                              distance_type = "spline",
                                                              model_distribution = "nb",
                                                              binned = TRUE,
                                                              df = 6,
                                                              Dmin = 0,
                                                              Dmax = 1.5e6, # recommended for HiChip data in manual
                                                              ssize = 0.01,
                                                              splineknotting = "uniform",
                                                              chrs = c("chr21","chr22")
                                                              )

# Check one chromosome
head(expanded_gi_list_with_valid_pairs_HiCDC[[1]])
# GInteractions object with 6 interactions and 4 metadata columns:
#   seqnames1   ranges1     seqnames2     ranges2 |         D    counts        gc       len
# <Rle> <IRanges>         <Rle>   <IRanges> | <integer> <numeric> <numeric> <numeric>
#   [1]     chr20    0-5000 ---     chr20      0-5000 |         0         0 -2.530504  -2.16879
# [2]     chr20    0-5000 ---     chr20  5000-10000 |      5000         0 -2.561440  -3.23152
# [3]     chr20    0-5000 ---     chr20 10000-15000 |     10000         0 -2.601628  -4.18493
# [4]     chr20    0-5000 ---     chr20 15000-20000 |     15000         0 -2.415055  -5.41408
# [5]     chr20    0-5000 ---     chr20 20000-25000 |     20000         0 -1.074367  -2.81809
# [6]     chr20    0-5000 ---     chr20 25000-30000 |     25000         0  0.131059  -7.14647
# -------

####################################   Visualisations   ##################################################

## Plot counts over interaction distances - takes too long
#plot(chr1_output[,1]$D, chr1_output[,2]$counts)

## Plot distribution of significance 
# hist(unique(chr1_output[,8]$qvalue), breaks = 100)
# 
# ### Filter to only include significant interactions
# filtered_ch1_output <- chr1_output[chr1_output$qvalue < 0.05]
# nrow(filtered_ch21_output)
# 
# hist(unique(filtered_ch21_output[,2]$counts), breaks = 100)
# 
# filtered_ch21_output[1,]
# 
# filtered_ch21_output_2 <- filtered_ch21_output[filtered_ch21_output$counts > 50]
# filtered_ch21_output_2

####################################   Write outputs   ##################################################

#write normalized counts (observed/expected) to a .hic file
hicdc2hic(expanded_gi_list_with_valid_pairs_HiCDC,
          hicfile=paste0(rds_path,'/HiCDC_output.hic'),
          mode='normcounts',
          gen_ver='galGal6')

#write results to a text file
gi_list_write(expanded_gi_list_with_valid_pairs_HiCDC,
              fname=paste0(rds_path,'/HiCDC_output.txt.gz'))
