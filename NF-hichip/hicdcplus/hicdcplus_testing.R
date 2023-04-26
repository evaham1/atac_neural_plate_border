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
    data_path = "./output/NF-hichip/HicDC/"
    
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
                   bin_type = "Bins-uniform", 
                   binsize = 5000 # resolution = 5kb
                   )
bintolen <- data.table::fread(paste0(rds_path,"test_bintolen.txt.gz"))
head(bintolen,20)

# Create gi_list from this bintolen object
gi_list<-generate_bintolen_gi_list(
  bintolen_path=paste0(rds_path,"/test_bintolen.txt.gz"),
  gen = "Ggallus", gen_ver = "galGal6")

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

#add .hic counts - need to edit valid pairs output to change '1' -> 'chr1' in column 2 and column 5
valid_pair_path = paste0(data_path, "NF_HiChip_r1_edited_v6.allValidPairs")
#valid_pair_file <- data.table::fread(valid_pair_path, sep = "\t", header = FALSE)
#head(valid_pair_file)

gi_list_with_valid_pairs <- add_hicpro_allvalidpairs_counts(gi_list, allvalidpairs_path = valid_pair_path)
gi_list_validate(gi_list_with_valid_pairs)
head(gi_list_with_valid_pairs)
length(gi_list_with_valid_pairs) # 42 - each element in list is a chromosome
gi_list_with_valid_pairs[[1]]
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

#expand features for modeling - adds 2D features in metadata handle? what does that mean?
expanded_gi_list_with_valid_pairs <- expand_1D_features(gi_list_with_valid_pairs)
expanded_gi_list_with_valid_pairs[[1]]
mcols(expanded_gi_list_with_valid_pairs[[1]])

#run HiC-DC+ on 2 cores
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
                                                              chrs = c('chr21', 'chr22')
                                                              )
head(expanded_gi_list_with_valid_pairs_HiCDC[[13]]) # chromosome that wasnt ran using HiCDC
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
#   regions: 2780 ranges and 2 metadata columns
# seqinfo: 1 sequence from an unspecified genome; no seqlengths
head(expanded_gi_list_with_valid_pairs_HiCDC[[14]]) # chromosome that WAS ran using HiCDC
# GInteractions object with 6 interactions and 8 metadata columns:
#   seqnames1   ranges1     seqnames2     ranges2 |         D    counts        gc       len        mu      sdev
# <Rle> <IRanges>         <Rle>   <IRanges> | <integer> <numeric> <numeric> <numeric> <numeric> <numeric>
#   [1]     chr21    0-5000 ---     chr21      0-5000 |         0         3  0.414755  -1.92886  10.35434   6.53781
# [2]     chr21    0-5000 ---     chr21  5000-10000 |      5000         2  0.500550  -3.06837   6.10078   4.16470
# [3]     chr21    0-5000 ---     chr21 10000-15000 |     10000         0 -0.670526  -5.92123   1.57612   1.52531
# [4]     chr21    0-5000 ---     chr21 15000-20000 |     15000         0 -0.600917  -5.92123   1.44458   1.44049
# [5]     chr21    0-5000 ---     chr21 25000-30000 |     25000         0 -0.874250  -3.21762   3.28106   2.55602
# [6]     chr21    0-5000 ---     chr21 30000-35000 |     30000         0  0.213983  -4.09739   2.46746   2.07527
# pvalue    qvalue
# <numeric> <numeric>
#   [1]  0.930362         1
# [2]  0.900999         1
# [3]  1.000000         1
# [4]  1.000000         1
# [5]  1.000000         1
# [6]  1.000000         1
# -------
#   regions: 1369 ranges and 2 metadata columns
# seqinfo: 1 sequence from an unspecified genome; no seqlengths


chr21_output <- expanded_gi_list_with_valid_pairs_HiCDC[[14]]
head(chr21_output)
hist(unique(chr21_output[,8]$qvalue))

filtered_ch21_output <- chr21_output[chr21_output$qvalue < 0.05]
nrow(filtered_ch21_output)

#write normalized counts (observed/expected) to a .hic file
hicdc2hic(expanded_gi_list_with_valid_pairs_HiCDC,
          hicfile=paste0(rds_path,'/Test_sample_combined_result_chr21_chr22.hic'),
          mode='normcounts',
          gen_ver='galGal6')

#write results to a text file
gi_list_write(expanded_gi_list_with_valid_pairs_HiCDC,
              fname=paste0(rds_path,'/Test_sample_combined_result.txt.gz'))

