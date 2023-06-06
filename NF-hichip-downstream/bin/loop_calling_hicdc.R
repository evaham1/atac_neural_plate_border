#!/usr/bin/env Rscript

### script to use HiCDCPlus to call significant loop from HiChip data
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
library(data.table)
setDTthreads(threads = 1) # to try to overcome data.table error in gi_list_write

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-b", "--binsize"), action = "store", type = "integer", help = "Size of uniform bins to make"),
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
    chrs = c("chr21", "chr22")
    opt$binsize = 5000
    opt$Dmin = 10000
    opt$Dmax = 1000000
    
    plot_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output/plots/"
    rds_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output/rds_files/"
    rds_path = "./output/NF-hichip-downstream/NF_HiChip_r3/HicDCPlus_output/rds_files/"
    
    data_path = "./output/NF-hichip-downstream/NF_HiChip_r1/edit_validpairs/" # for valid pairs NF replicate 1
    data_path = "./output/NF-hichip-downstream/NF_HiChip_r2/edit_validpairs/" # for valid pairs NF replicate 2
    data_path = "./output/NF-hichip-downstream/NF_HiChip_r3/edit_validpairs/" # for valid pairs NF replicate 3
    
    data_path = "./output/NF-hichip-downstream/bins/" # for bins
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    ncores = opt$cores
    chrs = NULL

    #### TEMP TO DEBUG:
    ncores = 8
    # chrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
    #             "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20")
    # chrs = c("chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30",
    #             "chr31", "chr32", "chr33", "chrZ", "chrW")
    # chrs = c("chr29")
    chrs = c("chr30", "chr31", "chr32", "chr33", "chrZ", "chrW")
    opt$binsize = 5000
    opt$Dmin = 10000
    opt$Dmax = 1000000
    ####
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

sessionInfo()

############################################################################################################
############################################# Call loops ###################################################
############################################################################################################

######################################   Read in bins   ####################################################

print("reading in bins...")

# Check bins object
bins <- data.table::fread(paste0(data_path, "rds_files/bins_bintolen.txt.gz"))
head(bins, 20)

# Create gi_list from this bintolen object ie make all possible interactions
gi_list<-generate_bintolen_gi_list(
  bintolen_path=paste0(data_path, "rds_files/bins_bintolen.txt.gz"),
  gen = "Ggallus", gen_ver = "galGal6")

# Check gi_list
gi_list_validate(gi_list) # passes without errors
#head(gi_list)
# GInteractions object with 1457234 interactions and 1 metadata column:
#   seqnames1           ranges1     seqnames2           ranges2 |         D
# <Rle>         <IRanges>         <Rle>         <IRanges> | <integer>
#   [1]     chr13            0-5000 ---     chr13            0-5000 |         0
# [2]     chr13            0-5000 ---     chr13        5000-10000 |      5000
# [3]     chr13            0-5000 ---     chr13       10000-15000 |     10000
# [4]     chr13            0-5000 ---     chr13       15000-20000 |     15000
# [5]     chr13            0-5000 ---     chr13       20000-25000 |     20000

print("bins read in!")

####################################   Read in ValidPairs   ###################################################

print("reading in HiChip counts...")

# Extract paths - should look like: 'WE_HiChip_r2_edited.allValidPairs'
valid_pair_path <- list.files(path = data_path, pattern = "*.allValidPairs", full.names = TRUE)
print(paste0("Valid pair path: ", valid_pair_path))

# Extract sample name
sample_name <- gsub(pattern = "_edited.allValidPairs", replacement = "", x = basename(valid_pair_path))
print(paste0("Detected sample name: ", sample_name))

# Add Validpairs to gi_list
gi_list_with_valid_pairs <- add_hicpro_allvalidpairs_counts(gi_list, allvalidpairs_path = valid_pair_path)

# Check file now
gi_list_validate(gi_list_with_valid_pairs)
# head(gi_list_with_valid_pairs)
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

print("HiChip counts read in!")

####################################   Run HiCDC+   ##################################################

print("running HiCDC+...")

# Expand features for modeling - 1D = genomic features for each bin, 2D = features belonging to an interaction.
expanded_gi_list_with_valid_pairs <- expand_1D_features(gi_list_with_valid_pairs)

# Run HiC-DC+
set.seed(1010) #HiC-DC downsamples rows for modeling
# finds significant interactions in HiC-DC readable matrix and expresses statistical significance of counts with p-val, q-val, FDR corrected pval (mu)
# ncore defaults to parallel::detectCores()-1

print("Number of cores detected:")
print(parallel::detectCores())

# expanded_gi_list_with_valid_pairs_HiCDC <- HiCDCPlus_parallel(expanded_gi_list_with_valid_pairs,
#                                                               # Covariates:
#                                                               covariates = NULL, # covariates to be considered in addition to D, defaults to all covariates besides D, counts, mu, sdev, pvalue, qvalue
#                                                               # Modelling params:
#                                                               distance_type = "spline", # distance covariate form, 'spline' or 'log'
#                                                               model_distribution = "nb", # 'nb' uses negative binomial model, 'nb_vardisp' uses nb model with distance specific dispersion inferred from data, 'nb_hurdle' uses legacy HiC-DC model
#                                                               df = 6, # degrees of freedom for the genomic distance spline function if distance_type='spline', defaults to 6 which corresponds to cubic splines
#                                                               ssize = 0.01, # distance stratified sampling size. increase recommended if model fails to converge, defaults to 0.01
#                                                               splineknotting = "uniform", # spline knotting strategy, either 'uniform' or 'count-based' (ie more closed spaces where counts are more dense)
#                                                               # Cores:
#                                                               ncore = 1, # number of cores to parallelize
#                                                               # Types of bins:
#                                                               binned = TRUE, # TRUE if uniform binning, FALSE if restriction enzyme fragment cutsites
#                                                               # Resulting loops params:
#                                                               Dmin = opt$Dmin, # minimum distance (included) to check for significant interactions, defaults to 0
#                                                               Dmax = opt$Dmax, # maximum distance (included) to check for significant interactions, 1.5e6 is recommended for HiChip data in manual
#                                                               # Which chroms (if interactive will be chr21 and ch22, if not will be all)
#                                                               chrs = chrs
#                                                               )

expanded_gi_list_with_valid_pairs_HiCDC <- HiCDCPlus(expanded_gi_list_with_valid_pairs,
                                                     # Covariates:
                                                     covariates = NULL, # covariates to be considered in addition to D, defaults to all covariates besides D, counts, mu, sdev, pvalue, qvalue
                                                     # Modelling params:
                                                     distance_type = "spline", # distance covariate form, 'spline' or 'log'
                                                     model_distribution = "nb", # 'nb' uses negative binomial model, 'nb_vardisp' uses nb model with distance specific dispersion inferred from data, 'nb_hurdle' uses legacy HiC-DC model
                                                     df = 6, # degrees of freedom for the genomic distance spline function if distance_type='spline', defaults to 6 which corresponds to cubic splines
                                                     ssize = 0.01, # distance stratified sampling size. increase recommended if model fails to converge, defaults to 0.01
                                                     splineknotting = "uniform", # spline knotting strategy, either 'uniform' or 'count-based' (ie more closed spaces where counts are more dense)
                                                     # Types of bins:
                                                     binned = TRUE, # TRUE if uniform binning, FALSE if restriction enzyme fragment cutsites
                                                     # Resulting loops params:
                                                     Dmin = opt$Dmin, # minimum distance (included) to check for significant interactions, defaults to 0
                                                     Dmax = opt$Dmax, # maximum distance (included) to check for significant interactions, 1.5e6 is recommended for HiChip data in manual
                                                     # Which chroms (if interactive will be chr21 and ch22, if not will be all)
                                                     chrs = chrs
                                                     )

# Check one chromosome
# head(expanded_gi_list_with_valid_pairs_HiCDC)
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

print("HiCDC+ run!")

# ####################################   Write outputs   ##################################################

# print("saving outputs...")

# # Save results as is for differential interaction testing
# gi_list_write(expanded_gi_list_with_valid_pairs_HiCDC, fname = paste0(rds_path, sample_name, '_HiCDC_output.txt.gz'))

# # Filter results by adjusted p value - done manually as can then remove NAs myself
# # There will be NAs for every interactions that falls outside DMin or DMax, so need to remove these from gi_list altogether
# filtered_list <- list()
# for (i in 1:length(expanded_gi_list_with_valid_pairs_HiCDC)){
#   name <- names(expanded_gi_list_with_valid_pairs_HiCDC)[i]
#   temp <- expanded_gi_list_with_valid_pairs_HiCDC[i][[1]]
#   temp_filtered <- temp[!is.na(temp$qvalue), ]
#   temp_filtered_1 <- temp_filtered[temp_filtered$qvalue < 0.05, ]
  
#   filtered_list[[name]] <- temp_filtered_1
# }

# gi_list_write(filtered_list,
#               fname = paste0(rds_path, sample_name, '_HiCDC_output_filtered.txt'),
#               rows = "all")

# print("outputs saved!")

# ####################################   Visualisations   ##################################################

# print("Visualisations...")

# # Read in bins bed object
# bins <- data.table::fread(paste0(data_path, "bed_files/bins.bed"))
# colnames(bins) <- c("chr", "start", "end", "bin_ID")
# head(bins)

# # Read in significant interactions object
# interactions <- data.table::fread(paste0(rds_path, sample_name, '_HiCDC_output_filtered.txt'))
# head(interactions)
# dim(interactions)

# #### Check bin sizes

# # calculate bin sizes
# bins <- as.data.frame(bins)
# bins <- bins %>% mutate(bin_width = as.numeric(bins$end) - as.numeric(bins$start))
# print("Bins:")
# head(bins)

# # how many different bin sizes are there?
# png(paste0(plot_path, "freq_of_different_bin_sizes.png"), width=40, height=20, units = 'cm', res = 200)
# plot(table(bins$bin_width))
# graphics.off()
# # most bins are 4999, but there are one of unique sizes

# # where are the bins that are not specified bin size?
# print("Bins which are not specified bin size:")
# print(bins[bins$bin_width != opt$binsize-1, ])
# # can see there is one bin which is not the specified bin size in each chromosome

# #### Check interaction distances

# # check the range of interactions in the sig interactions table
# print(paste0("Minimum interaction size: ", min(interactions$D)))
# print(paste0("Maximum interaction size: ", max(interactions$D)))

# # check which interactions are not divisible by bin size
# print("Interactions which distance is not divisible by bin size:")
# check.integer <- function(x) {x == round(x)}
# weird_length_interactions <- interactions[!check.integer(interactions$D / opt$binsize), ]
# print(head(weird_length_interactions))
# # these interactions are all interacting with the weirdly sized last bin

# # remove the weird sized interactions and extract counts vs distance
# normal_sized_interactions <- interactions[check.integer(interactions$D / opt$binsize), ]

# #### Check counts

# # plot a histogram of distribution of counts
# png(paste0(plot_path, "hist_of_counts.png"), width=40, height=20, units = 'cm', res = 200)
# hist(normal_sized_interactions$counts, breaks = 100)
# graphics.off()

# png(paste0(plot_path, "hist_of_log_counts.png"), width=40, height=20, units = 'cm', res = 200)
# hist(log(normal_sized_interactions$counts), breaks = 100)
# graphics.off()

# # plot summary stats of range of counts
# print("Counts summary:")
# print(summary(normal_sized_interactions$counts))

# #### Check interaction distances Vs counts

# # extract interactions distances and counts from normal sized bins
# data <- data.frame(distance = normal_sized_interactions$D,
#                    counts = normal_sized_interactions$counts)

# # plot means of counts per distance
# mean_data <- aggregate(counts ~ distance, data, mean)

# png(paste0(plot_path, "mean_distance_vs_counts.png"), width=40, height=20, units = 'cm', res = 200)
# ggplot(mean_data, aes(x=distance, y=counts)) +
#   geom_line() +
#   geom_point()
# graphics.off()

# png(paste0(plot_path, "log_mean_distance_vs_counts.png"), width=40, height=20, units = 'cm', res = 200)
# ggplot(mean_data, aes(x=log10(distance), y=log10(counts))) +
#   geom_line() +
#   geom_point()
# graphics.off()

# # plot medians of counts per distance
# median_data <- aggregate(counts ~ distance, data, median)

# png(paste0(plot_path, "median_distance_vs_counts.png"), width=40, height=20, units = 'cm', res = 200)
# ggplot(median_data, aes(x=distance, y=counts)) +
#   geom_line() +
#   geom_point()
# graphics.off()

# png(paste0(plot_path, "log_median_distance_vs_counts.png"), width=40, height=20, units = 'cm', res = 200)
# ggplot(median_data, aes(x=log10(distance), y=log10(counts))) +
#   geom_line() +
#   geom_point()
# graphics.off()

# #### Check significance values

# # Plot distribution of q values
# png(paste0(plot_path, "hist_q_val.png"), width=40, height=20, units = 'cm', res = 200)
# hist(unique(interactions$qvalue), breaks = 100)
# graphics.off()

# # Plot distribution of p values
# png(paste0(plot_path, "hist_p_val.png"), width=40, height=20, units = 'cm', res = 200)
# hist(unique(interactions$pvalue), breaks = 100)
# graphics.off()


# #### Check fragment sizes
# HiC_data <- fread(valid_pair_path)
# frag_sizes <- HiC_data$V8
# print(summary(frag_sizes))

# png(paste0(plot_path, "hist_frag_sizes.png"), width=40, height=20, units = 'cm', res = 200)
# hist(frag_sizes, breaks = 100)
# graphics.off()

# png(paste0(plot_path, "hist_log_frag_sizes.png"), width=40, height=20, units = 'cm', res = 200)
# hist(log(frag_sizes), breaks = 100)
# graphics.off()

# png(paste0(plot_path, "boxplot_log_frag_sizes.png"), width=40, height=20, units = 'cm', res = 200)
# BiocGenerics::boxplot(log(frag_sizes))
# graphics.off()