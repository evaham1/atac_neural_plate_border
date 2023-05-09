#!/usr/bin/env Rscript

### script to investigate the loops made by HiCDCPlus
print("script to investigate the loops made by HiCDCPlus")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(parallel)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(HiCDCPlus)
library(gridExtra)
library(grid)
library(VennDiagram)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

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
    
    plot_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output_investigating/plots/"
    rds_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output_investigating/rds_files/"

    data_path = "./output/NF-hichip-downstream/NF_HiChip_r1/HicDCPlus_output/rds_files/" # for HiCDCPlus output
    data_path = "./output/NF-hichip-downstream/bins/genes_intersect/" # for filtering by genes
    data_path = "./output/NF-hichip-downstream/bins/peaks_intersect/" # for filtering by peaks
    data_path = "./output/NF-hichip-downstream/bins/rds_files/" # for bins
    
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


######################################   Read in data   ####################################################

print("reading in data...")

# read in HiCDCPlus output
interactions <- data.table::fread(paste0(data_path, "rds_files/HiCDC_output_filtered.txt"))
head(interactions)
dim(interactions)

# read in genes intersected with bins
genes_bins <- data.table::fread(paste0(data_path, "tag_chroms_bins_intersected.bed"))
colnames(genes_bins) <- c("bin_chr", "bin_start", "bin_end", "bin_ID", "gene_chr", "gene_start", "gene_end", "gene_ID", "gene_name", "strand")
head(genes_bins)
dim(genes_bins)

# read in peaks intersected with bins
peaks_bins <- data.table::fread(paste0(data_path, "FullData_PeakSet_bins_intersected.bed"))
colnames(peaks_bins) <- c("bin_chr", "bin_start", "bin_end", "bin_ID", "peak_chr", "peak_start", "peak_end", "peak_ID", "dunno", "strand")
head(peaks_bins)
dim(peaks_bins)

# read in all bins
bins <- data.table::fread(paste0(data_path, "bins.bed"))
head(bins)

print("data read in!")

##############  Looking at anchors and how many times they appear    ###########################################

print("Investigating interaction anchor frequencies...")

# how many significant interactions are there
how_many_interactions <- nrow(interactions)
print(paste0("Number of sig interactions: ", how_many_interactions))

# extract all anchors in the significant interactions
anchors_I <- interactions[, 1:3]
colnames(anchors_I) <- c("chr", "start", "end")
anchors_J <- interactions[, 4:6]
colnames(anchors_J) <- c("chr", "start", "end")
all_anchors <- rbind(anchors_I, anchors_J)
all_anchors <- all_anchors %>% mutate(bin_ID = paste0(chr, "-", start+1, "-", end))
print("All anchors:")
print(head(all_anchors))

# frequency at which these anchors appear
print("Frequency at which anchors occur in significant interactions: ")
print(table(table(all_anchors$bin_ID)))

png(paste0(plot_path, "freq_of_anchors_in_interactions.png"), width=60, height=40, units = 'cm', res = 200)
plot(table(table(all_anchors$bin_ID)))
graphics.off()

# pull out the anchors that appear more than 100 times - can vary this threshold
highly_interacting_anchors <- table(all_anchors$bin_ID)[table(all_anchors$bin_ID) > 100]
print("highly interacting anchors: ")
print(highly_interacting_anchors)

##############  Looking at counts of interactions over distance    ###########################################

#### CHECKING BIN SIZES

# # check whether all bins are the same width
# bins <- as.data.frame(bins)
# first.word <- function(X){
#   unlist(strsplit(X, "-"))[1]
# }
# second.word <- function(X){
#   unlist(strsplit(X, "-"))[2]
# }
# third.word <- function(X){
#   unlist(strsplit(X, "-"))[3]
# }
# bins_ch22 <- bins_ch22 %>% 
#   mutate(chr = sapply(bins, first.word)) %>%
#   mutate(start = sapply(bins, second.word)) %>%
#   mutate(end = sapply(bins, third.word))
# bins_ch22 <- bins_ch22 %>%
#   mutate(bin_width = as.numeric(bins_ch22$end) - as.numeric(bins_ch22$start))

# # there are two bin lengths: 4999, 4461 (4461 is the last bin on the chromosome)
# unique(bins_ch22$bin_width)
# bins_ch22[bins_ch22$bin_width == 4461, ]
# tail(bins_ch22)

# #### CHECKING INTERACTION DISTANCES

# # check the range of interaction distances possible 
# summary(temp$D) # range of interaction distances is: min 0, max 2,000,000 ie 2 million ie 2 MB

# # check which interaction distances are divisible by 5000 (bin width)
# check.integer <- function(x) {
#   x == round(x)
# }
# length(temp$D) - sum(check.integer(temp$D / 5000))
# # 375 interactions do not have a distance divisible by 5000

# # check which bins these interaction sizes are connected to
# temp_weird_interactions <- temp[!check.integer(temp$D / 5000), ]
# temp_weird_interactions
# # these interactions are all interacting with the mis-sized last bin

# # remove these 'weird' interactions
# temp_normal_interactions <- temp[check.integer(temp$D / 5000), ]
# temp_normal_interactions

# # plot distance by counts
# data <- data.frame(distance = temp_normal_interactions$D,
#                    counts = temp_normal_interactions$counts)

# # using means per distance
# mean_data <- aggregate(counts ~ distance, data, mean)
# ggplot(mean_data, aes(x=log10(distance), y=log10(counts))) +
#   geom_line() +
#   geom_point()
# ggplot(mean_data, aes(x=distance, y=counts)) +
#   geom_line() +
#   geom_point()
# # looks like around 1 to 1.5mb the profile become jagged. 
# # currently the max loop length is set as 1.5mb, could maybe reduce it to 1mb

# # using medians per distance
# median_data <- aggregate(counts ~ distance, data, median)
# ggplot(median_data, aes(x=distance, y=counts)) +
#   geom_line() +
#   geom_point()
# ggplot(median_data, aes(x=log10(distance), y=log10(counts))) +
#   geom_line() +
#   geom_point()

# ## Plot distribution of significance 
# # hist(unique(temp$qvalue), breaks = 100)
# ########################


##############  Overlap of anchors with called peaks and genes    ###########################################

## how many bins in each type of data
bins_numbers <- data.frame(
  Total_Bins = length(unique(bins$bins)),
  Unique_Anchors = length(unique(all_anchors$bin_ID)),
  Gene_Bins = length(unique(genes_bins$bin_ID)),
  Peak_Bins = length(unique(peaks_bins$bin_ID))
)

png(paste0(plot_path, 'bin_counts_table.png'), height = 10, width = 20, units = 'cm', res = 400)
grid.arrange(tableGrob(bins_numbers, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# sanity check that all anchors are in the bins
if (sum(unique(all_anchors$bin_ID) %in% bins$bins)){
  print("All anchors are in the bins!")
}else{stop("PROBLEM: NOT ALL ANCHORS ARE IN THE BINS!")}

# Venn diagram of bins between interactions, peaks and genes
venn.diagram(
  x = list(unique(all_anchors$bin_ID), unique(genes_bins$bin_ID), unique(peaks_bins$bin_ID)),
  category.names = c("Interactions" , "Genes" , "Peaks"),
  filename = paste0(plot_path, 'bins_venn_diagram.png'),
  output=TRUE, disable.logging = TRUE,
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)


##############  Overlap of anchors with peaks and genes of interest    ###########################################