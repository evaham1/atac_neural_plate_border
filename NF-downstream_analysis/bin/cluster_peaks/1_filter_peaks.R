##### Input files: 1) summarised data .csv from SEACells computation (table of summarised peak counts across metacells)
##### and any 2) feature_metadata.csv from SEACells computation (table with all PeakSet information)
##### Normlises peak counts across seacells, then filters based on annotation and variability

# load libraries
library(getopt)
library(optparse)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
#library(hexbin)
library(gridExtra)
library(grid)
library(parallel)
library(data.table)
library(doSNOW)

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

    # combined SEACell outputs
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Combined_SEACell_outputs/"
    plot_path = "./output/NF-downstream_analysis/Processing/FullData/1_filter_peak_modules/plots/"
    rds_path = "./output/NF-downstream_analysis/Processing/FullData/1_filter_peak_modules/rds_files/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/csv_files/"
    
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Functions #######################################

# function to calculate variance and pick out top n features
calc_top_variable_features <- function(counts, n = 1000){
  # calculate variance
  variance <- apply(counts, 2, var)
  # rank and pick top features
  features <- names(variance[rank(-variance) < n+1])
  
  return(features)
}

############################## Read in SEACells summarised data + ArchR object #######################################

# Read in summarised peak counts across all seacells
SEACells_summarised <- as.matrix(fread(paste0(data_path, "Combined_summarised_by_metacells_counts.csv"), header = TRUE), rownames = 1)
#SEACells_summarised <- as.matrix(fread(paste0(data_path, "ss8_summarised_by_metacells_counts.csv"), header = TRUE), rownames = 1)

dim(SEACells_summarised)
SEACells_summarised[1:4, 1:4]
head(rownames(SEACells_summarised))

## Read in peak metadata
peak_metadata <- fread(paste0(data_path, "Feature_metadata.csv"), header = TRUE)
print(head(peak_metadata))

############################## 1) Normalise #######################################

## normalise each metacell by the total number of cut sites
normalised_counts <- t(apply(SEACells_summarised, 1, function(x) x/sum(x))) * 1000
normalised_counts[1:2, 1:2]
dim(normalised_counts)

## calculate new variance and plot 
variance <- apply(SEACells_summarised, 2, var)
print("Before normalising:")
print(summary(variance))

png(paste0(plot_path, "hist_variance_before_normalising.png"), width=60, height=40, units = 'cm', res = 200)
hist(variance, breaks = 1000)
graphics.off()

variance <- apply(normalised_counts, 2, var)
print("After normalising:")
print(summary(variance))

png(paste0(plot_path, "hist_variance_after_normalising.png"), width=60, height=40, units = 'cm', res = 200)
hist(variance, breaks = 1000)
graphics.off()


################ 2) Filter peaks by annotation #########################

# change this to directly use from feature set
peak_set <- getPeakSet(ArchR)
print(unique(peak_set$peakType))

included_peak_set <- peak_set[which(peak_set$peakType %in% c("Distal", "Intronic")), ]
print(paste0("Number of total peaks: ", length(peak_set$name)))
print(paste0("Number of peaks that are distal or intronic: ", length(included_peak_set$name)))

included_peaks <- included_peak_set$name

# filter summarised counts to only include these peaks
colnames(normalised_counts) <- gsub(':', '-', colnames(normalised_counts))
annot_filtered_matrix <- normalised_counts[, which(colnames(normalised_counts) %in% included_peaks)]
dim(annot_filtered_matrix)

############## 3) Filter peaks by variance #######################################

## pick peaks with top 10,000 variance (as test do top 100) determined by n param
top_features <- calc_top_variable_features(annot_filtered_matrix, n = n)
print(paste0("Number of top variable features: ", length(top_features)))

## filter matrix based on these features
variable_filtered_matrix <- annot_filtered_matrix[, top_features]
dim(variable_filtered_matrix)

variance <- apply(variable_filtered_matrix, 2, var)
print(summary(variance))

png(paste0(plot_path, "hist_variance_after_normalising_after_filtering.png"), width=60, height=40, units = 'cm', res = 200)
hist(variance, breaks = 1000)
graphics.off()


############################## Save filtered normalised summarised count data #######################################

write.csv(variable_filtered_matrix, paste0(rds_path, "Filtered_summarised_counts.csv"), row.names = TRUE)