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
library(gridExtra)
library(grid)
library(parallel)
library(data.table)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-n", "--nPeaks"), action = "store", type = "integer", help = "Number of top variable peaks to include", default = 100),
  make_option(c("", "--verbose"), action = "store", type = "logical", help = "Verbose", default = TRUE)
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
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/0_combining_outputs/csv_files/"
    plot_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/plots/"
    rds_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/"
    
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
SEACells_summarised <- fread(paste0(data_path, "Combined_summarised_by_metacells_counts.csv"), header = TRUE)
#SEACells_summarised <- as.matrix(fread(paste0(data_path, "ss8_summarised_by_metacells_counts.csv"), header = TRUE), rownames = 1)
print("Data read in!")

# Check input
print(dim(SEACells_summarised))
print("Preview of input df:")
print(SEACells_summarised[1:4, 1:4])

# Extract SEACell IDs from first column
SEACells_IDs <- as.vector(SEACells_summarised[,2])
print(head(SEACells_IDs))

# Clean up df
SEACells_summarised <- SEACells_summarised[,-1:-2]
print("Preview of input df after cleanup:")
print(SEACells_summarised[1:4, 1:4])

# Turn into numeric matrix for downstream processing
SEACells_summarised_numeric <- matrix(as.numeric(as.character(SEACells_summarised)), ncol = ncol(SEACells_summarised))

# Add SEACell IDs as rownames and peak IDs as colnames
rownames(SEACells_summarised_numeric) <- SEACells_IDs
colnames(SEACells_summarised_numeric) <- colnames(SEACells_summarised)

# Check resulting matrix
print(dim(SEACells_summarised_numeric))
print("Preview of summarised count df:")
print(SEACells_summarised_numeric[1:4, 1:4])

## Read in peak metadata
peak_metadata <- fread(paste0(data_path, "Feature_metadata.csv"), header = TRUE)
print("Preview of peakset df:")
print(head(peak_metadata))

############################## 1) Normalise #######################################
## normalise each metacell by the total number of cut sites

print("Normalising counts...")

normalised_counts <- t(apply(SEACells_summarised_numeric, 1, function(x) x/sum(x))) * 1000

print("Preview of normalised counts df:")
print(normalised_counts[1:2, 1:2])
dim(normalised_counts)

## calculate new variance and plot 
variance_before <- apply(SEACells_summarised_numeric, 2, var)
print("Before normalising:")
print(summary(variance_before))

png(paste0(plot_path, "hist_variance_before_normalising.png"), width=60, height=40, units = 'cm', res = 200)
hist(variance_before, breaks = 1000, xlim = c(0, 100))
graphics.off()

variance_after <- apply(normalised_counts, 2, var)
print("After normalising:")
print(summary(variance_after))

png(paste0(plot_path, "hist_variance_after_normalising.png"), width=60, height=40, units = 'cm', res = 200)
hist(variance_after, breaks = 1000, xlim = c(0, 0.0005))
graphics.off()

print("Counts normalised!")

################ 2) Filter peaks by annotation #########################
## only include intergenic and intronic peaks

print("Filtering by peak annotation...")

included_peak_set <- peak_metadata[which(peak_metadata$peakType %in% c("Distal", "Intronic")), ]
print(paste0("Number of total peaks: ", length(peak_metadata$name)))
print(paste0("Number of peaks that are distal or intronic: ", length(included_peak_set$name)))

included_peaks <- included_peak_set$name

# filter summarised counts to only include these peaks
colnames(normalised_counts) <- gsub(':', '-', colnames(normalised_counts))
annot_filtered_matrix <- normalised_counts[, which(colnames(normalised_counts) %in% included_peaks)]
dim(annot_filtered_matrix)

print("Peaks filtered!")

############## 3) Filter peaks by variance #######################################
## pick peaks with most variance across all cells

print("Filtering by variance...")

## pick peaks with top 10,000 variance (as test do top 100) determined by n param
top_features <- calc_top_variable_features(annot_filtered_matrix, n = opt$nPeaks)
print(paste0("Number of top variable features: ", length(top_features)))

## filter matrix based on these features
variable_filtered_matrix <- annot_filtered_matrix[, top_features]

# check new matrix
print("Preview of filtered normalised counts df:")
print(variable_filtered_matrix[1:2, 1:2])
dim(variable_filtered_matrix)

# Check new variance
variance <- apply(variable_filtered_matrix, 2, var)
print(summary(variance))

png(paste0(plot_path, "hist_variance_after_normalising_after_filtering.png"), width=60, height=40, units = 'cm', res = 200)
hist(variance, breaks = 1000)
graphics.off()

print("Peaks filtered!")

############################## Save filtered normalised summarised count data #######################################

write.csv(variable_filtered_matrix, paste0(rds_path, "Filtered_summarised_counts.csv"), row.names = TRUE, col.names = TRUE)