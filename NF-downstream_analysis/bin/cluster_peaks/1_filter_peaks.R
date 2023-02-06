##### Takes seacells summarised data .csv and filters peaks to reduce number for clustering

# load libraries
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(hexbin)
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
    addArchRThreads(threads = 1)

    # output from SEACells - summarised by metacells -> should change to output of seacells purity script
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/SEACells/4_exported_SEACells_data/rds_files/"
    label = "summarised_counts_1000.csv"
    # for transfer labels archr object -> should change to output of seacells purity script which will include extra cell metadata
    data_path = "./output/NF-downstream_analysis/Processing/TransferLabels/3_peak_call/rds_files/"
    n = 100
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    label = "AnnData_summarised_by_metacells_peak_counts.csv"
    n = 20000
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
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
SEACells_summarised <- as.matrix(fread(paste0(data_path, label), header = TRUE), rownames = 1)
dim(SEACells_summarised)

ArchR <- loadArchRProject(path = paste0(data_path, "TransferLabels_Save-ArchR"), force = FALSE, showLogo = TRUE)

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