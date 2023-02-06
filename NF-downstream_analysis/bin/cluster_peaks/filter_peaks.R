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
    label = "AnnData_summarised_by_metacells_peak_counts.csv"
    label = "summarised_counts_1000.csv"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


############################## Read in SEACells summarised data #######################################
data_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/SEACells/4_exported_SEACells_data/rds_files/"
SEACells_summarised <- read_csv(paste0(data_path, "AnnData_summarised_by_metacells_peak_counts.csv"))
dim(SEACells_summarised)

SEACells_summarised <- as.matrix(fread(paste0(data_path, "summarised_counts_1000.csv"), header = TRUE), rownames = 1)

hist(colSums(SEACells_summarised), breaks = 1000)
summary(colSums(SEACells_summarised))
dim(SEACells_summarised)

############################## 1) Normalise #######################################

## normalise each metacell by the total number of cut sites
normalised_counts <- t(apply(SEACells_summarised, 1, function(x) x/sum(x))) * 1000
normalised_counts[1:2, 1:2]
dim(normalised_counts)


## could add a minimum number of total cut sites here


# variance <- apply(SEACells_summarised, 2, var)
# hist(variance, breaks = 1000)
# summary(variance)


################ 2) Filter peaks by annotation #########################

peak_set <- getPeakSet(ArchR)
print(unique(peak_set$peakType))
included_peak_set <- peak_set[which(peak_set$peakType %in% c("Distal", "Intronic")), ]
print(length(included_peak_set$name))

# need to do this for the test matrix, for full matrix should return all peaks
included_peaks <- included_peak_set$name[included_peak_set$name %in% colnames(matrix)]
print(length(included_peaks))

# filter matrix to only include these peaks -> f1_matrix
f1_matrix <- matrix[, included_peaks]
print(paste0("Peak number before filtering: ", ncol(matrix)))
print(paste0("Peak number after filtering: ", ncol(f1_matrix)))
# 

############## 3) Filter peaks by variance #######################################

## pick peaks with top 10,000 variance
# as test do top 100

# function to calculate variance and pick out top n features
calc_top_variable_features <- function(counts, n = 1000){
  # calculate variance
  variance <- apply(counts, 2, var)
  # rank and pick top features
  features <- names(variance[rank(-variance) < n+1])
  
  return(features)
}

# variance <- apply(SEACells_summarised, 2, var)
# hist(variance, breaks = 1000)

top_features <- calc_top_variable_features(normalised_counts, n = 100)
length(top_features)
filtered_normalised_counts <- normalised_counts[, top_features]
dim(filtered_normalised_counts)
# variance <- apply(filtered_counts, 2, var)
# hist(variance, breaks = 1000)
# summary(variance)


############################## Subset summarised count data and save #######################################