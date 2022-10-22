#!/usr/bin/env Rscript

print("calculates peak modules")
# calculates peak modules by iterative correlation of chunks

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-v", "--verbose"), action = "store", type = "logical", help = "Verbose", default = FALSE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8

    # transfer_labels object
    data_path <- "./output/NF-downstream_analysis/Downstream_processing/transfer_labels/peak_call/rds_files/FullData_Save-ArchR/"
    ArchR <- loadArchRProject(path = data_path, force = FALSE, showLogo = TRUE)

    # read in test data
    matrix <- readMM("test_matrix.txt")
    colnames(matrix) <- scan("test_colnames.txt", character(), quote = "")
    rownames(matrix) <- scan("test_rownames.txt", character(), quote = "")

    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/rds_files/"
    rds_path = "./rds_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

###########     Functions
##########################################################
## function to create empty final correlation matrix
gen_corr_matrix <- function ( mat ){
  return(matrix(nrow = ncol(mat), ncol = ncol(mat),
                dimnames = list(colnames(mat), colnames(mat))))
}

## function to calculate feature pairs
calc_feature_pairs <- function ( mat ){
  return( t(combn(1:ncol(mat),2)) )
}

## function to generate subset coords - returns as table with cols 'start' and 'end'
calc_chunk_coords <- function( feature_pairs, k){
  start <- seq(1, nrow(feature_pairs), by = k)
  end <- c(start[-1] - 1, nrow(feature_pairs))
  return(data.frame(start = start, end = end))
}

## calculate correlation for the i-th chunk
chunk_cor <- function (mat, subset_feature_pairs) {
  # chunk the feature matrix to just get the columns corresponding to these peaks
  matrix_chunk <- mat[, unique(as.vector(subset_feature_pairs))]
  # run the correlation
  correlation_chunk <- cor(as.matrix(matrix_chunk), method = "spearman")
  # output the correlation chunk
  return(correlation_chunk)
}

## now we loop through all chunks and combine the result
fast_cor <- function (mat, chunk_size) {
  
  corr_mat <- gen_corr_matrix(mat)
  feature_pairs <- calc_feature_pairs(mat)
  coords <- calc_chunk_coords(feature_pairs, k = chunk_size)
  
  for (i in 1:nrow(coords)){
    chunk_coords <- coords[i, "start"]:coords[i, "end"]
    feature_pairs_chunk <- feature_pairs[chunk_coords, ]
    chunk_cor_mat <- chunk_cor(mat, feature_pairs_chunk)
    # populate the correct coords of the correlation_matrix
    corr_mat[rownames(chunk_cor_mat), colnames(chunk_cor_mat)] <- chunk_cor_mat
  }
  
  return(corr_mat)
}



##############################################################################

# extract peak matrix - takes a long time
peak_data <- getMatrixFromProject(ArchR, useMatrix = "PeakMatrix", threads = 1)
peak_matrix <- t(assays(peak_data)[[1]])
colnames(peak_matrix) <- rowData(peak_data)$name


### as extracting matrix takes a long time could write and read matrix
# # save peak matrix
# writeMM(peak_matrix, paste0(rds_path, "peak_matrix.txt"))
# write(colnames(peak_matrix), paste0(rds_path, "colnames.txt"))
# write(rownames(peak_matrix), paste0(rds_path, "rownames.txt"))

# # read in data
# matrix <- readMM("peak_matrix.txt")
# colnames(matrix) <- scan("colnames.txt", character(), quote = "")
# rownames(matrix) <- scan("rownames.txt", character(), quote = "")


# ######### Test matrix
# # create test matrix and save
# test_matrix <- matrix[1:1000, 1:10000]
# writeMM(test_matrix, "test_matrix.txt")
# write(colnames(test_matrix), "test_colnames.txt")
# write(rownames(test_matrix), "test_rownames.txt")


#########################################################################
### Filter peaks, try to cut down 260,000 -> 100,000
hist()
hist(colSums(matrix > 0), breaks = 1000)
summary(colSums(matrix > 0))

top_peaks <- sort(colSums(matrix > 0), decreasing = TRUE)[1:1000]
names(top_peaks)

filtered_matrix <- matrix[, names(top_peaks)]
dim(filtered_matrix)


############################################################################
### Cluster peaks

fast_cor(mat = matrix[1:100, 1:100], chunk_size = 5)





###################   OLD
#############   Test speed of correlating matrices of different sizes

# # cor with different sized matrices
# mini_matrix <- peak_matrix[, 1:1000]
# mini_matrix_dense <- as.matrix(mini_matrix)
# system.time(cor(mini_matrix_dense, method = "spearman"))
# #user  system elapsed 
# #0.690   0.000   0.692 

# mini_matrix <- peak_matrix[, 1:2000]
# mini_matrix_dense <- as.matrix(mini_matrix)
# system.time(cor(mini_matrix_dense, method = "spearman"))
# #user  system elapsed 
# #25.113   0.025  25.176 

# mini_matrix <- peak_matrix[, 1:5000]
# mini_matrix_dense <- as.matrix(mini_matrix)
# system.time(cor(mini_matrix_dense, method = "spearman"))
# #user  system elapsed 
# #150.106   0.891 151.164

# mini_matrix <- peak_matrix[, 1:10000]
# mini_matrix_dense <- as.matrix(mini_matrix)
# system.time(cor(mini_matrix_dense, method = "spearman"))
# #user  system elapsed 
# #580.238   0.931 582.086 

# # plot relative timings using standard cor function:
# plot_data <- data.frame(Number_of_peaks = c(1000, 2000, 5000, 10000),
#                         Seconds_to_compute = c(0.692, 25.176, 151.164, 582.086))
# ggplot(plot_data, aes(x = Number_of_peaks, y = Seconds_to_compute)) + geom_point() + geom_line()


# # see how long it takes to rank matrix, then should be able to just feed into pearsons test
# system.time(apply(mini_matrix_dense, 2, rank))
# #    user  system elapsed 
# #     1.156   0.000   1.157 
# system.time())
# frank(mini_matrix_dense, x = colnames(mini_matrix_dense))
# #    user  system elapsed 
# #   2.432   0.002   0.643 


# ranked_mini_matrix <- apply(mini_matrix_dense, 2, rank)
# ranked_mini_matrix[1:5, 1:5]
# dim(ranked_mini_matrix)
# max(ranked_mini_matrix[,4])