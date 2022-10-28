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
    
    # ss8 object
    data_path <- "./output/NF-downstream_analysis/Processing/ss8/PEAK_CALLING/peak_call/rds_files/ss8_Save-ArchR/"
    ArchR <- loadArchRProject(path = data_path, force = FALSE, showLogo = TRUE)

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

##############################################################################
############################### FUNCTIONS ####################################

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
    print(paste0("Chunk ", i))
    chunk_coords <- coords[i, "start"]:coords[i, "end"]
    feature_pairs_chunk <- feature_pairs[chunk_coords, ]
    chunk_cor_mat <- chunk_cor(mat, feature_pairs_chunk)
    # populate the correct coords of the correlation_matrix
    corr_mat[rownames(chunk_cor_mat), colnames(chunk_cor_mat)] <- chunk_cor_mat
  }
  
  return(corr_mat)
}

#################################################################################
############################## Read in data #####################################

# extract peak matrix - takes a long time
peak_data <- getMatrixFromProject(ArchR, useMatrix = "PeakMatrix", threads = 1)
matrix <- t(assays(peak_data)[[1]])
colnames(matrix) <- rowData(peak_data)$name

### as extracting matrix takes a long time could write and read matrix
# # save peak matrix
# writeMM(peak_matrix, paste0(rds_path, "peak_matrix.txt"))
# write(colnames(peak_matrix), paste0(rds_path, "colnames.txt"))
# write(rownames(peak_matrix), paste0(rds_path, "rownames.txt"))

# # read in data
# matrix <- readMM("peak_matrix.txt")
# colnames(matrix) <- scan("colnames.txt", character(), quote = "")
# rownames(matrix) <- scan("rownames.txt", character(), quote = "")


# ######### Generate Test matrix
# # create test matrix and save
# test_matrix <- matrix[1:1000, 1:10000]
# writeMM(test_matrix, "test_matrix.txt")
# write(colnames(test_matrix), "test_colnames.txt")
# write(rownames(test_matrix), "test_rownames.txt")


###########################################################################
################ Filter peaks based on annotation #########################

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

##############################################################################################
################ Filter peaks based on number of cells it is open in #########################

### Chose final number as 100,000 peaks

## Plot distribution of number of cells per peak
png(paste0(plot_path, 'Unfiltered_hist_cells_per_peak.png'), height = 20, width = 25, units = 'cm', res = 400)
hist(colSums(f1_matrix > 0), breaks = 1000)
graphics.off()

cells_per_peak <- summary(colSums(f1_matrix > 0))
table <- data.frame(Stats = names(cells_per_peak), Value = as.vector(cells_per_peak))
png(paste0(plot_path, 'Unfiltered_cells_per_peak_summary_table.png'), height = 30, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Number of cells peak is open in", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(table, rows=NULL, theme = ttheme_minimal()))
graphics.off()

## Filter peaks so take only ones open in many cells - top 100,000, change to 1,000 for test matrix
top_peaks <- sort(colSums(f1_matrix > 0), decreasing = TRUE)[1:100000]
print(length(names(top_peaks)))

f2_matrix <- f1_matrix[, names(top_peaks)]
print(paste0("Peak number before filtering: ", ncol(f1_matrix)))
print(paste0("Peak number after filtering: ", ncol(f2_matrix)))

## Redo plots to see where these filtered peaks lie in the distribution
png(paste0(plot_path, 'Filtered_hist_cells_per_peak.png'), height = 20, width = 25, units = 'cm', res = 400)
hist(colSums(f2_matrix > 0), breaks = 1000)
graphics.off()

cells_per_peak_filtered <- summary(colSums(f2_matrix > 0))
table <- data.frame(Stats = names(cells_per_peak_filtered), Value = as.vector(cells_per_peak_filtered))
png(paste0(plot_path, 'Filtered_cells_per_peak_summary_table.png'), height = 30, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Number of cells peak is open in", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(table, rows=NULL, theme = ttheme_minimal()))
graphics.off()


#################################################################################
############################## Cluster peaks #####################################

correlation_matrix <- fast_cor(mat = f2_matrix[1:100, 1:100], chunk_size = 1000)
save(correlation_matrix, file = paste0(rds_path, "correlation_matrix.RData"))



######################
#################################################################################
############################## Different clustering approach on just ss8  #####################################

peak_data <- getMatrixFromProject(ArchR, useMatrix = "PeakMatrix", threads = 1)
matrix <- t(assays(peak_data)[[1]])
colnames(matrix) <- rowData(peak_data)$name

#### extract LSI dimensions fo cell (30 for each cell)
dims <- getReducedDims(ArchRProj = ArchR)

## build a nearest neighbour graph from these dims
install.packages("FNN")
library(FNN)


iterative_NN_clustering <- function(dims, k, seed = 123){
  set.seed(seed)
  
  clusters <- data.frame()
  
  while(nrow(dims) > k){
    random <- sample(1:nrow(dims), 1)
    cell_IDs <- as.vector(get.knnx(data = dims, query = dims[random, , drop = FALSE], k = k)$nn.index)
    cell_names <- rownames(dims)[cell_IDs]
    
    clusters <- rbind(cell_names, clusters)
    
    dims <- dims[-cell_IDs, ]
  }
  
  return(clusters)
}

clusters <- iterative_NN_clustering(dims, k = 100)


# normalise so each peak has same range
# sum fragments per pseudocell, divide each peak count by this sum and then times by 1000

rownames(matrix)

# remove cells that are not clustered
filtered_matrix <- matrix[unlist(clusters), ]

# make cluster df long
rownames_to_column(clusters, var = "cluster_id")
pivot_longer(clusters, cols = colnames(clusters))

matrix %>% mutate()



