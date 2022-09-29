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

# Full data
data_path <- "./output/NF-downstream_analysis/ArchR_peak_exploration/transfer_labels/peak_call/rds_files/FullData_Save-ArchR/"
ArchR <- loadArchRProject(path = data_path, force = FALSE, showLogo = TRUE)
getAvailableMatrices(ArchR)

# just ss8
data_path <- "./output/NF-downstream_analysis/ArchR_preprocessing/FILTERING/ss8/postfiltering/peak_calling/rds_files/ss8_Save-ArchR/"
ArchR <- loadArchRProject(path = data_path, force = FALSE, showLogo = TRUE)
getAvailableMatrices(ArchR)

# extract peak matrix
peak_data <- getMatrixFromProject(ArchR, useMatrix = "PeakMatrix", threads = 1)
peak_matrix <- t(assays(peak_data)[[1]])
colnames(peak_matrix) <- rowData(peak_data)$idx # change this to 'names'

# make a test dataset


# cor with different sized matrices
mini_matrix <- peak_matrix[, 1:1000]
mini_matrix_dense <- as.matrix(mini_matrix)
system.time(cor(mini_matrix_dense, method = "spearman"))
#user  system elapsed 
#0.690   0.000   0.692 

mini_matrix <- peak_matrix[, 1:2000]
mini_matrix_dense <- as.matrix(mini_matrix)
system.time(cor(mini_matrix_dense, method = "spearman"))
#user  system elapsed 
#25.113   0.025  25.176 

mini_matrix <- peak_matrix[, 1:5000]
mini_matrix_dense <- as.matrix(mini_matrix)
system.time(cor(mini_matrix_dense, method = "spearman"))
#user  system elapsed 
#150.106   0.891 151.164

mini_matrix <- peak_matrix[, 1:10000]
mini_matrix_dense <- as.matrix(mini_matrix)
system.time(cor(mini_matrix_dense, method = "spearman"))
#user  system elapsed 
#580.238   0.931 582.086 

# plot relative timings using standard cor function:
plot_data <- data.frame(Number_of_peaks = c(1000, 2000, 5000, 10000),
                        Seconds_to_compute = c(0.692, 25.176, 151.164, 582.086))
ggplot(plot_data, aes(x = Number_of_peaks, y = Seconds_to_compute)) + geom_point() + geom_line()


# see how long it takes to rank matrix, then should be able to just feed into pearsons test
system.time(apply(mini_matrix_dense, 2, rank))
#    user  system elapsed 
#     1.156   0.000   1.157 
system.time())
frank(mini_matrix_dense, x = colnames(mini_matrix_dense))
#    user  system elapsed 
#   2.432   0.002   0.643 


ranked_mini_matrix <- apply(mini_matrix_dense, 2, rank)
ranked_mini_matrix[1:5, 1:5]
dim(ranked_mini_matrix)
max(ranked_mini_matrix[,4])




# fastcor with different sized matrices 
# - would have to first turn matrix into ranked one as this is pearson
install.packages("HiClimR") # needs a non-R package to install!
library(HiClimR)
fastCor(xt, nSplit = 1, upperTri = FALSE, optBLAS = FALSE, verbose = TRUE)