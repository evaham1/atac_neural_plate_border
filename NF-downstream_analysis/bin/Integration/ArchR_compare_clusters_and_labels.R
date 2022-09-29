#!/usr/bin/env Rscript

print("Compare the relationship of RNA labels with ATAC clusters")

############################## Load libraries #######################################
library(getopt)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)
library(presto)
library(Seurat)
library(plyr)
library(gtools)

############################## Set up script options #######################################
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    
    #data_path = "./output/NF-downstream_analysis/ArchR_integration//ss8/1_unconstrained_integration/rds_files/"
    data_path = "./output/NF-downstream_analysis/ArchR_integration/HH5/1_unconstrained_integration/rds_files/"
    #data_path = "./output/NF-downstream_analysis/ArchR_integration/FullData/1_unconstrained_integration/rds_files/"
    #plot_path = "./output/NF-downstream_analysis/ArchR_integration/FullData/2_identify_clusters/rds_files/"
    
    
    addArchRThreads(threads = 1) 
    
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

############################## Read in ArchR project  #######################################

# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

getAvailableMatrices(ArchR)

############################################################################################
############################## COLOURS #######################################

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

###### schelper cell type colours
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak', 
                              'PGC', 'BI', 'meso', 'endo')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                                "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                                "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                       'vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo')
# set colour palettes for UMAPs
atac_scHelper_old_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_old)]


###############################################################################################
############################## Rerun dim reduction with peaks #################################

plot_path = "./plots/UMAPs_with_different_dim_red"
dir.create(plot_path, recursive = T)

########################## Original dim reduction on tile matrix ###############################
ArchR_original <- ArchR

png(paste0(plot_path, 'Original_tile_matrix.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_original, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE)
graphics.off()

png(paste0(plot_path, 'Original_tile_matrix_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_original, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_old_cols)
graphics.off()

########################## Dim reduction on peak matrix: ALL peaks ###############################

ArchR_peaks <- addIterativeLSI(ArchR, useMatrix = "PeakMatrix", force = TRUE)
print("iterative LSI ran")
ArchR_peaks <- addUMAP(ArchR_peaks, force = TRUE)
print("UMAP added")
ArchR_peaks <- addClusters(ArchR_peaks, name = "clusters", resolution = 1, force = TRUE)
print("clustering ran")

png(paste0(plot_path, 'All_peaks.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_peaks, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE)
graphics.off()

png(paste0(plot_path, 'All_peaks_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_peaks, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_old_cols)
graphics.off()


########################## Dim reduction on peak matrix: DISTAL peaks ###############################

peak_matrix <- getMatrixFromProject(ArchR, useMatrix = "PeakMatrix", threads = 1)
peak_data <- getPeakSet(ArchR)
distal_indices <- peak_data[which(peak_data$peakType == "Distal"), ][, 12]

test <- peak_matrix[which(rowData(peak_matrix)$idx == "1"), ]
rowData(test)


# add code here to dim reduce, cluster and plot with this new peakset




