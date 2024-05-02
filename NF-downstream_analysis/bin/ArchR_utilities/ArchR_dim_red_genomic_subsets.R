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
    data_path = "./output/NF-downstream_analysis/ArchR_processing/HH5/INTEGRATING/1_unconstrained_integration/rds_files/"
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

set.seed(42)

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
                              'PGC', 'BI', 'meso', 'endo',
                              'Neural', 'Placodal', 'Non-neural', 'Contam')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                "#A5718D", "#3F918C", "#ed5e5f", "#9792A3")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                       'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo',
                                       'Neural', 'Placodal', 'Non-neural', 'Contam')
cols <- scHelper_cell_type_colours[unique(ArchR$transferred_scHelper_cell_type)]
cols_broad <- scHelper_cell_type_colours[unique(ArchR$transferred_scHelper_cell_type_broad)]


###############################################################################################
############################## Rerun dim reduction with peaks #################################

plot_path = "./plots/UMAPs_with_different_dim_red/"
dir.create(plot_path, recursive = T)

########################## Original dim reduction on tile matrix ###############################
ArchR_original <- ArchR

png(paste0(plot_path, 'Original_tile_matrix.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR,
              name = "transferred_scHelper_cell_type", pal = cols,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

png(paste0(plot_path, "Original_tile_matrix_broad.png"), width=30, height=40, units = 'cm', res = 200)
plotEmbedding(ArchR,
              name = "transferred_scHelper_cell_type_broad", pal = cols_broad,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

########################## Dim reduction on peak matrix: ALL peaks ###############################

ArchR_peaks <- addIterativeLSI(ArchR, useMatrix = "PeakMatrix", force = TRUE)
print("iterative LSI ran")
ArchR_peaks <- addUMAP(ArchR_peaks, force = TRUE)
print("UMAP added")
ArchR_peaks <- addClusters(ArchR_peaks, name = "clusters", resolution = 1, force = TRUE)
print("clustering ran")

png(paste0(plot_path, 'All_peaks.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_peaks,
              name = "transferred_scHelper_cell_type", pal = cols,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

png(paste0(plot_path, 'All_peaks_broad.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_peaks,
              name = "transferred_scHelper_cell_type_broad", pal = cols_broad,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()


########################## Dim reduction on peak matrix: DISTAL peaks ###############################

peak_set <- getPeakSet(ArchR)

# extract distal peaks
distal_peak_set <- peak_set[which(peak_set$peakType == "Distal"), ]

# overwrite peakset of ArchR object to only include distal peaks
ArchR_distal_peaks <- addPeakSet(ArchR, peakSet = distal_peak_set, force = TRUE)
getPeakSet(ArchR_distal_peaks)
ArchR_distal_peaks <- addPeakMatrix(ArchR_distal_peaks)

# re-run dim reduction with new peak matrix
ArchR_distal_peaks <- addIterativeLSI(ArchR_distal_peaks, useMatrix = "PeakMatrix", force = TRUE)
print("iterative LSI ran")
ArchR_distal_peaks <- addUMAP(ArchR_distal_peaks, force = TRUE)
print("UMAP added")
ArchR_distal_peaks <- addClusters(ArchR_distal_peaks, name = "clusters", resolution = 1, force = TRUE)
print("clustering ran")

png(paste0(plot_path, 'Distal_peaks.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_distal_peaks,
              name = "transferred_scHelper_cell_type", pal = cols,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

png(paste0(plot_path, 'Distal_peaks_broad.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_distal_peaks,
              name = "transferred_scHelper_cell_type_broad", pal = cols_broad,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

########################## Dim reduction on peak matrix: PROMOTER peaks ###############################

promoter_peak_set <- peak_set[which(peak_set$peakType == "Promoter"), ]

# overwrite peakset of ArchR object to only include distal peaks
ArchR_promoter_peaks <- addPeakSet(ArchR, peakSet = promoter_peak_set, force = TRUE)
getPeakSet(ArchR_promoter_peaks)
ArchR_promoter_peaks <- addPeakMatrix(ArchR_promoter_peaks)

# re-run dim reduction with new peak matrix
ArchR_promoter_peaks <- addIterativeLSI(ArchR_promoter_peaks, useMatrix = "PeakMatrix", force = TRUE)
print("iterative LSI ran")
ArchR_promoter_peaks <- addUMAP(ArchR_promoter_peaks, force = TRUE)
print("UMAP added")
ArchR_promoter_peaks <- addClusters(ArchR_promoter_peaks, name = "clusters", resolution = 1, force = TRUE)
print("clustering ran")

png(paste0(plot_path, 'Promoter_peaks.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_promoter_peaks,
              name = "transferred_scHelper_cell_type", pal = cols,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

png(paste0(plot_path, 'Promoter_peaks_broad.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_promoter_peaks,
              name = "transferred_scHelper_cell_type_broad", pal = cols_broad,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

########################## Dim reduction on peak matrix: PROMOTER+EXONIC peaks ###############################

genes_peak_set <- peak_set[which(peak_set$peakType %in% c("Exonic", "Promoter")), ]

# overwrite peakset of ArchR object to only include distal peaks
ArchR_genes_peaks <- addPeakSet(ArchR, peakSet = genes_peak_set, force = TRUE)
getPeakSet(ArchR_genes_peaks)
ArchR_genes_peaks <- addPeakMatrix(ArchR_genes_peaks)

# re-run dim reduction with new peak matrix
ArchR_genes_peaks <- addIterativeLSI(ArchR_genes_peaks, useMatrix = "PeakMatrix", force = TRUE)
print("iterative LSI ran")
ArchR_genes_peaks <- addUMAP(ArchR_genes_peaks, force = TRUE)
print("UMAP added")
ArchR_genes_peaks <- addClusters(ArchR_genes_peaks, name = "clusters", resolution = 1, force = TRUE)
print("clustering ran")

png(paste0(plot_path, 'Promoter_and_exonic_peaks.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_genes_peaks,
              name = "transferred_scHelper_cell_type", pal = cols,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

png(paste0(plot_path, 'Promoter_and_exonic_peaks_broad.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_genes_peaks,
              name = "transferred_scHelper_cell_type_broad", pal = cols_broad,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

########################## Dim reduction on peak matrix: EXONIC peaks ###############################

exonic_peak_set <- peak_set[which(peak_set$peakType %in% c("Exonic")), ]

# overwrite peakset of ArchR object to only include distal peaks
ArchR_exonic_peaks <- addPeakSet(ArchR, peakSet = exonic_peak_set, force = TRUE)
getPeakSet(ArchR_exonic_peaks)
ArchR_exonic_peaks <- addPeakMatrix(ArchR_exonic_peaks)

# re-run dim reduction with new peak matrix
ArchR_exonic_peaks <- addIterativeLSI(ArchR_exonic_peaks, useMatrix = "PeakMatrix", force = TRUE)
print("iterative LSI ran")
ArchR_exonic_peaks <- addUMAP(ArchR_exonic_peaks, force = TRUE)
print("UMAP added")
ArchR_exonic_peaks <- addClusters(ArchR_exonic_peaks, name = "clusters", resolution = 1, force = TRUE)
print("clustering ran")

png(paste0(plot_path, 'Exonic_peaks.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_exonic_peaks,
              name = "transferred_scHelper_cell_type", pal = cols,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

png(paste0(plot_path, 'Exonic_peaks_broad.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_exonic_peaks,
              name = "transferred_scHelper_cell_type_broad", pal = cols_broad,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

########################## Dim reduction on peak matrix: INTRONIC peaks ###############################

intronic_peak_set <- peak_set[which(peak_set$peakType %in% c("Intronic")), ]

# overwrite peakset of ArchR object to only include distal peaks
ArchR_intronic_peaks <- addPeakSet(ArchR, peakSet = intronic_peak_set, force = TRUE)
getPeakSet(ArchR_intronic_peaks)
ArchR_intronic_peaks <- addPeakMatrix(ArchR_intronic_peaks)

# re-run dim reduction with new peak matrix
ArchR_intronic_peaks <- addIterativeLSI(ArchR_intronic_peaks, useMatrix = "PeakMatrix", force = TRUE)
print("iterative LSI ran")
ArchR_intronic_peaks <- addUMAP(ArchR_intronic_peaks, force = TRUE)
print("UMAP added")
ArchR_intronic_peaks <- addClusters(ArchR_intronic_peaks, name = "clusters", resolution = 1, force = TRUE)
print("clustering ran")

png(paste0(plot_path, 'Intronic_peaks.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_intronic_peaks,
              name = "transferred_scHelper_cell_type", pal = cols,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

png(paste0(plot_path, 'Intronic_peaks_broad.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_intronic_peaks,
              name = "transferred_scHelper_cell_type_broad", pal = cols_broad,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

########################## Dim reduction on peak matrix: INTRONIC+DISTAL peaks ###############################

ex_peak_set <- peak_set[which(peak_set$peakType %in% c("Intronic", "Distal")), ]

# overwrite peakset of ArchR object to only include distal peaks
ArchR_ex_genes_peaks <- addPeakSet(ArchR, peakSet = ex_peak_set, force = TRUE)
getPeakSet(ArchR_ex_genes_peaks)
ArchR_ex_genes_peaks <- addPeakMatrix(ArchR_ex_genes_peaks)

# re-run dim reduction with new peak matrix
ArchR_ex_genes_peaks <- addIterativeLSI(ArchR_ex_genes_peaks, useMatrix = "PeakMatrix", force = TRUE)
print("iterative LSI ran")
ArchR_ex_genes_peaks <- addUMAP(ArchR_ex_genes_peaks, force = TRUE)
print("UMAP added")
ArchR_ex_genes_peaks <- addClusters(ArchR_ex_genes_peaks, name = "clusters", resolution = 1, force = TRUE)
print("clustering ran")

png(paste0(plot_path, 'Intronic_and_distal_peaks.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_ex_genes_peaks,
              name = "transferred_scHelper_cell_type", pal = cols,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

png(paste0(plot_path, 'Intronic_and_distal_peaks_broad.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR_ex_genes_peaks,
              name = "transferred_scHelper_cell_type_broad", pal = cols_broad,
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()