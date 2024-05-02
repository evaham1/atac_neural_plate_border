#!/usr/bin/env Rscript

print("Remove the HH4 stage from the full data transfer label RNA object")
# adapted from Alex's subset_cluster R script

############################## Load libraries #######################################
# Load packages
library(optparse)
library(Seurat)
library(future)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(scHelper)

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

  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/"
    rds_path = "./rds_files/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

set.seed(42)

############################## Load data #######################################

# input can be named anything but needs to be only one input file
label <- unique(list.files(data_path))
print(label) 

# Load seurat data
seurat_data <- readRDS(paste0(data_path, label))

############################## Remove HH4 #######################################

seurat_data <- subset(x = seurat_data, subset = stage == "HH4", invert = TRUE)
seurat_data

############################## Reprocess new subset #######################################

# Set RNA to default assay
DefaultAssay(seurat_data) <- "RNA"

# Re-run findvariablefeatures and scaling
seurat_data <- FindVariableFeatures(seurat_data, selection.method = "vst", nfeatures = 2000, assay = 'RNA')

# seurat_data <- ScaleData(seurat_data, features = rownames(seurat_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# Set Integrated to default assay
DefaultAssay(seurat_data) <- "integrated"

# Rescale data on integrated assay
seurat_data <- ScaleData(seurat_data, features = rownames(seurat_data), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

# PCA
seurat_data <- RunPCA(object = seurat_data, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=50, units = 'cm', res = 200)
DimHeatmap(seurat_data, dims = 1:30, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(seurat_data, return = 'plot')
graphics.off()

# automatically determine elbow
pc_cutoff <- ElbowCutoff(seurat_data)

# # if pc_cutoff is smaller than 7 then don't run with pc_cutoff-5 as too small to run UMAP
# cutoffs = ifelse(pc_cutoff < 7, c(pc_cutoff, pc_cutoff+5, pc_cutoff+10, pc_cutoff+15), c(pc_cutoff-5, pc_cutoff, pc_cutoff+5, pc_cutoff+10))
# png(paste0(plot_path, "UMAP_PCA_comparison.png"), width=40, height=30, units = 'cm', res = 200)
# PCALevelComparison(seurat_data, PCA_levels = cutoffs, cluster_res = opt$clustres)
# graphics.off()

seurat_data <- FindNeighbors(seurat_data, dims = 1:pc_cutoff, verbose = FALSE)
seurat_data <- RunUMAP(seurat_data, dims = 1:pc_cutoff, verbose = FALSE)

# # Find optimal cluster resolution
# png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
# ClustRes(seurat_object = seurat_data, by = 0.2, prefix = "integrated_snn_res.")
# graphics.off()

# # Cluster data
# seurat_data <- FindClusters(seurat_data, resolution = opt$clustres)

# # Plot UMAP for clusters and developmental stage
# png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
# ClustStagePlot(seurat_data, stage_col = "stage")
# graphics.off()

# Plot UMAP for developmental stage, clusters and integration (if subset contains more than one batch)
# schelper cell type colours
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                "#A5718D", "#3F918C", "#ed5e5f", "9792A3")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                       'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo',
                                       'Neural', 'Placodal', 'Non-neural', 'Contam')

# stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_cols = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_cols) <- stage_order

# set colour palettes for UMAPs
cols <- scHelper_cell_type_colours[as.character(unique(seurat_data$scHelper_cell_type))]

p1 <- DimPlot(seurat_data, pt.size = 1, reduction = "umap", group.by = "scHelper_cell_type", cols = cols, shuffle = TRUE)
p2 <- DimPlot(seurat_data, pt.size = 1, reduction = "umap", group.by = "stage", cols = stage_cols, shuffle = TRUE)

png(paste0(plot_path, 'UMAPs.png'), height = 10, width = 30, units = 'cm', res = 400)
p1 + p2
graphics.off()

############################## Save data #######################################

saveRDS(seurat_data, paste0(rds_path, "seurat_label_transfer_minus_HH4.RDS"), compress = FALSE)


