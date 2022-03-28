#!/usr/bin/env Rscript

# Load packages
library(getopt)
library(Seurat)
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

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores', 'c', 2, "integer",
  'split', 's', 2, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
plot_path = "./plots/"
rds_path = "./rds_files/"
data_path = "./input/"

dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

############################## UMAPs #######################################

# Plot UMAP for clusters and developmental stage
png(paste0(plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(seurat_data, stage_col = "stage")
graphics.off()

# Plot UMAP for developmental stage, clusters and integration (if subset contains more than one batch)
plots <- list()
if(length(unique(seurat_data$stage)) > 1){
    plots$stage_plot <- DimPlot(seurat_data, group.by = "stage") + ggtitle(paste("Developmental stage")) + theme(plot.title = element_text(hjust = 0.5))
}

plots$cluster_plot <- DimPlot(seurat_data, group.by = opt$meta_col) + ggtitle(paste(opt$meta_col)) + theme(plot.title = element_text(hjust = 0.5))

if(length(unique(seurat_data$run)) > 1){
    plots$integration_plot <- DimPlot(seurat_data, group.by = "run") + ggtitle(paste("Batches")) + theme(plot.title = element_text(hjust = 0.5))
}

png(paste0(plot_path, "UMAP.png"), width=20*length(plots), height=20, units = 'cm', res = 200)
do.call("grid.arrange", c(plots, nrow=1))
graphics.off()