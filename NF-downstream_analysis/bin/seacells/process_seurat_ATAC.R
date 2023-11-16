#!/usr/bin/env Rscript

print("Script to process seurat object made of ATAC metacells (gene score matrix used)")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(parallel)
library(Seurat)
library(dplyr)
library(tibble)
library(scHelper)
library(ggplot2)
library(future)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-i", "--input"), action = "store", type = "character", help = "Name of seurat input file to process", default = "seacells_seurat.RDS"),
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
    
    # Paths for testing ss8 locally
    data_path = "./local_test_data/convert_seacells_ATAC_to_seurat/rds_files/"
    rds_path = "./local_test_data/processed_seurat_ATAC/rds_files/"
    plot_path = "./local_test_data/processed_seurat_ATAC/plots/"
    
    # Paths for testing HH6 on NEMO
    #data_path = "./output/NF-downstream_analysis/Processing/HH6/SEACELLS_RNA_WF/3_SEACells_metadata_to_seurat/rds_files/"
    #rds_path = "./output/NF-downstream_analysis/Processing/HH6/SEACELLS_RNA_WF/4_Process_metacells/rds_files/"
    #plot_path = "./output/NF-downstream_analysis/Processing/HH6/SEACELLS_RNA_WF/4_Process_metacells/plots/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

###### STAGE COLOURS ####################
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

#####################################################################################
############################    Read in RDS object   #############################
#####################################################################################

seurat <- readRDS(paste0(data_path, opt$input))
print(seurat)
print(paste0("Number of genes in seurat object: ", length(rownames(seurat))))

DefaultAssay(object = seurat) <- "RNA"
DefaultAssay(object = seurat)

## Plot number of seacells and number of genes
df <- data.frame(dim(seurat))
rownames(df) <- c("Gene count: ", "SEACell cout: ")
png(paste0(plot_path, 'metacell_counts.png'), height = 5, width = 12, units = 'cm', res = 400)
grid.arrange(top=textGrob("Gene count and SEACell count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

########## Check for NA values
DefaultAssay(object = seurat) <- "RNA"
DefaultAssay(object = seurat)
sum(is.na(GetAssayData(object = seurat, slot = "counts"))) #0
sum(is.na(GetAssayData(object = seurat, slot = "data"))) #0
sum(is.na(GetAssayData(object = seurat, slot = "scale.data"))) #<0 x 0 matrix>

#####################################################################################
############################    Re-process 'RNA' slot   #############################
#####################################################################################
### Running everything on the 'RNA' slot only

################### Normalise counts #################

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

################### Calculate variable genes #################
#### Or should lift over variable genes from RNA metacells??
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, assay = 'RNA')

################### Dim reduction and clustering without regression #################
seurat <- ScaleData(seurat, features = rownames(seurat)) # scaling without regression

# PCA + UMAP
seurat <- RunPCA(object = seurat, verbose = FALSE)

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(seurat, return = 'plot')
graphics.off()

pc_cutoff <- ElbowCutoff(seurat)
seurat <- FindNeighbors(seurat, dims = 1:pc_cutoff, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:pc_cutoff, verbose = FALSE)

# UMAP of stage
seurat@meta.data$stage <- factor(seurat@meta.data$stage, levels = stage_order)
stage_cols <- stage_colours[levels(droplevels(seurat@meta.data$stage))]

png(paste0(plot_path, "stage_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'stage', label = TRUE, 
        label.size = 9, label.box = TRUE, repel = TRUE,
        pt.size = 10, cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = seurat, by = 0.2, prefix = "RNA_snn_res.")
graphics.off()

# Cluster using default clustering resolution
seurat <- FindClusters(seurat, resolution = 1.2)

# UMAP of clusters
png(paste0(plot_path, "clusters_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'seurat_clusters', label = TRUE, 
        label.size = ifelse(length(unique(seurat$stage)) == 1, 9, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(seurat$stage)) == 1, 6, 6), 
        shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

# Size of clusters
df <- as.data.frame(table(seurat@meta.data$seurat_clusters))
colnames(df) <- c("Cluster", "nCells")
png(paste0(plot_path, 'cluster_metacell_counts.png'), height = 20, width = 12, units = 'cm', res = 400)
grid.arrange(top=textGrob("", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

# Stage
png(paste0(plot_path, "stage_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'stage', label = TRUE, 
        label.size = ifelse(length(unique(seurat$stage)) == 1, 9, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(seurat$stage)) == 1, 6, 6), 
        shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

############################    Save seurat object   #############################

## save seacells seurat object
saveRDS(seurat, paste0(rds_path, "seacells_seurat_processed.RDS"), compress = FALSE)
