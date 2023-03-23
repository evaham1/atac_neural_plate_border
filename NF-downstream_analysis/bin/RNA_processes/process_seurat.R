#!/usr/bin/env Rscript

print("Script to process seurat object made of metacells")

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
    data_path = "./local_test_data/convert_seacells_to_seurat/rds_files/"
    rds_path = "./local_test_data/processed_seurat/rds_files/"
    plot_path = "./local_test_data/processed_seurat/plots/"
    
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

#####################################################################################
############################    Read in RDS object   #############################
#####################################################################################

seurat <- readRDS(paste0(data_path, opt$input))
print(seurat)

DefaultAssay(object = seurat) <- "RNA"
DefaultAssay(object = seurat)

########## Remove genes expressed in fewer than 5 cells
seurat <- DietSeurat(seurat, features = names(which(Matrix::rowSums(GetAssayData(seurat) > 0) >=5)))
seurat

########## Check for NA values
DefaultAssay(object = seurat) <- "RNA"
DefaultAssay(object = seurat)
sum(is.na(GetAssayData(object = seurat, slot = "counts"))) #0
sum(is.na(GetAssayData(object = seurat, slot = "data"))) #0
sum(is.na(GetAssayData(object = seurat, slot = "scale.data"))) #<0 x 0 matrix>

DefaultAssay(object = seurat) <- "integrated"
DefaultAssay(object = seurat)
sum(is.na(GetAssayData(object = seurat, slot = "counts"))) #0
sum(is.na(GetAssayData(object = seurat, slot = "data"))) #0
sum(is.na(GetAssayData(object = seurat, slot = "scale.data"))) #<0 x 0 matrix>

#####################################################################################
############################    Re-process 'RNA' slot   #############################
#####################################################################################

plot_path = "./plots/RNA_assay/"

########## Factors to regress out: MT percent, sex, cell cycle

DefaultAssay(object = seurat) <- "RNA"
DefaultAssay(object = seurat)

# 1) MT percentage: re-calculate using raw counts
seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
print("Added percent mt")
head(seurat@meta.data)

# 2) Sex: Have already added proportion of male to column 'sex'

# 3) Cell cycle: re-calculate using raw counts
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
print("Added cell cycle score")
head(seurat@meta.data)

## Normalising and Scaling
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, assay = 'RNA')
seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

print("raw counts re-processed")

seurat
head(seurat@assays$RNA@scale.data)

seurat_temp <- seurat

## Dim reduction
seurat <- RunPCA(object = seurat, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=65, units = 'cm', res = 200)
DimHeatmap(seurat, dims = 1:20, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(seurat, return = 'plot')
graphics.off()

pc_cutoff <- ElbowCutoff(seurat)

## Find neighbours and calculate UMAP
seurat <- FindNeighbors(seurat, dims = 1:pc_cutoff, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:pc_cutoff, verbose = FALSE)

print("dim reduction calculated on raw counts")

############################    Visualise on UMAPs   #############################

# Clusters
png(paste0(plot_path, "stage_UMAP.png"), width=40, height=20, units = 'cm', res = 200)
DimPlot(seurat, group.by = "stage", pt.size = 6)
graphics.off()

# QC metrics
png(paste0(plot_path, "percent.mt_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
FeaturePlot(object = seurat, features = "percent.mt", pt.size = 6)
graphics.off()

png(paste0(plot_path, "run_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
FeaturePlot(seurat, features = "run", pt.size = 6)
graphics.off()

png(paste0(plot_path, "sex_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
FeaturePlot(seurat, features = "sex", pt.size = 6)
graphics.off()

############################    Save seurat object   #############################

## save seacells seurat object
saveRDS(seurat, paste0(rds_path, "seacells_seurat_RNA_processed.RDS"), compress = FALSE)

############################################################################################
############################    Re-process 'Integrated' slot   #############################
############################################################################################

seurat <- seurat_temp

plot_path = "./plots/integrated_assay/"

DefaultAssay(object = seurat) <- "integrated"
DefaultAssay(object = seurat)

## Normalising and Scaling
seurat <- NormalizeData(seurat, normalization.method = "RC", scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, assay = "integrated")
seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

seurat
head(seurat@assays$integrated@scale.data)

print("integrated data re-processed")

## Dim reduction
seurat <- RunPCA(object = seurat, verbose = FALSE)

png(paste0(plot_path, "dimHM.png"), width=30, height=65, units = 'cm', res = 200)
DimHeatmap(seurat, dims = 1:20, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(seurat, return = 'plot')
graphics.off()

pc_cutoff <- ElbowCutoff(seurat)

## Find neighbours and calculate UMAP
seurat <- FindNeighbors(seurat, dims = 1:pc_cutoff, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:pc_cutoff, verbose = FALSE)

print("dim reduction calculated on integrated data")

############################    Visualise on UMAPs   #############################

# Clusters
png(paste0(plot_path, "stage_UMAP.png"), width=40, height=20, units = 'cm', res = 200)
DimPlot(seurat, group.by = "stage", pt.size = 6)
graphics.off()

# QC metrics
png(paste0(plot_path, "percent.mt_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
FeaturePlot(object = seurat, features = "percent.mt", pt.size = 6)
graphics.off()

png(paste0(plot_path, "run_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
FeaturePlot(object = seurat, features = "run", pt.size = 6)
graphics.off()

png(paste0(plot_path, "sex_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
FeaturePlot(seurat, features = "sex", pt.size = 6)
graphics.off()

############################    Save seurat object   #############################

## save seacells seurat object
saveRDS(seurat, paste0(rds_path, "seacells_seurat_integrated_processed.RDS"), compress = FALSE)
