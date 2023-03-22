#!/usr/bin/env Rscript

print("script to transfer seacell ids to seurat object AND create and process seurat object summarised by seacells")

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
  make_option(c("-m", "--metadata_file_name"), action = "store", type = "character", help = "Name of csv file which assigns cell ids to metacell ids", default = "Cell_metadata.csv"),
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
    data_path = "./local_test_data/test_inputs/test_input_seacells_meta_to_seurat/"
    rds_path = "./local_test_data/convert_seacells_to_seurat/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

################### Colours ########################

scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                                "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                                "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                       'vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo')

################### Functions ##########################

## function to aggregate matrix from seurat object and summarise by cell groupings
summarise_seurat_data <- function(seurat, data_slot = "counts", category = "SEACell"){
  
  # extract data into dataframe
  df <- GetAssayData(object = seurat, slot = data_slot)
  df <- as.data.frame(t(as.data.frame(df)))
  
  # convert cell ids to category ids
  category_ids <- select(seurat@meta.data, category)[,1]
  df <- df %>% mutate(category = category_ids)
  
  # aggregate df based on category
  df_summarised <- aggregate(. ~ category, df, sum)
  
  # format df so can be added back to seurat object
  df_summarised <- t(column_to_rownames(df_summarised, var = "category"))
  
  return(df_summarised)
}

#################### Read in data and add metadata to seurat object #########################

# Read in metadata (use metadata_file_name to find correct object, search only in /rds_files/ )
metacell_metadata <- read.csv(paste0(data_path, "/rds_files/", opt$metadata_file_name))
metacell_dictionary <- select(metacell_metadata, c("index", "SEACell"))

print("Metacell IDs read in")
print(head(metacell_dictionary))

# Read in seurat object (identify label first, file needs to be named as [LABEL]_clustered_data.RDS, only looks in input not in /rds_files/)
label <- setdiff(sub('_.*', '', list.files(data_path)), "rds")
print(label)
seurat <- readRDS(paste0(data_path, label, "_clustered_data.RDS"))

print("Seurat object read in")
print(seurat)

# Reorder seacells metadata to match cell order in seurat object
metacell_dictionary <- metacell_dictionary[match(rownames(seurat@meta.data), metacell_dictionary$index),]

print(head(metacell_dictionary))

# Add seacells metadata to seurat object
seurat <- AddMetaData(seurat, metacell_dictionary$SEACell, col.name = "SEACell")

print("SEACell IDs added to seurat object")

# Save seurat object
saveRDS(seurat, paste0(rds_path, "seurat.RDS"), compress = FALSE)

#################### Create new seurat object with only metacells #########################

# Filter seurat object to only include SEACells
seacells_seurat <- seurat[,colnames(seurat) %in% unique(metacell_dictionary$SEACell)]

print("new seurat object creates with only seacells")

#################### Add up counts across metacells #########################

###### RNA slot: 3 assays: counts (raw), data (normalised), scale.data -> only add up raw 'counts'
DefaultAssay(object = seacells_seurat) <- "RNA"
DefaultAssay(object = seacells_seurat)

summarised_RNA_counts <- summarise_seurat_data(seurat = seurat, data_slot = "counts", category = "SEACell")
seacells_seurat <- SetAssayData(object = seacells_seurat, slot = "counts", new.data = summarised_RNA_counts, assay = "RNA")

print("raw counts summarised")

###### Integrated slot: 2 assays: data, scale.data -> only add up 'data'
DefaultAssay(object = seacells_seurat) <- "integrated"
DefaultAssay(object = seacells_seurat)

summarised_integrated_data <- summarise_seurat_data(seurat = seurat, data_slot = "data", category = "SEACell")
seacells_seurat <- SetAssayData(object = seacells_seurat, slot = "data", new.data = summarised_integrated_data, assay = "integrated")

print("integrated counts summarised")

#####################################################################################
############################    Re-process 'RNA' slot   #############################
#####################################################################################

DefaultAssay(object = seacells_seurat) <- "RNA"

## Factors to regress out: MT percent, sex, cell cycle

# 1) MT percentage: re-calculate using raw counts
seacells_seurat <- PercentageFeatureSet(seacells_seurat, pattern = "^MT-", col.name = "percent.mt")

# 2) Sex: use the already calculated sex of the metacell ID OR recalculate a 'sex score' based on aggregated gene counts?

# 3) Cell cycle: re-calculate using raw counts
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seacells_seurat <- CellCycleScoring(seacells_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

## Dim reduction and scaling

seacells_seurat_processing <- NormalizeData(seacells_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seacells_seurat_processing <- FindVariableFeatures(seacells_seurat_processing, selection.method = "vst", nfeatures = 2000, assay = 'RNA')
seacells_seurat_processing <- ScaleData(seacells_seurat_processing, features = rownames(seacells_seurat_processing), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

print("raw counts re-processed")

############################################################################################
############################    Re-process 'Integrated' slot   #############################
############################################################################################

DefaultAssay(object = seacells_seurat) <- "integrated"

## Dim reduction and scaling

seacells_seurat_processing <- NormalizeData(seacells_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seacells_seurat_processing <- FindVariableFeatures(seacells_seurat_processing, selection.method = "vst", nfeatures = 2000, assay = 'RNA')
seacells_seurat_processing <- ScaleData(seacells_seurat_processing, features = rownames(seacells_seurat_processing), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))

png(paste0(plot_path, "dimHM.png"), width=30, height=65, units = 'cm', res = 200)
DimHeatmap(seacells_seurat_processing, dims = 1:20, balanced = TRUE, cells = 500)
graphics.off()

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(seacells_seurat_processing, return = 'plot')
graphics.off()

pc_cutoff <- ElbowCutoff(seacells_seurat_processing)

## Find neighbours and calculate UMAP

seacells_seurat_processing <- FindNeighbors(seacells_seurat_processing, dims = 1:pc_cutoff, verbose = FALSE)
seacells_seurat_processing <- RunUMAP(seacells_seurat_processing, dims = 1:pc_cutoff, verbose = FALSE)

## Visualise on UMAPs

png(paste0(plot_path, "UMAPs.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(seacells_seurat_processing)
graphics.off()

png(paste0(plot_path, "scHelper_cell_type_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
DimPlot(seacells_seurat_processing, group.by = "scHelper_cell_type", cols = scHelper_cell_type_colours, pt.size = 10)
graphics.off()

png(paste0(plot_path, "sex_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
DimPlot(seacells_seurat_processing, group.by = "sex", pt.size = 10)
graphics.off()

png(paste0(plot_path, "run_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
DimPlot(seacells_seurat_processing, group.by = "run", pt.size = 10)
graphics.off()

png(paste0(plot_path, "percent.mt_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
FeaturePlot(object = seacells_seurat_processing, features = "percent.mt", pt.size = 10)
graphics.off()

print("integrated counts re-processed")

## save seacells seurat object
saveRDS(seacells_seurat, paste0(rds_path, "seacells_seurat.RDS"), compress = FALSE)