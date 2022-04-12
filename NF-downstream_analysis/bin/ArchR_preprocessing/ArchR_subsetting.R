#!/usr/bin/env Rscript

print("ArchR_subset")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(gridExtra)
library(grid)
library(parallel)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
    make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
    make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
    make_option(c("-m", "--meta_col1"), action = "store", type = "character", help = "Name of first metadata column containing groups to subset", default = NULL),
    make_option(c("", "--meta_col2"), action = "store", type = "character", help = "Name of second metadata column containing groups to subset", default = NULL),
    make_option(c("-o", "--output"), action = "store", type = "character", help = "Name of output RDS file", default = 'seurat_subset'),
    make_option(c("-g", "--groups1"), action = "store", type = "character", help = "Classifications of cells (within meta_col1) to subset from dataset. \
    If multiple classifications are used to subest, must be provided as a comma separated list i.e. --groups celltype1,celltype2", default = NULL),
    make_option(c("", "--groups2"), action = "store", type = "character", help = "Classifications of cells (within meta_col2) to subset from dataset.", default = NULL),
    make_option(c("-i", "--invert1"), action = "store", type = "logical", help = "Boolean for whether to invert groups1 selection", default = FALSE),
    make_option(c("", "--invert2"), action = "store", type = "logical", help = "Boolean for whether to invert groups2 selection", default = FALSE),
    make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = FALSE))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

if(!is.null(opt$groups1)){
    opt$groups1 = strsplit(opt$groups1, ',')[[1]]
}
if(!is.null(opt$groups2)){
opt$groups2 = strsplit(opt$groups2, ',')[[1]]
}

if(is.na(opt$meta_col1)){
    stop("meta_col1 parameter must be provided. See script usage (--help)")
}

if(is.null(opt$groups1)){
    stop("groups1 parameter must be provided. See script usage (--help)")
}

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    
    plot_path = "../output/NF-downstream_analysis/4_ArchR_filter_clusters/plots/"
    rds_path = "../output/NF-downstream_analysis/4_ArchR_filter_clusters/rds_files/"
    data_path = "../output/NF-downstream_analysis/3_ArchR_clustering_prefiltering/rds_files/"

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


############################## Read in ArchR project #######################################
# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

############################## Subset ArchR object #######################################
# If invert is true, then subset the inverted groups from the seurat object
if(opt$invert1 == TRUE){
    opt$groups1 <- as.character(unique(seurat_data@meta.data[[opt$meta_col1]])[!unique(seurat_data@meta.data[[opt$meta_col1]]) %in% opt$groups1])
}
if(opt$invert2 == TRUE){
    opt$groups2 <- as.character(unique(seurat_data@meta.data[[opt$meta_col2]])[!unique(seurat_data@meta.data[[opt$meta_col2]]) %in% opt$groups2])
}

if(!is.null(opt$meta_col2)){
    seurat_subset <- subset(seurat_data, subset = !!as.symbol(opt$meta_col1) %in% opt$groups1 & !!as.symbol(opt$meta_col2) %in% opt$groups2)
    seurat_subset@meta.data[c(opt$meta_col1, opt$meta_col2)] <- lapply(seurat_subset@meta.data[c(opt$meta_col1, opt$meta_col2)], droplevels)
}else{
    seurat_subset <- subset(seurat_data, subset = !!as.symbol(opt$meta_col1) %in% opt$groups1)
    seurat_subset@meta.data[opt$meta_col1] <- lapply(seurat_subset@meta.data[opt$meta_col1], droplevels)
}

############################## Save new ArchR object #######################################
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)

############################## Plots #######################################
seurat_data@meta.data[c(opt$meta_col1, opt$meta_col2)] <- lapply(seurat_data@meta.data[c(opt$meta_col1, opt$meta_col2)], as.factor)

colours = ggPlotColours(length(levels(seurat_data@meta.data[[opt$meta_col1]])))

max_x = round(max(seurat_data@reductions$umap@cell.embeddings[,1]))
min_x = round(min(seurat_data@reductions$umap@cell.embeddings[,1]))
max_y = round(max(seurat_data@reductions$umap@cell.embeddings[,2]))
min_y = round(min(seurat_data@reductions$umap@cell.embeddings[,2]))

full_plot = DimPlot(seurat_data, group.by = opt$meta_col1, pt.size = 0.3, cols = colours, label = TRUE, label.size = 3, label.box = TRUE) +
                ylim(min_y, max_y) +
                xlim(min_x, max_x) +
                theme(legend.position = "none") +
                ggtitle("")

subset_plot = DimPlot(seurat_data, group.by = opt$meta_col1, pt.size = 0.3, cols = colours[levels(seurat_data@meta.data[[opt$meta_col1]]) %in% opt$groups1],
        cells = rownames(seurat_data@meta.data)[seurat_data@meta.data[[opt$meta_col1]] %in% opt$groups1]) +
                ylim(min_y, max_y) +
                xlim(min_x, max_x) +
                theme(legend.position = "none") +
                ggtitle("")

png(paste0(plot_path, 'cell_subset_umap.png'), width = 40, height = 20, units = 'cm', res = 200)
grid.arrange(full_plot, subset_plot, nrow = 1)
graphics.off()