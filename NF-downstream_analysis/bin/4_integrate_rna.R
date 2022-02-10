#!/usr/bin/env Rscript

############################## Load libraries #######################################
library(getopt)
library(Signac)
library(Seurat)
library(future)
library(tidyverse)
library(grid)
library(gridExtra)
library(clustree)
library(ggplot2)
library(dplyr)
library(scHelper)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
library(GenomicRanges)

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
    
    setwd("~/NF-downstream_analysis")
    ncores = 8

    plot_path = "../output/NF-downstream_analysis/integrate_rna/plots/"
    rds_path = "../output/NF-downstream_analysis/integrate_rna/rds_files/"
    data_path = "../output/NF-downstream_analysis/test_input/"
  }
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 305* 1024^3)
    plan()
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


############################## Read in Seurat RDS object and fragment files #######################################

seurat <- readRDS(paste0(data_path, "rds_files/seurat_GeneActivity.RDS"))
DefaultAssay(seurat) <- 'peaks'
print(seurat)

# read in fragment files
paths <- list.dirs(paste0(data_path, "cellranger_atac_output/"), recursive = FALSE, full.names = TRUE)
input <- data.frame(sample = sub('.*/', '', paths), 
                    matrix_path = paste0(paths, "/outs/filtered_peak_bc_matrix.h5"),
                    metadata_path = paste0(paths, "/outs/singlecell.csv"),
                    fragments_path = paste0(paths, "/outs/fragments.tsv.gz"))
new.paths <- as.list(input$fragments_path)
frags <- Fragments(seurat)  # get list of fragment objects
Fragments(seurat) <- NULL  # remove fragment information from assay

for (i in seq_along(frags)) {
  frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
}
Fragments(seurat) <- frags # assign updated list back to the object

# read in rna seurat object
seurat_rna <- readRDS(paste0(data_path, "seurat_label_transfer.RDS"))
seurat_rna

# UMAPs
plot1 <- DimPlot(
  object = seurat,
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot2 <- DimPlot(
  object = seurat_rna,
  group.by = 'scHelper_cell_type',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

png(paste0(plot_path, "UMAPs.png"), width=40, height=20, units = 'cm', res = 200)
plot1 + plot2
graphics.off()

# Integrate the RNA and ATAC data
transfer.anchors <- FindTransferAnchors(
  reference = seurat_rna,
  query = seurat,
  reduction = 'cca'
)
print("anchors calculated")

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = seurat_rna$scHelper_cell_type,
  weight.reduction = seurat[['lsi']],
  dims = 2:30
)
print("predicted labels")

seurat <- AddMetaData(object = seurat, metadata = predicted.labels)

# UMAPs
plot1 <- DimPlot(
  object = seurat,
  group.by = 'predicted.labels',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot2 <- DimPlot(
  object = seurat_rna,
  group.by = 'scHelper_cell_type',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

png(paste0(plot_path, "UMAPs_predicted_labels.png"), width=40, height=20, units = 'cm', res = 200)
plot1 + plot2
graphics.off()