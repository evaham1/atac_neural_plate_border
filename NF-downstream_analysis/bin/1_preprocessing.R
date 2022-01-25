#!/usr/bin/env Rscript

### Script to create seurat object, filter data, remove poor quality clusters and integrate across batches

############################## Load libraries #######################################
library(getopt)
library(Signac)
library(Seurat)
library(future)
library(tidyverse)
library(grid)
library(gridExtra)
library(clustree)
library(GenomeInfoDb)
library(ggplot2)
library(dplyr)
library(rtracklayer)


############################## Set up script options #######################################
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)
test = TRUE

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    setwd("~/NF-downstream_analysis")
    ncores = 8
    ref_path = "../output/NF-luslab_sc_multiomic/reference/"
    
    if(test == TRUE){
      plot_path = "../output/NF-downstream_analysis/1_preprocessing/TEST/plots/"
      rds_path = "../output/NF-downstream_analysis/1_preprocessing/TEST/rds_files/"
      data_path = "../output/NF-luslab_sc_multiomic/test/cellranger_atac_output/"
      }else{
      plot_path = "../output/NF-downstream_analysis/1_preprocessing/plots/"
      rds_path = "../output/NF-downstream_analysis/1_preprocessing/rds_files/"
      data_path = "../output/NF-luslab_sc_multiomic/full/cellranger_atac_output/"}
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/cellranger_atac_output/"
    ref_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 155* 1024^3)
    plan()
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Read in data and set up Signac object #######################################

# Make dataframe with stage and replicate info extracted from path
paths <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)

input <- data.frame(sample = sub('.*/', '', paths), 
                   matrix_path = paste0(paths, "/outs/filtered_peak_bc_matrix.h5"),
                   metadata_path = paste0(paths, "/outs/singlecell.csv"),
                   fragments_path = paste0(paths, "/outs/fragments.tsv.gz"))
print("path df made")

# Read in the 3 files needed in list format
counts_list <- apply(input, 1, function(x) Read10X_h5(filename = x[["matrix_path"]]))
metadata_list <- apply(input, 1, function(x) read.csv(file = x[["metadata_path"]], header = TRUE, row.names = 1))
fragments_list <- as.list(input$fragments_path)

print("data files read in")

# Build list of assays using these files
chrom_assays <- lapply(1:nrow(input), function(x) CreateChromatinAssay(
                                                    counts = counts_list[[x]],
                                                    sep = c(":", "-"),
                                                    fragments = fragments_list[[x]],
                                                    min.cells = 10,
                                                    min.features = 200))
print("ChromatinAssays made")

signac_datas <- lapply(1:nrow(input), function(x) CreateSeuratObject(
  counts = chrom_assays[[x]],
  project = input$sample[x],
  assay = "peaks",
  meta.data = metadata_list[[x]]))
print("seurat objects made")

# add annotations using chick gtf
gtf <- rtracklayer::import(paste0(ref_path, "genes.gtf.gz"))
gene.coords <- gtf[gtf$type == 'gene']

signac_list <- lapply(signac_datas, function(x) SetAssayData(x, slot = "annotation", new.data = gene.coords))
print("Annotations set")

# Init list of signac objects then merge
names(signac_list) <- input$sample
seurat_all <- merge(x = signac_list[[1]], y=signac_list[-1], add.cell.ids = names(signac_list), project = "chick.10x.atac")
print("seurat objects merged")

# Add metadata col for stage and flow cell
seurat_all@meta.data[["stage"]] <- substr(seurat_all@meta.data$orig.ident, 1, 3)
seurat_all@meta.data[["flow_cell"]] <- substr(seurat_all@meta.data$orig.ident, 5, 5)

# Convert metadata character cols to factors
seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)] <- lapply(seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)], as.factor)

# Calculate QC metrics
seurat_all <- NucleosomeSignal(object = seurat_all)
seurat_all <- TSSEnrichment(object = seurat_all, fast = FALSE)
seurat_all$pct_reads_in_peaks <- seurat_all$peak_region_fragments / seurat_all$passed_filters * 100

# save RDS
saveRDS(seurat_all, paste0(rds_path, "seurat_all.RDS"), compress = FALSE)