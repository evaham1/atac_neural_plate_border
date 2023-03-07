#!/usr/bin/env Rscript

print("converting from h5ad to seurat")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(parallel)
library(Seurat)
library(anndata)
library(SeuratDisk)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-i", "--input"), action = "store", type = "character", help = "Name of input file to convert"),
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
    data_path = "./local_test_data/SEACells_from_camp/SEACELLS_RNA_WF/SEACells_computation/rds_files/"
    opt$input = "AnnData_summarised_by_metacells.h5ad"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    label = "AnnData_summarised_by_metacells.h5ad"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


#################### R script to convert AnnData object to Seurat #########################

# https://satijalab.org/seurat/archive/v2.4/conversion_vignette.html



data_seurat@misc




# Use python anndata to import h5ad object
data_anndata <- read_h5ad(paste0(data_path, opt$input))

# Convert to seurat
data_seurat <- CreateSeuratObject(counts = t(as.matrix(data_anndata$X)), meta.data = data_anndata$obs)

pbmc3k <- Convert(data_anndata, to = "seurat")

head(data_anndata$varm)


Convert(paste0(data_path, opt$input), dest = "h5seurat", overwrite = TRUE)
pbmc3k <- LoadH5Seurat("pbmc3k_final.h5seurat")
pbmc3k
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 0 variable features)
#>  2 dimensional reductions calculated: pca, umap


# Some plots to check its worked
p1 <- TSNEPlot(seurat, group.by = "louvain", do.return = TRUE)
p2 <- VlnPlot(seurat, c("CST3", "NKG7", "PPBP"), group.by = "louvain", do.return = TRUE)
plot_grid(p1, p2)

# Save RDS
saveRDS(seurat, paste0(rds_path, "seurat.RDS"), compress = FALSE)

