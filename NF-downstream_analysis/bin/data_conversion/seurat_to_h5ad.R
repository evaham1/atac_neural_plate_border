#!/usr/bin/env Rscript

print("Script to convert seurat object to h5ad AnnData object for python processing")

# Load packages
library(Seurat)
library(SeuratDisk)
library(scHelper)
library(optparse) 

############################## Set up script options #######################################

# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-d", "--data_path"), action = "store", type = "character", help = "Name of data path that seurat object is in", default = './input/'),
  make_option(c("-i", "--input"), action = "store", type = "character", help = "Name of input seurat object if there is more than one in input"),
  make_option(c("-a", "--assay"), action = "store", type = "character", help = "Assay to export from seurat object ('integrated' or 'RNA')", default = 'integrated'),
  make_option(c("-o", "--outfile"), action = "store", type = "character", help = "Name of outfile"),
  make_option(c("-g", "--group_by"), action = "store", type = "character", help = "Name of metadata column containing groups to colour by", default = 'seurat_clusters'),
  make_option(c("", "--verbose"), action = "store_true", type = "logical", help = "Verbose", default = TRUE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    
    # interative path goes here
    data_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_RNA_WF/5_Classify_metacells/rds_files/"
    input = "Classified_metacells.RDS"
    plot_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_RNA_WF/6_RNA_Anndata_object_processed_classified/plots/"
    rds_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_RNA_WF/6_RNA_Anndata_object_processed_classified/rds_files/"
    
    # local interactive path
    data_path = "./local_test_data/state_classification_ATAC/rds_files/"
    input = "Classified_metacells.RDS"
    plot_path = "./local_test_data/ATAC_metacells_ready_to_integrate/plots/"
    rds_path = "./local_test_data/ATAC_metacells_ready_to_integrate/rds_files/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = opt$data_path
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

set.seed(42)

############################## Read in data #######################################

# If there is only one object in data_path use this, if there is more than one need 'input' arg to define which one to use

input_files <- list.files(data_path, full.names = TRUE, recursive = TRUE)
print(paste0("Input files: ", input_files))

if (length(input_files) > 1){
  print(paste0("Multiple input files detected, reading in only ", opt$input))
  if (is.null(opt$input)) {
    stop("ERROR: input arg not defined!")
  }
  seurat_object <- readRDS(paste0(data_path, opt$input))
} else {
  print("Only one input file detected, reading in now...")
  seurat_object <- readRDS(input_files)
}

print(seurat_object)
print("Seurat object read in successfully")

############################## Convert to H5ad #######################################

DefaultAssay(seurat_object) <- opt$assay
seurat_object <- DietSeurat(seurat_object, counts = TRUE, assays = opt$assay, dimreducs = c('pca', 'umap'))

# remove anything from misc slot as depending on the contents this can cause recursion errors
seurat_object@misc <- list()

# check whats in seurat metadata
print(colnames(seurat_object@meta.data))

# generate cell colours for group_by column
if(!opt[['group_by']] %in% colnames(seurat_object@meta.data)){
  stop('--group_by missing from seurat@meta.data')
}

if(opt$group_by == 'scHelper_cell_type'){
  colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
               "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
               "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3")
  
  cell_state <- factor(seurat_object@meta.data[['scHelper_cell_type']], levels = c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                                                                   'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                                                                   'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                                                                   'vFB', 'aNP', 'node', 'FB', 'pEpi'))
  names(colours) <- levels(cell_state)
  
} else {
  
  # Get ggcolours for cell states
  colours = ggPlotColours(length(unique(seurat_object@meta.data[[opt$group_by]])))
  
  if(class(seurat_object@meta.data[[opt$group_by]]) == 'factor'){
    cell_state <- droplevels(seurat_object@meta.data[[opt$group_by]])
  } else {
    cell_state <- as.factor(seurat_object@meta.data[[opt$group_by]])
  }
  names(colours) <- levels(cell_state)
}

# Add colours to new metadata colummn
seurat_object@meta.data[['cell_colours']] <- unname(colours[cell_state])
# If any na values are present in seurat identities, set colour to grey
seurat_object@meta.data[['cell_colours']][is.na(seurat_object@meta.data[['cell_colours']])] <- '#AAAAAA'

# Convert factor columns to character before converting to h5ad
i <- sapply(seurat_object@meta.data, is.factor)
seurat_object@meta.data[i] <- lapply(seurat_object@meta.data[i], as.character)

print("Converted to h5ad successfully")

############################## Save #######################################

# SaveH5Seurat sometimes encounters a recursion error. File is already written by this point so error can be ignored with try().
#try(SaveH5Seurat(seurat_object, filename = paste0(opt$outfile, '.h5Seurat')), silent = TRUE)
SaveH5Seurat(seurat_object, filename = paste0(opt$outfile, '.h5seurat'))
Convert(paste0(opt$outfile, '.h5seurat'), dest = "h5ad")

print("Saved h5ad successfully")

# Remove intermediate h5Seurat file
file.remove(paste0(opt$outfile, '.h5seurat'))

# Move .h5ad file to ./rds_files/
system(paste0("mv ", opt$outfile, ".h5ad ./rds_files/"))
