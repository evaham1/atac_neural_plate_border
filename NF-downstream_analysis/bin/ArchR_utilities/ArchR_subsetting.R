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
    make_option(c("-n", "--meta_col2"), action = "store", type = "character", help = "Name of second metadata column containing groups to subset", default = NULL),
    make_option(c("-g", "--groups1"), action = "store", type = "character", help = "Classifications of cells (within meta_col1) to subset from dataset. \
    If multiple classifications are used to subest, must be provided as a comma separated list i.e. --groups celltype1,celltype2", default = NULL),
    make_option(c("-p", "--groups2"), action = "store", type = "character", help = "Classifications of cells (within meta_col2) to subset from dataset.", default = NULL),
    make_option(c("-i", "--invert1"), action = "store", type = "logical", help = "Boolean for whether to invert groups1 selection", default = FALSE),
    make_option(c("-q", "--invert2"), action = "store", type = "logical", help = "Boolean for whether to invert groups2 selection", default = FALSE),
    make_option(c("-v", "--verbose"), action = "store", type = "logical", help = "Verbose", default = TRUE),
    make_option(c("", "--invert"), action = "store", type = "logical", help = "Invert subset", default = FALSE)
    )

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    addArchRThreads(threads = 1)
    
    data_path = "./output/NF-downstream_analysis/ArchR_integration/ss8/1_unconstrained_integration/rds_files/"
    opt$invert1 = TRUE
    opt$groups1 = "BI,PGC,meso,endo"
    opt$meta_col1 = "scHelper_cell_type_old"
    opt$meta_col2 = NULL
    
    data_path = "./output/NF-downstream_analysis/ArchR_peak_exploration/transfer_labels/rds_files/"
    opt$invert1 = FALSE
    opt$groups1 = "HH7_C4,HH7_C5,HH7_C6"
    opt$meta_col1 = "stage_clusters"
  
    
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

# Set up options
if(!is.null(opt$groups1)){
  opt$groups1 = strsplit(opt$groups1, ',')[[1]]}
if(!is.null(opt$groups2)){
  opt$groups2 = strsplit(opt$groups2, ',')[[1]]}

if(is.null(opt$meta_col1)){
  stop("meta_col1 parameter must be provided. See script usage (--help)")}
if(is.null(opt$groups1)){
  stop("groups1 parameter must be provided. See script usage (--help)")}

set.seed(42)

############################## Function to subset ArchR project #######################################
ArchR_Subset <- function(ArchR_object, meta_col1, meta_col2, groups1, groups2, invert1, invert2, invert = FALSE){
  
  print(paste0("Filtering on first meta_col: ", meta_col1))
  if (invert1 == FALSE){
    idxPass_1 <- which(ArchR_object@cellColData[,meta_col1] %in% groups1)
  } else { idxPass_1 <- which(!(ArchR_object@cellColData[,meta_col1] %in% groups1)) }
  idxPass <- idxPass_1
  
  if (is.null(meta_col2) == FALSE){
    print(paste0("Filtering on second meta_col: ", meta_col2))
    if (invert2 == FALSE){
      idxPass_2 <- which(ArchR_object@cellColData[,meta_col2] %in% groups2) 
    } else { idxPass_2 <- which(!(ArchR_object@cellColData[,meta_col2] %in% groups2)) }
    idxPass <- idxPass_1[(idxPass_1 %in% idxPass_2)] # take intersect of pass 1 and pass 2
  }
  
  cellsPass <- ArchR$cellNames[idxPass]
  ArchR_subset <- ArchR[cellsPass, ]
}

############################## Read in ArchR project #######################################

# If files are not in rds_files subdirectory look in input dir 
label <- unique(sub('_.*', '', list.files(data_path)))
print(label) 

if (length(label) == 0){
  data_path = "./input/"
  label <- sub('_.*', '', list.files(data_path))
  print(label)
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
} else {
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
}

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getPeakSet(ArchR)
getAvailableMatrices(ArchR)

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

#####################################################################################
############################## Subset ArchR object #######################################

ArchR_subset <- ArchR_Subset(ArchR, meta_col1 = opt$meta_col1, meta_col2 = opt$meta_col2, 
                             groups1 = opt$groups1, groups2 = opt$groups2,
                             invert1 = opt$invert1, invert2 = opt$invert2)
print("ArchR object subsetted")

saveArchRProject(ArchRProj = ArchR_subset, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)

#####################################################################################
############################## Visualisations #######################################

### Plot cell counts before and after subsetting per stage
unfiltered <- table(ArchR$stage)
filtered <- table(ArchR_subset$stage)
cell_counts <- as.data.frame(dplyr::bind_rows(unfiltered, filtered))
cell_counts[is.na(cell_counts)] = 0

png(paste0(plot_path, 'cell_counts_table_stages.png'), height = 10, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

### Plot cell counts before and after subsetting per cluster
unfiltered <- table(ArchR$clusters)
filtered <- table(ArchR_subset$clusters)
cell_counts <- as.data.frame(dplyr::bind_rows(unfiltered, filtered))
cell_counts[is.na(cell_counts)] = 0

png(paste0(plot_path, 'cell_counts_table_clusters.png'), height = 10, width = 30, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

### Plot cell counts before and after subsetting per stage_cluster
if ("stage_clusters" %in% colnames(ArchR@cellColData)){
  unfiltered <- table(ArchR$stage_clusters)
  filtered <- table(ArchR_subset$stage_clusters)
  cell_counts <- as.data.frame(dplyr::bind_rows(unfiltered, filtered))
  cell_counts[is.na(cell_counts)] = 0
  
  png(paste0(plot_path, 'cell_counts_table_stage_clusters.png'), height = 10, width = 60, units = 'cm', res = 400)
  grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
               tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
  graphics.off()
}

### Plot cell counts before and after subsetting per sc_helper_cell_type
if ("scHelper_cell_type_old" %in% colnames(ArchR@cellColData)){
  unfiltered <- table(ArchR$scHelper_cell_type_old)
  filtered <- table(ArchR_subset$scHelper_cell_type_old)
  cell_counts <- as.data.frame(dplyr::bind_rows(unfiltered, filtered))
  cell_counts[is.na(cell_counts)] = 0

  png(paste0(plot_path, 'cell_counts_table_scHelper_cell_type_old.png'), height = 10, width = 30, units = 'cm', res = 400)
  grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
              tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
  graphics.off()
}