#!/usr/bin/env Rscript

### Script to preprocess in ArchR
print("1.5_doublets_ArchR")

############################## Load libraries #######################################
library(getopt)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(parallel)
library(gridExtra)
library(grid)

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
    
    plot_path = "../output/NF-downstream_analysis/1.5_ArchR_doubelts/plots/"
    rds_path = "../output/NF-downstream_analysis/1.5_ArchR_doublets/rds_files/"
    data_path = "../output/NF-downstream_analysis/1_ArchR_preprocessing/rds_files/"

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

############################## Add doublet scores #######################################
ArchR <- addDoubletScores(
  input = ArchR,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1,
  outDir = plot_path, # output is a pdf with some plots and an rds with summary
  logFile = paste0(plot_path, "addDoubletScores"),
  force = TRUE
)
print("added doublet scores")

# check how much memory
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

############################## Plots #######################################
p <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "cellColData", 
  name = "DoubletScore", 
  embedding = "UMAP"
)
png(paste0(plot_path, "UMAP_DoubletScore.png"), width=20, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

############################## FILTERING #######################################
################################################################################

## look at how filtering doublets might affect cell counts
test <- filterDoublets(ArchR, filterRatio = 1)
print("filter ratio = 1")
test
test <- filterDoublets(ArchR, filterRatio = 2)
print("filter ratio = 2")
test
test <- filterDoublets(ArchR, filterRatio = 0.5)
print("filter ratio = 0.5")
test

## add something here to filter based on doublets?
ArchR_filtered <- ArchR

# save ArchR project
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)

############################ POST-FILTERING ####################################
################################################################################

unfiltered <- table(ArchR$stage)
filtered <- table(ArchR_filtered$stage)
cell_counts <- rbind(unfiltered, filtered)

png(paste0(plot_path, 'cell_counts_table.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()