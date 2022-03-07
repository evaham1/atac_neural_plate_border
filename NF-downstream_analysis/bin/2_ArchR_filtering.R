#!/usr/bin/env Rscript

### Script to preprocess in ArchR
print("2_filtering_ArchR")

############################## Load libraries #######################################
library(getopt)
library(future)
library(tidyverse)
library(grid)
library(gridExtra)
library(clustree)
library(GenomeInfoDb)
library(ggplot2)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(parallel)
library(ArchR)
library(GenomicFeatures)

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
    
    plot_path = "../output/NF-downstream_analysis/2_ArchR_filtering/plots/"
    rds_path = "../output/NF-downstream_analysis/2_ArchR_filtering/rds_files/"
    data_path = "../output/NF-downstream_analysis/1_ArchR_preprocessing/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 155* 1024^3)
    addArchRThreads(threads = ncores) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


############################## Read in ArchR project #######################################
ArchR <- loadArchRProject(path = paste0(data_path, "./rds_files/Save-ArchR"), force = FALSE, showLogo = TRUE)


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

############################## QC plots across samples #######################################
##############################################################################################

############################## Plot TSS Enrichment #######################################
p1 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
png(paste0(plot_path, 'TSS_enrichment_plot.png'), height = 15, width = 21, units = 'cm', res = 400)
print(p1)
graphics.off()

p2 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
png(paste0(plot_path, 'TSS_enrichment_vln.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p2)
graphics.off()

p2 <- plotTSSEnrichment(ArchRProj = ArchR)
png(paste0(plot_path, 'TSS_enrichment_plot.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p2)
graphics.off()

############################## Plot log10(Unique Fragments) #######################################
p3 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)
png(paste0(plot_path, 'fragment_count_ridge.png'), height = 15, width = 21, units = 'cm', res = 400)
print(p3)
graphics.off()

p4 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
png(paste0(plot_path, 'fragment_count_vln.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p4)
graphics.off()

############################## Plot nucleosome banding #######################################

p1 <- plotFragmentSizes(ArchRProj = ArchR)
png(paste0(plot_path, 'nucleosome_banding_plot.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p1)
graphics.off()

############# Plot log10(Unique Fragments) vs TSS enrichment score #######################
df <- getCellColData(ArchR, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

png(paste0(plot_path, 'fragments_vs_TSS.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p)
graphics.off()


############################## FILTERING #######################################
################################################################################

test <- filterDoublets(ArchR, filterRatio = 1)
print("filter ratio = 1")
test
test <- filterDoublets(ArchR, filterRatio = 2)
print("filter ratio = 2")
test
test <- filterDoublets(ArchR, filterRatio = 0.5)
print("filter ratio = 0.5")
test

# add something here to filter different samples to different thresholds?

# save ArchR project
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, "Save-ArchR"), load = FALSE)