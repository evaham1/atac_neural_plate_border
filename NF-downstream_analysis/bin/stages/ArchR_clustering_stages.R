#!/usr/bin/env Rscript

print("ArchR_clustering_stages")

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
library(hexbin)

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
    
    #plot_path = "../output/NF-downstream_analysis/3_ArchR_clustering/plots/"
    #data_path = "../output/NF-downstream_analysis/2_ArchR_filtering/"
    #rds_path = "../output/NF-downstream_analysis/3_ArchR_clustering/rds_files/"
    
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

addArchRThreads(threads = 8) 

############################## Read in ArchR project #######################################
# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

############################## Iterative LSI Dimensionality Reduction #################################
# need to look at setting seed for this
ArchR <- addIterativeLSI(
  ArchRProj = ArchR,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  seed = 1,
  force = TRUE
)
print("iterative LSI ran")

############################## Seurat graph-based clustering #################################
ArchR <- addClusters(
  input = ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "clusters",
  resolution = 2,
  force = TRUE
)
print("clustering ran")
print(table(ArchR$clusters))
# add a barchart here of cell numbers in each cluster?

############################## Which stage in which clusters (if have more than one stage) #################################
if (length(unique(ArchR$stage)) > 1){
  cM <- confusionMatrix(paste0(ArchR$clusters), paste0(ArchR$stage))
  cM <- cM / Matrix::rowSums(cM)
  p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
  )
  png(paste0(plot_path, "Cluster_stage_distribution.png"), width=25, height=20, units = 'cm', res = 200)
  print(p)
  graphics.off()
}

############################## Run and plot UMAP #################################
ArchR <- addUMAP(
  ArchRProj = ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)
print("UMAP added")

p1 <- plotEmbedding(ArchRProj = ArchR, colorBy = "cellColData", name = "stage", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = ArchR, colorBy = "cellColData", name = "clusters", embedding = "UMAP")

png(paste0(plot_path, "UMAPs.png"), width=60, height=40, units = 'cm', res = 200)
ggAlignPlots(p1, p2, type = "h")
graphics.off()

paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, "Save-ArchR"), load = FALSE)

############################## Plot QC metrics on UMAP #################################

p <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "cellColData", 
  name = "nFrags", 
  embedding = "UMAP"
)
png(paste0(plot_path, "UMAP_nFrags.png"), width=20, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

p <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "cellColData", 
  name = "DoubletScore", 
  embedding = "UMAP"
)
png(paste0(plot_path, "UMAP_DoubletScore.png"), width=20, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

p <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "cellColData", 
  name = "NucleosomeRatio", 
  embedding = "UMAP"
)
png(paste0(plot_path, "UMAP_NucleosomeRatio.png"), width=20, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

p <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "cellColData", 
  name = "TSSEnrichment", 
  embedding = "UMAP"
)
png(paste0(plot_path, "UMAP_TSSEnrichment.png"), width=20, height=20, units = 'cm', res = 200)
print(p)
graphics.off()