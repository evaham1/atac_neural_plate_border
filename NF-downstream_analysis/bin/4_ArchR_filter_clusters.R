#!/usr/bin/env Rscript

print("4_filter_clusters_ArchR")

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
    
    plot_path = "../output/NF-downstream_analysis/4_ArchR_filter_clusters/plots/"
    rds_path = "../output/NF-downstream_analysis/4_ArchR_filter_clusters/rds_files/"
    data_path = "../output/NF-downstream_analysis/3_ArchR_clustering/"
    
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

############################### FUNCTIONS #################################################
ArchR_IdentifyOutliers <- function(ArchR, group_by = 'Clusters', metrics, intersect_metrics = TRUE, quantiles){
  outlier <- list()
  if(!length(quantiles) == 2){
    stop('quantiles must be an array of length == 2')
  }
  for(metric in metrics){
    min = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1])
    max = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2])
    
    outlier[[metric]] <- as.tibble(getCellColData(ArchR)) %>%
      group_by((!!as.symbol(group_by))) %>%
      summarise(median = median((!!as.symbol(metric)))) %>%
      filter(median > max | median < min) %>%
      pull(!!as.symbol(group_by))
  }
  
  if(intersect_metrics){
    if(length(Reduce(intersect, outlier)) == 0){
      cat('No outliers detected!')
    } else {
      return(Reduce(intersect, outlier))
    }
  } else{
    if(length(as.character(unique(unlist(outlier)))) == 0){
      cat('No outliers detected!')
    } else {
      return(as.character(unique(unlist(outlier))))
    }
  }
}

############################## Read in ArchR project #######################################
ArchR <- loadArchRProject(path = paste0(data_path, "./rds_files/Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")


############################## nFrags #######################################
p1 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "clusters", 
  colorBy = "cellColData", 
  name = "nFrags",
  plotAs = "Violin"
)
png(paste0(plot_path, "VlnPlot_nFrags.png"), width=50, height=20, units = 'cm', res = 200)
p1
graphics.off()

# for now will not filter on this metric


############################## TSS Enrichment #######################################
p <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "clusters", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

png(paste0(plot_path, "VlnPlot_TSSEnrichment.png"), width=50, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

# violin with quantiles overlaid
metrics = "TSSEnrichment"
quantiles = c(0.2, 0.8)
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1]), linetype = "dashed", 
                   color = "red")
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2]), linetype = "dashed", 
                   color = "red")
png(paste0(plot_path, "VlnPlot_thresholds_TSSEnrichment.png"), width=50, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

# automatically identify outlier clusters using adapted scHelper function
outliers <- ArchR_IdentifyOutliers(ArchR, group_by = 'clusters', metrics = metrics, intersect_metrics = FALSE, quantiles = quantiles)

# highlight outlier clusters on UMAP
if (is.null(outliers) == FALSE){
  idxSample <- BiocGenerics::which(ArchR$clusters %in% outliers)
  cellsSample <- ArchR$cellNames[idxSample]
  p <- plotEmbedding(ArchR, colorBy = "cellColData", name = "clusters", embedding = "UMAP", highlightCells = cellsSample)
  png(paste0(plot_path, "UMAP_TSSEnrichment_outliers.png"), width=20, height=20, units = 'cm', res = 200)
  print(p)
  graphics.off()
}

############################## Nucleosome signal #######################################
p <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "clusters", 
  colorBy = "cellColData", 
  name = "NucleosomeRatio",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

png(paste0(plot_path, "VlnPlot_NucleosomeRatio.png"), width=50, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

# violin with quantiles overlaid
metrics = "NucleosomeRatio"
quantiles = c(0.2, 0.8)
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1]), linetype = "dashed", 
                   color = "red")
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2]), linetype = "dashed", 
                   color = "red")
png(paste0(plot_path, "VlnPlot_thresholds_NucleosomeRatio.png"), width=50, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

# automatically identify outlier clusters using adapted scHelper function
outliers <- ArchR_IdentifyOutliers(ArchR, group_by = 'clusters', metrics = metrics, intersect_metrics = FALSE, quantiles = quantiles)

# highlight outlier clusters on UMAP
if (is.null(outliers) == FALSE){
  idxSample <- BiocGenerics::which(ArchR$clusters %in% outliers)
  cellsSample <- ArchR$cellNames[idxSample]
  p <- plotEmbedding(ArchR, colorBy = "cellColData", name = "clusters", embedding = "UMAP", highlightCells = cellsSample)
  png(paste0(plot_path, "UMAP_NucleosomeRatio_outliers.png"), width=20, height=20, units = 'cm', res = 200)
  print(p)
  graphics.off()
}

##########################################################################################################################
############################## Filter on TSSEnrichment + NucleosomeRatio intersect #######################################
metrics = c("TSSEnrichment", "NucleosomeRatio")
quantiles = c(0.2, 0.8)
outliers <- ArchR_IdentifyOutliers(ArchR, group_by = 'clusters', metrics = metrics, intersect_metrics = TRUE, quantiles = quantiles)

# highlight intersect outlier clusters on UMAP
if (is.null(outliers) == FALSE){
  idxSample <- BiocGenerics::which(ArchR$clusters %in% outliers)
  cellsSample <- ArchR$cellNames[idxSample]
  p <- plotEmbedding(ArchR, colorBy = "cellColData", name = "clusters", embedding = "UMAP", highlightCells = cellsSample)
  png(paste0(plot_path, "UMAP_intersect_outliers.png"), width=20, height=20, units = 'cm', res = 200)
  print(p)
  graphics.off()
}

# filter ArchR object
if (is.null(outliers) == FALSE){
  ArchR <- addCellColData(ArchRProj = ArchR, data = rep("poor_quality", length(cellsSample)),
                            cells = cellsSample, name = "quality", force = TRUE)
  idxPass <- which(is.na(ArchR$quality) == TRUE)
  cellsPass <- ArchR$cellNames[idxPass]
  ArchR <- ArchR[cellsPass, ]
}

saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, "Save-ArchR"), load = FALSE)