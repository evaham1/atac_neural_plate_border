#!/usr/bin/env Rscript

print("clustering_ArchR")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)
library(clustree)

############################## Set up script options #######################################

addArchRThreads(threads = 1)
plot_path = "./output/NF-downstream_analysis/ArchR_QC_exploration/HH5/plots/"
dir.create(plot_path, recursive = T)

############################### FUNCTIONS - adapted from scHelper #################################################
ArchR_IdentifyOutliers <- function(ArchR, group_by = 'clusters', metrics, intersect_metrics = TRUE, quantiles){
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

############################## HH5 poor quality cells #######################################

data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/HH5/ArchR_clustering/rds_files/"
label <- sub('_.*', '', list.files(data_path))
print(label)
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

quantiles = c(0.2, 0.8)
metrics = "NucleosomeRatio"

# automatically identify outlier clusters using adapted scHelper function
outliers <- ArchR_IdentifyOutliers(ArchR, group_by = 'clusters', metrics = metrics, intersect_metrics = FALSE, quantiles = quantiles)

idxSample <- BiocGenerics::which(ArchR$clusters %in% outliers)
cellsSample <- ArchR$cellNames[idxSample]

### read in full data to map back to there

data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/NF-scATACseq_alignment_out/ArchR_clustering/rds_files/"
label <- sub('_.*', '', list.files(data_path))
print(label)
ArchR_FULL <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

png(paste0(plot_path, "UMAP_HH5_poor_quality_on_fulldata.png"), width=20, height=20, units = 'cm', res = 200)
print(plotEmbedding(ArchR_FULL, name = "clusters", highlightCells = cellsSample,
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 20, labelSize = 14, legendSize = 0, randomize = TRUE, labelAsFactors = FALSE))
graphics.off()


data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/HH6/ArchR_clustering/rds_files/"
label <- sub('_.*', '', list.files(data_path))
print(label)
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

quantiles = c(0.2, 0.8)
metrics = "NucleosomeRatio"

# automatically identify outlier clusters using adapted scHelper function
outliers <- ArchR_IdentifyOutliers(ArchR, group_by = 'clusters', metrics = metrics, intersect_metrics = FALSE, quantiles = quantiles)

idxSample <- BiocGenerics::which(ArchR$clusters %in% outliers)
cellsSample <- ArchR$cellNames[idxSample]
png(paste0(plot_path, "UMAP_HH6_poor_quality_on_fulldata.png"), width=20, height=20, units = 'cm', res = 200)
print(plotEmbedding(ArchR_FULL, name = "clusters", highlightCells = cellsSample,
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 20, labelSize = 14, legendSize = 0, randomize = TRUE, labelAsFactors = FALSE))
graphics.off()