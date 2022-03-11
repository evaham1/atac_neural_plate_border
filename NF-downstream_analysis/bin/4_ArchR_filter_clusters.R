#!/usr/bin/env Rscript

print("4_filter_clusters_ArchR")

############################## Load libraries #######################################
library(getopt)
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
    data_path = "../output/NF-downstream_analysis/3_ArchR_clustering_prefiltering/"

    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

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
  ArchR_filtered <- ArchR[cellsPass, ]
} else { 
  ArchR_filtered <- ArchR
}

# save filtered ArchR project
saveArchRProject(ArchRProj = ArchR_filtered, outputDirectory = paste0(rds_path, "Save-ArchR"), load = FALSE)

# plot cell counts before and after filtering
unfiltered <- table(ArchR$stage)
filtered <- table(ArchR_filtered$stage)
cell_counts <- rbind(unfiltered, filtered)

png(paste0(plot_path, 'cell_counts_table.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()