#!/usr/bin/env Rscript

print("ArchR_gene_scores")

############################## Load libraries #######################################
library(getopt)
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
library(presto)

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
    
    plot_path = "../output/NF-downstream_analysis/ArchR_split_stages_gene_scores/plots/"
    rds_path = "../output/NF-downstream_analysis/ArchR_split_stages_gene_scores/rds_files/"
    data_path = "../output/NF-downstream_analysis/ArchR_split_stages_clustering/"

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

############################### FUNCTIONS ####################################
# add a function here to extract top differentially expressed genes per cluster
# add a default function to run FeaturePlots and just input list of genes to plot?

############################## Read in ArchR project #######################################
ArchR <- loadArchRProject(path = paste0(data_path, "./rds_files/Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

#data_path = "./work/77/d57378441f7edc2ccab98f84a709da/"

############################## Calculate top gene markers and plot #################################

markers <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markers, cutOff = "Log2FC >= 0.5") # could make more stringent in future
top_markers <- tibble()
for (i in 1:length(markerList)){
  table <- as.tibble(markerList[[i]]) 
  table <- table %>% top_n(5, Log2FC) %>% mutate(cluster = i)
  top_markers <- rbind(top_markers, table)
}
if(nrow(top_markers) != 0){
  print("significant markers found")
  
  png(paste0(plot_path, 'top_genes.png'), height = 100, width = 30, units = 'cm', res = 400)
  grid.arrange(tableGrob(top_markers))
  dev.off()
  
  markerGenes <- top_markers$name
  heatmap <- markerHeatmap(
    seMarker = markers, 
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
    transpose = TRUE
  )
  png(paste0(plot_path, 'heatmap.png'), height = 30, width = 40, units = 'cm', res = 400)
  ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  graphics.off()
  
} else { print("No markers found that passed thresholds")}


############################## Feature plots of known marker genes #################################

# impute weights using MAGIC to plot better feature plots
ArchR <- addImputeWeights(ArchR)

# look for contaminating clusters
contaminating_markers <- c(
  'DAZL', #PGC
  'CDH5', 'TAL1', 'HBZ', # Blood island
  'CDX2', 'GATA6', 'ALX1', 'PITX2', 'TWIST1', 'TBXT', 'MESP1', #mesoderm
  'SOX17', 'CXCR4', 'FOXA2', 'NKX2-2', 'GATA6' #endoderm
)

p <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = contaminating_markers, 
  embedding = "UMAP"
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

png(paste0(plot_path, 'FeaturePlots_contaminating_markers.png'), height = 25, width = 25, units = 'cm', res = 400)
do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
graphics.off()

# look for late marker genes
late_markers <- c(
  "GATA3", "DLX5", "SIX1", "EYA2", #PPR
  "MSX1", "TFAP2A", "TFAP2B", #mix
  "PAX7", "CSRNP1", "SNAI2", "SOX10", #NC
  "SOX2", "SOX21" # neural
  )
p <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = late_markers, 
  embedding = "UMAP"
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
png(paste0(plot_path, 'FeaturePlots_late_markers.png'), height = 30, width = 25, units = 'cm', res = 400)
do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
graphics.off()

# look for ap marker genes
ap_markers <- c(
  "PAX2", "WNT4", "SIX3", "SHH" # no GBX2 in matrix
)
p <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = ap_markers, 
  embedding = "UMAP"
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
png(paste0(plot_path, 'FeaturePlots_ap_markers.png'), height = 15, width = 15, units = 'cm', res = 400)
do.call(cowplot::plot_grid, c(list(ncol = 2),p2))
graphics.off()

# look for early markers
early_markers <- c(
  "EPAS1", "BMP4", "YEATS4", "SOX3", "HOXB1", "ADMP", "EOMES")

p <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = early_markers, 
  embedding = "UMAP"
)
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
png(paste0(plot_path, 'FeaturePlots_early_markers.png'), height = 15, width = 15, units = 'cm', res = 400)
do.call(cowplot::plot_grid, c(list(ncol = 2),p2))
graphics.off()
