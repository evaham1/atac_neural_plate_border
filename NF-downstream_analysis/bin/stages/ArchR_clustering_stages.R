#!/usr/bin/env Rscript

print("ArchR_clustering_stages")

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
# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

#################################################################################
############################## PROCESSING #######################################

############## Iterative LSI Dimensionality Reduction ###########################
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

################## Seurat graph-based clustering #################################
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

# Plot number of cells in each cluster
cluster_cell_counts <- as.data.frame(table(ArchR$clusters))
colnames(cluster_cell_counts) <- c("Cluster ID", "Cell Count")

png(paste0(plot_path, 'cell_counts_table.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(tableGrob(cluster_cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

p<-ggplot(data=cluster_cell_counts, aes(x=`Cluster ID`, y=`Cell Count`)) +
  geom_bar(stat="identity")
png(paste0(plot_path, 'cell_counts_barchart.png'), height = 10, width = 10, units = 'cm', res = 400)
print(p)
graphics.off()

# Plot contribution of each stage to each cluster
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

#################################################################################
############################### QC PLOTS ########################################

############################ Plot QC metrics on UMAP #############################

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

######################## QC Vioin Plots #######################################

quantiles = c(0.2, 0.8)

##### nFrags
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

#### TSS Enrichment
p <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "clusters", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
metrics = "TSSEnrichment"
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

#### Nucleosome signal
p <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "clusters", 
  colorBy = "cellColData", 
  name = "NucleosomeRatio",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
metrics = "NucleosomeRatio"
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