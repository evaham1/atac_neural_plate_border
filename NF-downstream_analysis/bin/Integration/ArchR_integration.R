#!/usr/bin/env Rscript

print("ArchR integration")

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
library(Seurat)

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
    
    ncores = 8
    
    plot_path = "./output/NF-downstream_analysis/ArchR_integration/plots/"
    rds_path = "./output/NF-downstream_analysis/ArchR_integration/rds_files/"
    data_path_atac = "./output/NF-downstream_analysis/ArchR_clustering_postfiltering_twice/rds_files/"
    data_path_rna = "./output/NF-scRNAseq/minus_HH4_clustered/rds_files/"
    
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

############################## Read in ArchR project and seurat object #######################################
### SORT OUT INPUTS - in pairs RNA and ATAC objects

# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# load seurat object by reading in any rds object
rna_path <- list.files(path = data_path, pattern = "*.RDS", full.names = TRUE)
seurat <- readRDS(rna_path)

############################## UMAPs before integration #######################################
# UMAPs of RNA and ATAC data, with RNA coloured by cell state and ATAC by clusters
umap_rna <- DimPlot(seurat, group.by = "scHelper_cell_type")
umap_atac <- plotEmbedding(ArchR, embedding = "UMAP", colorBy = "cellColData", name = "clusters")

png(paste0(plot_path, 'UMAPs_before_integration.png'), height = 30, width = 20, units = 'cm', res = 400)
print(umap_rna + umap_atac)
graphics.off()

############################## Unconstrained integration #######################################
ArchR <- addGeneIntegrationMatrix(
  ArchRProj = ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seurat,
  addToArrow = TRUE,
  groupRNA = "scHelper_cell_type",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  force = TRUE
)
print("integration completed")

png(paste0(plot_path, 'UMAP_unconInt_scHelper_cell_type.png'), height = 20, width = 20, units = 'cm', res = 400)
print(plotEmbedding(ArchR, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un"))
graphics.off()

png(paste0(plot_path, 'UMAP_unconInt_Scores.png'), height = 20, width = 20, units = 'cm', res = 400)
print(plotEmbedding(ArchR, 
  colorBy = "cellColData", 
  name = "predictedScore_Un"))
graphics.off()

cM <- as.matrix(confusionMatrix(ArchR$clusters, ArchR$predictedGroup_Un))
scHelper_cell_types <- colnames(cM)[apply(cM, 1 , which.max)]
cluster_idents <- cbind(scHelper_cell_types, rownames(cM))

png(paste0(plot_path, 'cluster_labels_table.png'), height = 30, width = 10, units = 'cm', res = 400)
grid.arrange(tableGrob(cluster_idents, rows=NULL, theme = ttheme_minimal()))
graphics.off()

############################## Constrained integration #######################################
########  do I want to do this?? 

# # which seurat clusters correspond to which cell types?
# seurat_df <- data.frame(cluster_id = seurat$seurat_clusters,
#                         scHelper_cell_type = seurat$scHelper_cell_type)
# seurat_df <- seurat_df %>% remove_rownames(.) %>% distinct(., .keep_all = FALSE)
#Assign scRNA cell states to cluster numbers - automate for every cell state
#clust_FB_rna <- paste0(seurat_df$cluster_id[which(seurat_df$scHelper_cell_type == "FB")], collapse = "|")

# Identify scRNAseq cells that are in these clusters
# rna_FB_cells <- names(seurat$scHelper_cell_type)[which(seurat$scHelper_cell_type == "FB")]
# 
# #Assign scATAC to cell states to cluster numbers - automate for every cell state
# clust_FB_atac <- rownames(cM)[grep('FB', preClust)]
# atac_FB_cells <-ArchR$cellNames[ArchR$clusters %in% clust_FB_atac]
# 
# groupList <- SimpleList(
#   FB = SimpleList(
#     ATAC = atac_FB_cells,
#     RNA = rna_FB_cells
#   )
# )
# 
# # doesnt work - needs to include all ATAC cells
# # how would I do this, every cell type? then not going to differ..
# # just divide neural/NNE?? contam??
# ArchR <- addGeneIntegrationMatrix(
#   ArchRProj = ArchR, 
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "IterativeLSI",
#   seRNA = seurat,
#   addToArrow = FALSE, 
#   groupList = groupList,
#   groupRNA = "BioClassification",
#   nameCell = "predictedCell_Co",
#   nameGroup = "predictedGroup_Co",
#   nameScore = "predictedScore_Co"
# )
  

############################## Visualising pseudo-scRNA profiles #######################################
ArchR <- addImputeWeights(ArchR)
  
# look for late marker genes
late_markers <- c(
  "GATA3", "DLX5", "SIX1", "EYA2", #PPR
  "MSX1", "TFAP2A", "TFAP2B", #mix
  "PAX7", "CSRNP1", "SNAI2", "SOX10", #NC
  "SOX2", "SOX21" # neural
)
p1 <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = late_markers, 
  embedding = "UMAP"
)
p2 <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "GeneIntegrationMatrix", 
  name = late_markers, 
  embedding = "UMAP"
)

p1c <- lapply(p1, function(x){
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

p2c <- lapply(p2, function(x){
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

do.call(cowplot::plot_grid, c(list(ncol = 3), p1c))
do.call(cowplot::plot_grid, c(list(ncol = 3), p2c))




