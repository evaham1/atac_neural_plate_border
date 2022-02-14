#!/usr/bin/env Rscript

############################## Load libraries #######################################
library(getopt)
library(Signac)
library(Seurat)
library(future)
library(tidyverse)
library(grid)
library(gridExtra)
library(clustree)
library(ggplot2)
library(dplyr)
library(scHelper)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
library(GenomicRanges)

############################## Set up script options #######################################
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)
test = TRUE

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    setwd("~/NF-downstream_analysis")
    ncores = 8
    ref_path = "../output/NF-luslab_sc_multiomic/reference/"
    
    if(test == TRUE){
      plot_path = "../output/NF-downstream_analysis/gene_activity/TEST/plots/"
      rds_path = "../output/NF-downstream_analysis/gene_activity/TEST/rds_files/"
      data_path = "../output/NF-downstream_analysis/test_input/"
    }else{
      plot_path = "../output/NF-downstream_analysis/gene_activity/plots/"
      rds_path = "../output/NF-downstream_analysis/gene_activity/rds_files/"
      data_path = "../output/NF-downstream_analysis/test_input/"}
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 305* 1024^3)
    plan()
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


############################## Read in Seurat RDS object and fragment files #######################################

seurat <- readRDS(paste0(data_path, "rds_files/seurat_GeneActivity.RDS"))
print(seurat)
DefaultAssay(seurat) <- 'RNA'

# # read in fragment files
# paths <- list.dirs(paste0(data_path, "cellranger_atac_output/"), recursive = FALSE, full.names = TRUE)
# input <- data.frame(sample = sub('.*/', '', paths), 
#                     matrix_path = paste0(paths, "/outs/filtered_peak_bc_matrix.h5"),
#                     metadata_path = paste0(paths, "/outs/singlecell.csv"),
#                     fragments_path = paste0(paths, "/outs/fragments.tsv.gz"))
# new.paths <- as.list(input$fragments_path)
# frags <- Fragments(seurat)  # get list of fragment objects
# Fragments(seurat) <- NULL  # remove fragment information from assay

# for (i in seq_along(frags)) {
#   frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
# }
# Fragments(seurat) <- frags # assign updated list back to the object


######################################## SEX EFFECT #####################################################
sex_plot_path = paste0(plot_path, "sex_effect/")
dir.create(sex_plot_path, recursive = T)

png(paste0(sex_plot_path, 'FeaturePlot_sex_genes.png'), height = 20, width = 30, units = 'cm', res = 400)
FeaturePlot(
  object = seurat,
  features = c('W-Wpkci-7', 'W-NIPBLL', 'Z-SEMA6A', 'Z-ARID3C', 'W-Wpkci-7'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
graphics.off()

# Use W chromosome genes to K-means cluster the cells into male (zz) and female (zw)
W_genes <- as.matrix(seurat@assays$RNA[grepl("W-", rownames(seurat@assays$RNA)),])
k_clusters <- kmeans(t(W_genes), 2)
k_clusters <- data.frame(k_clusters$cluster)
seurat@meta.data$k_clusters <- k_clusters[match(colnames(seurat@assays$RNA), rownames(k_clusters)),]

# Get rownames for kmeans clusters 1 and 2
k_clus_1 <- rownames(seurat@meta.data[seurat@meta.data$k_clusters == 1,])
k_clus_2 <- rownames(seurat@meta.data[seurat@meta.data$k_clusters == 2,])

# K clustering identities are stochastic, so I need to identify which cluster is male and female
# Sum of W genes is order of magnitude greater in cluster 2 - these are the female cells
sumclus1 <- sum(W_genes[,k_clus_1])
sumclus2 <- sum(W_genes[,k_clus_2])

if(sumclus1 < sumclus2){
  k_male <- k_clus_1
  k_female <- k_clus_2
} else {
  k_female <- k_clus_1
  k_male <- k_clus_2
}

# Add sex data to meta.data
seurat@meta.data$sex <- unlist(lapply(rownames(seurat@meta.data), function(x)
  if(x %in% k_male){"male"} else if(x %in% k_female){"female"} else{stop("cell sex is not assigned")}))

png(paste0(sex_plot_path,"UMAP_sex.png"), height = 18, width = 18, units = "cm", res = 200)
DimPlot(seurat, group.by = "sex")
graphics.off()

# Make dataframe for mean Z expression in male cells
mean_Z_male <- data.frame(Z.mean = apply(seurat@assays$RNA[grepl("Z-", rownames(seurat@assays$RNA)), k_male], 1, mean))
# add 1 before log2 as log2(1) = 0
mean_Z_male <- log2(mean_Z_male + 1)

# Make dataframe for mean Z expression in female cells
mean_Z_female <- data.frame(Z.mean = apply(seurat@assays$RNA[grepl("Z-", rownames(seurat@assays$RNA)), k_female], 1, mean))
mean_Z_female <- log2(mean_Z_female + 1)

# Make dataframe for mean autosomal expression in male cells
mean_auto_male <- data.frame(auto.mean = apply(seurat@assays$RNA[!grepl("Z-", rownames(seurat@assays$RNA)) & !grepl("W-", rownames(seurat@assays$RNA)), k_male], 1, mean))
mean_auto_male <- log2(mean_auto_male + 1)

# Make dataframe for mean autosomal expression in male cells
mean_auto_female <- data.frame(auto.mean = apply(seurat@assays$RNA[!grepl("Z-", rownames(seurat@assays$RNA)) & !grepl("W-", rownames(seurat@assays$RNA)), k_female], 1, mean))
mean_auto_female <- log2(mean_auto_female + 1)

# Calculate FC by subtracting log2 expression from each other
FC <- list()
FC$Z <- mean_Z_male - mean_Z_female
FC$auto <-  mean_auto_male - mean_auto_female

# Plot boxplot of Z gene and autosomal expression in male vs female cells
png(paste0(sex_plot_path,"sex_kmeans_log2FC_boxplot.png"), height = 18, width = 18, units = "cm", res = 200)
boxplot(c(FC$Z, FC$auto),  ylab = "male - female log2 FC (mean normalised UMI +1)", names = c("Z chromosome genes", "autosomal genes"))
graphics.off()

############################## Top differentially expressed genes #######################################

# find marker genes for clusters and visualise them
seurat <- ScaleData(object = seurat, verbose = TRUE)
seurat

markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

png(paste0(plot_path, 'top_markers.png'), height = 40, width = 40, units = 'cm', res = 400)
grid.arrange(tableGrob(top10, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path, 'top_markers_heatmap.png'), height = 10, width = 20, units = 'cm', res = 400)
DoHeatmap(seurat, features = top10$gene) + NoLegend()
graphics.off()

###########################################################################################################


saveRDS(seurat, paste0(rds_path, "seurat_GeneActivity.RDS"))