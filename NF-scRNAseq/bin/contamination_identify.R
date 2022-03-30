#!/usr/bin/env Rscript

# Load packages
library(getopt)
library(Seurat)
library(future)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
library(scHelper)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    plot_path = "./output/NF-downstream_analysis/6_contamination_filt/plots/"
    rds_path = "./output/NF-downstream_analysis/6_contamination_filt/rds_files/"
    data_path = "./output/NF-downstream_analysis/5_cell_cycle/rds_files/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 32gb
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

#####################################################################################################
#                           Identify contamination (mesoderm and PGCs)                   #
#####################################################################################################

# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"

# Make gene list containing markers used to identify contamination clusters
genes <- list('PGC module' = 'DAZL',
              'Blood island module' = c('CDH5', 'TAL1', 'HBZ'),
              'Mesoderm module' = c('CDX2', 'GATA6', 'ALX1', 'PITX2', 'TWIST1', 'TBXT', 'MESP1'),
              'Endoderm module' = c('SOX17', 'CXCR4', 'FOXA2', 'NKX2-2', 'GATA6'))

# Calculate average module expression for contamination gene list
seurat_data <- AverageGeneModules(seurat_obj = seurat_data, gene_list = genes)

# Plot distribution of contamination gene modules
png(paste0(plot_path, "ContaminationClustersBoxPLot.png"), width = 40, height = 30, units = "cm", res = 200)
PlotCelltype(seurat_obj = seurat_data, gene_list = genes, quantiles = c(0.1, 0.90), ncol = 2)
graphics.off()

PGC_clusters <- IdentifyOutliers(seurat_obj = seurat_data, metrics = 'PGC module', quantiles = c(0.1, 0.90), intersect_metrics = FALSE)
BI_clusters <- IdentifyOutliers(seurat_obj = seurat_data, metrics = 'Blood island module', quantiles = c(0.1, 0.90), intersect_metrics = FALSE)
Mesoderm_clusters <- IdentifyOutliers(seurat_obj = seurat_data, metrics = 'Mesoderm module', quantiles = c(0.1, 0.90), intersect_metrics = FALSE)
Endoderm_clusters <- IdentifyOutliers(seurat_obj = seurat_data, metrics = 'Endoderm module', quantiles = c(0.1, 0.90), intersect_metrics = FALSE)

rownames(filter(seurat_data@meta.data, seurat_clusters %in% PGC_clusters))

which(rownames(seurat_data) %in% )

seurat_data@meta.data$seurat_clusters %in% BI_clusters

seurat_data@meta.data$old_cell_states <- 
  
seurat_data@meta.data %>% mutate(old_cell_states = 
                              )

seurat_data@meta.data <- seurat_data@meta.data %>%
  mutate(old_cell_states = case_when(seurat_clusters %in% PGC_clusters ~ "Contam: PGC",
         seurat_clusters %in% BI_clusters ~ "Contam: BI",
         seurat_clusters %in% Mesoderm_clusters ~ "Contam: Mesoderm",
         seurat_clusters %in% Endoderm_clusters ~ "Contam: Endoderm"))
         

unique(seurat_data@meta.data$old_cell_states)

png(paste0(plot_path, "IdentifiedContam_UMAP.png"), width=40, height=20, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = "old_cell_states")
graphics.off()


saveRDS(seurat_data, paste0(rds_path, "contamination_identified.RDS"), compress = FALSE)

#####################################################################################################
#                           Plots to confirm                   #
#####################################################################################################

# Plot UMAP for poor quality clusters
png(paste0(plot_path, "ContaminationClustersUMAP_PGC.png"), width=40, height=20, units = 'cm', res = 200)
ClusterDimplot(seurat_data, clusters = PGC_clusters, plot_title = 'Contamination PGC')
graphics.off()

png(paste0(plot_path, "ContaminationClustersUMAP_BI.png"), width=40, height=20, units = 'cm', res = 200)
ClusterDimplot(seurat_data, clusters = BI_clusters, plot_title = 'Contamination BI')
graphics.off()

png(paste0(plot_path, "ContaminationClustersUMAP_mesoderm.png"), width=40, height=20, units = 'cm', res = 200)
ClusterDimplot(seurat_data, clusters = Mesoderm_clusters, plot_title = 'Contamination Mesoderm')
graphics.off()

png(paste0(plot_path, "ContaminationClustersUMAP_endoderm.png"), width=40, height=20, units = 'cm', res = 200)
ClusterDimplot(seurat_data, clusters = Endoderm_clusters, plot_title = 'Contamination Endoderm')
graphics.off()

# Plot UMAPs for GOI
ncol = 4
png(paste0(plot_path, "UMAP_GOI.png"), width = ncol*10, height = 10*ceiling((length(unlist(genes))+1)/ncol), units = "cm", res = 200)
MultiFeaturePlot(seurat_obj = seurat_data, gene_list = unlist(genes), plot_clusters = T,
                 plot_stage = T, label = "", cluster_col = "integrated_snn_res.0.5", n_col = ncol)
graphics.off()

# Dotplot for identifying PGCs, Early mesoderm and Late mesoderm
png(paste0(plot_path, "dotplot_GOI.png"), width = 30, height = 12, units = "cm", res = 200)
DotPlot(seurat_data, features = unique(unlist(genes)))
graphics.off()
