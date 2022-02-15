#!/usr/bin/env Rscript

### Script to filter cells based on predicted GEX
print("4_gex_filtering")

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

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    setwd("~/NF-downstream_analysis")
    ncores = 8

    plot_path = "../output/NF-downstream_analysis/gex_filtering/plots/"
    rds_path = "../output/NF-downstream_analysis/gex_filtering/rds_files/"
    data_path = "../output/NF-downstream_analysis/test_input/"
    
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

print("data read in")

############################## Set stage colours #######################################
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order
stage_cols <- stage_colours[levels(droplevels(seurat@meta.data$stage))]

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

print("sex metadata added")

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
print(DoHeatmap(seurat, features = top10$gene) + NoLegend())
graphics.off()

print("top differentially expressed genes plotted")

#####################################################################################################
#                           Identify and remove contamination (mesoderm and PGCs)                   #
#####################################################################################################

contaminating_plot_path = paste0(plot_path, "contamination/")
dir.create(contaminating_plot_path, recursive = T)

# Make gene list containing markers used to identify contamination clusters
genes <- list('PGC module' = 'DAZL',
              'Blood island module' = c('CDH5', 'TAL1', 'HBZ'),
              'Mesoderm module' = c('CDX2', 'GATA6', 'ALX1', 'PITX2', 'TWIST1', 'TBXT', 'MESP1'),
              'Endoderm module' = c('SOX17', 'CXCR4', 'FOXA2', 'NKX2-2', 'GATA6'))


# Calculate average module expression for contamination gene list
seurat <- AverageGeneModules(seurat_obj = seurat, gene_list = genes)

# Plot distribution of contamination gene modules
png(paste0(contaminating_plot_path, "ContaminationClustersBoxPlot.png"), width = 40, height = 30, units = "cm", res = 200)
PlotCelltype(seurat_obj = seurat, gene_list = genes, quantiles = c(0.1, 0.90), ncol = 2)
graphics.off()

contaminating_clusters <- IdentifyOutliers(seurat_obj = seurat, metrics = names(genes), quantiles = c(0.1, 0.90), intersect_metrics = FALSE)

# Plot UMAP for poor contaminating clusters
png(paste0(contaminating_plot_path, "ContaminationClustersUMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClusterDimplot(seurat, clusters = contaminating_clusters, plot_title = 'Contamination')
graphics.off()

# Plot UMAPs for GOI
# ncol = 4
# png(paste0(contaminating_plot_path, "UMAP_GOI.png"), width = ncol*10, height = 10*ceiling((length(unlist(genes))+1)/ncol), units = "cm", res = 200)
# MultiFeaturePlot(seurat_obj = seurat, gene_list = unlist(genes), plot_clusters = T,
#                  plot_stage = T, label = "", n_col = ncol)
# graphics.off()

# Dotplot for identifying PGCs, Early mesoderm and Late mesoderm
png(paste0(contaminating_plot_path, "dotplot_GOI.png"), width = 30, height = 12, units = "cm", res = 200)
DotPlot(seurat, features = unique(unlist(genes)))
graphics.off()

print("contaminating populations plotted")

############################### Remove contaminating cells from clusters ########################################

filter_cells <- rownames(filter(seurat@meta.data, seurat_clusters %in% contaminating_clusters))
contamination_filt_data <- subset(seurat, cells = filter_cells, invert = T)

# Recluster after removing cells
DefaultAssay(seurat) <- 'peaks'
contamination_filt_data <- RunTFIDF(contamination_filt_data)
contamination_filt_data <- FindTopFeatures(contamination_filt_data, min.cutoff = 'q0')
contamination_filt_data <- RunSVD(contamination_filt_data)
contamination_filt_data <- RunUMAP(object = contamination_filt_data, reduction = 'lsi', dims = 2:30)
contamination_filt_data <- FindNeighbors(object = contamination_filt_data, reduction = 'lsi', dims = 2:30)
contamination_filt_data <- FindClusters(object = contamination_filt_data, verbose = FALSE, algorithm = 3)

png(paste0(contaminating_plot_path, "UMAP.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(object = contamination_filt_data, label = TRUE) + NoLegend()
graphics.off()

# UMAP for clusters and developmental stage
png(paste0(contaminating_plot_path, "ClustStagePlot_UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(contamination_filt_data)
graphics.off()

png(paste0(contaminating_plot_path, "stage_umap.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(contamination_filt_data, group.by = 'stage', label = TRUE, label.size = 12,
        label.box = TRUE, repel = TRUE,
        pt.size = 0.9, cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none",
                 plot.title = element_blank())
graphics.off()

print("contaminating populations filtered")

DefaultAssay(seurat) <- 'RNA'

# Find differentially expressed genes and plot heatmap of top DE genes for each cluster
markers <- FindAllMarkers(contamination_filt_data, only.pos = T, logfc.threshold = 0.25, assay = "RNA")
# get automated cluster order based on percentage of cells in adjacent stages
cluster_order = OrderCellClusters(seurat_object = contamination_filt_data, col_to_sort = 'seurat_clusters', sort_by = 'stage')
# Re-order genes in top15 based on desired cluster order in subsequent plot - this orders them in the heatmap in the correct order
top15 <- markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC) %>% arrange(factor(cluster, levels = cluster_order))
print(top15)

png(paste0(contaminating_plot_path, 'HM.top15.DE.contamination_filt_data.png'), height = 75, width = 100, units = 'cm', res = 500)
TenxPheatmap(data = contamination_filt_data, metadata = c("seurat_clusters", "stage"), custom_order_column = "seurat_clusters",
             custom_order = cluster_order, selected_genes = unique(top15$gene), gaps_col = "seurat_clusters", assay = 'RNA')
graphics.off()


# Plot feature plots for all variable genes
# Set RNA to default assay
DefaultAssay(contamination_filt_data) <- "RNA"

dir.create(paste0(contaminating_plot_path, 'feature_plots/'))
for(i in contamination_filt_data@assays$RNA@var.features){
    png(paste0(plot_path, 'feature_plots/', i, '.png'), height = 12, width = 12, units = 'cm', res = 100)
    print(
      FeaturePlot(contamination_filt_data, features = i, pt.size = 1.4) +
        theme_void() +
        theme(plot.title = element_blank(),
          legend.text = element_text(size=16),
          legend.key.size = unit(1, 'cm'))
        )
    graphics.off()
}

system(paste0("zip -rj ", contaminating_plot_path, "feature_plots.zip ", paste0(contaminating_plot_path, 'feature_plots/')))
unlink(paste0(contaminating_plot_path, 'feature_plots/'), recursive=TRUE, force=TRUE)

# Plot remaining cell counts
filt_counts <- contamination_filt_data@meta.data %>%
  group_by(orig.ident) %>%
  tally() %>%
  rename('filtered' := n)

png(paste0(contaminating_plot_path, 'remaining_cell_table.png'), height = 10, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(filt_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()




saveRDS(contamination_filt_data, paste0(rds_path, "contamination_filt_data.RDS"), compress = FALSE)