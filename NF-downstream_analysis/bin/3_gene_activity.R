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
    options(future.globals.maxSize = 395* 1024^3)
    plan()
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


############################## Read in Seurat RDS object and fragment files #######################################

seurat <- readRDS(paste0(data_path, "rds_files/seurat_all_filtered.RDS"))
print(seurat)

# read in fragment files
paths <- list.dirs(paste0(data_path, "cellranger_atac_output/"), recursive = FALSE, full.names = TRUE)
input <- data.frame(sample = sub('.*/', '', paths), 
                   matrix_path = paste0(paths, "/outs/filtered_peak_bc_matrix.h5"),
                   metadata_path = paste0(paths, "/outs/singlecell.csv"),
                   fragments_path = paste0(paths, "/outs/fragments.tsv.gz"))
new.paths <- as.list(input$fragments_path)
frags <- Fragments(seurat)  # get list of fragment objects
Fragments(seurat) <- NULL  # remove fragment information from assay

for (i in seq_along(frags)) {
  frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
}
Fragments(seurat) <- frags # assign updated list back to the object

############################## Nucleosome Banding Plot after filtering #######################################

# need to downsample first and then just using first 100bps of chromosome one as takes a long time
png(paste0(plot_path, "UMAP.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(object = seurat, label = TRUE) + NoLegend()
graphics.off()

print("about to subset")
seurat_small <- subset(seurat, downsample = 5)
print("seurat subsetted")
seurat_small

png(paste0(plot_path, "small_UMAP.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(object = seurat_small, label = TRUE) + NoLegend()
graphics.off()

png(paste0(plot_path, 'QC_Nucleosome_banding.png'), height = 15, width = 21, units = 'cm', res = 400)
FragmentHistogram(object = seurat_small, region = "chr1-1-100", group.by = 'stage')
graphics.off()

######################################## ESTIMATE GEX #####################################################

gene.activities <- GeneActivity(seurat)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
seurat[['RNA']] <- CreateAssayObject(counts = gene.activities)
seurat <- NormalizeData(
  object = seurat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat$nCount_RNA)
)

DefaultAssay(seurat) <- 'RNA'

png(paste0(plot_path, 'TestFeaturePlot.png'), height = 15, width = 34, units = 'cm', res = 400)
FeaturePlot(
  object = seurat,
  features = c('SIX1', 'EYA2', 'DLX5', 'SOX2', 'PAX7', 'TFAP2A'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
graphics.off()


###########################################################################################################
######################################## Top variable genes ###############################################

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
 
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat), 10)
 
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png(paste0(plot_path, 'MostVariable.png'), height = 15, width = 20, units = 'cm', res = 400)
plot2
graphics.off()

markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

png(paste0(plot_path, 'top_markers.png'), height = 10, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Top Cluster Markers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(top10, rows=NULL, theme = ttheme_minimal()))
graphics.off()

###########################################################################################################
######################################## SEX FILTERING #####################################################

png(paste0(plot_path, 'FeaturePlot_sex_genes.png'), height = 20, width = 30, units = 'cm', res = 400)
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

png(paste0(plot_path,"UMAP_sex.png"), height = 18, width = 18, units = "cm", res = 200)
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
png(paste0(plot_path,"sex_kmeans_log2FC_boxplot.png"), height = 18, width = 18, units = "cm", res = 200)
boxplot(c(FC$Z, FC$auto),  ylab = "male - female log2 FC (mean normalised UMI +1)", names = c("Z chromosome genes", "autosomal genes"))
graphics.off()

saveRDS(seurat, paste0(rds_path, "seurat_GeneActivity.RDS"))

# DefaultAssay(seurat) <- 'peaks'
# seurat_sex_filt <- ScaleData(seurat, vars.to.regress = c("sex"), verbose = TRUE)

# seurat_sex_filt <- RunSVD(seurat_sex_filt)
# seurat_sex_filt <- RunUMAP(object = seurat_sex_filt, reduction = 'lsi', dims = 2:30)
# seurat_sex_filt <- FindNeighbors(object = seurat_sex_filt, reduction = 'lsi', dims = 2:30)
# seurat_sex_filt <- FindClusters(object = seurat_sex_filt, verbose = FALSE, algorithm = 3)

# png(paste0(plot_path,"UMAP_sex_regressed.png"), height = 18, width = 18, units = "cm", res = 200)
# DimPlot(seurat_sex_filt, group.by = "sex")
# graphics.off()

# 


# ###########. TESTING
# seurat_small <- subset(x = seurat_sex_filt, subset = stage == "hh6")
# seurat_small <- subset(x = seurat_small, downsample = 1000)
# DimPlot(seurat_small)
# DimPlot(seurat_small, group.by = "sex")

# ### try to do dim reduction on atac data that has been scaled 
# # and regressed with 'sex' from metadata
# seurat_small <- RunSVD(seurat_small, scale.embeddings = FALSE)
# seurat_small <- RunUMAP(object = seurat_small, reduction = 'lsi', dims = 2:30)
# seurat_small <- FindNeighbors(object = seurat_small, reduction = 'lsi', dims = 2:30)
# seurat_small <- FindClusters(object = seurat_small, verbose = FALSE, algorithm = 3)

# DimPlot(seurat_small)
# DimPlot(seurat_small, group.by = "sex")
# DimPlot(seurat_small, group.by = "stage")

# ## dont seem to be able to run RunSVD on the scaled data rather than the TF-IDF, 
# # can only specify the assay (RNA or peaks)


# # try just running scaling and dim reduction on the imputed RNA
# DefaultAssay(seurat_small) <- 'RNA'
# seurat_small <- FindVariableFeatures(seurat_small, selection.method = "vst", nfeatures = 2000)
# seurat_small <- ScaleData(seurat_small, features = rownames(seurat_small), vars.to.regress = c("sex"), verbose = TRUE)

# seurat_small <- RunPCA(object = seurat_small, verbose = TRUE)
# print(ElbowCutoff(seurat_small, return = 'plot'))
# pc_cutoff <- ElbowCutoff(seurat_small)

# seurat_small <- FindNeighbors(seurat_small, dims = 1:pc_cutoff, verbose = TRUE)
# seurat_small <- FindClusters(seurat_small, resolution = 0.5, verbose = TRUE)

# seurat_small <- RunUMAP(seurat_small, dims = 1:pc_cutoff, verbose = TRUE)

# DimPlot(seurat_small)
# DimPlot(seurat_small, group.by = "sex")

# DefaultAssay(seurat_small) <- 'peaks'
# DimPlot(seurat_small)
# DimPlot(seurat_small, group.by = "sex")





###########################################################################################################
######################################## MT REGRESS #####################################################

# signac_rna_sex_mt_filt <- readRDS(paste0(plot_path, "signac_filtered_GeneActivity_sexfilt.RDS"))
# DefaultAssay(signac_rna_sex_mt_filt) <- 'RNA'

# markers <- FindAllMarkers(signac_rna_sex_mt_filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# names_mt_genes <- rownames(markers)[grepl("^MT", rownames(markers))]

# signac_rna_sex_mt_filt <- PercentageFeatureSet(signac_rna_sex_mt_filt, pattern = "^MT-", col.name = "percent.mt")
# signac_rna_sex_mt_filt <- ScaleData(signac_rna_sex_mt_filt, features = rownames(signac_rna_sex_mt_filt), vars.to.regress = c("percent.mt", "sex"))

# signac_rna_sex_mt_filt <- RunPCA(object = signac_rna_sex_mt_filt, verbose = FALSE)
# signac_rna_sex_mt_filt <- FindNeighbors(signac_rna_sex_mt_filt, dims = 1:10, verbose = FALSE)
# signac_rna_sex_mt_filt <- RunUMAP(signac_rna_sex_mt_filt, dims = 1:10, verbose = FALSE)
# signac_rna_sex_mt_filt <- FindClusters(signac_rna_sex_mt_filt, resolution = 0.5, verbose = FALSE)

# saveRDS(signac_rna_sex_mt_filt, paste0(plot_path, "signac_filtered_GeneActivity_sexmtfilt.RDS"))

# ###########################################################################################################
# ######################################## VISUALISATIONS #####################################################
# 
# signac <- readRDS(paste0(plot_path, "signac_filtered_GeneActivity_sexmtfilt.RDS"))
# DefaultAssay(signac) <- 'RNA'
# 
# DimPlot(signac_rna_sex_mt_filt, reduction = "umap")
# 
# markers <- FindAllMarkers(signac_rna_sex_mt_filt, only.pos = TRUE, logfc.threshold = 0.15)
# markers_0 <- markers[which(markers$cluster == 0),]
# markers_1 <- markers[which(markers$cluster == 1),]
# markers_2 <- markers[which(markers$cluster == 2),]
# markers_3 <- markers[which(markers$cluster == 3),]
# DoHeatmap(signac_rna_sex_mt_filt, features = markers_0$gene)
# markers %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC) -> top10
# DoHeatmap(signac_rna_sex_mt_filt, features = top10$gene) + NoLegend()
# 
# FeaturePlot(
#   object = signac_rna_sex_mt_filt,
#   features = c('W-Wpkci-7', 'W-NIPBLL', 'Z-SEMA6A', 'Z-ARID3C', 'W-Wpkci-7'),
#   pt.size = 0.1,
#   max.cutoff = 'q95',
#   ncol = 3
# )
# 
# top10_0 <- top10[1:10, ]
# top10_1 <- top10[11:12, ]
# top10_2 <- top10[13:18, ]
# top10_4 <- top10[19:28, ]
# 
# FeaturePlot(
#   object = signac_rna_sex_mt_filt,
#   features = top10_4$gene,
#   pt.size = 0.1,
#   max.cutoff = 'q95',
#   ncol = 3
# )
# 
# FeaturePlot(
#   object = signac_rna_sex_mt_filt,
#   features = c('EYA2', 'SIX1', 'PAX7', 'SOX2', 'SOX10'),
#   pt.size = 0.1,
#   max.cutoff = 'q95',
#   ncol = 3
# )