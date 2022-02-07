#!/usr/bin/env Rscript

### Script to create seurat object, filter data, remove poor quality clusters and integrate across batches

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
      data_path = "../output/NF-downstream_analysis/filtering/rds_files/"
      }else{
      plot_path = "../output/NF-downstream_analysis/gene_activity/plots/"
      rds_path = "../output/NF-downstream_analysis/gene_activity/rds_files/"
      data_path = "../output/NF-downstream_analysis/filtering/rds_files/"}
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 75* 1024^3)
    plan()
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

###########################################################################################################
## Function needed to calculate GEX
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}
###########################################################################################################


############################## Read filtered Seurat RDS object #######################################

seurat <- readRDS(paste0(data_path, "rds_files/seurat_all.RDS"))
print(seurat)

# read in fragment files
paths <- list.dirs(paste0(data_path, "cellranger_atac_output/"), recursive = FALSE, full.names = TRUE)
input <- data.frame(sample = sub('.*/', '', paths), 
                   matrix_path = paste0(paths, "/outs/filtered_peak_bc_matrix.h5"),
                   metadata_path = paste0(paths, "/outs/singlecell.csv"),
                   fragments_path = paste0(paths, "/outs/fragments.tsv.gz"))

fragments_list <- as.list(input$fragments_path)
print(fragments_list)

Fragments(seurat)


######################################## ESTIMATE GEX #####################################################

# test plot
#png(paste0(plot_path, "test_plot.png"), width=20, height=20, units = 'cm', res = 200)
#DimPlot(seurat)
#graphics.off()

####    WILL NEED TO COMBINE FRAGMENT FILES FOR THE SAMPLES USING A BASH SCRIPT AND READ THEM IN HERE?
# https://satijalab.org/signac/0.2/articles/merging.html

# 
# ### calculating gene activity doesnt work for chick data - last line:
# # Error in intI(i, n = x@Dim[1], dn[[1]], give.dn = FALSE) : 
# # 'NA' indices are not (yet?) supported for sparse Matrices
gene.activities <- GeneActivity(seurat)
# 
# ## extract code for GeneActivity function minus error catching to run myself
# annotation <- Annotation(object = signac_filtered)
# transcripts <- CollapseToLongestTranscript(ranges = annotation)
# transcripts <- Extend(x = transcripts, upstream = 2000, 
#                       downstream = 0)
# frags <- Fragments(object = signac_filtered[["peaks"]])
# cells <- colnames(x = signac_filtered[["peaks"]])
# counts <- FeatureMatrix(fragments = frags, features = transcripts, 
#                         cells = cells, verbose = TRUE)
# gene.key <- transcripts$gene_name
# names(x = gene.key) <- GRangesToString(grange = transcripts)
# rownames(x = counts) <- as.vector(x = gene.key[rownames(x = counts)])
# # this line doesnt run: counts <- counts[rownames(x = counts) != "", ]
# NAs <- grepl("^NA\\.", rownames(counts))
# counts_filtered <- counts[!NAs, ]
# counts_filtered <- counts_filtered[-1,]
# gene.activities <- counts_filtered
# 
# ################    NB lots of fragments are not assigned to a gene (NAs)
# #length(rownames(counts)) #24324
# #sum(grepl("^NA\\.", rownames(counts))) #10470
# #length(rownames(counts_pbmc)) #60666
# #sum(grepl("^NA\\.", rownames(counts_pbmc))) #0
# ################
# 
# # add the gene activity matrix to the Seurat object as a new assay and normalize it
# signac_filtered[['RNA']] <- CreateAssayObject(counts = gene.activities)
# 
# 
# ###########################################################################################################
# ######################################## GEX SCALING #####################################################
# 
# signac_rna <- signac_filtered
# 
# signac_rna <- NormalizeData(
#   object = signac_rna,
#   assay = 'RNA',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(signac_rna$nCount_RNA)
# )
# 
# DefaultAssay(signac_rna) <- 'RNA'
# signac_rna <- FindVariableFeatures(signac_rna, selection.method = "vst", nfeatures = 2000)
# 
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(signac_rna), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(signac_rna)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot1 + plot2
# 
# # scaling and clustering
# all.genes <- rownames(signac_rna)
# signac_rna <- ScaleData(signac_rna, features = all.genes)
# signac_rna <- RunPCA(signac_rna, features = VariableFeatures(object = signac_rna))
# signac_rna <- FindNeighbors(signac_rna, dims = 1:10)
# signac_rna <- FindClusters(signac_rna, resolution = 0.5)
# signac_rna <- RunUMAP(signac_rna, dims = 1:10)
# 
# saveRDS(signac_rna, paste0(plot_path, "signac_filtered_GeneActivity.RDS"))
# 
# ###########################################################################################################
# ######################################## SEX FILTERING #####################################################
# 
# signac_rna <- readRDS(paste0(plot_path, "signac_filtered_GeneActivity.RDS"))
# DefaultAssay(signac_rna) <- 'RNA'
# 
# DimPlot(signac_rna, reduction = "umap")
# 
# markers <- FindAllMarkers(signac_rna, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# markers %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC) -> top10
# DoHeatmap(signac_rna, features = top10$gene) + NoLegend()
# 
# 
# #### need to remove MT and sex genes (will later regress properly)
# names_sex_genes_W <- rownames(markers)[grepl("^W", rownames(markers))]
# names_sex_genes_Z <- rownames(markers)[grepl("^Z", rownames(markers))]
# 
# sex_genes <- grepl("^W|^Z", rownames(markers))
# markers_filt <- markers[!sex_genes,]
# markers_filt %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC) -> top10
# DoHeatmap(signac_rna, features = top10$gene) + NoLegend()
# 
# FeaturePlot(
#   object = signac_rna,
#   features = c('W-Wpkci-7', 'W-NIPBLL', 'Z-SEMA6A', 'Z-ARID3C', 'W-Wpkci-7'),
#   pt.size = 0.1,
#   max.cutoff = 'q95',
#   ncol = 3
# )
# 
# # Use W chromosome genes to K-means cluster the cells into male (zz) and female (zw)
# W_genes <- as.matrix(signac_rna@assays$RNA[grepl("W-", rownames(signac_rna@assays$RNA)),])
# k_clusters <- kmeans(t(W_genes), 2)
# k_clusters <- data.frame(k_clusters$cluster)
# signac_rna@meta.data$k_clusters <- k_clusters[match(colnames(signac_rna@assays$RNA), rownames(k_clusters)),]
# 
# # Get rownames for kmeans clusters 1 and 2
# k_clus_1 <- rownames(signac_rna@meta.data[signac_rna@meta.data$k_clusters == 1,])
# k_clus_2 <- rownames(signac_rna@meta.data[signac_rna@meta.data$k_clusters == 2,])
# 
# # K clustering identities are stochastic, so I need to identify which cluster is male and female
# # Sum of W genes is order of magnitude greater in cluster 2 - these are the female cells
# sumclus1 <- sum(W_genes[,k_clus_1])
# sumclus2 <- sum(W_genes[,k_clus_2])
# 
# if(sumclus1 < sumclus2){
#   k_male <- k_clus_1
#   k_female <- k_clus_2
# } else {
#   k_female <- k_clus_1
#   k_male <- k_clus_2
# }
# 
# # Add sex data to meta.data
# signac_rna@meta.data$sex <- unlist(lapply(rownames(signac_rna@meta.data), function(x)
#   if(x %in% k_male){"male"} else if(x %in% k_female){"female"} else{stop("cell sex is not assigned")}))
# 
# 
# # Make dataframe for mean Z expression in male cells
# mean_Z_male <- data.frame(Z.mean = apply(signac_rna@assays$RNA[grepl("Z-", rownames(signac_rna@assays$RNA)), k_male], 1, mean))
# # add 1 before log2 as log2(1) = 0
# mean_Z_male <- log2(mean_Z_male + 1)
# 
# # Make dataframe for mean Z expression in female cells
# mean_Z_female <- data.frame(Z.mean = apply(signac_rna@assays$RNA[grepl("Z-", rownames(signac_rna@assays$RNA)), k_female], 1, mean))
# mean_Z_female <- log2(mean_Z_female + 1)
# 
# # Make dataframe for mean autosomal expression in male cells
# mean_auto_male <- data.frame(auto.mean = apply(signac_rna@assays$RNA[!grepl("Z-", rownames(signac_rna@assays$RNA)) & !grepl("W-", rownames(signac_rna@assays$RNA)), k_male], 1, mean))
# mean_auto_male <- log2(mean_auto_male + 1)
# 
# # Make dataframe for mean autosomal expression in male cells
# mean_auto_female <- data.frame(auto.mean = apply(signac_rna@assays$RNA[!grepl("Z-", rownames(signac_rna@assays$RNA)) & !grepl("W-", rownames(signac_rna@assays$RNA)), k_female], 1, mean))
# mean_auto_female <- log2(mean_auto_female + 1)
# 
# # Calculate FC by subtracting log2 expression from each other
# FC <- list()
# FC$Z <- mean_Z_male - mean_Z_female
# FC$auto <-  mean_auto_male - mean_auto_female
# 
# # Plot boxplot of Z gene and autosomal expression in male vs female cells
# png(paste0(plot_path,"sex_kmeans_log2FC_boxplot.png"), height = 18, width = 18, units = "cm", res = 200)
# boxplot(c(FC$Z, FC$auto),  ylab = "male - female log2 FC (mean normalised UMI +1)", names = c("Z chromosome genes", "autosomal genes"))
# graphics.off()
# 
# signac_rna_sex_filt <- signac_rna
# 
# DefaultAssay(signac_rna_sex_filt) <- "RNA"
# 
# signac_rna_sex_filt <- ScaleData(signac_rna_sex_filt, features = rownames(signac_rna_sex_filt), vars.to.regress = c("sex"), verbose = TRUE)
# 
# signac_rna_sex_filt <- RunPCA(object = signac_rna_sex_filt, verbose = FALSE)
# signac_rna_sex_filt <- FindNeighbors(signac_rna_sex_filt, dims = 1:10, verbose = FALSE)
# signac_rna_sex_filt <- RunUMAP(signac_rna_sex_filt, dims = 1:10, verbose = FALSE)
# signac_rna_sex_filt <- FindClusters(signac_rna_sex_filt, resolution = 0.5, verbose = FALSE)
# 
# saveRDS(signac_rna_sex_filt, paste0(plot_path, "signac_filtered_GeneActivity_sexfilt.RDS"))
# 
# ###########################################################################################################
# ######################################## MT REGRESS #####################################################
# 
# signac_rna_sex_mt_filt <- readRDS(paste0(plot_path, "signac_filtered_GeneActivity_sexfilt.RDS"))
# DefaultAssay(signac_rna_sex_mt_filt) <- 'RNA'
# 
# markers <- FindAllMarkers(signac_rna_sex_mt_filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# names_mt_genes <- rownames(markers)[grepl("^MT", rownames(markers))]
# 
# signac_rna_sex_mt_filt <- PercentageFeatureSet(signac_rna_sex_mt_filt, pattern = "^MT-", col.name = "percent.mt")
# signac_rna_sex_mt_filt <- ScaleData(signac_rna_sex_mt_filt, features = rownames(signac_rna_sex_mt_filt), vars.to.regress = c("percent.mt", "sex"))
# 
# signac_rna_sex_mt_filt <- RunPCA(object = signac_rna_sex_mt_filt, verbose = FALSE)
# signac_rna_sex_mt_filt <- FindNeighbors(signac_rna_sex_mt_filt, dims = 1:10, verbose = FALSE)
# signac_rna_sex_mt_filt <- RunUMAP(signac_rna_sex_mt_filt, dims = 1:10, verbose = FALSE)
# signac_rna_sex_mt_filt <- FindClusters(signac_rna_sex_mt_filt, resolution = 0.5, verbose = FALSE)
# 
# saveRDS(signac_rna_sex_mt_filt, paste0(plot_path, "signac_filtered_GeneActivity_sexmtfilt.RDS"))
# 
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