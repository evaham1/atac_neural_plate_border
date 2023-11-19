#!/usr/bin/env Rscript

print("Script to process seurat object made of metacells")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(parallel)
library(Seurat)
library(dplyr)
library(tibble)
library(scHelper)
library(ggplot2)
library(future)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-i", "--input"), action = "store", type = "character", help = "Name of seurat input file to process", default = "seacells_seurat.RDS"),
  make_option(c("", "--verbose"), action = "store", type = "logical", help = "Verbose", default = TRUE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    
    # Paths for testing ss8 locally
    #data_path = "./local_test_data/convert_seacells_to_seurat/rds_files/"
    #rds_path = "./local_test_data/processed_seurat/rds_files/"
    #plot_path = "./local_test_data/processed_seurat/plots/"
    
    # Paths for testing HH6 on NEMO
    data_path = "./output/NF-downstream_analysis/Processing/HH6/SEACELLS_RNA_WF/3_SEACells_metadata_to_seurat/rds_files/"
    rds_path = "./output/NF-downstream_analysis/Processing/HH6/SEACELLS_RNA_WF/4_Process_metacells/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/HH6/SEACELLS_RNA_WF/4_Process_metacells/plots/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

########################       CELL STATE COLOURS    ########################################
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                        'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                        'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                        'aNP', 'FB', 'vFB', 'node', 'streak', 
                        'PGC', 'BI', 'meso', 'endo')

scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                          "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                          "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                          "#786D73", "#581845", "#9792A3", "#BBB3CB")

names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                 'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                 'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                 'vFB', 'aNP', 'node', 'FB', 'pEpi',
                                 'PGC', 'BI', 'meso', 'endo')

stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order
############################################################################################

#####################################################################################
############################    Read in RDS object   #############################
#####################################################################################

seurat <- readRDS(paste0(data_path, opt$input))
print(seurat)
print(paste0("Number of genes in seurat object: ", length(rownames(seurat))))

DefaultAssay(object = seurat) <- "RNA"
DefaultAssay(object = seurat)

## Plot number of seacells and number of genes
df <- data.frame(dim(seurat))
rownames(df) <- c("Gene count: ", "SEACell cout: ")
png(paste0(plot_path, 'metacell_counts.png'), height = 5, width = 12, units = 'cm', res = 400)
grid.arrange(top=textGrob("Gene count and SEACell count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

########## Check for NA values
DefaultAssay(object = seurat) <- "RNA"
DefaultAssay(object = seurat)
sum(is.na(GetAssayData(object = seurat, slot = "counts"))) #0
sum(is.na(GetAssayData(object = seurat, slot = "data"))) #0
sum(is.na(GetAssayData(object = seurat, slot = "scale.data"))) #<0 x 0 matrix>

DefaultAssay(object = seurat) <- "integrated"
DefaultAssay(object = seurat)
sum(is.na(GetAssayData(object = seurat, slot = "counts"))) #0
sum(is.na(GetAssayData(object = seurat, slot = "data"))) #0
sum(is.na(GetAssayData(object = seurat, slot = "scale.data"))) #<0 x 0 matrix>

#####################################################################################
############################    Re-process 'RNA' slot   #############################
#####################################################################################
### Running everything on the 'RNA' slot only

################### Add potential factors to regress out: MT percent, sex, cell cycle, run #################

DefaultAssay(object = seurat) <- "RNA"
DefaultAssay(object = seurat)

# MT percentage: re-calculate using raw counts
seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt")
print("Added percent mt")
head(seurat@meta.data)

# Cell cycle: re-calculate using raw counts - "Error: Insufficient data values to produce 24 bins."
# as there was no cell cycle effect at ss8 anyway (or any other variable other than run) I'm just going to comment out for now
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# print("Added cell cycle score")
# head(seurat@meta.data)

# Run: proportion of metacells from run 1 already in metadata

# Sex: proportion of metacells from male already in metadata

################### Normalise and find variable features #################

seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, assay = 'RNA')

################### Dim reduction and clustering without regression #################

seurat <- ScaleData(seurat, features = rownames(seurat)) # scaling without regression

# PCA + UMAP
seurat <- RunPCA(object = seurat, verbose = FALSE)

png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
ElbowCutoff(seurat, return = 'plot')
graphics.off()

pc_cutoff <- ElbowCutoff(seurat)
seurat <- FindNeighbors(seurat, dims = 1:pc_cutoff, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:pc_cutoff, verbose = FALSE)

# Cluster using default clustering resolution
seurat <- FindClusters(seurat, resolution = 1)

################### Check if need to regress for sex #################

plot_path = "./plots/investigating_what_to_regress/"
dir.create(plot_path, recursive = T)

# # Scale data and regress
# sex_data <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c("sex"))

# # PCA + UMAP
# sex_data <- RunPCA(object = sex_data, verbose = FALSE)
# pc_cutoff <- ElbowCutoff(sex_data)
# sex_data <- FindNeighbors(sex_data, dims = 1:pc_cutoff, verbose = FALSE)
# sex_data <- RunUMAP(sex_data, dims = 1:pc_cutoff, verbose = FALSE)

# # UMAP before and after regressing out
# png(paste0(plot_path, "sex.png"), width=40, height=20, units = 'cm', res = 200)
# print(gridExtra::grid.arrange(FeaturePlot(seurat, features = "sex", pt.size = 6),
#                               FeaturePlot(sex_data, features = "sex", pt.size = 6),
#                               ncol = 2))
# graphics.off()

### No difference to regress out sex effect - don't do it!

################### Check if need to regress for percent.mt #################

# # Scale data and regress
# mt_data <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c("percent.mt"))

# # PCA + UMAP
# mt_data <- RunPCA(object = mt_data, verbose = FALSE)
# pc_cutoff <- ElbowCutoff(mt_data)
# mt_data <- FindNeighbors(mt_data, dims = 1:pc_cutoff, verbose = FALSE)
# mt_data <- RunUMAP(mt_data, dims = 1:pc_cutoff, verbose = FALSE)

# # UMAP before and after regressing out
# png(paste0(plot_path, "percent.mt.png"), width=40, height=20, units = 'cm', res = 200)
# print(gridExtra::grid.arrange(FeaturePlot(seurat, features = "percent.mt", pt.size = 6),
#                               FeaturePlot(mt_data, features = "percent.mt", pt.size = 6),
#                               ncol = 2))
# graphics.off()

### No difference to regress out percent.mt - don't do it!

################### Check if need to regress for cell cycle #################

# # Scale data and regress
# cc_data <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c("S.Score", "G2M.Score"))

# # PCA + UMAP
# cc_data <- RunPCA(object = cc_data, verbose = FALSE)
# pc_cutoff <- ElbowCutoff(cc_data)
# cc_data <- FindNeighbors(cc_data, dims = 1:pc_cutoff, verbose = FALSE)
# cc_data <- RunUMAP(cc_data, dims = 1:pc_cutoff, verbose = FALSE)

# # UMAP before and after regressing out
# png(paste0(plot_path, "S.Score.png"), width=40, height=20, units = 'cm', res = 200)
# print(gridExtra::grid.arrange(FeaturePlot(seurat, features = "S.Score", pt.size = 6),
#                               FeaturePlot(cc_data, features = "S.Score", pt.size = 6),
#                               ncol = 2))
# graphics.off()

# png(paste0(plot_path, "G2M.Score"), width=40, height=20, units = 'cm', res = 200)
# print(gridExtra::grid.arrange(FeaturePlot(seurat, features = "G2M.Score", pt.size = 6),
#                               FeaturePlot(cc_data, features = "G2M.Score", pt.size = 6),
#                               ncol = 2))
# graphics.off()

### No difference to regress out cell cycle effect - don't do it!

################### Check if need to regress for run #################

# Scale data and regress
run_data <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c("run"))

# PCA + UMAP
run_data <- RunPCA(object = run_data, verbose = FALSE)
pc_cutoff <- ElbowCutoff(run_data)
run_data <- FindNeighbors(run_data, dims = 1:pc_cutoff, verbose = FALSE)
run_data <- RunUMAP(run_data, dims = 1:pc_cutoff, verbose = FALSE)

# UMAP before and after regressing out
png(paste0(plot_path, "run.png"), width=40, height=20, units = 'cm', res = 200)
print(gridExtra::grid.arrange(FeaturePlot(seurat, features = "run", pt.size = 6),
                              FeaturePlot(run_data, features = "run", pt.size = 6),
                              ncol = 2))
graphics.off()

### There is a run effect - cluster to see how bad it is
# Cluster using default clustering resolution
run_data <- FindClusters(run_data, resolution = 1)

# Clusters before and after regressing out
png(paste0(plot_path, "run_clusterst.png"), width=40, height=20, units = 'cm', res = 200)
print(gridExtra::grid.arrange(DimPlot(seurat, group.by = "seurat_clusters", pt.size = 6),
                              DimPlot(run_data, group.by = "seurat_clusters", pt.size = 6),
                              ncol = 2))
graphics.off()

# Plot QC for each cluster - not regressed
png(paste0(plot_path, "run_QCPlot_not_regressed.png"), width=28, height=28, units = 'cm', res = 200)
QCPlot(seurat, plot_quantiles = TRUE, y_elements = "run")
graphics.off()

# Plot QC for each cluster - regressed
png(paste0(plot_path, "run_QCPlot_regressed.png"), width=28, height=28, units = 'cm', res = 200)
QCPlot(run_data, plot_quantiles = TRUE, y_elements = "run")
graphics.off()

# Run does impact clustering so need to regress for this


######################################################################################################
############################    Visualise final seurat object on UMAPs   #############################

plot_path = "./plots/"
dir.create(plot_path, recursive = T)

final_seurat <- run_data

# Find optimal cluster resolution
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = final_seurat, by = 0.2, prefix = "RNA_snn_res.")
graphics.off()

# Decide clustering resolution for annotations here!
final_seurat <- FindClusters(final_seurat, resolution = 1.4)

# Clusters
png(paste0(plot_path, "clusters_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(final_seurat, group.by = 'seurat_clusters', label = TRUE, 
        label.size = ifelse(length(unique(final_seurat$stage)) == 1, 9, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(final_seurat$stage)) == 1, 6, 6), 
        shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

# Size of clusters
df <- as.data.frame(table(final_seurat@meta.data$seurat_clusters))
colnames(df) <- c("Cluster", "nCells")
png(paste0(plot_path, 'cluster_metacell_counts.png'), height = 20, width = 12, units = 'cm', res = 400)
grid.arrange(top=textGrob("", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

# UMAP of stage
final_seurat@meta.data$stage <- factor(final_seurat@meta.data$stage, levels = stage_order)
stage_cols <- stage_colours[levels(droplevels(final_seurat@meta.data$stage))]

png(paste0(plot_path, "stage_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(final_seurat, group.by = 'stage', label = TRUE, 
        label.size = 9, label.box = TRUE, repel = TRUE,
        pt.size = 10, cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

# QC metrics
png(paste0(plot_path, "percent.mt_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
FeaturePlot(object = final_seurat, features = "percent.mt", pt.size = 6)
graphics.off()

png(paste0(plot_path, "run_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
FeaturePlot(final_seurat, features = "run", pt.size = 6)
graphics.off()

png(paste0(plot_path, "sex_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
FeaturePlot(final_seurat, features = "sex", pt.size = 6)
graphics.off()

# png(paste0(plot_path, "cell_cycle_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
# FeaturePlot(final_seurat, features = "S.Score", pt.size = 6)
# graphics.off()

# # QC for each cluster
# png(paste0(plot_path, "cluster_QCPlot.png"), width=28, height=28, units = 'cm', res = 200)
# QCPlot(final_seurat, plot_quantiles = TRUE, y_elements = c("run", "sex", "S.Score", "percent.mt"))
# graphics.off()

# QC for each cluster
png(paste0(plot_path, "cluster_QCPlot.png"), width=28, height=28, units = 'cm', res = 200)
QCPlot(final_seurat, plot_quantiles = TRUE, y_elements = c("run", "sex", "percent.mt"))
graphics.off()

# Proportion-based cell type assignments
#scHelper_cols <- scHelper_cell_type_colours[levels(droplevels(final_seurat@meta.data$scHelper_cell_type_by_proportion))]
scHelper_cols <- scHelper_cell_type_colours[unique(final_seurat@meta.data$scHelper_cell_type_by_proportion)]
png(paste0(plot_path, "proportion_based_cell_types_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(final_seurat, group.by = 'scHelper_cell_type_by_proportion', label = TRUE, 
        label.size = ifelse(length(unique(final_seurat$stage)) == 1, 9, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(final_seurat$stage)) == 1, 6, 6), 
        cols = scHelper_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()


############################    Save seurat object   #############################

## save seacells seurat object
saveRDS(final_seurat, paste0(rds_path, "seacells_seurat_processed.RDS"), compress = FALSE)

# ############################################################################################
# ############################    Re-process 'Integrated' slot   #############################
# ############################################################################################
# 
# seurat <- seurat_temp
# 
# plot_path = "./plots/integrated_assay/"
# dir.create(plot_path, recursive = T)
# 
# DefaultAssay(object = seurat) <- "integrated"
# DefaultAssay(object = seurat)
# 
# ## Normalising and Scaling
# seurat <- NormalizeData(seurat, normalization.method = "RC", scale.factor = 10000)
# seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, assay = "integrated")
# seurat <- ScaleData(seurat, features = rownames(seurat), vars.to.regress = c("percent.mt", "sex", "S.Score", "G2M.Score"))
# 
# seurat
# head(seurat@assays$integrated@scale.data)
# 
# print("integrated data re-processed")
# 
# ## Dim reduction
# seurat <- RunPCA(object = seurat, verbose = FALSE)
# 
# png(paste0(plot_path, "dimHM.png"), width=30, height=65, units = 'cm', res = 200)
# DimHeatmap(seurat, dims = 1:20, balanced = TRUE, cells = 500)
# graphics.off()
# 
# png(paste0(plot_path, "ElbowCutoff.png"), width=30, height=20, units = 'cm', res = 200)
# ElbowCutoff(seurat, return = 'plot')
# graphics.off()
# 
# pc_cutoff <- ElbowCutoff(seurat)
# 
# ## Find neighbours and calculate UMAP
# seurat <- FindNeighbors(seurat, dims = 1:pc_cutoff, verbose = FALSE)
# seurat <- RunUMAP(seurat, dims = 1:pc_cutoff, verbose = FALSE)
# 
# print("dim reduction calculated on integrated data")
# 
# ############################    Visualise on UMAPs   #############################
# 
# # Clusters
# png(paste0(plot_path, "stage_UMAP.png"), width=40, height=20, units = 'cm', res = 200)
# DimPlot(seurat, group.by = "stage", pt.size = 6)
# graphics.off()
# 
# # QC metrics
# png(paste0(plot_path, "percent.mt_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
# FeaturePlot(object = seurat, features = "percent.mt", pt.size = 6)
# graphics.off()
# 
# png(paste0(plot_path, "run_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
# FeaturePlot(object = seurat, features = "run", pt.size = 6)
# graphics.off()
# 
# png(paste0(plot_path, "sex_UMAP.png"), width=10, height=10, units = 'cm', res = 200)
# FeaturePlot(seurat, features = "sex", pt.size = 6)
# graphics.off()
# 
# ############################    Save seurat object   #############################
# 
# ## save seacells seurat object
# saveRDS(seurat, paste0(rds_path, "seacells_seurat_integrated_processed.RDS"), compress = FALSE)
