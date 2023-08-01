#!/usr/bin/env Rscript

print("Integrate RNA and ATAC seurat objects using harmony and paired nearest neighbour")
# transfer labels new

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(gridExtra)
library(grid)
library(ArchR)
library(Seurat)
library(Signac)
library(scMEGA)
library(harmony)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("", "--verbose"), action = "store", type = "logical", help = "Verbose", default = FALSE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    addArchRThreads(threads = 1)
    
    ## ss8 (for faster testing)
    data_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA/rds_files/"
    rds_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA_integrated/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA_integrated/plots/"
    
    ## full data (real thing)
    # data_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/scMEGA/rds_files/"
    # rds_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/scMEGA_integrated/rds_files/"
    # plot_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/scMEGA_integrated/plots/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/rds_files/"
    rds_path = "./rds_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

set.seed(42)

############################## Read in seurat objects #######################################

print("reading in data...")

obj.rna <- readRDS(paste0(data_path, "./seurat_label_transfer.RDS"))
obj.atac <- readRDS(paste0(data_path, "./ATAC_seurat.RDS"))
gene.activity <- readRDS(paste0(data_path, "./gene_score_matrix.RDS"))

print("data read in!")

############################## Plot UMAPs #######################################

print("plotting UMAPs..")

## set cols
###### schelper cell type colours
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

# set colour palettes for UMAPs
rna_cols <- scHelper_cell_type_colours[as.character(unique(obj.rna$scHelper_cell_type))]
atac_cols <- scHelper_cell_type_colours[as.character(unique(obj.atac$scHelper_cell_type))]

## UMAPs
p1 <- DimPlot(obj.rna, pt.size = 1, reduction = "umap", group.by = "scHelper_cell_type", cols = rna_cols, shuffle = TRUE) +
  ggtitle("scRNA-seq")
p2 <- DimPlot(obj.atac, pt.size = 1, reduction = "umap_iLSI", group.by = "scHelper_cell_type", cols = atac_cols, shuffle = TRUE) +
  ggtitle("snATAC-seq")
png(paste0(plot_path, 'UMAPs_scHelper_cell_type.png'), height = 10, width = 30, units = 'cm', res = 400)
p1 + p2
graphics.off()

stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_cols = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_cols) <- stage_order

p1 <- DimPlot(obj.rna, pt.size = 1, reduction = "umap", group.by = "stage", cols = stage_cols, shuffle = TRUE) +
  ggtitle("scRNA-seq")
p2 <- DimPlot(obj.atac, pt.size = 1, reduction = "umap_iLSI", group.by = "stage", cols = stage_cols, shuffle = TRUE) +
  ggtitle("snATAC-seq")
png(paste0(plot_path, 'UMAPs_stage.png'), height = 10, width = 24, units = 'cm', res = 400)
p1 + p2
graphics.off()

print("UMAPs plotted!")

############################## Coembed RNA and ATAC #######################################

print("Coembedding modalities...")

# set pca dim reduction of rna seurat as 'dr' slot so slots are named the same for atac and rna
obj.rna[["dr"]] <- obj.rna[["pca"]]

## coembedding
obj.coembed <- CoembedData(
  obj.rna,
  obj.atac, 
  gene.activity, 
  weight.reduction = "dr", # 'dr' slot holds the iterative LSI for atac and pca for rna
  verbose = FALSE
)

p1 <- DimPlot(obj.coembed, shuffle = TRUE, label = TRUE, reduction = "umap",
              group.by = "tech", )
p2 <- DimPlot(obj.coembed, shuffle = TRUE, label = TRUE, reduction = "umap",
              group.by = "scHelper_cell_type", cols = atac_cols)
p3 <- DimPlot(obj.coembed, shuffle = TRUE, label = TRUE, reduction = "umap",
              group.by = "stage", cols = stage_cols)

png(paste0(plot_path, 'UMAP_coembed_pre_integration.png'), height = 10, width = 32, units = 'cm', res = 400)
p1 + p2 + p3
graphics.off()

############################## Run harmony and re-plot #######################################

print("running harmony...")

## run batch correction to integrate atac and rna
obj.coembed <- RunHarmony(
  obj.coembed,
  group.by.vars = c("tech"),
  reduction = "pca",
  max.iter.harmony = 30,
  dims.use = 1:30,
  project.dim = FALSE,
  plot_convergence = FALSE
)

## coembedding after batch correction
obj.coembed <- RunUMAP(
  obj.coembed,
  dims = 1:30,
  reduction = 'harmony',
  reduction.name = "umap_harmony",
  reduction.ke = 'umapharmony_',
  verbose = FALSE,
  min.dist = 0.4
)

p1 <- DimPlot(obj.coembed, shuffle = TRUE, label = TRUE, reduction = "umap_harmony",
              group.by = "tech", )
p2 <- DimPlot(obj.coembed, shuffle = TRUE, label = TRUE, reduction = "umap_harmony",
              group.by = "scHelper_cell_type", cols = atac_cols)
p3 <- DimPlot(obj.coembed, shuffle = TRUE, label = TRUE, reduction = "umap_harmony",
              group.by = "stage", cols = stage_cols)

png(paste0(plot_path, 'UMAP_coembed_post_integration.png'), height = 10, width = 32, units = 'cm', res = 400)
p1 + p2 + p3
graphics.off()

p <- DimPlot(obj.coembed, group.by = "scHelper_cell_type", label = FALSE,
             reduction = "umap_harmony", shuffle = TRUE, cols = atac_cols)
png(paste0(plot_path, 'UMAP_coembed_post_integration_cell_type.png'), height = 10, width = 12, units = 'cm', res = 400)
p
graphics.off()

p <- DimPlot(obj.coembed, group.by = "stage", label = FALSE,
             reduction = "umap_harmony", shuffle = TRUE, cols = stage_cols)
png(paste0(plot_path, 'UMAP_coembed_post_integration_stage.png'), height = 10, width = 12, units = 'cm', res = 400)
p
graphics.off()

print("integration run!")

############################## Clustering #######################################

print("clustering...")

## clustering
obj.coembed <- FindNeighbors(obj.coembed, reduction = "harmony", dims = 1:30)
obj.coembed <- FindClusters(obj.coembed, resolution = 0.1, verbose = FALSE)

p <- DimPlot(obj.coembed, group.by = "RNA_snn_res.0.1", label = TRUE,
             reduction = "umap_harmony", shuffle = TRUE) +
  xlab("UMAP1") + ylab("UMAP2")
png(paste0(plot_path, 'UMAP_coembed_post_integration_clustered.png'), height = 10, width = 12, units = 'cm', res = 400)
p
graphics.off()


## proportion of cells in each cluster
p1 <- CellPropPlot(obj.coembed, 
                   group.by = "tech", 
                   prop.in = "RNA_snn_res.0.1")
p2 <- CellPropPlot(obj.coembed, 
                   group.by = "scHelper_cell_type", 
                   prop.in = "RNA_snn_res.0.1",
                   cols = atac_cols)
p3 <- CellPropPlot(obj.coembed, 
                   group.by = "stage", 
                   prop.in = "RNA_snn_res.0.1",
                   cols = stage_cols)
png(paste0(plot_path, 'cluster_proportions.png'), height = 10, width = 32, units = 'cm', res = 400)
p1 + p2 + p3
graphics.off()

# ## find gene markers per cluster and plot as dotplot
# all.markers <- FindAllMarkers(obj.coembed, 
#                               only.pos = TRUE, 
#                               min.pct = 0.5, logfc.threshold = 0.5)
# df <- all.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 3, order_by = avg_log2FC)

# p <- DotPlot(obj.coembed, features = unique(df$gene)) + RotatedAxis()
# png(paste0(plot_path, 'cluster_DotPlot.png'), height = 13, width = 32, units = 'cm', res = 400)
# print(p)
# graphics.off()

# ## UMAP split by modality
# p <- DimPlot(obj.coembed, group.by = "RNA_snn_res.0.1", label = TRUE,
#              reduction = "umap_harmony", shuffle = TRUE, split.by = "tech") +
#   xlab("UMAP1") + ylab("UMAP2")
# png(paste0(plot_path, 'UMAPs_post_integration_clustered.png'), height = 13, width = 22, units = 'cm', res = 400)
# p
# graphics.off()

print("clustering run!")

############################## Create fake multimodal data #######################################

print("pairing cells...")

# pair cells between modalities
df.pair <- PairCells(object = obj.coembed, reduction = "harmony",
                     pair.by = "tech", ident1 = "ATAC", ident2 = "RNA")

# save the cell pairings
write.csv(df.pair, file = paste0(rds_path, "Cell_pairings.csv"), row.names = FALSE)

# only keep paired cells in the seurat object
sel_cells <- c(df.pair$ATAC, df.pair$RNA)
coembed.sub2 <- obj.coembed[, sel_cells]

# see how many cells are left after filtering
cell_counts <- data.frame(dim(obj.coembed)[2], dim(obj.atac)[2], dim(obj.rna)[2], dim(coembed.sub2)[2])
colnames(cell_counts) <- c("Before pairing total", "Before pairing ATAC", "Before pairing RNA", "After pairing total")

png(paste0(plot_path, 'cell_counts_after_pairing.png'), height = 10, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# plot UMAP split by tech
options(repr.plot.height = 5, repr.plot.width = 10)
png(paste0(plot_path, 'UMAPs_post_integration_clustered_split_by_tech.png'), height = 13, width = 22, units = 'cm', res = 400)
DimPlot(coembed.sub2, reduction = "umap_harmony", 
        split.by = "tech")
graphics.off()

## create paired object
obj.pair <- CreatePairedObject(df.pair = df.pair, 
                               object = coembed.sub2,
                               use.assay1 = "RNA", 
                               use.assay2 = "ATAC")

print("cells paired!")

############################## Save data #######################################

## save paired object
saveRDS(obj.pair, paste0(rds_path, "paired_object.RDS"), compress = FALSE)