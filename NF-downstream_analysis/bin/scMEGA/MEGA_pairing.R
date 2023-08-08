#!/usr/bin/env Rscript

print("Pair RNA and ATAC cells in integrated object using previously calculated cell pairings from ArchR integration")
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
options(Seurat.object.assay.version = 'v5')

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
    
    ## ss8
    data_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA/ArchR-to_seurat/" # for ATAC and gene score
    data_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA/ArchR-to_seurat/"
    rds_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA_integrated/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA_integrated/plots/"
    
    ## full data
    # data_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/scMEGA/rds_files/"
    # rds_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/scMEGA_integrated/rds_files/"
    # plot_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/scMEGA_integrated/plots/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/"
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
getArchRThreads()

############################## Read in seurat objects #######################################

print("reading in data...")

# if reading in integrated object
obj.coembed <- readRDS(paste0(data_path, "./rds_files/TransferLabel_integrated_object.RDS"))
print(obj.coembed)

# read in diffusion object
obj.diff <- readRDS(paste0(data_path, "./rds_files/TransferLabel_diffusion_object.RDS"))
print(obj.diff)

# read in cell pairings from archr integration - in input/rds_files folder
df.pair <- read.csv(paste0(data_path, "./rds_files/archr_cell_pairings.csv"))
head(df.pair)

print("data read in!")

############################## Try to get UMAPs to work here! #######################################

# set colours

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
stage_cols = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_cols) <- stage_order
cols <- scHelper_cell_type_colours[as.character(unique(obj.coembed$scHelper_cell_type))]

## comembed 

p1 <- DimPlot(obj.coembed, pt.size = 1, shuffle = TRUE, label = TRUE, reduction = "umap", group.by = "tech")
p2 <- DimPlot(obj.coembed, pt.size = 1, shuffle = TRUE, label = TRUE, reduction = "umap", group.by = "scHelper_cell_type", cols = cols)
p3 <- DimPlot(obj.coembed, pt.size = 1, shuffle = TRUE, label = TRUE, reduction = "umap", group.by = "stage", cols = stage_cols)

png(paste0(plot_path, '1_UMAP_coembed_pre_integration.png'), height = 10, width = 32, units = 'cm', res = 400)
print(p1 + p2 + p3)
graphics.off()

p1 <- DimPlot(obj.coembed, pt.size = 1, shuffle = TRUE, label = FALSE, reduction = "umap", group.by = "tech")
png(paste0(plot_path, '1_UMAP_coembed_pre_integration_tech.png'), height = 10, width = 12, units = 'cm', res = 400)
print(p1)
graphics.off()

p1 <- DimPlot(obj.coembed, pt.size = 1, shuffle = TRUE, label = FALSE, reduction = "umap", group.by = "scHelper_cell_type", cols = cols)
png(paste0(plot_path, '1_UMAP_coembed_pre_integration_celltype.png'), height = 10, width = 12, units = 'cm', res = 400)
print(p1)
graphics.off()

p1 <- DimPlot(obj.coembed, pt.size = 1, shuffle = TRUE, label = FALSE, reduction = "umap", group.by = "stage", cols = stage_cols)
png(paste0(plot_path, '1_UMAP_coembed_pre_integration_stage.png'), height = 10, width = 12, units = 'cm', res = 400)
print(p1)
graphics.off()

png(paste0(plot_path, '1_UMAP_coembed_pre_integration_celltype_split_by_tech.png'), height = 13, width = 22, units = 'cm', res = 400)
DimPlot(obj.coembed, pt.size = 1, reduction = "umap", group.by = "scHelper_cell_type", cols = cols, split.by = "tech")
graphics.off()

png(paste0(plot_path, '1_UMAP_coembed_pre_integration_stage_split_by_tech.png'), height = 13, width = 22, units = 'cm', res = 400)
DimPlot(obj.coembed, pt.size = 1, reduction = "umap", group.by = "stage", cols = stage_cols, split.by = "tech")
graphics.off()

## harmony 

p1 <- DimPlot(obj.coembed, pt.size = 1, shuffle = TRUE, label = TRUE, reduction = "umap_harmony", group.by = "tech")
p2 <- DimPlot(obj.coembed, pt.size = 1, shuffle = TRUE, label = TRUE, reduction = "umap_harmony", group.by = "scHelper_cell_type", cols = cols)
p3 <- DimPlot(obj.coembed, pt.size = 1, shuffle = TRUE, label = TRUE, reduction = "umap_harmony", group.by = "stage", cols = stage_cols)

png(paste0(plot_path, '2_UMAP_coembed_post_integration.png'), height = 10, width = 32, units = 'cm', res = 400)
print(p1 + p2 + p3)
graphics.off()

p <- DimPlot(obj.coembed, pt.size = 1, group.by = "tech", label = FALSE, reduction = "umap_harmony", shuffle = TRUE)
png(paste0(plot_path, '2_UMAP_coembed_post_integration_tech.png'), height = 10, width = 12, units = 'cm', res = 400)
print(p)
graphics.off()

p <- DimPlot(obj.coembed, pt.size = 1, group.by = "scHelper_cell_type", label = FALSE, reduction = "umap_harmony", shuffle = TRUE, cols = cols)
png(paste0(plot_path, '2_UMAP_coembed_post_integration_celltype.png'), height = 10, width = 12, units = 'cm', res = 400)
print(p)
graphics.off()

p <- DimPlot(obj.coembed, pt.size = 1, group.by = "stage", label = FALSE, reduction = "umap_harmony", shuffle = TRUE, cols = stage_cols)
png(paste0(plot_path, '2_UMAP_coembed_post_integration_stage.png'), height = 10, width = 12, units = 'cm', res = 400)
print(p)
graphics.off()

png(paste0(plot_path, '2_UMAP_coembed_post_integration_celltype_split_by_tech.png'), height = 13, width = 22, units = 'cm', res = 400)
DimPlot(obj.coembed, pt.size = 1, reduction = "umap_harmony", group.by = "scHelper_cell_type", cols = cols, split.by = "tech")
graphics.off()

png(paste0(plot_path, '2_UMAP_coembed_post_integration_stage_split_by_tech.png'), height = 13, width = 22, units = 'cm', res = 400)
DimPlot(obj.coembed, pt.size = 1, reduction = "umap_harmony", group.by = "stage", cols = stage_cols, split.by = "tech")
graphics.off()

## diffusion

p1 <- DimPlot(obj.diff, pt.size = 1, shuffle = TRUE, label = TRUE, reduction = "umap_harmony", group.by = "tech")
p2 <- DimPlot(obj.diff, pt.size = 1, shuffle = TRUE, label = TRUE, reduction = "umap_harmony", group.by = "scHelper_cell_type", cols = cols)
p3 <- DimPlot(obj.diff, pt.size = 1, shuffle = TRUE, label = TRUE, reduction = "umap_harmony", group.by = "stage", cols = stage_cols)

png(paste0(plot_path, '3_UMAP_coembed_post_diffusion.png'), height = 10, width = 32, units = 'cm', res = 400)
print(p1 + p2 + p3)
graphics.off()

p <- DimPlot(obj.diff, pt.size = 1, group.by = "tech", label = FALSE, reduction = "dm", shuffle = TRUE)
png(paste0(plot_path, '3_UMAP_coembed_post_diffusion_tech.png'), height = 10, width = 12, units = 'cm', res = 400)
print(p)
graphics.off()

p <- DimPlot(obj.diff, pt.size = 1, group.by = "scHelper_cell_type", label = FALSE, reduction = "dm", shuffle = TRUE, cols = cols)
png(paste0(plot_path, '3_UMAP_coembed_post_diffusion_celltype.png'), height = 10, width = 12, units = 'cm', res = 400)
print(p)
graphics.off()

p <- DimPlot(obj.diff, pt.size = 1, group.by = "stage", label = FALSE, reduction = "dm", shuffle = TRUE, cols = stage_cols)
png(paste0(plot_path, '3_UMAP_coembed_post_diffusion_stage.png'), height = 10, width = 12, units = 'cm', res = 400)
print(p)
graphics.off()

png(paste0(plot_path, '3_UMAP_coembed_post_diffusion_celltype_split_by_tech.png'), height = 13, width = 22, units = 'cm', res = 400)
DimPlot(obj.diff, pt.size = 1, reduction = "dm", group.by = "scHelper_cell_type", cols = cols, split.by = "tech")
graphics.off()

png(paste0(plot_path, '3_UMAP_coembed_post_diffusion_stage_split_by_tech.png'), height = 13, width = 22, units = 'cm', res = 400)
DimPlot(obj.diff, pt.size = 1, reduction = "dm", group.by = "stage", cols = stage_cols, split.by = "tech")
graphics.off()

############################## Clean up cell pairings #######################################

print("pairing cells...")

#### THIS IS HOW scMEGA DOES IT: but it seems to run forever, so instead use pre-computed pairings from ArchR integration
# # pair cells between modalities
# df.pair <- PairCells(object = obj.coembed, reduction = "harmony",
#                      pair.by = "tech", ident1 = "ATAC", ident2 = "RNA")
# 
# # save the cell pairings
# write.csv(df.pair, file = paste0(rds_path, "Cell_pairings.csv"), row.names = FALSE)

# how many unique ATAC and RNA cells left in paired object
print("cell numbers:")
length(unique(df.pair$ATAC)) # 86217
length(unique(df.pair$RNA)) # 1895
dim(df.pair) # 86217     3

# the RNA transfer labels object doesn't include contamination, so 191 RNA cells which are in the df.pair are missing
# so first lets remove these cells from the df.pair object, so there should be no contam cells in the final paired seurat object
RNA_cells_to_remove <- setdiff(unique(df.pair$RNA), Cells(obj.coembed))
length(RNA_cells_to_remove) # 191
filtered_df_pair <- df.pair %>% filter(!RNA %in% RNA_cells_to_remove)
head(filtered_df_pair)
dim(filtered_df_pair)

# ############################## Pair integrated object #######################################

# # how many unique ATAC and RNA cells left in paired object
# print("RNA data")
# length(unique(filtered_df_pair$RNA))
# sum(unique(filtered_df_pair$RNA) %in% Cells(obj.coembed))
# print("ATAC data")
# length(unique(filtered_df_pair$ATAC))
# sum(unique(filtered_df_pair$ATAC) %in% Cells(obj.coembed))

# # only keep paired cells in the seurat object
# sel_cells <- c(filtered_df_pair$ATAC, filtered_df_pair$RNA)
# coembed.sub2 <- obj.coembed[, sel_cells]
# print(coembed.sub2)

# ## create paired object
# obj.pair <- CreatePairedObject(df.pair = filtered_df_pair,
#                                object = coembed.sub2,
#                                use.assay1 = "RNA",
#                                use.assay2 = "ATAC")

# # see how many cells are left in the paired object
# cell_counts <- data.frame(dim(obj.coembed)[2], dim(coembed.sub2)[2], dim(obj.pair)[2])
# colnames(cell_counts) <- c("Before removed contam", "Before paired obj", "After paired obj")

# png(paste0(plot_path, 'Harmony_cell_counts_after_creating_paired_object.png'), height = 10, width = 20, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# print("cells paired!")

# # UMAP
# ###### schelper cell type colours
# scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
#                                 "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
#                                 "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
#                                 "#786D73", "#581845", "#9792A3", "#BBB3CB")
# names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
#                                        'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
#                                        'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
#                                        'vFB', 'aNP', 'node', 'FB', 'pEpi',
#                                        'PGC', 'BI', 'meso', 'endo')
# cols <- scHelper_cell_type_colours[as.character(unique(obj.pair$scHelper_cell_type))]

# p1 <- DimPlot(obj.pair, group.by = "scHelper_cell_type", shuffle = TRUE, label = TRUE, reduction = "umap_harmony", cols = cols)
# png(paste0(plot_path, 'Harmony_UMAP_paired_cell_type.png'), height = 10, width = 14, units = 'cm', res = 400)
# p1
# graphics.off()

# ## save paired object
# label = "TransferLabel"
# saveRDS(obj.pair, paste0(rds_path, label, "integrated_paired_object.RDS"), compress = FALSE)