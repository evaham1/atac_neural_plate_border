#!/usr/bin/env Rscript

print("Combine the RNA and ATAC cells using the single cell mappings from ArchR integration")

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

# if reading in transfer labels object
obj.rna <- readRDS(paste0(data_path, "seurat_label_transfer_minus_HH4.RDS"))

# read in atac data - in input/rds_files folder
obj.atac <- readRDS(paste0(data_path, "rds_files/ATAC_seurat.RDS"))

# read in cell pairings from archr integration - in input/rds_files folder
df.pair <- read.csv(paste0(data_path, "./rds_files/archr_cell_pairings.csv"))
head(df.pair)

# read in gene activity matrix - in input/rds_files folder
gene.activity <- readRDS(paste0(data_path, "rds_files/gene_score_matrix.RDS"))

print("data read in!")

############################## Plot UMAPs #######################################

print("plotting UMAPs..")

## set cols
# schelper cell type colours
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                "#A5718D", "#3F918C", "#ed5e5f", "9792A3")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                       'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo',
                                       'Neural', 'Placodal', 'Non-neural', 'Contam')

# stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_cols = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_cols) <- stage_order

# set colour palettes for UMAPs
rna_cols <- scHelper_cell_type_colours[as.character(unique(obj.rna$scHelper_cell_type))]
atac_cols <- scHelper_cell_type_colours[as.character(unique(obj.atac$scHelper_cell_type))]

## UMAPs
p1 <- DimPlot(obj.rna, pt.size = 1, reduction = "umap", group.by = "scHelper_cell_type", cols = rna_cols, shuffle = TRUE) +
  ggtitle("scRNA-seq")
p2 <- DimPlot(obj.atac, pt.size = 1, reduction = "umap_iLSI", group.by = "scHelper_cell_type", cols = atac_cols, shuffle = TRUE) +
  ggtitle("snATAC-seq")
png(paste0(plot_path, '0_UMAPs_scHelper_cell_type.png'), height = 10, width = 30, units = 'cm', res = 400)
p1 + p2
graphics.off()

p1 <- DimPlot(obj.rna, pt.size = 1, reduction = "umap", group.by = "stage", cols = stage_cols, shuffle = TRUE) +
  ggtitle("scRNA-seq")
p2 <- DimPlot(obj.atac, pt.size = 1, reduction = "umap_iLSI", group.by = "stage", cols = stage_cols, shuffle = TRUE) +
  ggtitle("snATAC-seq")
png(paste0(plot_path, '0_UMAPs_stage.png'), height = 10, width = 24, units = 'cm', res = 400)
p1 + p2
graphics.off()

print("UMAPs plotted!")

############################## Combine RNA and ATAC #######################################

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

print("Coembedding complete!")
print(obj.coembed)

## save paired object
saveRDS(obj.coembed, paste0(rds_path, "TransferLabel_comembed_object.RDS"), compress = FALSE)

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

############################## Pair object #######################################

# how many unique ATAC and RNA cells left in paired object
print("RNA data")
length(unique(filtered_df_pair$RNA))
sum(unique(filtered_df_pair$RNA) %in% Cells(obj.coembed))
print("ATAC data")
length(unique(filtered_df_pair$ATAC))
sum(unique(filtered_df_pair$ATAC) %in% Cells(obj.coembed))

# only keep paired cells in the seurat object
sel_cells <- c(filtered_df_pair$ATAC, filtered_df_pair$RNA)
coembed.sub2 <- obj.coembed[, sel_cells]
print(coembed.sub2)

############################## Create multiome object #######################################

## create paired object
obj.pair <- CreatePairedObject(df.pair = df.pair,
                               object = filtered_df_pair,
                               use.assay1 = "RNA",
                               use.assay2 = "ATAC")

# see how many cells are left in the paired object
cell_counts <- data.frame(dim(obj.coembed)[2], dim(coembed.sub2)[2], dim(obj.pair)[2])
colnames(cell_counts) <- c("Before removed contam", "Before paired obj", "After paired obj")

png(paste0(plot_path, 'Cell_counts_after_creating_paired_object.png'), height = 10, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

print("cells paired!")

############################## Go back to iLSI dim red #######################################

# change cell names to ATAC cell ids
obj.pair <- RenameCells(obj.pair, new.names = df.pair$ATAC, old.names = df.pair$cell_name)

# lift over dim reduction from atac object onto paired object
obj.pair[["iLSI"]] <- obj.atac[["iLSI"]]

# run UMAP on iLSI
obj.pair <- RunUMAP(obj.pair, 
                    dims = 1:30, 
                    reduction = 'iLSI',
                    reduction.name = "umap_iLSI",
                    reduction.ke = 'umapiLSI_',
                    verbose = FALSE,
                    min.dist = 0.4)

print("UMAP recalculated on iLSI!")

############################## UMAPs #######################################

###### plot scHelepr cell type and stage
p1 <- DimPlot(obj.pair, pt.size = 2, group.by = "scHelper_cell_type", shuffle = TRUE, label = FALSE, reduction = "umap_iLSI", cols = cols)
p2 <- DimPlot(obj.pair, pt.size = 2, group.by = "stage", shuffle = TRUE, label = FALSE, reduction = "umap_iLSI", cols = stage_cols)
png(paste0(plot_path, 'Paired_UMAPs.png'), height = 10, width = 24, units = 'cm', res = 400)
p1 + p2
graphics.off()

###### add broad cell types and plot
scHelper_cell_types <- as.data.frame(obj.pair$scHelper_cell_type)

broad <- scHelper_cell_types %>% mutate(broad = mapvalues(obj.pair$scHelper_cell_type, 
                                                          from=c("NP", "aNP", "iNP", "pNP", "eN", "vFB", "FB", "MB", "HB", "eCN", "eN",
                                                                 'PPR', 'aPPR', 'pPPR',
                                                                 'eNPB', 'NPB', 'aNPB', 'pNPB',
                                                                 'NC', 'dNC',
                                                                 'NNE', 'pEpi',
                                                                 'EE', 'meso', 'endo', 'BI', 'PGC'),
                                                          to=c(
                                                            rep("Neural", 11),
                                                            rep("Placodal", 3),
                                                            rep("NPB", 4),
                                                            rep("NC", 2),
                                                            rep("Non-neural", 2),
                                                            rep("Contam", 5)
                                                          )
))
obj.pair$scHelper_cell_type_broad <- broad$broad
print("Broad scHelper cell type labels added")

cols <- scHelper_cell_type_colours[as.character(unique(obj.pair$scHelper_cell_type_broad))]

p1 <- DimPlot(obj.pair, pt.size = 2, group.by = "scHelper_cell_type_broad", shuffle = TRUE, label = FALSE, reduction = "umap_iLSI", cols = cols)
p2 <- DimPlot(obj.pair, pt.size = 2, group.by = "stage", shuffle = TRUE, label = FALSE, reduction = "umap_iLSI", cols = stage_cols)
png(paste0(plot_path, 'Paired_UMAPs_broad.png'), height = 10, width = 24, units = 'cm', res = 400)
p1 + p2
graphics.off()

############################## Save #######################################

## save paired object
saveRDS(obj.pair, paste0(rds_path, "TransferLabel_paired_object.RDS"), compress = FALSE)