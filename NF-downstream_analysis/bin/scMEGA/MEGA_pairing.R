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
obj.coembed <- readRDS(paste0(data_path, "./rds_files/TEST_OBJECT.RDS"))
print(obj.coembed)

# read in cell pairings from archr integration - in input/rds_files folder
df.pair <- read.csv(paste0(data_path, "./rds_files/archr_cell_pairings.csv"))
head(df.pair)

print("data read in!")

############################## Create fake multimodal data #######################################

print("pairing cells...")

#### THIS IS HOW scMEGA DOES IT: but it seems to run forever, so instead use pre-computed pairings from ArchR integration
# # pair cells between modalities
# df.pair <- PairCells(object = obj.coembed, reduction = "harmony",
#                      pair.by = "tech", ident1 = "ATAC", ident2 = "RNA")
# 
# # save the cell pairings
# write.csv(df.pair, file = paste0(rds_path, "Cell_pairings.csv"), row.names = FALSE)

# only keep paired cells in the seurat object
sel_cells <- c(df.pair$ATAC, df.pair$RNA)
coembed.sub2 <- obj.coembed[, sel_cells]
print(coembed.sub2)

# how many unique ATAC and RNA cells left in paired object
print("cell numbers:")
length(unique(df.pair$ATAC)) # 86217
length(unique(df.pair$RNA)) # 1895
dim(df.pair) # 86217     3

# # see how many cells are left after filtering
# cell_counts <- data.frame(dim(obj.atac)[2], dim(obj.rna)[2], length(unique(df.pair$RNA)), length(unique(df.pair$ATAC)))
# colnames(cell_counts) <- c("Before pairing ATAC", "Before pairing RNA", "After pairing RNA", "After pairing ATAC")

# png(paste0(plot_path, 'cell_counts_after_pairing.png'), height = 10, width = 20, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# plot UMAP split by tech
options(repr.plot.height = 5, repr.plot.width = 10)
png(paste0(plot_path, 'UMAPs_post_integration_clustered_split_by_tech.png'), height = 13, width = 22, units = 'cm', res = 400)
DimPlot(coembed.sub2, reduction = "umap_harmony", 
        split.by = "tech")
graphics.off()

# the RNA transfer labels object doesn't include contamination, so 191 RNA cells which are in the df.pair are missing
# so first lets remove these cells from the df.pair object, so there should be no contam cells in the final paired seurat object
RNA_cells_to_remove <- setdiff(unique(df.pair$RNA), Cells(coembed.sub2))
length(RNA_cells_to_remove) # 191
filtered_df_pair <- df.pair %>% filter(!RNA %in% RNA_cells_to_remove)
head(filtered_df_pair)
dim(filtered_df_pair)

## debug createpairedobject
object = coembed.sub2
use.assay1 = "RNA"
use.assay2 = "ATAC"

print("RNA data")
length(unique(filtered_df_pair$RNA))
sum(unique(filtered_df_pair$RNA) %in% Cells(object))
print("ATAC data")
length(unique(filtered_df_pair$ATAC))
sum(unique(filtered_df_pair$ATAC) %in% Cells(object))


print("Debugging:")
print("line 1")
rna.counts <- GetAssayData(object, assay = use.assay1, slot = "counts")[, df.pair$RNA]
print("line 2")
atac.counts <- GetAssayData(object, assay = use.assay2, slot = "counts")[, df.pair$ATAC]
print("line 3")
colnames(rna.counts) <- df.pair$cell_name
print("line 4")
colnames(atac.counts) <- df.pair$cell_name
print("line 5")
obj.pair <- CreateSeuratObject(counts = rna.counts, assay = use.assay1)
print("line 6")
obj.pair[[use.assay2]] <- CreateChromatinAssay(counts = atac.counts, sep = sep, min.cells = 10)
print("line 7")
for (reduction in names(object@reductions)) {
        embedding <- Embeddings(object, reduction = reduction)[df.pair$RNA, 
            ]
        rownames(embedding) <- df.pair$cell_name
        obj.pair[[reduction]] <- CreateDimReducObject(embeddings = embedding, 
            assay = DefaultAssay(obj.pair))
    }
meta.data <- object@meta.data[df.pair$RNA, ]
rownames(meta.data) <- df.pair$cell_name
obj.pair <- AddMetaData(obj.pair, metadata = meta.data)

# ## create paired object
# obj.pair <- CreatePairedObject(df.pair = df.pair,
#                                object = coembed.sub2,
#                                use.assay1 = "RNA",
#                                use.assay2 = "ATAC")

print("cells paired!")

# UMAP
p1 <- DimPlot(obj.pair, group.by = "scHelper_cell_type", shuffle = TRUE, label = TRUE, reduction = "umap_harmony", cols = atac_cols)

png(paste0(plot_path, 'UMAP_paired_cell_type.png'), height = 10, width = 14, units = 'cm', res = 400)
p1
graphics.off()

# see how many cells are left in the paired object
cell_counts <- data.frame(dim(coembed.sub2)[2], dim(obj.pair)[2])
colnames(cell_counts) <- c("Before paired obj", "After paired obj")

png(paste0(plot_path, 'cell_counts_after_creating_paired_object.png'), height = 10, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

############################## Save data #######################################

## save paired object
saveRDS(obj.pair, paste0(rds_path, label, "_paired_object.RDS"), compress = FALSE)