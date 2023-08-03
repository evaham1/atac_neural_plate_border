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

############################## Read in seurat objects #######################################

# read in atac data - in input/rds_files folder
obj.coembed <- readRDS(paste0(data_path, "rds_files/TEST_OBJECT.RDS"))

############################## Create fake multimodal data #######################################

print("pairing cells...")

print(sessionInfo())

# try this if the paircells fails again
print(obj.coembed)

print("debugging pairing:")
object <- obj.coembed
pair.by <- "tech"
ident1 = "ATAC"
ident2 = "RNA"
reduction = "harmony"

print("line 1:")
obj.1 <- object[, object@meta.data[[pair.by]] == ident1]
print("line 2:")
obj.2 <- object[, object@meta.data[[pair.by]] == ident2]
print("line 3:")
embedding.atac <- Seurat::Embeddings(object = obj.1, reduction = reduction)
print("line 4:")
embedding.rna <- Seurat::Embeddings(object = obj.2, reduction = reduction)
print("line 5:")
embedding <- rbind(embedding.atac, embedding.rna)
print("line 6:")
n.cells <- dim(embedding)[1]

saveRDS(obj.coembed, paste0(rds_path, "TEST_OBJECT.RDS"), compress = FALSE)

# pair cells between modalities
df.pair <- PairCells(object = obj.coembed, reduction = "harmony",
                     pair.by = "tech", ident1 = "ATAC", ident2 = "RNA")

# # save the cell pairings
# write.csv(df.pair, file = paste0(rds_path, "Cell_pairings.csv"), row.names = FALSE)

# # only keep paired cells in the seurat object
# sel_cells <- c(df.pair$ATAC, df.pair$RNA)
# coembed.sub2 <- obj.coembed[, sel_cells]

# # see how many cells are left after filtering
# cell_counts <- data.frame(dim(obj.coembed)[2], dim(obj.atac)[2], dim(obj.rna)[2], dim(coembed.sub2)[2])
# colnames(cell_counts) <- c("Before pairing total", "Before pairing ATAC", "Before pairing RNA", "After pairing total")

# png(paste0(plot_path, 'cell_counts_after_pairing.png'), height = 10, width = 20, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# # plot UMAP split by tech
# options(repr.plot.height = 5, repr.plot.width = 10)
# png(paste0(plot_path, 'UMAPs_post_integration_clustered_split_by_tech.png'), height = 13, width = 22, units = 'cm', res = 400)
# DimPlot(coembed.sub2, reduction = "umap_harmony", 
#         split.by = "tech")
# graphics.off()

# ## create paired object
# obj.pair <- CreatePairedObject(df.pair = df.pair, 
#                                object = coembed.sub2,
#                                use.assay1 = "RNA", 
#                                use.assay2 = "ATAC")

# print("cells paired!")

# ############################## Save data #######################################

# ## save paired object
# saveRDS(obj.pair, paste0(rds_path, label, "_paired_object.RDS"), compress = FALSE)