#!/usr/bin/env Rscript

print("Convert from ArchR object into seurat object")
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
    
    #data_path = "./output/NF-downstream_analysis/Processing/ss8/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/"
    #rds_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA/rds_files/"
    data_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/rds_files/TransferLabels_Save-ArchR/"
    rds_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/scMEGA/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/scMEGA/plots/"
    
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
getArchRThreads()

############################## Read in ArchR project #######################################

# If files are not in rds_files subdirectory look in input dir 
label <- unique(sub('_.*', '', list.files(data_path)))
print(label) 

if (length(label) == 0){
  data_path = "./input/"
  label <- sub('_.*', '', list.files(data_path))
  print(label)
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
} else {
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
}

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getPeakSet(ArchR)
getAvailableMatrices(ArchR)

############################## Extract data #######################################

# extract peak count matrix
peak_counts <- getMatrixFromProject(ArchR, useMatrix = "PeakMatrix")
count_data <- assay(peak_counts)
ranges <- rowRanges(peak_counts)$name
rownames(count_data) <- ranges
count_data[1:3, 1:3]

# extract cell metadata
cell_metadata <- as.data.frame(ArchR@cellColData)
head(cell_metadata)

############################## Make Seurat obj #######################################

print("Making chromatin assay...")

# generate chromatin assay object
chrom_assay <- CreateChromatinAssay(
  counts = count_data,
  ranges = rowRanges(peak_counts),
  sep = c("_", "_"),
  min.cells = 100
)

print("Making seurat object...")

# generate ATAC seurat object
obj.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC",
  meta.data = cell_metadata,
  names.field = 1, 
  names.delim = "#")

print("Seurat object made!")

obj.atac

print("Adding dimension reduced matrix...")

# add dimension reduced matrix - iterative LSI
reduced_dims <- ArchR@reducedDims$IterativeLSI$matSVD
print(head(reduced_dims))

obj.atac[["iLSI"]] <- CreateDimReducObject(embeddings = reduced_dims,
                                              assay = DefaultAssay(obj.atac),
                                              key = "iLSI_")
obj.atac[["dr"]] <- CreateDimReducObject(embeddings = reduced_dims,
                                           assay = DefaultAssay(obj.atac),
                                           key = "dr_")


############################## Run UMAPs #######################################

print("Running UMAP...")

# run UMAP
obj.atac <- RunUMAP(obj.atac, 
                    dims = 1:30, 
                    reduction = 'iLSI',
                    reduction.name = "umap_iLSI",
                    reduction.ke = 'umapiLSI_',
                    verbose = FALSE,
                    min.dist = 0.4)

###### schelper cell type colours
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                                "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                                "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                       'vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo')
cols <- scHelper_cell_type_colours[as.character(unique(obj.atac$scHelper_cell_type))]

# plot UMAP
p1 <- DimPlot(obj.atac, group.by = "clusters", pt.size = 1,
              reduction = "umap_iLSI")
p2 <- DimPlot(obj.atac, group.by = "scHelper_cell_type", pt.size = 1,
              reduction = "umap_iLSI", cols = cols)
png(paste0(plot_path, 'UMAPs_stages.png'), height = 10, width = 20, units = 'cm', res = 400)
p1 + p2
graphics.off()

############################## Save seurat object #######################################

print("Saving seurat object...")

# save seurat object
saveRDS(obj.atac, paste0(rds_path, "ATAC_seurat.RDS"), compress = FALSE)

############################## Extract and save gene score matrix #######################################

print("Extracting gene score matrix...") 

### Now need to extract gene score matrix and save as rds
gene_counts <- getMatrixFromProject(ArchR, useMatrix = "GeneScoreMatrix")
gene_df <- assay(gene_counts)
gene_names <- rowData(gene_counts)$name
rownames(gene_df) <- gene_names
gene_df[1:3, 1:3]

saveRDS(gene_df, paste0(rds_path, "gene_score_matrix.RDS"), compress = FALSE)

############################## Extract and save cell pairings #######################################

print("Extracting ArchR cell pairings...")

# extract cell ids of paired ATAC and RNA cells from the ArchR single cell integration
archr_cell_pairings <- as.data.frame(getCellColData(ArchR, select = c("predictedCell", "predictedScore")))
archr_cell_pairings <- rownames_to_column(archr_cell_pairings)
# optionally filter here based on score?
archr_cell_pairings <- archr_cell_pairings[,-3]
colnames(archr_cell_pairings) <- c("ATAC", "RNA")
# add new cell ids
archr_cell_pairings <- archr_cell_pairings %>%
  mutate(cell_name = paste0("cell_", row_number()))

head(archr_cell_pairings)

write.csv(archr_cell_pairings, paste0(rds_path, "archr_cell_pairings.csv"), row.names = FALSE, col.names = TRUE)