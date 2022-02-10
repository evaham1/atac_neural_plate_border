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

    plot_path = "../output/NF-downstream_analysis/integrate_rna/plots/"
    rds_path = "../output/NF-downstream_analysis/integrate_rna/rds_files/"
    data_path = "../output/NF-downstream_analysis/test_input/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("future::multisession", workers = ncores)
    options(future.globals.maxSize = 305* 1024^3)
    plan()
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Read in RDS objects and fragment files #######################################

seurat <- readRDS(paste0(data_path, "rds_files/seurat_GeneActivity.RDS"))
DefaultAssay(seurat) <- 'peaks'
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

# read in rna seurat object
seurat_rna <- readRDS(paste0(data_path, "seurat_label_transfer.RDS"))
seurat_rna

############################## Set colours - WILL NEED TO CHANGE TO HH over hh #######################################
stage_order <- c("hh5", "hh6", "hh7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

stage_cols <- stage_colours[levels(droplevels(seurat@meta.data$stage))]

scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                                "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                                "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                       'vFB', 'aNP', 'node', 'FB', 'pEpi')

scHelper_cell_type_cols <- scHelper_cell_type_colours[levels(droplevels(seurat_rna@meta.data$scHelper_cell_type))]

############################## INTEGRATION #######################################

# Stage UMAPs
plot1 <- DimPlot(seurat, group.by = 'stage', label = TRUE, label.size = 12,
        label.box = TRUE, repel = TRUE,
        pt.size = 0.9, cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none") +
  ggtitle('scATAC-seq')
plot2 <- DimPlot(seurat_rna, group.by = 'stage', label = TRUE, label.size = 12,
                 label.box = TRUE, repel = TRUE,
                 pt.size = 0.9, cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none") +
  ggtitle('scRNA-seq')

png(paste0(plot_path, "stage_UMAPs.png"), width=40, height=20, units = 'cm', res = 200)
plot1 + plot2
graphics.off()

# Cluster UMAPs
plot1 <- DimPlot(seurat, group.by = 'seurat_clusters', label = TRUE, label.size = 12,
                 label.box = TRUE, repel = TRUE,
                 pt.size = 0.9, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none") +
  ggtitle('scATAC-seq')
plot2 <- DimPlot(seurat_rna, group.by = 'scHelper_cell_type', label = TRUE, label.size = 12,
                 label.box = TRUE, repel = TRUE,
                 pt.size = 0.9, cols = scHelper_cell_type_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none") +
  ggtitle('scRNA-seq')

png(paste0(plot_path, "cluster_UMAPs.png"), width=40, height=20, units = 'cm', res = 200)
plot1 + plot2
graphics.off()

# Integrate the RNA and ATAC data
DefaultAssay(seurat) <- 'RNA'

transfer.anchors <- FindTransferAnchors(
  reference = seurat_rna,
  query = seurat,
  reduction = 'cca'
)
print("anchors calculated")

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = seurat_rna$scHelper_cell_type,
  weight.reduction = seurat[['lsi']],
  dims = 2:30
)
print("predicted labels")

seurat <- AddMetaData(object = seurat, metadata = predicted.labels)

# Cluster UMAPs
plot1 <- DimPlot(seurat, group.by = 'predicted.id', label = TRUE, label.size = 12,
                 label.box = TRUE, repel = TRUE,
                 pt.size = 0.9, cols = scHelper_cell_type_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none") +
  ggtitle('scATAC-seq')
plot2 <- DimPlot(seurat_rna, group.by = 'scHelper_cell_type', label = TRUE, label.size = 12,
                 label.box = TRUE, repel = TRUE,
                 pt.size = 0.9, cols = scHelper_cell_type_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none") +
  ggtitle('scRNA-seq')

png(paste0(plot_path, "cluster_UMAPs_predicted_ids.png"), width=40, height=20, units = 'cm', res = 200)
plot1 + plot2
graphics.off()

saveRDS(seurat, paste0(rds_path, "seurat_atac_transfer_labels.RDS"))


############################## CO-EMBEDDING #######################################
### only for visualisation, not sure if great this was shown only with multiomic data

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(seurat_rna)
refdata <- GetAssayData(seurat_rna, assay = "RNA", slot = "data")[genes.use, ]


##### whats going on here? how is this different to GeneActivity() ?
# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = seurat[["lsi"]],
                           dims = 2:30)
seurat[["imputation"]] <- imputation

coembed <- merge(x = seurat_rna, y = seurat)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

png(paste0(plot_path, "coembedded_UMAPs.png"), width=40, height=20, units = 'cm', res = 200)
DimPlot(coembed, group.by = c("stage", "scHelper_cell_type"))
graphics.off()


