#!/usr/bin/env Rscript

print("ArchR integration")

############################## Load libraries #######################################
library(getopt)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)
library(presto)
library(Seurat)

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
    
    ncores = 8
    
    plot_path = "./output/NF-downstream_analysis/ArchR_integration/plots/"
    rds_path = "./output/NF-downstream_analysis/ArchR_integration/rds_files/"
    #data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/7_ArchR_clustering_postfiltering_twice/rds_files/"
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ss8/ArchR_clustering//rds_files/"
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Read in ArchR project and seurat object #######################################

# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label[1], "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# load seurat object by reading in any rds object
rna_path <- list.files(path = data_path, pattern = "*.RDS", full.names = TRUE)
seurat_data <- readRDS(rna_path)


############################################################################################
############################## Pre-Integration Plots #######################################

#### Prepare RNA labels and colours #######

#### combine contamination and old scHelper_cell_state labels
contam <- seurat_data@meta.data$contamination
original <- seurat_data@meta.data$scHelper_cell_type_original
old_labels <- contam
old_labels[is.na(contam)] <- as.character(original[is.na(contam)])
seurat_data@meta.data$old_labels <- old_labels

###### stage colours
stage_order <- c("HH4", "HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#E78AC3", "#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order
stage_cols <- stage_colours[levels(droplevels(seurat_data@meta.data$stage))]

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
seurat_data@meta.data$old_labels <- factor(seurat_data@meta.data$old_labels, levels = scHelper_cell_type_order)
scHelper_new_cols <- scHelper_cell_type_colours[levels(droplevels(seurat_data@meta.data$scHelper_cell_type))]
scHelper_old_cols <- scHelper_cell_type_colours[levels(droplevels(seurat_data@meta.data$old_labels))]

###### UMAPs
umap_rna_new <- DimPlot(seurat_data, group.by = 'scHelper_cell_type', label = TRUE, 
                    label.size = ifelse(length(unique(seurat_data$stage)) == 1, 9, 3),
                    label.box = TRUE, repel = TRUE,
                    pt.size = ifelse(length(unique(seurat_data$stage)) == 1, 1.2, 1), 
                    cols = scHelper_new_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
umap_rna_old <- DimPlot(seurat_data, group.by = 'old_labels', label = TRUE, 
                    label.size = ifelse(length(unique(seurat_data$stage)) == 1, 9, 3),
                    label.box = TRUE, repel = TRUE,
                    pt.size = ifelse(length(unique(seurat_data$stage)) == 1, 1.2, 1), 
                    cols = scHelper_old_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())

png(paste0(plot_path, 'UMAPs_old_new_rna.png'), height = 20, width = 40, units = 'cm', res = 400)
print(umap_rna_new + umap_rna_old)
graphics.off()

############################## UMAPs before integration #######################################
# UMAPs of RNA and ATAC data, with RNA coloured by cell state and ATAC by clusters
umap_atac <- plotEmbedding(ArchR, name = "clusters", plotAs = "points", size = 1.8, baseSize = 0, labelSize = 8, legendSize = 0)

png(paste0(plot_path, 'UMAPs_before_integration_new_scHelper_cell_states.png'), height = 20, width = 40, units = 'cm', res = 400)
print(umap_rna_new + umap_atac)
graphics.off()

png(paste0(plot_path, 'UMAPs_before_integration_old_scHelper_cell_states.png'), height = 20, width = 40, units = 'cm', res = 400)
print(umap_rna_old + umap_atac)
graphics.off()


################################################################################################
############################## Unconstrained integration #######################################

ArchR <- addGeneIntegrationMatrix(
  ArchRProj = ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seurat_data,
  addToArrow = TRUE,
  groupRNA = "scHelper_cell_type",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  force = TRUE,
  threads = 1
)
print("integration completed")

# use matched RNA cells to add new and old labels to ATAC cells
extracted_rna_labels <- seurat_data@meta.data[ArchR$predictedCell_Un, c("scHelper_cell_type", "old_labels")]
ArchR$scHelper_cell_type_new <- extracted_rna_labels[, "scHelper_cell_type"]
ArchR$scHelper_cell_type_old <- extracted_rna_labels[, "old_labels"]

# save integrated ArchR project
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label[1], "_Save-ArchR"), load = FALSE)


#############################################################################################
############################## Post-Integration Plots #######################################

# set colour palettes for UMAPs
atac_scHelper_new_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_new)]
atac_scHelper_old_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_old)]

############################## RNA cell labels on ATAC data #######################################

### New labels
png(paste0(plot_path, 'UMAP_unconInt_scHelper_cell_type_new.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type_new", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_new_cols, labelAsFactors = FALSE)
graphics.off()

png(paste0(plot_path, 'UMAP_unconInt_scHelper_cell_type_new_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type_new", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_new_cols)
graphics.off()

### Old labels 
png(paste0(plot_path, 'UMAP_unconInt_scHelper_cell_type_old.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE)
graphics.off()

png(paste0(plot_path, 'UMAP_unconInt_scHelper_cell_type_old_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_old_cols)
graphics.off()

############################## Assign cluster labels to ATAC data #######################################

### New labels
cM <- confusionMatrix(paste0(ArchR$clusters), paste0(ArchR$scHelper_cell_type_new))
  cM <- cM / Matrix::rowSums(cM)
  p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
  )
png(paste0(plot_path, "Cluster_new_labels_distribution.png"), width=25, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

cM <- as.matrix(confusionMatrix(ArchR$clusters, ArchR$scHelper_cell_type_new))
scHelper_cell_types <- colnames(cM)[apply(cM, 1 , which.max)]
cluster_idents <- cbind(scHelper_cell_types, rownames(cM))

png(paste0(plot_path, 'scHelper_cell_type_new_cluster_idents_table.png'), height = 20, width = 10, units = 'cm', res = 400)
grid.arrange(tableGrob(cluster_idents, rows=NULL, theme = ttheme_minimal()))
graphics.off()

ArchR$cluster_new_labels <- mapLabels(ArchR$cluster_new_labels, newLabels = cluster_idents, oldLabels = ArchR$clusters)

### Old labels
cM <- confusionMatrix(paste0(ArchR$clusters), paste0(ArchR$scHelper_cell_type_old))
  cM <- cM / Matrix::rowSums(cM)
  p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
  )
png(paste0(plot_path, "Cluster_old_labels_distribution.png"), width=25, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

cM <- as.matrix(confusionMatrix(ArchR$clusters, ArchR$scHelper_cell_type_old))
scHelper_cell_types <- colnames(cM)[apply(cM, 1 , which.max)]
cluster_idents <- cbind(scHelper_cell_types, rownames(cM))

png(paste0(plot_path, 'scHelper_cell_type_old_cluster_labels_table.png'), height = 20, width = 10, units = 'cm', res = 400)
grid.arrange(tableGrob(cluster_idents, rows=NULL, theme = ttheme_minimal()))
graphics.off()

ArchR$cluster_old_labels <- mapLabels(ArchR$cluster_old_labels, newLabels = cluster_idents, oldLabels = ArchR$clusters)

############################## Integration scores plots #######################################

png(paste0(plot_path, 'UMAP_unconInt_Scores.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "predictedScore_Un", plotAs = "points", size = 1.8, baseSize = 0, 
              legendSize = 10)
graphics.off()

png(paste0(plot_path, "VlnPlot_unconInt_Scores.png"), width=50, height=20, units = 'cm', res = 200)
plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", 
  name = "predictedScore_Un", plotAs = "Violin")
graphics.off()

############################## Visualising pseudo-scRNA profiles #######################################
ArchR <- addImputeWeights(ArchR)
  
# look for late marker genes
late_markers <- c(
  "GATA3", "DLX5", "SIX1", "EYA2", #PPR
  "MSX1", "TFAP2A", "TFAP2B", #mix
  "PAX7", "CSRNP1", "SNAI2", "SOX10", #NC
  "SOX2", "SOX21" # neural
)
p1 <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = late_markers, 
  embedding = "UMAP"
)
p2 <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "GeneIntegrationMatrix", 
  name = late_markers, 
  embedding = "UMAP"
)

p1c <- lapply(p1, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

p2c <- lapply(p2, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

png(paste0(plot_path, 'late_markers_GeneScoreMatrix.png'), height = 40, width = 25, units = 'cm', res = 400)
do.call(cowplot::plot_grid, c(list(ncol = 3), p1c))
graphics.off()

png(paste0(plot_path, 'late_markers_GeneIntegrationMatrix.png'), height = 40, width = 25, units = 'cm', res = 400)
do.call(cowplot::plot_grid, c(list(ncol = 3), p2c))
graphics.off()

############################## Constrained integration #######################################
########  do I want to do this?? 

# # which seurat clusters correspond to which cell types?
# seurat_df <- data.frame(cluster_id = seurat$seurat_clusters,
#                         scHelper_cell_type = seurat$scHelper_cell_type)
# seurat_df <- seurat_df %>% remove_rownames(.) %>% distinct(., .keep_all = FALSE)
#Assign scRNA cell states to cluster numbers - automate for every cell state
#clust_FB_rna <- paste0(seurat_df$cluster_id[which(seurat_df$scHelper_cell_type == "FB")], collapse = "|")

# Identify scRNAseq cells that are in these clusters
# rna_FB_cells <- names(seurat$scHelper_cell_type)[which(seurat$scHelper_cell_type == "FB")]
# 
# #Assign scATAC to cell states to cluster numbers - automate for every cell state
# clust_FB_atac <- rownames(cM)[grep('FB', preClust)]
# atac_FB_cells <-ArchR$cellNames[ArchR$clusters %in% clust_FB_atac]
# 
# groupList <- SimpleList(
#   FB = SimpleList(
#     ATAC = atac_FB_cells,
#     RNA = rna_FB_cells
#   )
# )
# 
# # doesnt work - needs to include all ATAC cells
# # how would I do this, every cell type? then not going to differ..
# # just divide neural/NNE?? contam??
# ArchR <- addGeneIntegrationMatrix(
#   ArchRProj = ArchR, 
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "IterativeLSI",
#   seRNA = seurat,
#   addToArrow = FALSE, 
#   groupList = groupList,
#   groupRNA = "BioClassification",
#   nameCell = "predictedCell_Co",
#   nameGroup = "predictedGroup_Co",
#   nameScore = "predictedScore_Co"
# )