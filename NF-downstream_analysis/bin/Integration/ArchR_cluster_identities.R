#!/usr/bin/env Rscript

print("ArchR cell state plots")

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
library(plyr)

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
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/HH5/ArchR_clustering//rds_files/"
    
    # already integrated
    data_path = "./output/NF-downstream_analysis/ArchR_integration//ss8/1_unconstrained_integration/rds_files/"
    data_path = "./output/NF-downstream_analysis/ArchR_integration/HH5/1_unconstrained_integration/rds_files/"
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


############################### FUNCTIONS ####################################

# function to count how many cells from each cluster/sample are assigned the same label/cluster
cell_counts <- function(ArchR = ArchR, group1 = "clusters", group2 = "Sample") {
  group1_data <- getCellColData(ArchR, select = group1)[,1]
  group1_cell_counts <- as.data.frame(table(group1_data))
  colnames(group1_cell_counts) <- c("ID", "Total_count")
  
  group2_cell_counts <- data.frame()
  group2_data <- getCellColData(ArchR, select = group2)[,1]
  data_group1 <- getCellColData(ArchR, select = group1)[,1]
  for (i in unique(group1_data)) {
    cells <- ArchR$cellNames[BiocGenerics::which(data_group1 == i)]
    if (length(cells) > 1){
      ArchR_subset <- ArchR[cells, ]
      data_group2 <- getCellColData(ArchR_subset, select = group2)[,1]
      group2_cell_counts_i <- as.data.frame(table(data_group2)) %>%
        pivot_wider(names_from = data_group2, values_from = Freq) %>% 
        add_column(ID = !!i)
      group2_cell_counts <- rbind.fill(group2_cell_counts, group2_cell_counts_i)
    }
  }
  
  cell_counts <- merge(group1_cell_counts, group2_cell_counts)
  
  
  # Ordering rows and columns to better visualise
  if (group1 == "clusters"){
    cell_counts <- cell_counts %>% 
      mutate(ID = as.numeric(gsub('^.', '', ID))) %>%
      arrange(ID)
    }
  
  if (group2 == "clusters"){
    new_names <- as.numeric(gsub('^.', '', colnames(cell_counts)[3:length(colnames(cell_counts))]))
    colnames(cell_counts)[3:length(colnames(cell_counts))] <- new_names
    cell_counts <- cell_counts[, c("ID", "Total_count", 1:max(new_names))]
  }
  
  grid.arrange(tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
}

# function to make heatmap showing contribution of cell groups to other cell groups
cell_counts_heatmap <- function(ArchR = ArchR, group1 = "scHelper_cell_type_new", group2 = "clusters") {
  group1_data <- getCellColData(ArchR, select = group1)[,1]
  group2_data <- getCellColData(ArchR, select = group2)[,1]
  cM <- confusionMatrix(paste0(group2_data), paste0(group1_data))
  cM <- cM / Matrix::rowSums(cM)
  
  p <- pheatmap::pheatmap(
    mat = cM,
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
  )
}

############################## Read in ArchR project and seurat object #######################################

# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")


############################################################################################
############################## COLOURS #######################################

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
# set colour palettes for UMAPs
atac_scHelper_new_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_new)]
atac_scHelper_old_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_old)]

#############################################################################################
############################## Assign cluster labels #######################################

########### New labels
plot_path = "./plots/new_labels/"
dir.create(plot_path, recursive = T)

# visualise distribution across clusters
png(paste0(plot_path, 'counts_by_cluster_table.png'), height = 25, width = 40, units = 'cm', res = 400)
cell_counts(ArchR = ArchR, group1 = "scHelper_cell_type_new", group2 = "clusters")
graphics.off()

png(paste0(plot_path, "cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
cell_counts_heatmap(ArchR = ArchR, group1 = "scHelper_cell_type_new", group2 = "clusters")
graphics.off()

# assign cluster labels
cM <- as.matrix(confusionMatrix(ArchR$clusters, ArchR$scHelper_cell_type_new))
scHelper_cell_types <- colnames(cM)[apply(cM, 1 , which.max)]
cluster_idents <- cbind(scHelper_cell_types, rownames(cM))

png(paste0(plot_path, 'assigned_cluster_idents_table.png'), height = 20, width = 10, units = 'cm', res = 400)
grid.arrange(tableGrob(cluster_idents, rows=NULL, theme = ttheme_minimal()))
graphics.off()

new_labels <- cluster_idents[,1]
names(new_labels) <- cluster_idents[,2]
ArchR$cluster_new_labels <- mapLabels(ArchR$clusters, newLabels = new_labels)

p1 <- plotEmbedding(ArchR, name = "cluster_new_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_new_cols, labelAsFactors = FALSE)
p2 <- plotEmbedding(ArchR, name = "scHelper_cell_type_new", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_new_cols, labelAsFactors = FALSE)

png(paste0(plot_path, 'assigned_cluster_idents_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
print(p1)
graphics.off()

png(paste0(plot_path, 'assigned_cluster_idents_UMAP_comparison.png'), height = 20, width = 40, units = 'cm', res = 400)
print(p1 + p2)
graphics.off()

########### Old labels
plot_path = "./plots/old_labels/"
dir.create(plot_path, recursive = T)

# visualise distribution across clusters
png(paste0(plot_path, 'counts_by_cluster_table.png'), height = 25, width = 40, units = 'cm', res = 400)
cell_counts(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "clusters")
graphics.off()

png(paste0(plot_path, "cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
cell_counts_heatmap(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "clusters")
graphics.off()

# assign cluster labels
cM <- as.matrix(confusionMatrix(ArchR$clusters, ArchR$scHelper_cell_type_old))
scHelper_cell_types <- colnames(cM)[apply(cM, 1 , which.max)]
cluster_idents <- cbind(scHelper_cell_types, rownames(cM))

png(paste0(plot_path, 'assigned_cluster_idents_table.png'), height = 20, width = 10, units = 'cm', res = 400)
grid.arrange(tableGrob(cluster_idents, rows=NULL, theme = ttheme_minimal()))
graphics.off()

new_labels <- cluster_idents[,1]
names(new_labels) <- cluster_idents[,2]
ArchR$cluster_old_labels <- mapLabels(ArchR$clusters, newLabels = new_labels)

p1 <- plotEmbedding(ArchR, name = "cluster_old_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_new_cols, labelAsFactors = FALSE)
p2 <- plotEmbedding(ArchR, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_new_cols, labelAsFactors = FALSE)

png(paste0(plot_path, 'assigned_cluster_idents_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
print(p1)
graphics.off()

png(paste0(plot_path, 'assigned_cluster_idents_UMAP_comparison.png'), height = 20, width = 40, units = 'cm', res = 400)
print(p1 + p2)
graphics.off()