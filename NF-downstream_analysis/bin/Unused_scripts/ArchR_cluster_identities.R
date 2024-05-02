#!/usr/bin/env Rscript

print("ArchR cell state plots and assign identities based on label proportions")

############################## Load libraries #######################################
library(getopt)
library(optparse)
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
library(gtools)

############################## Set up script options #######################################

option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-m", "--min_threshold"), action = "store", type = "integer", help = "minimum percentage of cells a label must contain", default = 40),
  make_option(c("-l", "--max_label"), action = "store", type = "integer", help = "maximum number of labels a cluster can have", default = 3),
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

    data_path = "./output/NF-downstream_analysis/Processing/ss8/ARCHR_INTEGRATING_WF/Single_cell_integration/rds_files/"
    data_path = "./output/NF-downstream_analysis/Processing/HH5/ARCHR_INTEGRATING_WF/Single_cell_integration/rds_files/"
    
    data_path = "./output/NF-downstream_analysis/Processing/HH5/INTEGRATING/Single_cell_integration/rds_files/"
    
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

############################## Read in ArchR project  #######################################

# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
print('data read in')
print(ArchR)

# check that gene score matrix and gene integration matrix are available
getAvailableMatrices(ArchRProj = ArchR)

############################################################################################
############################## COLOURS #######################################

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

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
atac_scHelper_old_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type)]


############################## Gene scores plots #######################################
#### compare gene scores with integrated gene exp values

# plot_path = "./plots/gene_scores_vs_integrated_gex/"
# dir.create(plot_path, recursive = T)

# ArchR <- addImputeWeights(ArchR)
# print("Impute weights added")

# # look for late marker genes
# late_markers <- c(
#   "GATA3", "DLX5", "SIX1", "EYA2", #PPR
#   "MSX1", "TFAP2A", "TFAP2B", #mix
#   "PAX7", "CSRNP1", "SNAI2", "SOX10", #NC
#   "SOX2", "SOX21" # neural
# )

# png(paste0(plot_path, 'late_markers_GeneScoreMatrix.png'), height = 40, width = 25, units = 'cm', res = 400)
# scHelper::ArchR_FeaturePlotGrid(ArchR, matrix = "GeneScoreMatrix", late_markers)
# graphics.off()

# png(paste0(plot_path, 'late_markers_GeneIntegrationMatrix.png'), height = 40, width = 25, units = 'cm', res = 400)
# scHelper::ArchR_FeaturePlotGrid(ArchR, matrix = "GeneIntegrationMatrix", late_markers)
# graphics.off()


##################### Distribution of labels across clusters ##################################

plot_path = "./plots/label_by_cluster_distribution/"
dir.create(plot_path, recursive = T)

# visualise distribution across clusters: table of cell counts
png(paste0(plot_path, 'label_by_cluster_cell_number_table.png'), height = 25, width = 40, units = 'cm', res = 400)
scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters", print_table = TRUE, scHelper_cell_type_order = scHelper_cell_type_order)
graphics.off()

# visualise distribution across clusters: confusion matrix
png(paste0(plot_path, "label_by_cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
scHelper::ArchRCellCountsHeatmap(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters")
graphics.off()

# visualise distribution across clusters: table of cell percentages
cell_counts <- scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
percentage_counts <- as.data.frame(lapply(cell_counts, function(x) (x / sum(x))*100))
rownames(percentage_counts) <- rownames(cell_counts)

png(paste0(plot_path, 'label_by_cluster_cell_percentage_table.png'), height = 25, width = 40, units = 'cm', res = 400)
grid.arrange(tableGrob(round(percentage_counts, 2), theme = ttheme_minimal()))
graphics.off()

# visualise distribution across clusters: piecharts
counts <- scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
png(paste0(plot_path, "label_by_cluster_piecharts.png"), width=50, height=40, units = 'cm', res = 200)
scHelper::CellLabelPieCharts(counts, col = scHelper_cell_type_colours)
graphics.off()

##################### Label clusters based on thresholds ##################################

plot_path = "./plots/label_by_cluster_distribution/assigned_cluster_labels/"
dir.create(plot_path, recursive = T)

identities <- c()

for(i in 1:ncol(percentage_counts)) {       # for-loop over columns (ie clusters)
  column <- percentage_counts[ , i]
  names(column) <- rownames(percentage_counts)
  
  identity <- ifelse(
    # do any labels pass the minimum threshold? 
    # if NO labels label more than the minimum threshold, OR MORE THAN max_label labels pass the minimum threshold -> mixed
    any(sum(column > opt$min_threshold) == 0 | sum(column > opt$min_threshold) > opt$max_label), "MIXED", # condition - mixed
                     ifelse(
                       sum(column > opt$min_threshold) == 1, names(column[column > opt$min_threshold]), # condition 2 - monolabel
                            paste(names( sort(column[column > opt$min_threshold], decreasing = TRUE) ), collapse='/') # condition 3 - multilabel
                     ) 
  )
  identities <- c(identities, identity)
  
}

# assign cluster labels
cluster_idents <- cbind(colnames(percentage_counts), identities)

png(paste0(plot_path, 'assigned_cluster_idents_table.png'), height = 20, width = 10, units = 'cm', res = 400)
grid.arrange(tableGrob(cluster_idents, rows=NULL, theme = ttheme_minimal()))
graphics.off()

new_labels <- cluster_idents[,2]
names(new_labels) <- cluster_idents[,1]
ArchR$cluster_labels <- mapLabels(ArchR$clusters, newLabels = new_labels)

# save ArchR object with cluster labels
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label[1], "_Save-ArchR"), load = FALSE)

# plot cluster labels on UMAPs
p1 <- plotEmbedding(ArchR, name = "cluster_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, labelAsFactors = FALSE)
p2 <- plotEmbedding(ArchR, name = "scHelper_cell_type", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE)

png(paste0(plot_path, 'assigned_cluster_idents_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
print(p1)
graphics.off()

png(paste0(plot_path, 'assigned_cluster_idents_UMAP_comparison.png'), height = 20, width = 40, units = 'cm', res = 400)
print(p1 + p2)
graphics.off()