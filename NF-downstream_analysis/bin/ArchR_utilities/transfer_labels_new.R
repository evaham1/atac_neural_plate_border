#!/usr/bin/env Rscript

print("ArchR transfer cell metadata from stage data to full data")
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
library(parallel)

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
    
    data_path = "./output/NF-downstream_analysis/test_inputs/test_input_transfer_labels/"
    
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

############################## Read in ArchR projects #######################################
# Read in all data
files <- list.files(data_path, full.names = TRUE)
print(paste0("All input files: ", files))

stages_data <- grep("FullData", files, invert = T, value = TRUE) # source data from which labels are extracted
print(paste0("Stages data: ", stages_data))

full_data <- grep("FullData", files, invert = F, value = TRUE) # destination data where labels will be transfered onto
print(paste0("Destination data: ", full_data))
ArchR_full <- loadArchRProject(path = full_data, force = FALSE, showLogo = TRUE)
print(paste0("Number of cells in fulldata: ", length(rownames(ArchR_full@cellColData))))

################## Transfer cluster info from stages onto full data #######################################
######## stages: 'clusters' -> TransferLabels: 'stage_clusters'

# extract cell IDs from each of the stages datasets
all_cell_ids <- c()
all_cluster_ids <- c()
for (i in 1:5) {
  stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
  print(paste0("Iteration: ", i, ", Stage: ", stage))
  
  ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
  cell_ids <- ArchR$cellNames
  print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
  all_cell_ids <- c(all_cell_ids, cell_ids)
  
  stage <- unique(ArchR$stage)
  cluster_ids <- paste0(stage, "_", ArchR$clusters)
  print(length(cluster_ids))
  all_cluster_ids <- c(all_cluster_ids, cluster_ids)
  
}
print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
print(paste0("Length of all stage cluster ids: ", length(all_cluster_ids)))

# check that the stage cell ids are all found in the fulldata object
intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
if (intersect_cell_id_length != length(all_cell_ids)){
  stop("Error! Not all stage cell ids match with cell ids in Full data")
}

# add the stage_clusters to the full dataset
ArchR_full <- addCellColData(ArchRProj = ArchR_full, 
                        data = all_cluster_ids,
                        cells = all_cell_ids, 
                        name = "stage_clusters",
                        force = TRUE)
print(table(ArchR_full$stage_clusters))

# filter out cells that have NA in the stage_clusters (ie have been removed from stages because they are contam)
print(paste0("number of NA cell ids: ", sum(is.na(ArchR_full$stage_clusters))))

idxSample <- BiocGenerics::which(!is.na(ArchR_full$stage_clusters))
cellsSample <- ArchR_full$cellNames[idxSample]
ArchR_filtered <- ArchR_full[cellsSample, ]

################## Transfer cluster identity info from stages onto full data #######################################
######## stages: 'scHelper_cell_type_old' -> TransferLabels: 'stage_scHelper_cell_type_old'

# extract cell IDs from each of the stages datasets
all_cell_ids <- c()
all_cluster_ids <- c()
for (i in 1:5) {
  stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
  print(paste0("Iteration: ", i, ", Stage: ", stage))
  
  ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
  cell_ids <- ArchR$cellNames
  print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
  all_cell_ids <- c(all_cell_ids, cell_ids)
  
  stage <- unique(ArchR$stage)
  cluster_ids <- paste0(stage, "_", ArchR$scHelper_cell_type_old)
  print(length(cluster_ids))
  all_cluster_ids <- c(all_cluster_ids, cluster_ids)
  
}
print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
print(paste0("Length of all stage cluster labels: ", length(all_cluster_ids)))

# check that the stage cell ids are all found in the fulldata object
intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
if (intersect_cell_id_length != length(all_cell_ids)){
  stop("Error! Not all stage cell ids match with cell ids in Full data")
}

# add the stage_clusters to the full dataset
ArchR_filtered <- addCellColData(ArchRProj = ArchR_filtered, 
                        data = all_cluster_ids,
                        cells = all_cell_ids, 
                        name = "stage_scHelper_cell_type_old",
                        force = TRUE)
print(table(ArchR_filtered$stage_scHelper_cell_type_old))

################## Transfer cluster identity info from stages onto full data #######################################
######## stages: 'cluster_labels' -> TransferLabels: 'stage_cluster_labels'

# extract cell IDs from each of the stages datasets
all_cell_ids <- c()
all_cluster_ids <- c()
for (i in 1:5) {
  stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
  print(paste0("Iteration: ", i, ", Stage: ", stage))
  
  ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
  cell_ids <- ArchR$cellNames
  print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
  all_cell_ids <- c(all_cell_ids, cell_ids)
  
  stage <- unique(ArchR$stage)
  cluster_ids <- paste0(stage, "_", ArchR$cluster_labels)
  print(length(cluster_ids))
  all_cluster_ids <- c(all_cluster_ids, cluster_ids)
  
}
print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
print(paste0("Length of all stage cluster labels: ", length(all_cluster_ids)))

# check that the stage cell ids are all found in the fulldata object
intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
if (intersect_cell_id_length != length(all_cell_ids)){
  stop("Error! Not all stage cell ids match with cell ids in Full data")
}

# add the stage_clusters to the full dataset
ArchR_filtered <- addCellColData(ArchRProj = ArchR_filtered, 
                        data = all_cluster_ids,
                        cells = all_cell_ids, 
                        name = "stage_cluster_labels",
                        force = TRUE)
print(table(ArchR_filtered$stage_cluster_labels))

# save transfer_labels data
saveArchRProject(ArchRProj = ArchR_filtered, outputDirectory = paste0(rds_path, "TransferLabels_Save-ArchR"), load = FALSE)


#####################################################################################
############################## Visualisations #######################################

### Plot cell counts before and after subsetting per stage
unfiltered <- table(ArchR$stage)
filtered <- table(ArchR_filtered$stage)
cell_counts <- as_tibble(rbind(unfiltered, filtered))
cell_counts <- cbind(cell_counts, Total = rowSums(cell_counts))

png(paste0(plot_path, 'cell_counts_table_stages.png'), height = 10, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

#### Stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

### Plot stages and stage clusters
p1 <- plotEmbedding(ArchR_filtered, 
                    name = "stage",
                    plotAs = "points", size = ifelse(length(unique(ArchR_filtered$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    pal = stage_colours, randomize = TRUE)
p2 <- plotEmbedding(ArchR_filtered, 
                    name = "stage_clusters",
                    plotAs = "points", size = ifelse(length(unique(ArchR_filtered$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0,
                    randomize = TRUE)

png(paste0(plot_path, "cluster_UMAPs.png"), width=60, height=40, units = 'cm', res = 200)
ggAlignPlots(p1, p2, type = "h")
graphics.off()

### Plot stage_scHelper_cell_type and stage_cluster_labels
p1 <- plotEmbedding(ArchR_filtered, 
                    name = "stage_scHelper_cell_type_old",
                    plotAs = "points", size = ifelse(length(unique(ArchR_filtered$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    pal = stage_colours, randomize = TRUE)
p2 <- plotEmbedding(ArchR_filtered, 
                    name = "stage_cluster_labels",
                    plotAs = "points", size = ifelse(length(unique(ArchR_filtered$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0,
                    randomize = TRUE)

png(paste0(plot_path, "stage_scHelper_cell_type_old_cell_type_UMAPs.png"), width=60, height=40, units = 'cm', res = 200)
ggAlignPlots(p1, p2, type = "h")
graphics.off()