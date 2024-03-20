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
  make_option(c("-t", "--target_name"), action = "store", type = "character", help = "File name of the ArchR full data object on which to transfer labels", default = "FullData"),
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
    addArchRThreads(threads = 1) 
    
    rds_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/plots/"
    
    ## working interactively:
    stages_data <- c("./output/NF-downstream_analysis/Processing/HH5/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH5_Save-ArchR",
                     "./output/NF-downstream_analysis/Processing/HH6/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH6_Save-ArchR",
                     "./output/NF-downstream_analysis/Processing/HH7/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/HH7_Save-ArchR",
                     "./output/NF-downstream_analysis/Processing/ss4/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss4_Save-ArchR",
                     "./output/NF-downstream_analysis/Processing/ss8/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss8_Save-ArchR")
    full_data <- "./output/NF-downstream_analysis/Processing/FullData/Peak_call/rds_files/FullData_Save-ArchR/"
    
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

############################## Read in ArchR projects #######################################
# Read in all data
files <- list.files(data_path, full.names = TRUE)
print(paste0("All input files: ", files))

stages_data  <- grep(opt$target_name, files, invert = T, value = TRUE) # source data from which labels are extracted
print(paste0("Stages data: ", stages_data))

full_data <- grep(opt$target_name, files, invert = F, value = TRUE) # destination data where labels will be transfered onto

print(paste0("Destination data: ", full_data))
ArchR_full <- loadArchRProject(path = full_data, force = FALSE, showLogo = TRUE)
print(paste0("Number of cells in fulldata: ", length(rownames(ArchR_full@cellColData))))

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR_full)
getAvailableMatrices(ArchR_full)

################## Transfer metacell info from stages onto full data #######################################

print("Transferring SEACell_ID...")

# extract cell IDs from each of the stages datasets
all_cell_ids <- c()
all_group_ids <- c()
for (i in 1:5) {
  stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
  print(paste0("Iteration: ", i, ", Stage: ", stage))
  
  ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
  cell_ids <- ArchR$cellNames
  print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
  all_cell_ids <- c(all_cell_ids, cell_ids)
  
  stage <- unique(ArchR$stage)

  group_ids <- ArchR$SEACell_ID
  print(length(group_ids))
  all_group_ids <- c(all_group_ids, group_ids)
  
}
print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
print(paste0("Length of all stage cluster ids: ", length(all_group_ids)))

# check that the stage cell ids are all found in the fulldata object
intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
if (intersect_cell_id_length != length(all_cell_ids)){
  stop("Error! Not all stage cell ids match with cell ids in Full data")
}

# add the stage_clusters to the full dataset
ArchR_full <- addCellColData(ArchRProj = ArchR_full, 
                        data = all_group_ids,
                        cells = all_cell_ids, 
                        name = "SEACell_ID",
                        force = TRUE)

##############
print("Transferring RNA_SEACell_ID...")

# extract cell IDs from each of the stages datasets
all_cell_ids <- c()
all_group_ids <- c()
for (i in 1:5) {
  stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
  print(paste0("Iteration: ", i, ", Stage: ", stage))
  
  ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
  cell_ids <- ArchR$cellNames
  print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
  all_cell_ids <- c(all_cell_ids, cell_ids)
  
  stage <- unique(ArchR$stage)

  group_ids <- ArchR$RNA_SEACell_ID
  print(length(group_ids))
  all_group_ids <- c(all_group_ids, group_ids)
  
}
print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
print(paste0("Length of all stage cluster ids: ", length(all_group_ids)))

# check that the stage cell ids are all found in the fulldata object
intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
if (intersect_cell_id_length != length(all_cell_ids)){
  stop("Error! Not all stage cell ids match with cell ids in Full data")
}

# add the stage_clusters to the full dataset
ArchR_full <- addCellColData(ArchRProj = ArchR_full, 
                        data = all_group_ids,
                        cells = all_cell_ids, 
                        name = "RNA_SEACell_ID",
                        force = TRUE)

##############
print("Transferring SEACell_scHelper_cell_type...")

# extract cell IDs from each of the stages datasets
all_cell_ids <- c()
all_group_ids <- c()
for (i in 1:5) {
  stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
  print(paste0("Iteration: ", i, ", Stage: ", stage))
  
  ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
  cell_ids <- ArchR$cellNames
  print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
  all_cell_ids <- c(all_cell_ids, cell_ids)
  
  stage <- unique(ArchR$stage)

  group_ids <- ArchR$SEACell_scHelper_cell_type
  print(length(group_ids))
  all_group_ids <- c(all_group_ids, group_ids)
  
}
print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
print(paste0("Length of all stage cluster ids: ", length(all_group_ids)))

# check that the stage cell ids are all found in the fulldata object
intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
if (intersect_cell_id_length != length(all_cell_ids)){
  stop("Error! Not all stage cell ids match with cell ids in Full data")
}

# add the stage_clusters to the full dataset
ArchR_full <- addCellColData(ArchRProj = ArchR_full, 
                        data = all_group_ids,
                        cells = all_cell_ids, 
                        name = "SEACell_scHelper_cell_type",
                        force = TRUE)

##############
print("Transferring SEACell_scHelper_cell_type_broad...")

# extract cell IDs from each of the stages datasets
all_cell_ids <- c()
all_group_ids <- c()
for (i in 1:5) {
  stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
  print(paste0("Iteration: ", i, ", Stage: ", stage))
  
  ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
  cell_ids <- ArchR$cellNames
  print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
  all_cell_ids <- c(all_cell_ids, cell_ids)
  
  stage <- unique(ArchR$stage)

  group_ids <- ArchR$SEACell_scHelper_cell_type_broad
  print(length(group_ids))
  all_group_ids <- c(all_group_ids, group_ids)
  
}
print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
print(paste0("Length of all stage cluster ids: ", length(all_group_ids)))

# check that the stage cell ids are all found in the fulldata object
intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
if (intersect_cell_id_length != length(all_cell_ids)){
  stop("Error! Not all stage cell ids match with cell ids in Full data")
}

# add the stage_clusters to the full dataset
ArchR_full <- addCellColData(ArchRProj = ArchR_full, 
                        data = all_group_ids,
                        cells = all_cell_ids, 
                        name = "SEACell_scHelper_cell_type_broad",
                        force = TRUE)


# ################## Transfer cluster info from stages onto full data #######################################
# ######## stages: 'clusters' -> TransferLabels: 'stage_clusters'

# print("Transferring cluster IDs...")

# # extract cell IDs from each of the stages datasets
# all_cell_ids <- c()
# all_cluster_ids <- c()
# for (i in 1:5) {
#   stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
#   print(paste0("Iteration: ", i, ", Stage: ", stage))
  
#   ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
#   cell_ids <- ArchR$cellNames
#   print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
#   all_cell_ids <- c(all_cell_ids, cell_ids)
  
#   stage <- unique(ArchR$stage)
#   cluster_ids <- paste0(stage, "_", ArchR$clusters)
#   print(length(cluster_ids))
#   all_cluster_ids <- c(all_cluster_ids, cluster_ids)
  
# }
# print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
# print(paste0("Length of all stage cluster ids: ", length(all_cluster_ids)))

# # check that the stage cell ids are all found in the fulldata object
# intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
# print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
# if (intersect_cell_id_length != length(all_cell_ids)){
#   stop("Error! Not all stage cell ids match with cell ids in Full data")
# }

# # add the stage_clusters to the full dataset
# ArchR_full <- addCellColData(ArchRProj = ArchR_full, 
#                         data = all_cluster_ids,
#                         cells = all_cell_ids, 
#                         name = "stage_clusters",
#                         force = TRUE)
# print(table(ArchR_full$stage_clusters))

# # filter out cells that have NA in the stage_clusters (ie have been removed from stages because they are contam)
# ### update shouldn't filter out any cells but run anyway
# print(paste0("number of NA cell ids: ", sum(is.na(ArchR_full$stage_clusters))))

# idxSample <- BiocGenerics::which(!is.na(ArchR_full$stage_clusters))
# cellsSample <- ArchR_full$cellNames[idxSample]
# ArchR_filtered <- ArchR_full[cellsSample, ]

# ################## Transfer cluster identity info from stages onto full data #######################################
# ######## stages: 'scHelper_cell_type' -> TransferLabels: 'stage_scHelper_cell_type'

# print("Transferring scHelper cell type...")

# # extract cell IDs from each of the stages datasets
# all_cell_ids <- c()
# all_cluster_ids <- c()
# for (i in 1:5) {
#   stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
#   print(paste0("Iteration: ", i, ", Stage: ", stage))
  
#   ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
#   cell_ids <- ArchR$cellNames
#   print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
#   all_cell_ids <- c(all_cell_ids, cell_ids)
  
#   stage <- unique(ArchR$stage)
#   cluster_ids <- paste0(stage, "_", ArchR$scHelper_cell_type)
#   print(length(cluster_ids))
#   all_cluster_ids <- c(all_cluster_ids, cluster_ids)
  
# }
# print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
# print(paste0("Length of all stage cluster labels: ", length(all_cluster_ids)))

# # check that the stage cell ids are all found in the fulldata object
# intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
# print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
# if (intersect_cell_id_length != length(all_cell_ids)){
#   stop("Error! Not all stage cell ids match with cell ids in Full data")
# }

# # add the labels to the full dataset
# ArchR_filtered <- addCellColData(ArchRProj = ArchR_filtered, 
#                         data = all_cluster_ids,
#                         cells = all_cell_ids, 
#                         name = "stage_scHelper_cell_type",
#                         force = TRUE)
# print(table(ArchR_filtered$stage_scHelper_cell_type))

# # remove the stage prefix for another slot
# ArchR_filtered <- addCellColData(ArchRProj = ArchR_filtered, 
#                                  data = substring(all_cluster_ids, 5),
#                                  cells = all_cell_ids, 
#                                  name = "scHelper_cell_type",
#                                  force = TRUE)
# print(table(ArchR_filtered$scHelper_cell_type))

# ################## Transfer cluster identity info from stages onto full data #######################################
# ######## stages: 'scHelper_cell_type_broad' -> TransferLabels: 'stage_scHelper_cell_type_broad'

# print("Transferring scHelper cell type broad...")

# # extract cell IDs from each of the stages datasets
# all_cell_ids <- c()
# all_cluster_ids <- c()
# for (i in 1:5) {
#   stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
#   print(paste0("Iteration: ", i, ", Stage: ", stage))
  
#   ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
#   cell_ids <- ArchR$cellNames
#   print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
#   all_cell_ids <- c(all_cell_ids, cell_ids)
  
#   stage <- unique(ArchR$stage)
#   cluster_ids <- paste0(stage, "_", ArchR$scHelper_cell_type_broad)
#   print(length(cluster_ids))
#   all_cluster_ids <- c(all_cluster_ids, cluster_ids)
  
# }
# print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
# print(paste0("Length of all stage cluster labels: ", length(all_cluster_ids)))

# # check that the stage cell ids are all found in the fulldata object
# intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
# print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
# if (intersect_cell_id_length != length(all_cell_ids)){
#   stop("Error! Not all stage cell ids match with cell ids in Full data")
# }

# # add the stage_clusters to the full dataset
# ArchR_filtered <- addCellColData(ArchRProj = ArchR_filtered, 
#                                  data = all_cluster_ids,
#                                  cells = all_cell_ids, 
#                                  name = "stage_scHelper_cell_type_broad",
#                                  force = TRUE)
# print(table(ArchR_filtered$stage_scHelper_cell_type_broad))

# # remove the stage prefix for another slot
# ArchR_filtered <- addCellColData(ArchRProj = ArchR_filtered, 
#                                  data = substring(all_cluster_ids, 5),
#                                  cells = all_cell_ids, 
#                                  name = "scHelper_cell_type_broad",
#                                  force = TRUE)
# print(table(ArchR_filtered$scHelper_cell_type_broad))

# ################## Transfer rna cell label info from stages onto full data #######################################
# ######## stages: 'predictedCell_Un' -> TransferLabels: 'predictedCell_Un'

# print("Transferring predictedCell_Un...")

# # extract cell IDs from each of the stages datasets
# all_cell_ids <- c()
# all_transfer_ids <- c()
# for (i in 1:5) {
#   stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
#   print(paste0("Iteration: ", i, ", Stage: ", stage))
  
#   ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
#   cell_ids <- ArchR$cellNames
#   print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
#   all_cell_ids <- c(all_cell_ids, cell_ids)
  
#   transfer_ids <- paste0(ArchR$predictedCell_Un)
#   print(length(transfer_ids))
#   all_transfer_ids <- c(all_transfer_ids, transfer_ids)
  
# }
# print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
# print(paste0("Length of all stage transfer ids: ", length(all_transfer_ids)))

# # check that the stage cell ids are all found in the fulldata object
# intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
# print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
# if (intersect_cell_id_length != length(all_cell_ids)){
#   stop("Error! Not all stage cell ids match with cell ids in Full data")
# }

# # add the stage_clusters to the full dataset
# ArchR_filtered <- addCellColData(ArchRProj = ArchR_filtered, 
#                                  data = all_transfer_ids,
#                                  cells = all_cell_ids, 
#                                  name = "predictedCell_Un",
#                                  force = TRUE)
# print(table(ArchR_filtered$predictedCell_Un))

# ################## Transfer rna cell label info from stages onto full data #######################################
# ######## stages: 'predictedScore_Un' -> TransferLabels: 'predictedScore_Un'

# print("Transferring predictedScore_Un...")

# # extract cell IDs from each of the stages datasets
# all_cell_ids <- c()
# all_transfer_ids <- c()
# for (i in 1:5) {
#   stage <- substr(sub('.*\\/', '', stages_data[i]), 1, 3)
#   print(paste0("Iteration: ", i, ", Stage: ", stage))
  
#   ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  
#   cell_ids <- ArchR$cellNames
#   print(paste0("length of cell_ids in ", stage, ": ", length(cell_ids)))
#   all_cell_ids <- c(all_cell_ids, cell_ids)
  
#   transfer_ids <- paste0(ArchR$predictedScore_Un)
#   print(length(transfer_ids))
#   all_transfer_ids <- c(all_transfer_ids, transfer_ids)
  
# }
# print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))
# print(paste0("Length of all stage transfer ids: ", length(all_transfer_ids)))

# # check that the stage cell ids are all found in the fulldata object
# intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
# print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
# if (intersect_cell_id_length != length(all_cell_ids)){
#   stop("Error! Not all stage cell ids match with cell ids in Full data")
# }

# # add the stage_clusters to the full dataset
# ArchR_filtered <- addCellColData(ArchRProj = ArchR_filtered, 
#                                  data = all_transfer_ids,
#                                  cells = all_cell_ids, 
#                                  name = "predictedScore_Un",
#                                  force = TRUE)
# print(table(ArchR_filtered$predictedScore_Un))

#####################################################################################
############################## Save data #######################################

# save transfer_labels data
saveArchRProject(ArchRProj = ArchR_full, outputDirectory = paste0(rds_path, "TransferLabels_Save-ArchR"), load = FALSE)

#####################################################################################
############################## Visualisations #######################################


#### Stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

###### schelper cell type colours ~ 29 cell states
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak', 
                              'PGC', 'BI', 'meso', 'endo', 'MIXED', 'Unmapped',
                              'Neural', 'Placodal', 'Non-neural', 'Contam')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                "#A5718D", "#3F918C", "#ed5e5f", "#9792A3",
                                "#7C8483", "#EAEAEA")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                       'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo',
                                       'Neural', 'Placodal', 'Non-neural', 'Contam',
                                       'MIXED', 'Unmapped')

p1 <- plotEmbedding(ArchR_full, 
                    name = "SEACell_ID",
                    plotAs = "points", size = ifelse(length(unique(ArchR_full$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    randomize = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
png(paste0(plot_path, "SEACell_IDs_UMAPs.png"), width=60, height=100, units = 'cm', res = 200)
print(p1)
graphics.off()

png(paste0(plot_path, "SEACell_scHelper_cell_type.png"), width=30, height=40, units = 'cm', res = 200)
plotEmbedding(ArchR_full,
              name = "SEACell_scHelper_cell_type",
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE,
              pal = scHelper_cell_type_colours) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

png(paste0(plot_path, "SEACell_scHelper_cell_type_broad.png"), width=30, height=40, units = 'cm', res = 200)
plotEmbedding(ArchR_full,
              name = "SEACell_scHelper_cell_type_broad",
              plotAs = "points", size = 2.5,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE,
              pal = scHelper_cell_type_colours) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()