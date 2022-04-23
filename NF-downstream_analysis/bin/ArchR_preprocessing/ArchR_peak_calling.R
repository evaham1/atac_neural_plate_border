#!/usr/bin/env Rscript

print("peak_calling_ArchR")

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
library(clustree)
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
    
    plot_path = "./output/NF-downstream_analysis/ArchR_peak_calling/FullData/plots/"
    rds_path = "./output/NF-downstream_analysis/ArchR_peak_calling/FullData/rds_files/"
    
    #data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ss8/ArchR_clustering/rds_files/"
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/7_ArchR_clustering_postfiltering_twice/rds_files/"
    
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
cell_counts <- function(ArchR = ArchR, group1 = "clusters", group2 = "Sample") {
  group1_data <- getCellColData(ArchR, select = group1)[,1]
  group1_cell_counts <- as.data.frame(table(group1_data))
  colnames(group1_cell_counts) <- c("ID", "Total_count")
  
  group2_cell_counts <- data.frame()
  group2_data <- getCellColData(ArchR, select = group2)[,1]
  for (i in unique(group1_data)) {
    data_group1 <- getCellColData(ArchR, select = group1)[,1]
    cells <- ArchR$cellNames[BiocGenerics::which(data_group1 == i)]
    ArchR_subset <- ArchR[cells, ]
    data_group2 <- getCellColData(ArchR_subset, select = group2)[,1]
    group2_cell_counts_i <- as.data.frame(table(data_group2)) %>%
      pivot_wider(names_from = data_group2, values_from = Freq) %>% 
      add_column(ID = !!i)
    group2_cell_counts <- rbind.fill(group2_cell_counts, group2_cell_counts_i)
  }
  
  cell_counts <- merge(group1_cell_counts, group2_cell_counts)
  cell_counts[is.na(cell_counts)] <- 0
  
  ## Ordering rows and columns to better visualise
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


pseudoreplicate_counts <- function(ArchR = ArchR, pseudo_replicates) {
  unlisted <- unlist(pseudo_replicates, recursive=FALSE)
  
  pseudoreplicate_cell_counts_samples <- data.frame()
  for (i in c(1:length(unlisted))){
    group_name <- names(unlisted[i])
    cell_IDs <- unlisted[i][[1]]
    ArchR_pseudo_replicate <- ArchR[cell_IDs, ]
    sample_cell_counts <- as.data.frame(table(ArchR_pseudo_replicate$Sample)) %>%
      pivot_wider(names_from = Var1, values_from = Freq) %>% 
      add_column(pseudo_replicate_ID = group_name)
    pseudoreplicate_cell_counts_samples <- rbind.fill(pseudoreplicate_cell_counts_samples, sample_cell_counts)
  }

  pseudoreplicate_cell_counts_samples <- pseudoreplicate_cell_counts_samples %>%
    mutate(Cluster_ID = gsub("\\..*","", pseudo_replicate_ID)) %>%
    mutate(Cluster_ID = as.numeric(parse_number(Cluster_ID))) %>%
    arrange(Cluster_ID) %>% 
    mutate(sum_contributing_samples = rowSums(is.na(pseudoreplicate_cell_counts_samples[-which(names(pseudoreplicate_cell_counts_samples) %in% c("Cluster_ID", "pseudo_replicate_ID"))]) == FALSE))
  
  pseudoreplicate_cell_counts_samples <- pseudoreplicate_cell_counts_samples[,order(colnames(pseudoreplicate_cell_counts_samples))]
  pseudoreplicate_cell_counts_samples[is.na(pseudoreplicate_cell_counts_samples)] <- 0
  
  grid.arrange(tableGrob(pseudoreplicate_cell_counts_samples, rows=NULL, theme = ttheme_minimal()))
}


############################## Read in ArchR project #######################################

# If files are not in rds_files subdirectory look in input dir
label <- sub('_.*', '', list.files(data_path))
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

###### stage colours
#stage_order <- c("HH4", "HH5", "HH6", "HH7", "ss4", "ss8")
#stage_colours = c("#E78AC3", "#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
#names(stage_colours) <- stage_order


##############################################################################################
############################## Generating pseudo-replicates ##################################

plot_path_1 <- paste0(plot_path, "pseudoreplicates/")
dir.create(plot_path_1, recursive = T)

# Plot number of cells in each cluster that come from each sample
png(paste0(plot_path_1, 'cell_counts_by_sample_table.png'), height = 25, width = 30, units = 'cm', res = 400)
cell_counts(ArchR = ArchR, group1 = "clusters", group2 = "Sample")
graphics.off()

# Make pseudo replicates and see which samples these cells come from
pseudo_replicates <- addGroupCoverages(ArchR, groupBy = "clusters", returnGroups = TRUE)

png(paste0(plot_path_1, 'cell_counts_by_pseudoreplicate_table.png'), height = 80, width = 30, units = 'cm', res = 400)
pseudoreplicate_counts(ArchR, pseudo_replicates)
graphics.off()


#####  Make actual pseudo-replicates for peak calling:
# merge cells within each designated cell group for the generation of pseudo-bulk replicates 
# and then merge these replicates into a single insertion coverage file.
ArchR_pseudo_replicates <- addGroupCoverages(ArchR, groupBy = "clusters", threads = 1, returnGroups = FALSE)
print("pseudo replicates created")

##############################################################################################
############################## Call peaks on pseudo-replicates ###############################

# idxSample <- BiocGenerics::which(ArchR$clusters %in% c("C14", "C15"))
# cellsSample <- ArchR$cellNames[idxSample]
# ArchR <- ArchR[cellsSample, ]

#pathToMacs2 <- findMacs2()

#ArchR_pseudo_replicates <- addReproduciblePeakSet(
#  ArchRProj = ArchR_pseudo_replicates, 
#  groupBy = "clusters", 
#  pathToMacs2 = pathToMacs2
#)
#print("peaks called using Macs2")
#getPeakSet(ArchR_pseudo_replicates)

#ArchR <- addReproduciblePeakSet(
#  ArchRProj = ArchR_pseudo_replicates, 
#  groupBy = "clusters",
#  peakMethod = "Tiles",
#  method = "p",
#  threads = 1
#)
