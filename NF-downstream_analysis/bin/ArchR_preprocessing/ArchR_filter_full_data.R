#!/usr/bin/env Rscript

print("ArchR filter full data based on stage data")
# read in filtered 5 stage data sets, extract Cell IDs from them and use that to filter full data

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(gridExtra)
library(grid)
library(parallel)

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
    
    #plot_path = "./output/NF-downstream_analysis/4_ArchR_filter_clusters/plots/"
    #rds_path = "./output/NF-downstream_analysis/4_ArchR_filter_clusters/rds_files/"
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/test_input/"

    addArchRThreads(threads = 1)
    
    opt$invert1 = TRUE
    opt$groups1 = "BI,PGC,meso,endo"
    opt$meta_col1 = "scHelper_cell_type_old"
  
    
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

set.seed(42)

############################## Read in ArchR projects #######################################
# Read in all data
files <- list.files(data_path, full.names = TRUE)
print(files)
stages_data <- grep("FullData", files, invert = T, value = TRUE) # source data from which labels are extracted
print(paste0("Stages data: ", stages_data))
full_data <- grep("FullData", files, invert = F, value = TRUE) # destination data where labels will be transfered onto
print(paste0("Destination data: ", full_data))

# load full data
ArchR_full <- loadArchRProject(path = full_data, force = FALSE, showLogo = TRUE)
print(paste0("Number of cells in fulldata: ", length(rownames(ArchR_full@cellColData))))

# extract cell IDs from each of the stages datasets
all_cell_ids <- c()
for (i in 1:5) {
  print(i)
  ArchR <- loadArchRProject(path = stages_data[i], force = TRUE, showLogo = FALSE)
  cell_ids <- ArchR$cellNames
  print(length(cell_ids))
  all_cell_ids <- c(all_cell_ids, cell_ids)
}
print(paste0("Length of all stage cell ids: ", length(all_cell_ids)))

# check cell counts and IDs match between full data and stages
intersect_cell_id_length <- sum(all_cell_ids %in% rownames(ArchR_full@cellColData))
print(paste0("Total number of stage cell ids found in full data: ", intersect_cell_id_length))
if (intersect_cell_id_length != length(all_cell_ids)){
  stop("Error! Not all stage cell ids match with cell ids in Full data")
}

# filter full data based on cell IDs from filtered stages
ArchR_filtered <- ArchR_full[all_cell_ids, ]

# save filtered data
saveArchRProject(ArchRProj = ArchR_filtered, outputDirectory = paste0(rds_path, "FullData_Save-ArchR"), load = FALSE)

#####################################################################################
############################## Visualisations #######################################

### Plot cell counts before and after subsetting per stage
unfiltered <- table(ArchR_full$stage)
filtered <- table(ArchR_filtered$stage)
cell_counts <- as_tibble(rbind(unfiltered, filtered))
cell_counts <- cbind(cell_counts, Total = rowSums(cell_counts))

png(paste0(plot_path, 'cell_counts_table_stages.png'), height = 10, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()