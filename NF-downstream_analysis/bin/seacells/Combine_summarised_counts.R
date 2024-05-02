#!/usr/bin/env Rscript

print("take .csv files from SEACell computation for all stages (HH5_cell_metadata.csv, HH5_feature_metadata.csv and HH5_summarised_by_metacell_counts.csv)")
print("and combine them into one dataframe")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(tidyverse)
library(plyr)
library(dplyr)
library(parallel)
library(data.table)

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
    
    # interactive paths
    data_path = "./output/temp_renamed_seacell_outputs/" # made this folder and copied across inputs manually
    data_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_ATAC_WF/2_SEACells_computed_renamed/csv_files/" # just for renamed ATAC SEACell outputs
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Integrated_SEACells_label_transfer/rds_files/" # just for the integrated outputs
    rds_path = "./output/NF-downstream_analysis/Processing/FullData/Combined_SEACell_outputs/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    rds_path = "./csv_files/"
    data_path = "./input/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(rds_path, recursive = T)
}



############################## FUNCTIONS #######################################

# check if list of dfs are identical
all_identical <- function(df_list) {
  len <- length(df_list)
  for(i in 2:len){
    if(!identical(df_list[[i]], df_list[[i-1]])){
      return(FALSE)
    }
  }
  return(TRUE)
}

############################## Identify .csv files and list of stages #######################################
input_files <- list.files(path = data_path, pattern = "*.csv", full.names = TRUE)
print(paste0("Input paths detected: ", input_files))


############################## Check feature metadata files are the same and write out one of them if they are #######################################

print("Checking feature metadata...")

# Extract paths
feature_metadata_paths <- list.files(path = data_path, pattern = "*_feature_metadata.csv", full.names = TRUE)
print(paste0("Feature metadata paths: ", feature_metadata_paths))

# Check all feature files are identical
feature_metadata_files <- lapply(feature_metadata_paths, read.csv)
ifelse (all_identical(feature_metadata_files), print("All feature metadata feature files identical!"), stop("Not all metadata feature files are identical!"))

# Extract one feature file and format it
feature_metadata <- feature_metadata_files[[1]]
feature_metadata <- feature_metadata[,-1]
print("Checking if any peak IDs are duplicated:")
table(duplicated(feature_metadata$X))

# Check
print("Preview of feature metadata:")
feature_metadata[1:5, 1:5]

# Save
write.csv(feature_metadata, paste0(rds_path, 'Feature_metadata.csv'))

print("Feature metadata checked and saved!")

############################## Combine cell metadata files #######################################

print("Checking cell metadata...")

# Extract paths
cell_metadata_paths <- list.files(path = data_path, pattern = "*_cell_metadata.csv", full.names = TRUE)
print(paste0("Cell metadata paths: ", cell_metadata_paths))

# Read in 
cell_metadata_files <- lapply(cell_metadata_paths, read.csv)

# Add another column to cell metadata with seacell ID-stage eg SEACell-84-ss8
cell_metadata_files <- lapply(cell_metadata_files, function(x) mutate(x, SEACell_stage = paste0(SEACell, "-", stage)))

# Combine all cell metacell csvs into one
combined_cell_metacell <- ldply(cell_metadata_files, data.frame)
combined_cell_metacell <- combined_cell_metacell[,-1]

# Check
print("Preview of cell metadata:")
combined_cell_metacell[1:5,]

# Save
write.csv(combined_cell_metacell, paste0(rds_path, 'Combined_cell_metadata.csv'))

print("Cell metadata files combined!")

############################## Combine summarised peak counts #######################################

print("Checking summarised peak counts...")

# Extract paths
count_paths <- list.files(path = data_path, pattern = "*_summarised_by_metacells_counts.csv", full.names = TRUE)
print(paste0("Summarised count paths: ", count_paths))

count_files <- lapply(count_paths, fread)

# Extract stages in same order than count files are in (hopefully they don't swap...)
count_file_names <- sub(".*/", "", count_paths)
stages <- lapply(count_file_names, function(x) substr(x, 1, 3))
print(stages)

## Edit seacell ID to be seacell ID-stage eg SEACell-84-ss8
for (i in 1:length(stages)){
  stage <- stages[i]
  data <- count_files[[i]][,-1]
  
  print("preview of data before formatting:")
  print(data[1:4, 1:4])
  
  data <- data %>% mutate(SEACell = paste0(V1, "-", stage))
  data <- data[, -1]
  data <- data %>% select(SEACell, everything())
  
  print("preview of data after formatting:")
  print(data[1:4, 1:4])
  
  assign(paste0(stage, "_data"), data)
}

# Should all have the same number of columns - can have variable row numbers
print("Checking dimensions of resulting dfs: ")
print(dim(HH5_data))
print(dim(HH6_data))
print(dim(HH7_data))
print(dim(ss4_data))
print(dim(ss8_data))

# Should all have the same number of columns - can have variable row numbers
print("Checking format of resulting dfs: ")
print(HH5_data[1:2, 1:2])
print(HH6_data[1:2, 1:2])
print(HH7_data[1:2, 1:2])
print(ss4_data[1:2, 1:2])
print(ss8_data[1:2, 1:2])

# Combine all into one
combined_df <- rbindlist(list(HH5_data, HH6_data, HH7_data, ss4_data, ss8_data))

## write out new csv
write.csv(combined_df, paste0(rds_path, 'Combined_summarised_by_metacells_counts.csv'))

print("Summarised peak counts combined!")

############################## Combine SEACell integrated metadata files #######################################

print("Checking SEACell integrated metadata...")

# Extract paths
SEACell_integrated_metadata_paths <- list.files(path = data_path, pattern = "*_SEACells_integration_map.csv", full.names = TRUE)
print(paste0("SEACell integrated metadata paths: ", SEACell_integrated_metadata_paths))

SEACell_integrated_metadata_files <- list()

for (i in 1:length(SEACell_integrated_metadata_paths)){
  
  print(i)
  
  # extract correct path
  path <- SEACell_integrated_metadata_paths[i]
  
  # detect stage from file name
  filename <- sub(".*/", "", path)
  stage <- sub("_.*", "", filename)
  print(paste0("Stage detected from file name: ", stage))
  
  # add stage to SEACell IDs as otherwise they wont be unique when you combine the stages together
  file <- read.csv(path)
  file <- file %>% 
    mutate(ATAC = paste0(ATAC, "-", stage)) %>%
    mutate(RNA = paste0(RNA, "-", stage))
  
  # remove first column
  file <- file[,-1]
  
  # add altered df to list of dfs
  SEACell_integrated_metadata_files[[i]] <- file
  
}

## check new integrated metadata
print(length(SEACell_integrated_metadata_files))
print("Preview of each df in list of integrated SEACell metadata dfs:")
print(SEACell_integrated_metadata_files[[1]][1:5, ])
print(SEACell_integrated_metadata_files[[2]][1:5, ])
print(SEACell_integrated_metadata_files[[3]][1:5, ])
print(SEACell_integrated_metadata_files[[4]][1:5, ])
print(SEACell_integrated_metadata_files[[5]][1:5, ])

## combine all cell metacell csvs into one
combined_SEACell_metacell <- ldply(SEACell_integrated_metadata_files, data.frame)
combined_SEACell_metacell <- combined_SEACell_metacell[,-1]
dim(combined_SEACell_metacell)
print(combined_SEACell_metacell[1:5, ])

## write out new csv
write.csv(combined_SEACell_metacell, paste0(rds_path, 'Combined_SEACell_integrated_metadata.csv'))

print("Cell integrated metadata files combined!")