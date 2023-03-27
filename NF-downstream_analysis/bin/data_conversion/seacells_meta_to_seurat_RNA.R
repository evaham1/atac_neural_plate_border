#!/usr/bin/env Rscript

print("script to do two things:")
print("1) transfer seacell ids to seurat object")
print("2) Create a seurat object summarised by seacells with summarised metadata")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(parallel)
library(Seurat)
library(dplyr)
library(tibble)
library(scHelper)
library(ggplot2)
library(future)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-m", "--metadata_file_name"), action = "store", type = "character", help = "Name of csv file which assigns cell ids to metacell ids", default = "exported_data/Cell_metadata.csv"),
  make_option(c("-t", "--categories"), action = "store", type = "character", help = "Which categories to use to check for purity", default = "run,sex,stage,seurat_clusters,scHelper_cell_type"),
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
    data_path = "./local_test_data/test_inputs/test_input_seacells_meta_to_seurat/"
    rds_path = "./local_test_data/convert_seacells_to_seurat/"
    plot_path = "./local_test_data/convert_seacells_to_seurat//plots/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


################### Functions ##########################

## function to aggregate matrix from seurat object and summarise by cell groupings
summarise_seurat_data <- function(seurat, data_slot = "counts", category = "SEACell"){
  
  # extract data into dataframe
  df <- GetAssayData(object = seurat, slot = data_slot)
  df <- as.data.frame(t(as.data.frame(df)))
  
  # convert cell ids to category ids
  category_ids <- select(seurat@meta.data, category)[,1]
  df <- df %>% mutate(category = category_ids)
  
  # aggregate df based on category
  df_summarised <- aggregate(. ~ category, df, sum)
  
  # format df so can be added back to seurat object
  df_summarised <- t(column_to_rownames(df_summarised, var = "category"))
  
  return(df_summarised)
}

# function to take seurat object and a category and make a freq table of how frequently metacells found in each category
calculate_metacell_frequencies <- function(seurat, metacell_slot = "SEACell", category = "stage"){
  
  df <- data.frame(FetchData(object = seurat, vars = c(metacell_slot, category)))
  colnames(df) <- c("Metacell", "Category")
  
  freq <- df %>%
    group_by(Metacell, Category) %>% 
    dplyr::summarize(count = n()) %>%
    mutate(Metacell = str_split(Metacell, "-", simplify = TRUE)[ , 2]) %>%
    mutate(Metacell = as.numeric(Metacell)) %>% 
    arrange(Metacell) %>% 
    mutate(count = as.numeric(count))
  
  print(paste0("Number of metacells: ", length(unique(freq$Metacell))))
  print(paste0("Number of categories: ", length(unique(freq$Category))))
  
  return(freq)
}

## Function to take freq table and turn it into proportions per metacell
calculate_metacell_proportions <- function(freq_table){
  
  # calculate total cell counts per metacell
  totals_df <- aggregate(count ~ Metacell, data=freq_table, FUN=sum)
  totals <- totals_df$count
  print(length(totals))
  
  # calculate proportions per metacell
  prop_table <- freq_table %>%
    mutate(prop = count/totals[Metacell+1])
  
  return(prop_table)
}

## Function to plot piechart of how many metacells pass threshold for proportion of cells coming from one label
piechart_proportion_threshold <- function(prop_table, threshold = 0.5){
  
  # filter cells to only include those that pass threshold
  passed_cells <- prop_table %>% filter(prop > threshold)
  # number of cells that pass threshold:
  passed_cells <- length(unique(passed_cells$Metacell))
  # number of cells that didn't pass threshold:
  failed_cells <- length(unique(prop_table$Metacell)) - passed_cells
  # plot piechart
  slices <- c(passed_cells, failed_cells)
  lbls <- c(paste0("Passed: ", passed_cells), paste0("Didn't pass: ", failed_cells))
  
  return(pie(slices, labels = lbls, main = paste0("Number of cells that passed threshold (", threshold, ") of label proportions")))
  
}


#######################################################################################
######################    Add SEACells IDs to seurat object   #########################
#######################################################################################

# Read in metadata (use metadata_file_name to find correct object)
metacell_metadata <- read.csv(paste0(data_path, opt$metadata_file_name))
metacell_dictionary <- select(metacell_metadata, c("index", "SEACell"))

print("Metacell IDs read in")
print(head(metacell_dictionary))

# Read in seurat object (identify label first, file needs to be named as [LABEL]_clustered_data.RDS, only looks in input not in /rds_files/)
label <- setdiff(sub('_.*', '', list.files(data_path)), c("rds", "exported", "plots"))
print(label)
seurat <- readRDS(paste0(data_path, label, "_clustered_data.RDS"))

print("Seurat object read in")
print(seurat)

# Reorder seacells metadata to match cell order in seurat object
metacell_dictionary <- metacell_dictionary[match(rownames(seurat@meta.data), metacell_dictionary$index),]

print(head(metacell_dictionary))

# Add seacells metadata to seurat object
seurat <- AddMetaData(seurat, metacell_dictionary$SEACell, col.name = "SEACell")

# Save seurat object
saveRDS(seurat, paste0(rds_path, "seurat.RDS"), compress = FALSE)

#####################################################################################
######################    Create summarised seurat object   #########################
#####################################################################################

############################## Explore individual seacell purity #######################################

categories <- strsplit(opt$categories, ",")[[1]]

## loop through categories and check each of them for purity in metacells
for (cat in categories) {
  
  print(paste0("Checking purity of ", cat))
  if (!(cat %in% colnames(seurat@meta.data))){
    stop("Category not in seurat cell col data!")}
  
  plot_path_temp = paste0(plot_path, cat, "/")
  dir.create(plot_path_temp, recursive = T)
  
  # calculate frequencies
  freq_table <- calculate_metacell_frequencies(seurat, metacell_slot = "SEACell", category = cat)
  
  # calculate proportions of labels in metacells
  prop_table <- calculate_metacell_proportions(freq_table)
  
  # plot the relative proportions of labels in each metacell
  png(paste0(plot_path_temp, "Hist_all_proportions.png"), width=25, height=20, units = 'cm', res = 200)
  hist(prop_table$prop)
  graphics.off()
  
  # plot the max relative proportions of labels in each metacell
  max_prop_table <- prop_table %>% group_by(Metacell) %>% dplyr::summarise(prop = max(prop))
  png(paste0(plot_path_temp, "Hist_max_proportions_per_metacell.png"), width=25, height=20, units = 'cm', res = 200)
  hist(max_prop_table$prop)
  graphics.off()
  
  ## how many metacells have > 40% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.4.png"), width=25, height=20, units = 'cm', res = 200)
  piechart_proportion_threshold(prop_table, threshold = 0.4)
  graphics.off()
  
  ## how many metacells have > 50% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.5.png"), width=25, height=20, units = 'cm', res = 200)
  piechart_proportion_threshold(prop_table, threshold = 0.5)
  graphics.off()
  
  ## how many metacells have > 75% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.75.png"), width=25, height=20, units = 'cm', res = 200)
  piechart_proportion_threshold(prop_table, threshold = 0.75)
  graphics.off()
  
  ## how many metacells have > 90% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.9.png"), width=25, height=20, units = 'cm', res = 200)
  piechart_proportion_threshold(prop_table, threshold = 0.9)
  graphics.off()
  
}

#################### Add up counts across metacells #########################

###### RNA slot: 3 assays: counts (raw), data (normalised), scale.data -> only add up raw 'counts'
DefaultAssay(object = seurat) <- "RNA"
DefaultAssay(object = seurat)
summarised_RNA_counts <- summarise_seurat_data(seurat = seurat, data_slot = "counts", category = "SEACell")

dim(summarised_RNA_counts) # 18683   215
sum(is.na(summarised_RNA_counts)) # 0 NA values

print("raw counts summarised")

###### Integrated slot: 2 assays: data, scale.data -> only add up 'data'
DefaultAssay(object = seurat) <- "integrated"
DefaultAssay(object = seurat)
summarised_integrated_data <- summarise_seurat_data(seurat = seurat, data_slot = "data", category = "SEACell")

dim(summarised_integrated_data) # 18683   215
sum(is.na(summarised_integrated_data)) # 0 NA values

print("integrated counts summarised")

#################### Create cell metadata for SEACells #########################

### For 'sex' and 'run' there are only two options, so for seacells can express these metadata values as a proportion #
### For 'stage' this can just be carried over from any cell as should always be 100% for each metacell

seacell_ids_stripped <- sort(as.numeric(str_split(colnames(summarised_RNA_counts), "-", simplify = TRUE)[ , 2]))
temp_metadata <- data.frame(seacell_ids_stripped)
colnames(temp_metadata) <- "Metacell"

# Sex (proportion cells in metacell which are MALE)
sex_freq_table <- calculate_metacell_frequencies(seurat, metacell_slot = "SEACell", category = "sex")
sex_prop_table <- calculate_metacell_proportions(sex_freq_table)
sex_male_prop <- sex_prop_table %>%
  subset(Category == "male") %>%
  select(c(Metacell, prop))

temp_metadata_sex <- merge(data.frame(sex_male_prop), temp_metadata, by = "Metacell", all = TRUE)
table(is.na(temp_metadata_sex))
temp_metadata_sex <- temp_metadata_sex %>% replace(is.na(.), 0)
colnames(temp_metadata_sex) <- c("Metacell", "sex")

# Run (proportion cells in metacell which are from RUN 1)
run_freq_table <- calculate_metacell_frequencies(seurat, metacell_slot = "SEACell", category = "run")
run_prop_table <- calculate_metacell_proportions(run_freq_table)
run_1_prop <- run_prop_table %>%
  subset(Category == "1") %>%
  select(c(Metacell, prop))

temp_metadata_run <- merge(data.frame(run_1_prop), temp_metadata_sex, by = "Metacell", all = TRUE)
table(is.na(temp_metadata_run))
temp_metadata_run <- temp_metadata_run %>% replace(is.na(.), 0)
colnames(temp_metadata_run) <- gsub("prop", "run", colnames(temp_metadata_run))

# Stage (just carry over)
stage_freq_table <- calculate_metacell_frequencies(seurat, metacell_slot = "SEACell", category = "stage")
stage_prop_table <- calculate_metacell_proportions(stage_freq_table)
stage_prop_table <- stage_prop_table %>%
  select(c(Metacell, Category))

temp_metadata_stage <- merge(data.frame(stage_prop_table), temp_metadata_run, by = "Metacell", all = TRUE)
colnames(temp_metadata_stage) <- gsub("Category", "stage", colnames(temp_metadata_stage))

# Turn stripped IDs back into normal SEACell IDs
seacells_seurat_metadata <- temp_metadata_stage
seacells_seurat_metadata$Metacell <- sub("^", "SEACell-", seacells_seurat_metadata$Metacell)
seacells_seurat_metadata <- column_to_rownames(seacells_seurat_metadata, var = "Metacell")

head(seacells_seurat_metadata)

#################### Create new seurat object #########################

# Create object using summarised RNA counts and newly created cell metadata
seacells_seurat <- CreateSeuratObject(
  project = "seacells_seurat",
  counts = summarised_RNA_counts,
  assay = "RNA",
  names.field = 2,
  names.delim = "-",
  meta.data = seacells_seurat_metadata
)

# Add summarised integrated counts to seurat object
seacells_seurat[["integrated"]] <- CreateAssayObject(counts = summarised_integrated_data)

# Validate that the object now contains multiple assays
Assays(seacells_seurat)

# Check metadata
head(seacells_seurat@meta.data)

# Print object
print(seacells_seurat)

print("created summarised seurat object!")

#################### Save new seurat object #########################

## save seacells seurat object
saveRDS(seacells_seurat, paste0(rds_path, "seacells_seurat.RDS"), compress = FALSE)