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
    data_path = "./local_test_data/test_input_seacells_meta_to_seurat/"
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
dir.create("/rds_files_full/", recursive = T)
saveRDS(seurat, paste0("/rds_files_full/", "seurat.RDS"), compress = FALSE)

## Plot number of cells and number of genes in original object
df <- data.frame(dim(seurat))
rownames(df) <- c("Gene count: ", "Cell count: ")
png(paste0(plot_path, 'original_cell_counts.png'), height = 5, width = 12, units = 'cm', res = 400)
grid.arrange(top=textGrob("Gene count and cell count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

#####################################################################################
######################    Create summarised seurat object   #########################
#####################################################################################

############################## Explore individual seacell purity #######################################

categories <- strsplit(opt$categories, ",")[[1]]

print("All seurat metadata:")
print(colnames(seurat@meta.data))

## loop through categories and check each of them for purity in metacells
for (cat in categories) {
  
  print(paste0("Checking purity of ", cat))
  if (!(cat %in% colnames(seurat@meta.data))){
    stop("Category not in seurat cell col data!")}
  
  plot_path_temp = paste0(plot_path, cat, "/")
  dir.create(plot_path_temp, recursive = T)
  
  # calculate frequencies
  prop_table <- SEACells_MetacellFrequencies(seurat, metacell_slot = "SEACell", category = cat, calc_proportions = TRUE)
  
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
  SEACells_PiechartProportionThreshold(prop_table, threshold = 0.4)
  graphics.off()
  
  ## how many metacells have > 50% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.5.png"), width=25, height=20, units = 'cm', res = 200)
  SEACells_PiechartProportionThreshold(prop_table, threshold = 0.5)
  graphics.off()
  
  ## how many metacells have > 75% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.75.png"), width=25, height=20, units = 'cm', res = 200)
  SEACells_PiechartProportionThreshold(prop_table, threshold = 0.75)
  graphics.off()
  
  ## how many metacells have > 90% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.9.png"), width=25, height=20, units = 'cm', res = 200)
  SEACells_PiechartProportionThreshold(prop_table, threshold = 0.9)
  graphics.off()
  
}

#################### Add up counts across metacells #########################

###### RNA slot: 3 assays: counts (raw), data (normalised), scale.data -> only add up raw 'counts'
DefaultAssay(object = seurat) <- "RNA"
DefaultAssay(object = seurat)
summarised_RNA_counts <- SEACells_SummariseSeuratData(seurat = seurat, data_slot = "counts", category = "SEACell")

dim(summarised_RNA_counts) # 18683   215
sum(is.na(summarised_RNA_counts)) # 0 NA values

print("raw counts summarised")

###### Integrated slot: 2 assays: data, scale.data -> only add up 'data'
DefaultAssay(object = seurat) <- "integrated"
DefaultAssay(object = seurat)
summarised_integrated_data <- SEACells_SummariseSeuratData(seurat = seurat, data_slot = "data", category = "SEACell")

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
sex_prop_table <- SEACells_MetacellFrequencies(seurat, metacell_slot = "SEACell", category = "sex", calc_proportions = TRUE)
sex_male_prop <- sex_prop_table %>%
  subset(Category == "male") %>%
  select(c(Metacell, prop))

temp_metadata_sex <- merge(data.frame(sex_male_prop), temp_metadata, by = "Metacell", all = TRUE)
table(is.na(temp_metadata_sex))
temp_metadata_sex <- temp_metadata_sex %>% replace(is.na(.), 0)
colnames(temp_metadata_sex) <- c("Metacell", "sex")

# Run (proportion cells in metacell which are from RUN 1)
run_prop_table <- SEACells_MetacellFrequencies(seurat, metacell_slot = "SEACell", category = "run", calc_proportions = TRUE)
run_1_prop <- run_prop_table %>%
  subset(Category == "1") %>%
  select(c(Metacell, prop))

temp_metadata_run <- merge(data.frame(run_1_prop), temp_metadata_sex, by = "Metacell", all = TRUE)
table(is.na(temp_metadata_run))
temp_metadata_run <- temp_metadata_run %>% replace(is.na(.), 0)
colnames(temp_metadata_run) <- gsub("prop", "run", colnames(temp_metadata_run))

# Stage (just carry over)
stage_prop_table <- SEACells_MetacellFrequencies(seurat, metacell_slot = "SEACell", category = "stage", calc_proportions = TRUE)
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

## Plot number of cells and number of genes in summarised metacells object
df <- data.frame(dim(seacells_seurat))
rownames(df) <- c("Gene count: ", "SEACell count: ")
png(paste0(plot_path, 'SEACell_counts.png'), height = 5, width = 12, units = 'cm', res = 400)
grid.arrange(top=textGrob("Gene count and SEAcell count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

#################### Save new seurat object #########################

## save seacells seurat object
saveRDS(seacells_seurat, paste0(rds_path, "seacells_seurat.RDS"), compress = FALSE)