#!/usr/bin/env Rscript

print("inputs: 1) SEACells ATAC seurat object, 2) csv mapping single cell IDs to SEACell IDs, 3) csv mapping SEACell IDs to integrate cell type label")
print("script uses csvs to 1) make consensus RNA-ATAC SEACell mapping -> csv, 2) visualise these new SEACell labels on ATAC SEACell seurat object and 3) maps this SEACell map for ATAC onto ATAC single cell data -> csv")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(parallel)
library(Seurat)

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
    
    ### Interactively will have to read from a few different places
    data_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_ATAC_WF/5_Classify_metacells/" # for seurat object
    data_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_ATAC_WF/2_SEACells_computation/exported_data/" # for single cell to seacell mapping
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Integrated_SEACells/" # integrated seacell mapping

    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    rds_path = "./rds_files/"
    plot_path = "./plots/"
    data_path = "./input/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(rds_path, recursive = T)
  dir.create(plot_path, recursive = T)
}

########################       CELL STATE COLOURS    ########################################
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

############################## Read data #######################################

# List all files in input
input_paths <- list.files(data_path)
print(input_paths)

# 1) Read in seurat object
seurat_path <- list.files(path = paste0(data_path, "rds_files/"), pattern = ".RDS", full.names = TRUE)
print(paste0("Seurat path: ", seurat_path))
seurat <- readRDS(seurat_path)
print(seurat)

# 2) Read in ATAC single cell to SEACell mapping
SEACell_map_path <- list.files(path = paste0(data_path, "csv_files/"), pattern = "*_cell_metadata.csv", full.names = TRUE)
#SEACell_map_path <- list.files(path = paste0(data_path, ""), pattern = "*Cell_metadata.csv", full.names = TRUE) #interactive
print(paste0("SEACell map path: ", SEACell_map_path))

SEACell_map <- read.csv(SEACell_map_path)
print(head(SEACell_map))
print(dim(SEACell_map))

# 3) Read in ATAC SEACell to RNA SEACell mapping
Integration_map_path <- list.files(path = paste0(data_path, "integrated_output/"), pattern = "*_mappings_cell_type", full.names = TRUE)
print(paste0("Integration map path: ", Integration_map_path))

Integration_maps <- lapply(Integration_map_path, read.csv)
print(paste0("Number of integration maps (ie k values run): ", length(Integration_maps)))

############################## 1) Generate consensus SEACell integration map RNA-ATAC #######################################

# Initiate combined integration map with mapping from k=1
combined_integration_map <- Integration_maps[[1]]

# Check that all ATAC IDs in this table are unique
if (length(unique(combined_integration_map$ATAC)) == dim(combined_integration_map)[1]){
  print("All ATAC SEACell IDs are unique in k=1 mapping!") }else{
    stop("Not all ATAC SEACell IDs are unique in k=1 mapping!")
  }

# Clean up df
combined_integration_map <- combined_integration_map[, -1]
combined_integration_map <- combined_integration_map %>% mutate(k = 1)
print("Preview of integration map of k=1:")
head(combined_integration_map)

for (i in 2:length(Integration_maps)){
  print(i)
  
  # Get mapping from this k value and clean up df
  map <- Integration_maps[[i]]
  map <- map[, -1]
  map <- map %>% mutate(k = i)
  
  # remove any matches to ATAC cells which have already been matched
  unique_map <- map %>% 
    filter(!(ATAC %in% combined_integration_map$ATAC))
  print(head(unique_map))
  print(dim(unique_map))
  
  # add any new matches to combined df
  combined_integration_map <- rbind(combined_integration_map, unique_map)
  
}

dim(combined_integration_map)
head(combined_integration_map)
mapped <- unique(combined_integration_map$ATAC) %in% unique(SEACell_map$SEACell)
print(paste0("Number of ATAC SEACell IDs mapped to an RNA SEACell ID: ", length(mapped)))
print(paste0("Number of ATAC SEACell IDs NOT mapped to an RNA SEACell ID: ", length(unique(SEACell_map$SEACell)) - length(mapped)))

## dealing with duplicate ATAC IDs by:
# identify which ones map to > 1 scHelper cell type, remove these ones (v few)
# for the rest just remove duplicates as the scHelper_cell_type will be the same even if the RNA SEACell ID is different
duplicated_ATAC_IDs <- combined_integration_map$ATAC[duplicated(combined_integration_map$ATAC)]
print(paste0("Number of duplicated ATAC SEACell IDs in integration maps: ", length(duplicated_ATAC_IDs)))

duplicates_map <- combined_integration_map %>% filter(ATAC %in% duplicated_ATAC_IDs) %>%
  arrange(ATAC) %>%
  group_by(ATAC, scHelper_cell_type) %>% 
  dplyr::mutate(duplicated_cell_type = n()>1)
SEACells_to_remove <- duplicates_map[which(duplicates_map$duplicated_cell_type == FALSE), ]$ATAC
  
print(paste0("Number of ATAC SEACell duplicates that map to more than one cell type: ", length(SEACells_to_remove)))
print("Removing these SEACells from mapping!")

filtered_integration_map <- combined_integration_map %>% 
  dplyr::filter(!ATAC %in% SEACells_to_remove) %>% 
  distinct(ATAC, .keep_all = TRUE)

dim(filtered_integration_map)
head(filtered_integration_map)

## Extract stage from cell metadata to use for naming file
stage <- substr(SEACell_map$index[1], 8, 10)
print(paste0("Stage detected: ", stage))

write.csv(filtered_integration_map, paste0(rds_path, stage, '_filtered_SEACells_integration_map.csv'))

############################## 2) Visualise integration result on ATAC SEACell seurat #######################################

## add new transferred labels to seurat object
map1 <- filtered_integration_map %>% arrange(ATAC)
colnames(map1)[colnames(map1) == 'scHelper_cell_type'] <- 'scHelper_cell_type_integration'
map2 <- rownames_to_column(seurat@meta.data, var = "ATAC")
map <- merge(map1, map2, by = "ATAC", all = TRUE)
metadata <- column_to_rownames(map, var = "ATAC")
head(metadata)

seurat <- AddMetaData(seurat, metadata = metadata$RNA, col.name = "Integrated_RNA_SEACell_ID")
seurat <- AddMetaData(seurat, metadata = metadata$scHelper_cell_type_integration, col.name = "scHelper_cell_type_from_integration")
seurat <- AddMetaData(seurat, metadata = metadata$k, col.name = "Mapping_k")

## save seacells seurat object with new metadata
saveRDS(seurat, paste0(rds_path, stage, "_seacells_seurat_integrated.RDS"), compress = FALSE)

## plot new metadata on SEACell UMAPs
png(paste0(plot_path, "stage_UMAP.png"), width=25, height=20, units = 'cm', res = 200)
DimPlot(seurat, group.by = "stage", pt.size = 10)
graphics.off()

scHelper_cols <- scHelper_cell_type_colours[levels(droplevels(seurat_data@meta.data$scHelper_cell_type_from_integration))]

png(paste0(plot_path, "scHelper_cell_type_from_integration_UMAP.png"), width=25, height=20, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = 'scHelper_cell_type_from_integration', label = TRUE, 
        label.size = 9, label.box = TRUE, repel = TRUE,
        pt.size = 10,
        cols = scHelper_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

png(paste0(plot_path, "mapping_k_UMAP.png"), width=25, height=20, units = 'cm', res = 200)
FeaturePlot(seurat, features = "Mapping_k", pt.size = 10)
graphics.off()


############################## 3) Generate ATAC SEACell to ATAC single cell map #######################################

df_new <- merge(SEACell_map, filtered_integration_map, by.x = "SEACell", by.y = "ATAC", all.x = TRUE)
single_cell_integrated_map <- df_new %>% select(c("SEACell", "index", "RNA", "scHelper_cell_type"))

head(single_cell_integrated_map)

write.csv(single_cell_integrated_map, paste0(rds_path, stage, '_ATAC_singlecell_integration_map.csv'))
