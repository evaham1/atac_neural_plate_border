#!/usr/bin/env Rscript 

print("inputs: 1) csv mapping SEACell IDs to integrate cell type label, 2) csv mapping single cell ATAC IDs to SEACEll IDs, 3) ATAC SEACell seurat object")
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
library(gridExtra)
library(grid)
library(tidyverse)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-k", "--k_cutoff"), action = "store", type = "integer", help = "integration param k cutoff to use", default = 3),
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
    
    ### Interactively will have to read from a few different places
    # 1) ATAC SEACell to RNA SEACell mapping
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Integrated_SEACells/"
    # 2) ATAC SEACell to single cell mapping
    data_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_ATAC_WF/2_SEACells_computation/exported_data/"
    # 3) ATAC SEACell seurat object
    data_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_ATAC_WF/5_Classify_metacells/"
    
    
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

set.seed(42)

########################       CELL STATE COLOURS    ########################################
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

###### STAGE COLOURS ####################
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

###### K-MAPPING COLOURS ####################
# k_colors <- c("#f21111", "#ef3e2a", "#ea573f", "#e56a54", "#de7c69", "#d58c7f", "#c99b94", "#bbaaab", "a8b8c1", "#EAEAEA")
# names(k_colors) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "Unmapped")

k_colors <- c("#f21111", "#f63535", "#f95959", "#fd7d7d", "#ffa1a1", "#ffc5c5", "#FFDADA", "#ffe0e0", "#ffe5e5", "#ffebec", "#fff0f0", "#fff5f5", "#EAEAEA")
names(k_colors) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "Unmapped")

############################## Read data #######################################

# List all files in input
input_paths <- list.files(data_path)
print(input_paths)

# 1) Read in ATAC SEACell to RNA SEACell mapping
Integration_map_path <- list.files(path = paste0(data_path, "integrated_output/"), pattern = "*_mappings_cell_type", full.names = TRUE)
print(paste0("Integration map path: ", Integration_map_path))

Integration_maps <- lapply(Integration_map_path, read.csv)
print(paste0("Number of integration maps (ie k values run): ", length(Integration_maps)))

# 2) Read in ATAC single cell to SEACell mapping
SEACell_map_path <- list.files(path = paste0(data_path, "csv_files/"), pattern = "*_cell_metadata.csv", full.names = TRUE)
#SEACell_map_path <- list.files(path = paste0(data_path, ""), pattern = "*Cell_metadata.csv", full.names = TRUE) #interactive
print(paste0("SEACell map path: ", SEACell_map_path))

SEACell_map <- read.csv(SEACell_map_path)
print(head(SEACell_map))
print(dim(SEACell_map))

# 3) Read in SEACell ATAC seurat object
seurat_path <- list.files(path = paste0(data_path, "rds_files/"), pattern = ".RDS", full.names = TRUE)
print(paste0("Seurat path: ", seurat_path))
seurat <- readRDS(seurat_path)
print(seurat)


#############################################################################################################################
#####################            1) Generate consensus SEACell integration map RNA-ATAC            ##########################
#############################################################################################################################

plot_path = "./plots/generating_consenus_integration_map/"
dir.create(plot_path, recursive = T)

############### 1.1) Clean up data ##############

# How many ATAC SEACells are there 
all_SEACells <- unique(SEACell_map$SEACell)
print(paste0("How many SEACells total: ", length(all_SEACells)))

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

############### 1.2) Make combined map by iterating through each k value ##############

# Initiate count of unique ATAC IDs mapped
labelled_cell_count <- c(dim(combined_integration_map)[1])

# Initiate count of duplicates in maps
duplicated_ATAC_IDs <- combined_integration_map$ATAC[duplicated(combined_integration_map$ATAC)]
duplicated_cell_count <- c(length(duplicated_ATAC_IDs))

# Iterate through each k value and add new SEACell mappings to combined_integration_map
for (i in 2:length(Integration_maps)){
  print(i)
  
  # Get mapping from this k value and clean up df
  map <- Integration_maps[[i]]
  map <- map[, -1]
  map <- map %>% mutate(k = i)
  
  # remove any matches to ATAC cells which have already been matched
  unique_map <- map %>% 
    filter(!(ATAC %in% combined_integration_map$ATAC))
  
  # count how many unique ATAC SEAcell IDs this iteration adds
  print(paste0("Number of unique ATAC SEACells additionally labelled with k=", i, ": ", length(unique(unique_map$ATAC))))
  labelled_cell_count <- c(labelled_cell_count, length(unique(unique_map$ATAC)))
  
  # count how many duplicated ATAC SEACell IDs have in this iteration
  duplicated_ATAC_IDs <- unique_map$ATAC[duplicated(unique_map$ATAC)]
  print(paste0("Number of duplicated ATAC SEACell IDs mapping with k=", i, ": ", length(duplicated_ATAC_IDs)))
  duplicated_cell_count <- c(duplicated_cell_count, length(duplicated_ATAC_IDs))
  
  # add any new matches to combined df
  combined_integration_map <- rbind(combined_integration_map, unique_map)
  
}

print("Preview of combined integration map:")
head(combined_integration_map)
dim(combined_integration_map)

############### 1.3) Set cut-off K value ##############

# use plots to see when using less stringent k value doesn't add signficiantly more maps
# png(paste0(plot_path, "13_k_values_plot.png"), width = 12, height = 12, units = 'cm', res = 200)
# plot(labelled_cell_count, ylim = c(0, length(all_SEACells) + 2), cex.axis = 1.5, cex.lab=2)
# abline(h = length(all_SEACells), col = "blue", cex = 1.5)
# abline(v = opt$k_cutoff, col = "red", cex = 1.5)
# graphics.off()

png(paste0(plot_path, "13_k_values_plot.png"), width = 12, height = 12, units = 'cm', res = 200)
plot(labelled_cell_count, cex.axis = 1.5, cex.lab=2)
abline(v = opt$k_cutoff, col = "red", cex = 1.5)
graphics.off()

df <- data.frame(k = c(1:length(labelled_cell_count)),
                 SEACell_count = labelled_cell_count)
png(paste0(plot_path, '13_k_values_table.png'), height = 18, width = 6, units = 'cm', res = 400)
grid.arrange(top=textGrob(" ", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

# use user-defined k value to remove data
cutoff_integration_map <- combined_integration_map %>% filter(k <= opt$k_cutoff)
dim(cutoff_integration_map)

# plot how many ATAC SEACell IDs have been removed from this filtering step
df <- data.frame(before_filter = length(unique(combined_integration_map$ATAC)),
                 after_filter = length(unique(cutoff_integration_map$ATAC)),
                 nFiltered = length(unique(combined_integration_map$ATAC)) - length(unique(cutoff_integration_map$ATAC))
                 )
png(paste0(plot_path, '13_k_values_filtering.png'), height = 5, width = 12, units = 'cm', res = 400)
grid.arrange(top=textGrob(" ", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

############### 1.4) Deal with single ATAC SEACells mapping to multiple RNA SEACells ##############

write.csv(cutoff_integration_map, "cutoff_integration_map.csv")

# initial number of ATAC seacells that are mapped
print(paste0("Number of ATAC SEACells that have been mapped: ", length(unique(cutoff_integration_map$ATAC))))

# which ATAC metacells map to more than one RNA metacell = multimapped_ATAC_IDs
multimapped_ATAC_IDs <- unique(cutoff_integration_map$ATAC[duplicated(cutoff_integration_map$ATAC)])
duplicates_map <- cutoff_integration_map %>% dplyr::filter(ATAC %in% multimapped_ATAC_IDs) %>%
  dplyr::arrange(ATAC) %>%
  dplyr::group_by(ATAC) %>%
  dplyr::mutate(duplicated_cell_type = n_distinct(scHelper_cell_type_by_proportion) == 1)
print(paste0("Number of ATAC SEACells that map to more than one RNA metacell: ", length(multimapped_ATAC_IDs)))

# which ATAC metacells map to more than one RNA metacell cell type = multimapped_cell_type_ATAC_IDs
multimapped_cell_type_ATAC_IDs <- unique(duplicates_map[which(duplicates_map$duplicated_cell_type == FALSE), ]$ATAC)
print(paste0("Number of ATAC SEACells that map to more than one cell type: ", length(multimapped_cell_type_ATAC_IDs)))

# which ATAC metacells map to more than one BROAD RNA metacell cell type = multimapped_broad_cell_type_ATAC_IDs
mapping <- c(
  "NP" = "Neural", "aNP" = "Neural", "iNP" = "Neural", "pNP" = "Neural", 
  "eN" = "Neural", "vFB" = "Neural", "FB" = "Neural", 
  "MB" = "Neural", "HB" = "Neural", "eCN" = "Neural", "eN" = "Neural",
  'PPR' = "Placodal", 'aPPR' = "Placodal", 'pPPR' = "Placodal",
  'eNPB' = "NPB", 'NPB' = "NPB", 'aNPB' = "NPB", 'pNPB' = "NPB",
  'NC' = "NC", 'dNC' = "NC",
  'NNE' = "Non-neural", 'pEpi' = "Non-neural",
  'EE' = "Contam", 'meso' = "Contam", 'endo' = "Contam", 'BI' = "Contam", 'PGC' = "Contam", 'streak' = 'Contam',
  'MIXED' = 'MIXED'
)
duplicates_map_broad <- duplicates_map %>% 
  dplyr::mutate(broad = case_when(
    scHelper_cell_type_by_proportion %in% names(mapping) ~ mapping[scHelper_cell_type_by_proportion],
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(ATAC %in% multimapped_cell_type_ATAC_IDs) %>%
  dplyr::group_by(ATAC) %>%
  dplyr::mutate(duplicated_cell_type = n_distinct(broad) == 1)

multimapped_broad_cell_type_ATAC_IDs <- unique(duplicates_map_broad[which(duplicates_map_broad$duplicated_cell_type == FALSE), ]$ATAC)
print(paste0("Number of ATAC SEACells that map to more than one BROAD cell type: ", length(multimapped_broad_cell_type_ATAC_IDs)))

# 1) resolve multimaps that dont matter as they are different RNA metacells but to the same cell state
filtered_integration_map <- cutoff_integration_map %>% 
  dplyr::filter(!ATAC %in% multimapped_cell_type_ATAC_IDs) %>% 
  dplyr::distinct(ATAC, .keep_all = TRUE)

if (length(unique(filtered_integration_map$ATAC)) ==  length(unique(cutoff_integration_map$ATAC)) - length(multimapped_cell_type_ATAC_IDs)){
  print("Length of filtered df correct!") } else {
    stop("ERROR: Length of filtered df is incorrect!")
  }

if (sum(duplicated(filtered_integration_map$ATAC)) == 0) {
  print("Filtered integration map has no duplicates!") } else {
    stop("ERROR: Filtered integration map has duplicates!")
  }

# 2) add in the multimaps which map to the same BROAD cell type
broad_labelling <- duplicates_map_broad %>%
  dplyr::filter(duplicated_cell_type == TRUE) %>%
  dplyr::mutate(scHelper_cell_type_by_proportion = broad) %>%
  dplyr::select(RNA, ATAC, scHelper_cell_type_by_proportion, k) %>% 
  dplyr::distinct(ATAC, .keep_all = TRUE)
length(unique(broad_labelling$ATAC))
filtered_integration_map <- base::rbind(filtered_integration_map, broad_labelling)

print(length(unique(cutoff_integration_map$ATAC)))
print(length(unique(cutoff_integration_map$ATAC)) - length(multimapped_broad_cell_type_ATAC_IDs))
print(length(unique(filtered_integration_map$ATAC)))

if (length(unique(filtered_integration_map$ATAC)) ==  length(unique(cutoff_integration_map$ATAC)) - length(multimapped_broad_cell_type_ATAC_IDs)){
  print("Length of filtered df correct!") } else {
    stop("ERROR: Length of filtered df is incorrect!")
  }

if (sum(duplicated(filtered_integration_map$ATAC)) == 0) {
  print("Filtered integration map has no duplicates!") } else {
    stop("ERROR: Filtered integration map has duplicates!")
  }

# plot how many seacells are multimapped to different levels
df <- data.frame(n_multimapped_SEACells = length(multimapped_ATAC_IDs),
                 n_multimapped_cell_type_SEACells = length(multimapped_cell_type_ATAC_IDs),
                 n_multimapped_cell_type_broad_SEACells = length(multimapped_broad_cell_type_ATAC_IDs)
)
png(paste0(plot_path, 'multimapped_numbers.png'), height = 5, width = 30, units = 'cm', res = 400)
grid.arrange(top=textGrob(" ", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

# plot how many ATAC SEACell IDs have been removed from this filtering step
df <- data.frame(before_filter = length(unique(cutoff_integration_map$ATAC)),
                 after_filter = length(unique(filtered_integration_map$ATAC)),
                 nFiltered = length(unique(cutoff_integration_map$ATAC)) - length(unique(filtered_integration_map$ATAC))
)
png(paste0(plot_path, '14_filtering.png'), height = 5, width = 12, units = 'cm', res = 400)
grid.arrange(top=textGrob(" ", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

############### 1.5) Plot some final stats ##############

# Plot how many ATAC SEACell IDs mapped at the end
mapped_SEACells <- all_SEACells[(all_SEACells %in% filtered_integration_map$ATAC)]
unmapped_SEACells <- all_SEACells[!(all_SEACells %in% filtered_integration_map$ATAC)]

df <- as.data.frame(c(length(mapped_SEACells), length(unmapped_SEACells), length(all_SEACells)))
rownames(df) <- c("ATAC IDs mapped to an RNA ID:", "ATAC IDs NOT mapped to an RNA ID:", "Total ATAC IDs")
colnames(df) <- "SEACell counts"
png(paste0(plot_path, '15_Filtered_how_many_metacells_mapped.png'), height = 8, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob(" ", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

# Plot how many ATAC SEACell IDs in final table were mapped using different k values
table <- as.data.frame(table(filtered_integration_map$k))
colnames(table) <- c("k value", "How many SEACells mapped")
png(paste0(plot_path, '15_Filtered_how_many_metacells_mapped_from_each_k.png'), height = 15, width = 15, units = 'cm', res = 400)
grid.arrange(top=textGrob(" ", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(table, theme = ttheme_minimal()))
graphics.off()

# Plot how many of each cell type there is
df <- as.data.frame(table(filtered_integration_map$scHelper_cell_type_by_proportion))
colnames(df) <- c("Identity", "Frequency")
png(paste0(plot_path, '15_Filtered_how_many_metacells_mapped_to_each_cell_state.png'), height = 15, width = 15, units = 'cm', res = 400)
grid.arrange(top=textGrob(" ", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

############### 1.6) Add rows with 'NA' for ATAC SEACell IDs which are not mapped ##############

## Add ATAC IDs to df which are NOT mapped
unmapped_df <- data.frame (RNA  = rep("Unmapped", length(unmapped_SEACells)),
                           ATAC = unmapped_SEACells,
                           scHelper_cell_type_by_proportion  = rep("Unmapped", length(unmapped_SEACells)),
                           k  = rep("Unmapped", length(unmapped_SEACells))
)
full_integration_map <- base::rbind(filtered_integration_map, unmapped_df)

# check final dataframe includes all SEACell ATAC IDs
if (length(all_SEACells) == sum(all_SEACells %in% full_integration_map$ATAC)) {
  print("Final integration map contains all SEACell ATAC IDs!") } else {
    stop("ERROR: Final integration map does not include all SEACell ATAC IDs!")
  }

# check final dataframe has no duplicates
if (sum(duplicated(full_integration_map$ATAC)) == 0) {
  print("Final integration map has no duplicates!") } else {
    stop("ERROR: Final integration map has duplicates!")
  }

############### 1.7) Add stage data to dataframe ##############

## Detect stage from cell metadata + add to df + use to name file
stage <- base::substr(SEACell_map$index[1], 8, 10)
print(paste0("Stage detected: ", stage))
full_integration_map <- full_integration_map %>% dplyr::mutate(Stage = stage)

############### 1.8) Save processed integration map ##############

write.csv(full_integration_map, paste0(rds_path, stage, '_SEACells_integration_map.csv'))


#############################################################################################################################
#####################            2) Visualise integration result on ATAC SEACell seurat           ##########################
#############################################################################################################################

print("Visualising label transfer results...")
plot_path = "./plots/seurat_visualise/"
dir.create(plot_path, recursive = T)

## add new transferred labels to seurat object
map1 <- full_integration_map %>% dplyr::arrange(ATAC)
colnames(map1)[colnames(map1) == 'scHelper_cell_type_by_proportion'] <- 'scHelper_cell_type_by_proportion'
map2 <- tibble::rownames_to_column(seurat@meta.data, var = "ATAC")
map <- base::merge(map1, map2, by = "ATAC", all = TRUE)
metadata <- tibble::column_to_rownames(map, var = "ATAC")
head(metadata)
print(sum(is.na(metadata)))

# add broad cell state labels
mapping <- c(
  "NP" = "Neural", "aNP" = "Neural", "iNP" = "Neural", "pNP" = "Neural", 
  "eN" = "Neural", "vFB" = "Neural", "FB" = "Neural", 
  "MB" = "Neural", "HB" = "Neural", "eCN" = "Neural", "eN" = "Neural",
  'PPR' = "Placodal", 'aPPR' = "Placodal", 'pPPR' = "Placodal",
  'eNPB' = "NPB", 'NPB' = "NPB", 'aNPB' = "NPB", 'pNPB' = "NPB",
  'NC' = "NC", 'dNC' = "NC",
  'NNE' = "Non-neural", 'pEpi' = "Non-neural",
  'EE' = "Contam", 'meso' = "Contam", 'endo' = "Contam", 'BI' = "Contam", 'PGC' = "Contam", 'streak' = 'Contam',
  'MIXED' = 'MIXED', 'Unmapped' = 'Unmapped',
  'Neural' = 'Neural', 'Placodal' = 'Placodal', 'NC' = 'NC', 'NPB' = 'NPB', 'Contam' = 'Contam'
)
metadata <- metadata %>% 
  dplyr::mutate(broad = case_when(
    scHelper_cell_type_by_proportion %in% names(mapping) ~ mapping[scHelper_cell_type_by_proportion],
    TRUE ~ NA_character_
  ))

# add to metadata of seurat object
seurat <- AddMetaData(seurat, metadata = metadata$RNA, col.name = "Integrated_RNA_SEACell_ID")
seurat <- AddMetaData(seurat, metadata = metadata$scHelper_cell_type_by_proportion, col.name = "scHelper_cell_type_by_proportion")
seurat <- AddMetaData(seurat, metadata = metadata$broad, col.name = "scHelper_cell_type_broad_by_proportion")
seurat <- AddMetaData(seurat, metadata = metadata$k, col.name = "Mapping_k")

#   Set levels
seurat@meta.data$scHelper_cell_type_by_proportion <- base::factor(seurat@meta.data$scHelper_cell_type_by_proportion, levels = scHelper_cell_type_order)
seurat@meta.data$scHelper_cell_type_broad_by_proportion <- base::factor(seurat@meta.data$scHelper_cell_type_broad_by_proportion, levels = scHelper_cell_type_order)
seurat@meta.data$Mapping_k <- base::factor(seurat@meta.data$Mapping_k, levels = names(k_colors))

head(seurat@meta.data)
sum(is.na((seurat@meta.data)))

## save seacells seurat object with new metadata
saveRDS(seurat, paste0(rds_path, stage, "_seacells_seurat_integrated.RDS"), compress = FALSE)

## plot new metadata on SEACell UMAPs

# Plot stage
seurat@meta.data$stage <- base::factor(seurat@meta.data$stage, levels = stage_order)
stage_cols <- stage_colours[base::levels(base::droplevels(seurat@meta.data$stage))]

png(paste0(plot_path, "2_stage_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'stage', label = TRUE, 
        label.size = 9, label.box = TRUE, repel = TRUE,
        pt.size = 10,
        cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

# Plot transferred cell type labels
scHelper_cols <- scHelper_cell_type_colours[base::levels(base::droplevels(seurat@meta.data$scHelper_cell_type_by_proportion))]

png(paste0(plot_path, "2_scHelper_cell_type_by_proportion_from_integration_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'scHelper_cell_type_by_proportion', label = TRUE, 
        label.size = ifelse(length(unique(seurat$stage)) == 1, 9, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(seurat$stage)) == 1, 6, 6), 
        cols = scHelper_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

png(paste0(plot_path, "2_scHelper_cell_type_by_proportion_from_integration_UMAP_no_label.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'scHelper_cell_type_by_proportion', label = FALSE, 
        label.box = FALSE,
        pt.size = ifelse(length(unique(seurat$stage)) == 1, 6, 6), 
        cols = scHelper_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

# Plot transferred broad cell type labels
scHelper_cols <- scHelper_cell_type_colours[base::levels(base::droplevels(seurat@meta.data$scHelper_cell_type_broad_by_proportion))]

png(paste0(plot_path, "2_scHelper_cell_type_broad_by_proportion_from_integration_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'scHelper_cell_type_broad_by_proportion', label = TRUE, 
        label.size = ifelse(length(unique(seurat$stage)) == 1, 9, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(seurat$stage)) == 1, 6, 6), 
        cols = scHelper_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

png(paste0(plot_path, "2_scHelper_cell_type_broad_by_proportion_from_integration_UMAP_no_label.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'scHelper_cell_type_broad_by_proportion', label = FALSE, 
        label.box = FALSE,
        pt.size = ifelse(length(unique(seurat$stage)) == 1, 6, 6), 
        cols = scHelper_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

# Plot k values of transferred labels
k_cols <- k_colors[base::levels(base::droplevels(seurat@meta.data$Mapping_k))]

png(paste0(plot_path, "2_mapping_k_UMAP.png"), width=25, height=20, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'Mapping_k', label = FALSE, 
        pt.size = 10,
        cols = k_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(plot.title = element_blank())
graphics.off()

#############################################################################################################################
#########################            3) Generate ATAC SEACell to ATAC single cell map           #############################
#############################################################################################################################

print("Projecting metacell labels back to single cells...")

# extract metacell map from seurat object
print("Metacell map: ")
metacell_map <- tibble::rownames_to_column(seurat@meta.data, var = "ATAC")
print(head(metacell_map))

# project back onto single cells
print("Single cell map: ")
singlecell_map <- base::merge(SEACell_map, metacell_map, by.x = "SEACell", by.y = "ATAC", all.x = TRUE)
print(head(singlecell_map))

# select which columns want to keep
print("Single cell map with columns selected:")
singlecell_map <- singlecell_map %>% dplyr::select(c("SEACell", "index", "Integrated_RNA_SEACell_ID", "scHelper_cell_type_by_proportion", "scHelper_cell_type_broad_by_proportion", "Mapping_k"))
print(head(singlecell_map))

# save
write.csv(singlecell_map, paste0(rds_path, stage, '_ATAC_singlecell_integration_map.csv'))