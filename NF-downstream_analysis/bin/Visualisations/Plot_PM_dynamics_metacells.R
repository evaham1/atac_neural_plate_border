#!/usr/bin/env Rscript

print("Calculate the dynamics of PM accessiblity across metacells by latent time weighted by lineage probs (GAMs)")

############################## Load libraries #######################################
library(optparse)
library(future)
library(pheatmap)
library(tidyverse)
library(Antler)
library(RColorBrewer)
library(scHelper)
library(ComplexHeatmap) # Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. DOI: 10.1093/bioinformatics/btw313
library(data.table)
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
    
    # data paths for the different inputs
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/" # normalised count matrix
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/rds_files/FullData/" # the full data peal modules
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Metacell_metadata_latent_time/" # latent time on metacells metadata 
    # output paths:
    rds_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/4_PM_Dynamics/FullData/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/4_PM_Dynamics/FullData/plots/"
    
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

########################       CELL STATE COLOURS    ########################################
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'Non-neural',
                              'PPR', 'aPPR', 'pPPR', 'Placodal',
                              'eNPB', 'NPB', 'aNPB', 'pNPB',
                              'NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'Neural',
                              'node', 'streak', 'PGC', 'BI', 'meso', 'endo', 'Contam',
                              'MIXED', 'Unmapped')
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
########################       STAGE COLOURS     ###########################################
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
names(stage_colours) <- stage_order
############################################################################################
lineage_colours = c('placodal' = '#3F918C', 'NC' = '#DE4D00', 'neural' = '#8000FF')

########################################################################################################
#                                 Read in data and clean up                               #
########################################################################################################

########## COMBINED SEACELL METADATA ############# - with average latent time and lineage probabilities 

metadata <- read.csv(paste0(data_path, "rds_files/Combined_SEACell_integrated_metadata_latent_time.csv"), row.names = 'SEACell_ID')

# Add stage to metadata using SEACell IDs
substrRight <- function(x, n){
  sapply(x, function(xx)
    substr(xx, (nchar(xx)-n+1), nchar(xx))
  )
}
metadata <- metadata %>% mutate(stage = substrRight(rownames(metadata), 3))
metadata <- metadata[,-1]

# Change cell names to match matrix
rownames(metadata) <- gsub('-', '_', rownames(metadata))

# Check metadata
print(head(metadata))

print("Metadata read in!")

# Set rownames as a column and write the metadata out again
metadata_to_save <- rownames_to_column(metadata, var = "Rownames")
write_csv(metadata_to_save, paste0(rds_path, "Combined_SEACell_integrated_metadata_latent_time.csv"))

########## NORMALISED COUNTS MATRIX ############# - are these scaled?

print("Reading in normalised count data...")

# read in SEACells data
SEACells_normalised_summarised <- fread(paste0(data_path, "Filtered_normalised_summarised_counts.csv"), header = TRUE)
print("Normalised data read in!")

# Extract SEACell IDs from first column
SEACells_IDs <- SEACells_normalised_summarised$V1
length(SEACells_IDs)

# Check for duplicates
table(duplicated(SEACells_IDs))

# Clean up df
SEACells_normalised_summarised <- SEACells_normalised_summarised[,-1]

# Turn into numeric matrix for downstream processing
SEACells_normalised_summarised_numeric <- as.matrix(sapply(SEACells_normalised_summarised, as.numeric))

# Add SEACell IDs as rownames
rownames(SEACells_normalised_summarised_numeric) <- SEACells_IDs

# change cell names for Antler
rownames(SEACells_normalised_summarised_numeric) <- gsub('-', '_', rownames(SEACells_normalised_summarised_numeric))

# Check resulting matrix
print(dim(SEACells_normalised_summarised_numeric))
print("Preview of summarised count df:")
print(SEACells_normalised_summarised_numeric[1:4, 1:4])

# Overwrite cleaned data
SEACells_normalised_summarised <- SEACells_normalised_summarised_numeric

print("Normalised data read in!")


########## PEAK MODULES ############# 
antler_data <- readRDS(paste0(data_path, 'antler.RDS'))

pms <- antler_data$gene_modules$lists$unbiasedPMs$content
names(pms)

print("Antler data read in!")


#### check cell ids in metadata and accessibility data match
if (nrow(SEACells_normalised_summarised_numeric) == nrow(metadata) &
    nrow(metadata) == sum(rownames(SEACells_normalised_summarised_numeric) %in% rownames(metadata))){
  print("Metacell IDs match") } else {stop("Problem! Metacell IDs of accessibility data and metacell data dont match!!")}


########################################################################################################
#                                 Plot GAMs                               #
########################################################################################################

#### !!! edit the below code which was alex's to use for seurat RNA, instead for matrix directly ATAC

# DefaultAssay(seurat_data) <- "RNA"

# # Iteratively get expression data for each gene module and bind to tidy dataframe
# plot_data <- data.frame()
for(module in names(pms)){
  print(module)
  
  # select peaks from that peak module
  peaks <- unique(unlist(pms[[module]]))
  print(paste0("Number of peaks in PMs ", length(peaks)))
  
  # subset accessibility matrix for peaks in that peak module
  if ( length(as.vector(peaks)) == length(as.vector(peaks[peaks %in% colnames(SEACells_normalised_summarised)])) ){
    filtered_normalised_matrix <- SEACells_normalised_summarised[, as.vector(peaks)]
  } else {stop("ERROR: PM peaks are not found in the filtered peak matrix!")}

  # prep plotting data
  temp <- merge(filtered_normalised_matrix, metadata[,c('rna_latent_time', 'rna_lineage_NC_probability', 'rna_lineage_neural_probability', 'rna_lineage_placodal_probability'), drop=FALSE], by=0)
  plot_data <- temp %>%
    column_to_rownames('Row.names') %>%
    pivot_longer(!c(rna_latent_time, rna_lineage_NC_probability, rna_lineage_neural_probability, rna_lineage_placodal_probability)) %>%
    dplyr::rename(scaled_accessibility = value) %>%
    dplyr::rename(peak = name) %>%
    pivot_longer(cols = !c(rna_latent_time, peak, scaled_accessibility)) %>%
    dplyr::rename(lineage_probability = value) %>%
    dplyr::rename(lineage = name) %>%
    dplyr::group_by(lineage) %>%
    dplyr::mutate(lineage = unlist(strsplit(lineage, '_'))[3]) %>%
    dplyr::bind_rows() %>%
    dplyr::ungroup()
  
  # make gam plot
  plot = ggplot(plot_data, aes(x = rna_latent_time, y = scaled_accessibility)) +
    geom_smooth(method="gam", formula = y ~ s(x, bs = "cr", k = 7), se=FALSE, mapping = aes(weight = lineage_probability, color = lineage, group=lineage)) +
    xlab("Latent time") + ylab("Scaled accessibility") +
    theme_classic() +
    scale_colour_manual(values=lineage_colours) +
    ggtitle(module)
  
  # print and save plot
  png(paste0(plot_path, module, '_k7.png'), width = 18, height = 12, res = 200, units = 'cm')
  print(plot)
  graphics.off()
  
}


########################################################################################################
#                                 Plot distribution of stages across latent time                       #
########################################################################################################

latent_times <- metadata %>%
  dplyr::select(c("rna_latent_time", "stage"))
# HH5_latent_times <- latent_times %>%
#   dplyr::filter(stage == "HH5")
# HH6_latent_times <- latent_times %>%
#   dplyr::filter(stage == "HH6")
# HH7_latent_times <- latent_times %>%
#   dplyr::filter(stage == "HH7")
# ss4_latent_times <- latent_times %>%
#   dplyr::filter(stage == "ss4")
# ss8_latent_times <- latent_times %>%
#   dplyr::filter(stage == "ss8")

stage_cols = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")

plot <- ggplot(latent_times, aes(x = rna_latent_time, fill = stage)) +
  geom_density(alpha = 0.7, color = NA) + 
  scale_fill_manual(values = stage_cols) +
  theme_minimal() +
  xlab("Latent time") + ylab("Metacell Density") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 18))

png(paste0(plot_path, 'Stage_distribution_across_latent_time.png'), width = 18, height = 12, res = 400, units = 'cm')
print(plot)
graphics.off()



########################################################################################################
#                                 Plot more heatmaps of subsets of PMs                       #
########################################################################################################

# filter out the temporal PMs PM5, PM8 and PM9
PMs_to_plot <- subset(antler_data$gene_modules$lists$unbiasedPMs$content, 
                      !(names(antler_data$gene_modules$lists$unbiasedPMs$content) %in% c("FullData_PM5", "FullData_PM8", "FullData_PM9")))

# subset matrix to only include peaks that are in PMs
peaks <- unique(unlist(PMs_to_plot))
print("Number of peaks in Full Data PMs:")
length(peaks)
if ( length(as.vector(peaks)) == length(as.vector(peaks[peaks %in% colnames(SEACells_normalised_summarised)])) ){
  filtered_normalised_matrix <- SEACells_normalised_summarised[, as.vector(peaks)]
} else {stop("ERROR: PM peaks are not found in the filtered peak matrix!")}

# filter out unmapped and contaminating cell states
seacell_filtered_metadata <- metadata %>% filter(!scHelper_cell_type_by_proportion %in% c("Unmapped", "streak", "meso", "endo", "BI", "pEpi", "Contam", "MIXED"))
seacell_filtered_normalised_matrix <- filtered_normalised_matrix[rownames(seacell_filtered_metadata), ]

# prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
order <- scHelper_cell_type_order[scHelper_cell_type_order %in% seacell_filtered_metadata$scHelper_cell_type_by_proportion]
scHelper_cell_type_cols <- scHelper_cell_type_colours[order]

# Prepare plot data - ordering by scHelper cell type and then by hclust
plot_data <- PrepPeakModuleHeatmap(seacell_filtered_normalised_matrix, seacell_filtered_metadata, col_order = c('stage', 'scHelper_cell_type_by_proportion'),
                                   custom_order_column = "scHelper_cell_type_by_proportion", custom_order = order,
                                   hclust_SEACells = TRUE, hclust_SEACells_within_groups = TRUE,
                                   peak_modules = PMs_to_plot, peak_row_annotation = TRUE,
                                   log_path = NULL, scale_data = TRUE)

# Plot heatmap
plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                bottom_annotation = CreateCellTypeAnnotation(plot_data, scHelper_cell_type_cols),
                top_annotation = CreateStageAnnotation(plot_data, stage_colours),
                col = PurpleAndYellow())

png(paste0(plot_path, 'All_peak_modules_not_temporal_mapped_not_contam_SEACells_ordered_by_cell_type.png'), width = 60, height = 40, res = 400, units = 'cm')
print(plot)
graphics.off()


# only keep placodal modules PM1, PM2, PM3, PM4
PMs_to_plot <- subset(antler_data$gene_modules$lists$unbiasedPMs$content, 
                      (names(antler_data$gene_modules$lists$unbiasedPMs$content) %in% c("FullData_PM1", "FullData_PM2", "FullData_PM3", "FullData_PM4")))

# subset matrix to only include peaks that are in PMs
peaks <- unique(unlist(PMs_to_plot))
print("Number of peaks in Full Data PMs:")
length(peaks)
if ( length(as.vector(peaks)) == length(as.vector(peaks[peaks %in% colnames(SEACells_normalised_summarised)])) ){
  filtered_normalised_matrix <- SEACells_normalised_summarised[, as.vector(peaks)]
} else {stop("ERROR: PM peaks are not found in the filtered peak matrix!")}

# filter out unmapped and contaminating cell states
seacell_filtered_metadata <- metadata %>% filter(!scHelper_cell_type_by_proportion %in% c("Unmapped", "streak", "meso", "endo", "BI", "pEpi", "Contam", "MIXED"))
seacell_filtered_normalised_matrix <- filtered_normalised_matrix[rownames(seacell_filtered_metadata), ]

# prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
order <- scHelper_cell_type_order[scHelper_cell_type_order %in% seacell_filtered_metadata$scHelper_cell_type_by_proportion]
scHelper_cell_type_cols <- scHelper_cell_type_colours[order]

# Prepare plot data - ordering by scHelper cell type and then by hclust
plot_data <- PrepPeakModuleHeatmap(seacell_filtered_normalised_matrix, seacell_filtered_metadata, col_order = c('stage', 'scHelper_cell_type_by_proportion'),
                                   custom_order_column = "scHelper_cell_type_by_proportion", custom_order = order,
                                   hclust_SEACells = TRUE, hclust_SEACells_within_groups = TRUE,
                                   peak_modules = PMs_to_plot, peak_row_annotation = TRUE,
                                   log_path = NULL, scale_data = TRUE)

# Plot heatmap
plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                bottom_annotation = CreateCellTypeAnnotation(plot_data, scHelper_cell_type_cols),
                top_annotation = CreateStageAnnotation(plot_data, stage_colours),
                col = PurpleAndYellow())

png(paste0(plot_path, 'Placodal_peak_modules_mapped_not_contam_SEACells_ordered_by_cell_type.png'), width = 60, height = 40, res = 400, units = 'cm')
print(plot)
graphics.off()

########################################################################################################
#                                 Add PMs as features to seurat object                               #
########################################################################################################

## for each peak module, calculate the average scaled accessibility score per metacell
module_scores <- data.frame(matrix(nrow = 2156, ncol = 0))
rownames(module_scores) <- rownames(SEACells_normalised_summarised)

for(module in names(pms)){
  module_mat <- SEACells_normalised_summarised[, as.vector(unique(unlist(pms[[module]])))]
  module_scores <- module_scores %>% mutate(!!module := rowMeans(module_mat))
}

module_scores <- tibble::rownames_to_column(module_scores, var = "SEACell_ID")

write_csv(module_scores, paste0(rds_path, "PM_avg_scores.csv"))