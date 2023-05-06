#!/usr/bin/env Rscript

print("Calculate Peak modules using Antler")

## need to make sure writes PMs out properly
## need to go over heatmap plots with newly re-processed data
## also need to go over params for calculating PMs, maybe make less stringent

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
    
    # interactive local paths
    #data_path = "./local_test_data/peak_count_matrix_to_cluster/"
    #plot_path = "./local_test_data/clustered_peaks/plots"
    #rds_path = "./clustered_peaks/rds_files"
    
    # interactive NEMO paths
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/"
    # output paths:
    rds_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/plots/"
    PMs_path = ""
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    PMs_path = "./PMs/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  dir.create(PMs_path, recursive = T)
}

########################       CELL STATE COLOURS    ########################################
scHelper_cell_type_order <- c('node', 'streak', 'PGC', 'BI', 'meso', 'endo',
                              'EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'Unmapped')

scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                                "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                                "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB", "#EAEAEA")

names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                       'vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo', 'Unmapped')
########################       STAGE COLOURS     ###########################################
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
names(stage_colours) <- stage_order
############################################################################################

########################       FUNCTIONS    ########################################

## Function to write peak modules
export_antler_modules <- function(antler_object, publish_dir, names_list){
  for(gm_list in names_list){
    mods = antler_object$gene_modules$lists[[gm_list]]$content
    for (i in seq(length(mods))) {
      modname = base::names(mods)[i]
      if (is.null(modname)) {
        modname = paste0("PM: ", i)
      }
      write(paste0(modname, "; ", paste0(mods[[i]], collapse = ", ")), file = paste0(publish_dir, '/', gm_list, '.txt'), append = TRUE)
    }
  }
}

## Function to generate plot data for Complex Heatmap with ordering of data by cell metadata and peak modules
# peak_normalised_matrix = normalised count matrix with seacells as rows and peaks as columns
# cell_metadata = each row a metacell, columns correspond to metadata eg scHelper_cell_type, stage, etc
# col_order = which columns of cell_metadata to use to order cells, if you specify more than one will respect ordering (e.g. if c('stage', 'cell_type'), cells will be ordered first by stage and then by cell type)
# custom_order_column = for one columns of the cell_metadata you can specify a custom order of variables by which to order cells, this param is to select which column you use (eg 'cell_type')
# custom_order = for one columns of the cell_metadata you can specify a custom order of variables by which to order cells, this param is to input the custom ordering (eg c('NC', 'PPR', 'NPB'))
# order_SEACells = whether or not to hclust SEACells within their grouping and order by that
PrepPeakModuleHeatmap <- function (peak_normalised_matrix, cell_metadata, 
                                   col_order, custom_order_column = metadata[1], custom_order = NULL, order_SEACells = FALSE,
                                   peak_modules, peak_row_annotation = TRUE,
                                   scale_data = TRUE) 
{
  
  ### Cell-level ordering and annotations ###
  
  # Initiate column anndata
  col_ann <- cell_metadata %>% mutate_if(is.character, as.factor)
  
  # If 'custom_order' is set use this to reorder cells
  if (!is.null(custom_order)) {
    if (!setequal(custom_order, unique(col_ann[[custom_order_column]]))) {
      stop("custom_order factors missing from custom_order_column \n\n")
    }
    col_ann[[custom_order_column]] <- factor(col_ann[[custom_order_column]], levels = custom_order)
    col_ann <- col_ann[order(col_ann[[custom_order_column]]),]
  }
  
  # If 'col_order' is use these columns to order cells
  if (!is.null(col_order)) {
    col_ann <- col_ann[do.call("order", c(col_ann[col_order],
                                          list(decreasing = FALSE))), , drop = FALSE]
  }
  
  # Optionally hclust cells - either within cell group if 'col_order' has been specified, or just on all cells
  if (order_SEACells == TRUE) {
    if (is.null(col_order)) {
      dist_mat <- dist(peak_normalised_matrix, method = "euclidean")
      hclust_avg <- hclust(dist_mat, method = "average")
      ordered_SEACells <- hclust_avg$labels[c(hclust_avg$order)]
      col_ann <- col_ann[order(match(rownames(col_ann), ordered_SEACells)), , drop = FALSE]
    } else {
      cell_groups <- split(col_ann, col_ann[[tail(col_order, n=1)]])
      CellGroups_ordered_SEACells <- c()
      for (i in names(cell_groups)) {
        mat <- peak_normalised_matrix[rownames(cell_groups[[i]]), ]
        dist_mat <- dist(mat, method = "euclidean")
        hclust_avg <- hclust(dist_mat, method = "average")
        ordered_SEACells <- hclust_avg$labels[c(hclust_avg$order)]
        CellGroups_ordered_SEACells[[i]] <- ordered_SEACells
      }
      col_ann <- col_ann[order(match(rownames(col_ann), unlist(CellGroups_ordered_SEACells))), , drop = FALSE]
    }
  }
  
  ### Peak-level ordering and annotations ###
  
  # Optionally annotate peaks by their modules
  if (peak_row_annotation == TRUE) {
    row_ann <- stack(peak_modules) %>% dplyr::rename(`Peak Modules` = ind) %>%
      column_to_rownames("values")
  } else {
    row_ann <- NA
  }
  
  ### Prepare data for plotting ###
  
  # Order matrix by row and column annotation orders
  plot_data <- t(peak_normalised_matrix)[unlist(peak_modules), rownames(col_ann)]
  
  # Optionally scale
  if (scale_data) {
    cat("Scaling data \n")
    plot_data <- t(scale(t(plot_data)))
    plot_data <- replace(plot_data, plot_data >= 2, 2)
    plot_data <- replace(plot_data, plot_data <= -2, -2)
  }
  
  ### Output plotting data and annotations ###
  
  output <- list(plot_data = plot_data,
                 row_ann = row_ann,
                 col_ann = col_ann)
  return(output)
  
}


## Function to create bottom annotation for Complex Heatmap
# plot_data has to be generated by `PrepPeakModuleHeatmap` function and include $col_ann, scHelper_cell_type_colors should be named and ordered vector of colours
create_scHelper_cell_type_bottom_annotation <- function(plot_data, scHelper_cell_type_colors){
  return(
    HeatmapAnnotation(scHelper_cell_type = anno_simple(x = as.character(plot_data$col_ann$scHelper_cell_type),
                                                       col = scHelper_cell_type_colors, height = unit(0.5, "cm")), show_annotation_name = FALSE,
                      labels = anno_mark(at = cumsum(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths) - floor(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths/2),
                                         labels = rle(as.character(plot_data$col_ann$scHelper_cell_type))$values,
                                         which = "column", side = 'bottom',
                                         labels_gp = gpar(fontsize = 10), lines_gp = gpar(lwd=2)))
  )
}

## Function to create stage top annotation from ComplexHeatmap
# Needs plot_data generated from PrepPeakModuleHeatmap and named ordered vector of stage colours
create_stage_top_annotation <- function(plot_data, stage_colours){
  return(
    HeatmapAnnotation(stage = anno_block(gp = gpar(fill = stage_colours),
                                         labels = levels(plot_data$col_ann$stage),
                                         labels_gp = gpar(col = "white", fontsize = 20, fontface='bold')),
                      simple_anno_size = unit(1, "cm"),
                      annotation_label = "stage", gp = gpar(fontsize = 20))
  )
}


########################################################################################################
#                                 Read in data and clean up                               #
########################################################################################################

########## RAW COUNTS MATRIX -> FOR ANTLER #############

print("Reading in raw count data...")

# read in SEACells data
SEACells_summarised <- fread(paste0(data_path, "Filtered_summarised_counts.csv"), header = TRUE)
print("Raw data read in!")

# Check input
print("Preview of input df:")
print(SEACells_summarised[1:4, 1:4])
print(dim(SEACells_summarised))

# Extract SEACell IDs from first column
SEACells_IDs <- SEACells_summarised$V1
print(head(SEACells_IDs))
length(SEACells_IDs)

# Clean up df
SEACells_summarised <- SEACells_summarised[,-1]
print("Preview of input df after cleanup:")
print(SEACells_summarised[1:4, 1:4])
dim(SEACells_summarised)

# Turn into numeric matrix for downstream processing
SEACells_summarised_numeric <- as.matrix(sapply(SEACells_summarised, as.numeric))  

# Add SEACell IDs as rownames
rownames(SEACells_summarised_numeric) <- SEACells_IDs

# change cell names for Antler
rownames(SEACells_summarised_numeric) <- gsub('-', '_', rownames(SEACells_summarised_numeric))

# Check resulting matrix
print(dim(SEACells_summarised_numeric))
print("Preview of summarised count df:")
print(SEACells_summarised_numeric[1:4, 1:4])

# Overwrite cleaned data
SEACells_summarised <- SEACells_summarised_numeric

print("Raw data read in!")

########## NORMALISED COUNTS MATRIX -> FOR PLOTTING #############

print("Reading in normalised count data...")

# read in SEACells data
SEACells_normalised_summarised <- fread(paste0(data_path, "Filtered_normalised_summarised_counts.csv"), header = TRUE)
print("Normalised data read in!")

# Extract SEACell IDs from first column
SEACells_IDs <- SEACells_normalised_summarised$V1
length(SEACells_IDs)

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

########## COMBINED SEACELL METADATA #############

print("Reading in SEACell metadata...")

metadata <- read.csv(paste0(data_path, "Combined_SEACell_integrated_metadata.csv"), row.names = 'ATAC')

# Add stage to metadata using SEACell IDs
substrRight <- function(x, n){
  sapply(x, function(xx)
    substr(xx, (nchar(xx)-n+1), nchar(xx))
  )
}
metadata <- metadata %>% mutate(stage = substrRight(rownames(metadata), 3))

# Change cell names to match matrix
rownames(metadata) <- gsub('-', '_', rownames(metadata))

# Check metadata
print(head(metadata))

print("Metadata read in!")

########################################################################################################
#                                 Calculate peak modules on full data                                  #
########################################################################################################
################## use raw count matrix for this as Antler has its own normalisation step ##############

print("Calculating peak modules for full data...")

# generate fake metadata needed for Antler
pheno_data <- data.frame(row.names = rownames(SEACells_summarised),
                         "timepoint" = rep(1, nrow(SEACells_summarised)),
                         "treatment" = rep("null", nrow(SEACells_summarised)),
                         "replicate_id" = rep(1, nrow(SEACells_summarised))
)

# create antler folder
antler_path = "./antler_FullData/"
dir.create(antler_path)

# save pheno data
write.table(pheno_data, file = paste0(antler_path, "phenoData.csv"), row.names = T, sep = "\t", col.names = T)

# save count data
write.table(t(SEACells_summarised), file = paste0(antler_path, "assayData.csv"), row.names = T, sep = "\t", col.names = T, quote = F)

# Create Antler object
antler_data <- Antler$new(output_folder = paste0(plot_path, "FullData/"), num_cores = ncores)
antler_data$load_dataset(folder_path = antler_path)

# Normalize data
antler_data$normalize(method = 'CPM')

# Calculate GMs unbiasedly

# antler default params:
# antler$gene_modules$identify(
#     name                  = "unbiasedGMs",
#     corr_t                = 0.3,  # the Spearman correlation treshold
#     corr_min              = 3,    # min. number of genes a gene must correlate with
#     mod_min_cell          = 10,   # min. number of cells expressing the module
#     mod_consistency_thres = 0.3,  # ratio of expressed genes among "positive" cells
#     process_plots         = TRUE) # plot optimal module number heuristics

antler_data$gene_modules$identify(
  name                  = "unbiasedPMs",
  corr_t                = 0.3,
  corr_min              = 3,
  mod_min_cell          = 5,   # have reduced this because using SEACells not individual cells
  mod_consistency_thres = 0.3,
  process_plots         = TRUE)

# rename peak modules
names(antler_data$gene_modules$lists$unbiasedPMs$content) <- paste0("PM", 1:length(antler_data$gene_modules$lists$unbiasedPMs$content))

# how many peak modules were generated
print(paste0("Number of peak modules made: ", length(antler_data$gene_modules$lists$unbiasedPMs$content)))

# export peak modules
temp_path = paste0(PMs_path, "FullData/")
dir.create(temp_path)
export_antler_modules(antler_data, publish_dir = temp_path, names_list = "unbiasedPMs")

# save antler object
temp_path = paste0(rds_path, "FullData/")
dir.create(temp_path)
saveRDS(antler_data, paste0(temp_path, 'antler.RDS'))

print("Peak modules calculated for full data!")

# ########################################################################################################
# #                             Calculate peak modules on each stage                                     #
# ########################################################################################################
# ################## use raw count matrix for this as Antler has its own normalisation step ##############

print("Calculating peak modules for each stage...")

number_of_PMs_calculated <- c(length(antler_data$gene_modules$lists$unbiasedPMs$content))
corr_t_range <- c(0.3, 0.25, 0.4, 0.35, 0.4) # have adjusted these so you get between 10-25 PMs per stage

for (i in seq(1:length(stage_order))){
  
  # Extract stage
  stage <- stage_order[i]
  print(paste0("Calculating peak modules for stage: ", stage))
  
  # subset data matrix
  SEACells_summarised_temp <- SEACells_summarised[grep(stage, rownames(SEACells_summarised)), ]
  dim(SEACells_summarised_temp)
  
  # generate fake metadata needed for Antler
  pheno_data <- data.frame(row.names = rownames(SEACells_summarised_temp),
                           "timepoint" = rep(1, nrow(SEACells_summarised_temp)),
                           "treatment" = rep("null", nrow(SEACells_summarised_temp)),
                           "replicate_id" = rep(1, nrow(SEACells_summarised_temp))
  )
  
  # create antler folder
  antler_path = paste0("./antler_", stage, "/")
  dir.create(antler_path)
  
  # save pheno data
  write.table(pheno_data, file = paste0(antler_path, "phenoData.csv"), row.names = T, sep = "\t", col.names = T)
  
  # save count data
  write.table(t(SEACells_summarised_temp), file = paste0(antler_path, "assayData.csv"), row.names = T, sep = "\t", col.names = T, quote = F)
  
  # Create Antler object
  antler_temp <- Antler$new(output_folder = paste0(plot_path, stage, "/"), num_cores = ncores)
  antler_temp$load_dataset(folder_path = antler_path)
  
  # Normalize data
  antler_temp$normalize(method = 'CPM')
  
  # Calculate GMs unbiasedly
  antler_temp$gene_modules$identify(
    name                  = "unbiasedPMs",
    corr_t                = corr_t_range[i],  # the Spearman correlation threshold
    corr_min              = 3,    # min. number of genes a gene must correlate with
    mod_min_cell          = 5,   # min. number of cells expressing the module
    mod_consistency_thres = 0.3,  # ratio of expressed genes among "positive" cells
    process_plots         = TRUE)
  
  # rename peak modules
  names(antler_temp$gene_modules$lists$unbiasedPMs$content) <- paste0("PM", 1:length(antler_temp$gene_modules$lists$unbiasedPMs$content))
  
  # how many peak modules were generated
  print(paste0("Number of peak modules made: ", length(antler_temp$gene_modules$lists$unbiasedPMs$content)))
  number_of_PMs_calculated <- c(number_of_PMs_calculated, length(antler_temp$gene_modules$lists$unbiasedPMs$content))
  
  # export peak modules
  temp_path = paste0(PMs_path, stage, "/")
  dir.create(temp_path)
  export_antler_modules(antler_temp, publish_dir = temp_path, names_list = "unbiasedPMs")
  
  # save antler object
  temp_path = paste0(rds_path, stage, "/")
  dir.create(temp_path)
  saveRDS(antler_temp, paste0(temp_path, 'antler.RDS'))
  
  print(paste0("Peak modules calculated for stage: ", stage))
  
}

print("All peak modules calculated!")

########################################################################################################
#                     Visualise peak modules calculated on full data                                  #
########################################################################################################

print("Plotting peak modules on full data...")

# set plot path
temp_plot_path = paste0(plot_path, "FullData/")
dir.create(temp_plot_path, recursive = T)

# read in antler RDS file:
antler_data <- readRDS(paste0(rds_path, 'FullData/', 'antler.RDS'))

# subset matrix to only include peaks that are in PMs
peaks <- unlist(antler_data$gene_modules$lists$unbiasedPMs$content)
length(peaks)
filtered_normalised_matrix <- SEACells_normalised_summarised[, peaks]

# prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
order <- scHelper_cell_type_order[scHelper_cell_type_order %in% metadata$scHelper_cell_type]
scHelper_cell_type_colours <- scHelper_cell_type_colours[order]

########  Plot all peak modules ordered by stage and then by cell type ########

# Prepare plot data - ordering by stage and then within that by scHelper_cell_type with custom order
plot_data <- PrepPeakModuleHeatmap(filtered_normalised_matrix, metadata, col_order = c('stage', 'scHelper_cell_type'),
                                   custom_order_column = "scHelper_cell_type", custom_order = order,
                                   order_SEACells = TRUE,
                                   peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE)

# Plot heatmap
plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                bottom_annotation = create_scHelper_cell_type_bottom_annotation(plot_data, scHelper_cell_type_colours),
                top_annotation = create_stage_top_annotation(plot_data, stage_colours), 
                col = PurpleAndYellow())

png(paste0(temp_plot_path, 'All_peak_modules.png'), width = 60, height = 80, res = 400, units = 'cm')
plot
graphics.off()

########  Plot all peak modules ordered celltypes on cells which are mapped and not contam ########

# subset matrix to only include SEACells that mapped and not contam
seacell_filtered_metadata <- metadata %>% filter(!scHelper_cell_type %in% c("Unmapped", "streak", "meso", "endo", "BI"))
seacell_filtered_normalised_matrix <- filtered_normalised_matrix[rownames(seacell_filtered_metadata), ]

# prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
order <- scHelper_cell_type_order[scHelper_cell_type_order %in% seacell_filtered_metadata$scHelper_cell_type]
scHelper_cell_type_cols <- scHelper_cell_type_colours[order]

# Prepare plot data - ordering by scHelper cell type and then by hclust
plot_data <- PrepPeakModuleHeatmap(seacell_filtered_normalised_matrix, seacell_filtered_metadata, col_order = c('scHelper_cell_type'),
                                   custom_order_column = "scHelper_cell_type", custom_order = order, order_SEACells = TRUE,
                                   peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE)

# Plot heatmap
plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                bottom_annotation = create_scHelper_cell_type_bottom_annotation(plot_data, scHelper_cell_type_cols),
                top_annotation = create_stage_top_annotation(plot_data, stage_colours), 
                col = PurpleAndYellow())

png(paste0(temp_plot_path, 'All_peak_modules_mapped_not_contam_SEACells_ordered_by_cell_type.png.png'), width = 60, height = 130, res = 400, units = 'cm')
print(plot)
graphics.off()

print("Full data peak modules plotted!")

########################################################################################################
#                     Visualise peak modules calculated on stages                                        #
########################################################################################################

print("Plotting stage peak modules...")

for (i in seq(1:length(stage_order))){
  
  # Extract stage
  temp_stage <- stage_order[i]
  print(paste0("Plotting peak modules for stage: ", temp_stage))
  
  # set plot path
  temp_plot_path = paste0(plot_path, temp_stage, "/")
  dir.create(temp_plot_path, recursive = T)
  
  # read in antler RDS file:
  antler_data <- readRDS(paste0(rds_path, temp_stage, '/', 'antler.RDS'))
  
  # subset matrix to only include peaks that are in PMs
  peaks <- unlist(antler_data$gene_modules$lists$unbiasedPMs$content)
  length(peaks)
  filtered_normalised_matrix <- SEACells_normalised_summarised[, peaks]
  
  # subset matrix to only include SEACells that are in this stage
  filtered_normalised_matrix <- filtered_normalised_matrix[grep(temp_stage, rownames(filtered_normalised_matrix)), ]
  
  # prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
  stage_metadata <- metadata %>% filter(stage == temp_stage)
  order <- scHelper_cell_type_order[scHelper_cell_type_order %in% stage_metadata$scHelper_cell_type]
  scHelper_cell_type_cols <- scHelper_cell_type_colours[order]
  
  ########  Plot all peak modules ordered cell type ########
  
  # Prepare plot data - ordering by scHelper cell type and then by hclust
  plot_data <- PrepPeakModuleHeatmap(filtered_normalised_matrix, stage_metadata, col_order = c('scHelper_cell_type'),
                                     custom_order_column = "scHelper_cell_type", custom_order = order, order_SEACells = TRUE,
                                     peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE)
  
  # Plot heatmap
  plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                  show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                  row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                  bottom_annotation = create_scHelper_cell_type_bottom_annotation(plot_data, scHelper_cell_type_cols),
                  top_annotation = create_stage_top_annotation(plot_data, stage_colours[i]), 
                  col = PurpleAndYellow())
  
  png(paste0(temp_plot_path, 'All_peak_modules_ordered_by_celltype.png'), width = 60, height = 130, res = 400, units = 'cm')
  print(plot)
  graphics.off()
  
  ########  Plot all peak modules ordered hclust ########
  
  # Prepare plot data - ordering by hclust only
  plot_data <- PrepPeakModuleHeatmap(filtered_normalised_matrix, stage_metadata, col_order = NULL,
                                     order_SEACells = TRUE,
                                     peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE)
  
  # Plot heatmap
  plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                  show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                  row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                  bottom_annotation = create_scHelper_cell_type_bottom_annotation(plot_data, scHelper_cell_type_cols),
                  top_annotation = create_stage_top_annotation(plot_data, stage_colours[i]), 
                  col = PurpleAndYellow())
  
  png(paste0(temp_plot_path, 'All_peak_modules_ordered_by_hclust.png'), width = 60, height = 130, res = 400, units = 'cm')
  print(plot)
  graphics.off()
  
  ########  Plot all peak modules ordered celltypes on cells which are mapped and not contam ########
  
  # subset matrix to only include SEACells that mapped and not contam
  seacell_filtered_metadata <- stage_metadata %>% filter(!scHelper_cell_type %in% c("Unmapped", "streak", "meso", "endo", "BI"))
  seacell_filtered_normalised_matrix <- filtered_normalised_matrix[rownames(seacell_filtered_metadata), ]
  
  # prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
  order <- scHelper_cell_type_order[scHelper_cell_type_order %in% seacell_filtered_metadata$scHelper_cell_type]
  scHelper_cell_type_cols <- scHelper_cell_type_colours[order]
  
  # Prepare plot data - ordering by scHelper cell type and then by hclust
  plot_data <- PrepPeakModuleHeatmap(seacell_filtered_normalised_matrix, seacell_filtered_metadata, col_order = c('scHelper_cell_type'),
                                     custom_order_column = "scHelper_cell_type", custom_order = order, order_SEACells = TRUE,
                                     peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE)
  
  # Plot heatmap
  plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                  show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                  row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                  bottom_annotation = create_scHelper_cell_type_bottom_annotation(plot_data, scHelper_cell_type_cols),
                  top_annotation = create_stage_top_annotation(plot_data, stage_colours[i]), 
                  col = PurpleAndYellow())
  
  png(paste0(temp_plot_path, 'All_peak_modules_mapped_not_contam_SEACells_ordered_by_cell_type.png.png'), width = 60, height = 130, res = 400, units = 'cm')
  print(plot)
  graphics.off()
  
  
  
  ## Finished plots
  print(paste0(temp_stage, " peak modules plotted!"))
  
}

print("Stage peak modules plotted!")