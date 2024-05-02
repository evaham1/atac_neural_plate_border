#!/usr/bin/env Rscript

print("Calculate Peak modules using Antler")

## may need to go over params for calculating PMs, maybe make less stringent

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
    PMs_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/PMs/"
    
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

########## NEED TO UPDATE FUNCTION IN SCHELPER:
## function to export an antler object either as a list of peak ids or as a bed file
ExportAntlerModules <- function (antler_object, publish_dir, names_list = "unbiasedPMs") {
  # extract data based on name of slot in antler object
  for (gm_list in names_list) {
    mods = antler_object$gene_modules$lists[[gm_list]]$content
    
    # loop through each module
    df_full = data.frame(matrix(nrow = 0, ncol = 4))
    for (i in seq(length(mods))) {
      # for each module in that slot name it by numbering system
      modname = base::names(mods)[i]
      if (is.null(modname)) {
        modname = paste0("PM: ", i)
      }
      # write out peaks as a txt file
      write(paste0(modname, "; ", paste0(mods[[i]], collapse = ", ") ), 
            file = paste0(publish_dir, "/", gm_list, ".txt"), 
            append = TRUE)
      # write out peaks in a df to turn into bed file
      df = data.frame(matrix(nrow = 0, ncol = 4))
      for (j in seq(length(mods[[i]])) ) {
        chr <- strsplit(mods[[i]][j], "-")[[1]][1]
        start <- strsplit(mods[[i]][j], "-")[[1]][2]
        end <- strsplit(mods[[i]][j], "-")[[1]][3]
        df[j, ] <- c(chr, start, end, paste0(modname, "-", j))
      }
      df_full <- rbind(df_full, df)
    }
    write.table(df_full, paste0(publish_dir, "/", gm_list, "_with_chr.bed"), sep="\t", 
                row.names=FALSE, col.names = FALSE, quote = FALSE)
    df_full$X1 <- substring(df_full$X1, 4)
    write.table(df_full, paste0(publish_dir, "/", gm_list, "_without_chr.bed"), sep="\t", 
                row.names=FALSE, col.names = FALSE, quote = FALSE)
  }
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


########################################################################################################
#                                 Read in data and clean up                               #
########################################################################################################

########## RAW COUNTS MATRIX -> FOR ANTLER #############

print("Reading in raw count data...")

# read in SEACells data
SEACells_summarised <- fread(paste0(data_path, "Filtered_raw_summarised_counts.csv"), header = TRUE)
print("Raw data read in!")

# Check input
print("Preview of input df:")
print(SEACells_summarised[1:4, 1:4])
print(dim(SEACells_summarised))

# Extract SEACell IDs from first column
SEACells_IDs <- SEACells_summarised$V1
print(head(SEACells_IDs))
length(SEACells_IDs)

# Check for duplicates
table(duplicated(SEACells_IDs))

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
metadata <- metadata[,-1]

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
names(antler_data$gene_modules$lists$unbiasedPMs$content) <- paste0("FullData_PM", 1:length(antler_data$gene_modules$lists$unbiasedPMs$content))

# how many peak modules were generated
print(paste0("Number of peak modules made: ", length(antler_data$gene_modules$lists$unbiasedPMs$content)))

# export peak modules as list of peaks
temp_path = paste0(PMs_path, "FullData/")
dir.create(temp_path)
ExportAntlerModules(antler_data, publish_dir = temp_path, names_list = "unbiasedPMs")

# export peak modules as a bed file for motif analysus
# antler_data

# save antler object
temp_path = paste0(rds_path, "FullData/")
dir.create(temp_path)
saveRDS(antler_data, paste0(temp_path, 'antler.RDS'))

print("Peak modules calculated for full data!")

########################################################################################################
#                             Calculate peak modules on each stage                                     #
########################################################################################################
################## use raw count matrix for this as Antler has its own normalisation step ##############

print("Calculating peak modules for each stage...")

number_of_PMs_calculated <- c(length(antler_data$gene_modules$lists$unbiasedPMs$content))
corr_t_range <- c(0.3, 0.4, 0.6, 0.6, 0.6) # have adjusted these so you get between 10-25 PMs per stage

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
  names(antler_temp$gene_modules$lists$unbiasedPMs$content) <- paste0(stage, "_PM", 1:length(antler_temp$gene_modules$lists$unbiasedPMs$content))
  
  # how many peak modules were generated
  print(paste0("Number of peak modules made: ", length(antler_temp$gene_modules$lists$unbiasedPMs$content)))
  number_of_PMs_calculated <- c(number_of_PMs_calculated, length(antler_temp$gene_modules$lists$unbiasedPMs$content))
  
  # export peak modules
  temp_path = paste0(PMs_path, stage, "/")
  dir.create(temp_path)
  ExportAntlerModules(antler_temp, publish_dir = temp_path, names_list = "unbiasedPMs")
  
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
peaks <- unique(unlist(antler_data$gene_modules$lists$unbiasedPMs$content))
print("Number of peaks in Full Data PMs:")
length(peaks)
if ( length(as.vector(peaks)) == length(as.vector(peaks[peaks %in% colnames(SEACells_normalised_summarised)])) ){
  filtered_normalised_matrix <- SEACells_normalised_summarised[, as.vector(peaks)]
} else {stop("ERROR: PM peaks are not found in the filtered peak matrix!")}

# prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
order <- scHelper_cell_type_order[scHelper_cell_type_order %in% metadata$scHelper_cell_type_by_proportion]
scHelper_cell_type_colours <- scHelper_cell_type_colours[order]

########  Plot all peak modules ordered by stage and then by cell type ########

# Prepare plot data - ordering by stage and then within that by scHelper_cell_type with custom order
plot_data <- PrepPeakModuleHeatmap(filtered_normalised_matrix, metadata, 
                                     col_order = c('stage', 'scHelper_cell_type_by_proportion'), custom_order_column = "scHelper_cell_type_by_proportion", custom_order = order, 
                                     hclust_SEACells = TRUE, hclust_SEACells_within_groups = TRUE,
                                     peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE,
                                     log_path = NULL, scale_data = TRUE)

# Plot heatmap
plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                bottom_annotation = CreateCellTypeAnnotation(plot_data, scHelper_cell_type_colours),
                top_annotation = CreateStageAnnotation(plot_data, stage_colours),
                col = PurpleAndYellow())

png(paste0(temp_plot_path, 'All_peak_modules.png'), width = 60, height = 80, res = 400, units = 'cm')
plot
graphics.off()

png(paste0(temp_plot_path, 'All_peak_modules_shorter.png'), width = 60, height = 40, res = 400, units = 'cm')
plot
graphics.off()

########  Plot all peak modules ordered celltypes on cells which are mapped and not contam ########

# subset matrix to only include SEACells that mapped and not contam
seacell_filtered_metadata <- metadata %>% filter(!scHelper_cell_type_by_proportion %in% c("Unmapped", "streak", "meso", "endo", "BI", "pEpi", "Contam", "MIXED"))
seacell_filtered_normalised_matrix <- filtered_normalised_matrix[rownames(seacell_filtered_metadata), ]

# prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
order <- scHelper_cell_type_order[scHelper_cell_type_order %in% seacell_filtered_metadata$scHelper_cell_type_by_proportion]
scHelper_cell_type_cols <- scHelper_cell_type_colours[order]

# Prepare plot data - ordering by scHelper cell type and then by hclust
plot_data <- PrepPeakModuleHeatmap(seacell_filtered_normalised_matrix, seacell_filtered_metadata, col_order = c('stage', 'scHelper_cell_type_by_proportion'),
                                   custom_order_column = "scHelper_cell_type_by_proportion", custom_order = order,
                                   hclust_SEACells = TRUE, hclust_SEACells_within_groups = TRUE,
                                   peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE,
                                   log_path = NULL, scale_data = TRUE)

# Plot heatmap
plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                bottom_annotation = CreateCellTypeAnnotation(plot_data, scHelper_cell_type_cols),
                top_annotation = CreateStageAnnotation(plot_data, stage_colours),
                col = PurpleAndYellow())

png(paste0(temp_plot_path, 'All_peak_modules_mapped_not_contam_SEACells_ordered_by_cell_type.png'), width = 60, height = 130, res = 400, units = 'cm')
print(plot)
graphics.off()

png(paste0(temp_plot_path, 'All_peak_modules_mapped_not_contam_SEACells_ordered_by_cell_type_shorter.png'), width = 60, height = 40, res = 400, units = 'cm')
print(plot)
graphics.off()

########  Plot all peak modules on NC subset ########

# subset matrix to only include SEACells that mapped and not contam
seacell_filtered_metadata <- metadata %>% filter(scHelper_cell_type_by_proportion %in% c("NC", "dNC"))
seacell_filtered_normalised_matrix <- filtered_normalised_matrix[rownames(seacell_filtered_metadata), ]

# prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
order <- scHelper_cell_type_order[scHelper_cell_type_order %in% seacell_filtered_metadata$scHelper_cell_type_by_proportion]
scHelper_cell_type_cols <- scHelper_cell_type_colours[order]

# Prepare plot data - ordering by scHelper cell type and then by hclust
plot_data <- PrepPeakModuleHeatmap(seacell_filtered_normalised_matrix, seacell_filtered_metadata, col_order = c('stage', 'scHelper_cell_type_by_proportion'),
                                   custom_order_column = "scHelper_cell_type_by_proportion", custom_order = order,
                                   hclust_SEACells = TRUE, hclust_SEACells_within_groups = TRUE,
                                   peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE,
                                   log_path = NULL, scale_data = TRUE)

# Plot heatmap
plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                bottom_annotation = CreateCellTypeAnnotation(plot_data, scHelper_cell_type_cols),
                top_annotation = CreateStageAnnotation(plot_data, stage_colours),
                col = PurpleAndYellow())

png(paste0(temp_plot_path, 'All_peak_modules_mapped_NC_SEACells_ordered_by_cell_type.png'), width = 10, height = 20, res = 400, units = 'cm')
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
  peaks <- unique(unlist(antler_data$gene_modules$lists$unbiasedPMs$content))
  print("Number of peaks in PMs:")
  length(peaks)
  if ( length(as.vector(peaks)) == length(as.vector(peaks[peaks %in% colnames(SEACells_normalised_summarised)])) ){
    filtered_normalised_matrix <- SEACells_normalised_summarised[, as.vector(peaks)]
  } else {stop("ERROR: PM peaks are not found in the filtered peak matrix!")}
  
  # subset matrix to only include SEACells that are in this stage
  filtered_normalised_matrix <- filtered_normalised_matrix[grep(temp_stage, rownames(filtered_normalised_matrix)), ]
  
  # prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
  stage_metadata <- metadata %>% filter(stage == temp_stage)
  order <- scHelper_cell_type_order[scHelper_cell_type_order %in% stage_metadata$scHelper_cell_type_by_proportion]
  scHelper_cell_type_cols <- scHelper_cell_type_colours[order]
  
  ########  Plot all peak modules ordered cell type, within each cell type hclust ########
  
  # Prepare plot data - ordering by scHelper cell type and then by hclust
  plot_data <- PrepPeakModuleHeatmap(filtered_normalised_matrix, stage_metadata, 
                                     col_order = c('scHelper_cell_type_by_proportion'), custom_order_column = "scHelper_cell_type_by_proportion", custom_order = order, 
                                     hclust_SEACells = TRUE, hclust_SEACells_within_groups = TRUE,
                                     peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE,
                                     log_path = NULL, scale_data = TRUE)
  
  # Plot heatmap
  plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                  show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                  row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                  bottom_annotation = CreateCellTypeAnnotation(plot_data, scHelper_cell_type_cols),
                  top_annotation = CreateStageAnnotation(plot_data, stage_colours[i]), 
                  col = PurpleAndYellow())
  
  png(paste0(temp_plot_path, 'All_peak_modules_ordered_by_celltype.png'), width = 60, height = 130, res = 400, units = 'cm')
  print(plot)
  graphics.off()

  png(paste0(temp_plot_path, 'All_peak_modules_ordered_by_celltype_shorter.png'), width = 20, height = 30, res = 400, units = 'cm')
  print(plot)
  graphics.off()
  
  ########  Plot all peak modules ordered hclust ########
  # 
  # # Prepare plot data - ordering by hclust only
  # plot_data <- PrepPeakModuleHeatmap(filtered_normalised_matrix, stage_metadata, 
  #                                    col_order = c('scHelper_cell_type'), custom_order_column = "scHelper_cell_type", custom_order = order, 
  #                                    hclust_SEACells = TRUE, hclust_SEACells_within_groups = FALSE,
  #                                    peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE,
  #                                    log_path = paste0(temp_plot_path, "logs/all_cells_hclust/"), scale_data = TRUE)
  # # Plot heatmap
  # plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
  #                 show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
  #                 row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
  #                 bottom_annotation = CreateCellTypeAnnotation(plot_data, scHelper_cell_type_cols),
  #                 top_annotation = CreateStageAnnotation(plot_data, stage_colours[i]), 
  #                 col = PurpleAndYellow())
  # 
  # png(paste0(temp_plot_path, 'All_peak_modules_ordered_by_hclust.png'), width = 60, height = 130, res = 400, units = 'cm')
  # print(plot)
  # graphics.off()
  
  ########  Plot all peak modules ordered celltypes on cells which are mapped and not contam ########
  
  # subset matrix to only include SEACells that mapped and not contam
  seacell_filtered_metadata <- stage_metadata %>% filter(!scHelper_cell_type_by_proportion %in% c("Unmapped", "streak", "meso", "endo", "BI", "pEpi", "Contam", "MIXED"))
  seacell_filtered_normalised_matrix <- filtered_normalised_matrix[rownames(seacell_filtered_metadata), ]
  
  # prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
  order <- scHelper_cell_type_order[scHelper_cell_type_order %in% seacell_filtered_metadata$scHelper_cell_type_by_proportion]
  scHelper_cell_type_cols <- scHelper_cell_type_colours[order]
  
  # Prepare plot data - ordering by scHelper cell type and then by hclust
  plot_data <- PrepPeakModuleHeatmap(seacell_filtered_normalised_matrix, seacell_filtered_metadata, 
                                     col_order = c('scHelper_cell_type_by_proportion'), custom_order_column = "scHelper_cell_type_by_proportion", custom_order = order, 
                                     hclust_SEACells = TRUE, hclust_SEACells_within_groups = TRUE,
                                     peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE,
                                     log_path = NULL, scale_data = TRUE)
  
  # Plot heatmap
  plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                  show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                  row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                  bottom_annotation = CreateCellTypeAnnotation(plot_data, scHelper_cell_type_cols),
                  top_annotation = CreateStageAnnotation(plot_data, stage_colours[i]), 
                  col = PurpleAndYellow())
  
  png(paste0(temp_plot_path, 'All_peak_modules_mapped_not_contam_SEACells_ordered_by_cell_type.png'), width = 60, height = 130, res = 400, units = 'cm')
  print(plot)
  graphics.off()

  png(paste0(temp_plot_path, 'All_peak_modules_mapped_not_contam_SEACells_ordered_by_cell_type_shorter.png'), width = 20, height = 30, res = 400, units = 'cm')
  print(plot)
  graphics.off()

  if (c("NC", "dNC") %in% stage_metadata$scHelper_cell_type_by_proportion){

    ########  Plot all peak modules on NC subset ########
  
    # subset matrix to only include SEACells that mapped and not contam
    seacell_filtered_metadata <- stage_metadata %>% filter(scHelper_cell_type_by_proportion %in% c("NC", "dNC"))
    seacell_filtered_normalised_matrix <- filtered_normalised_matrix[rownames(seacell_filtered_metadata), ]
    
    # prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
    order <- scHelper_cell_type_order[scHelper_cell_type_order %in% seacell_filtered_metadata$scHelper_cell_type_by_proportion]
    scHelper_cell_type_cols <- scHelper_cell_type_colours[order]
    
    # Prepare plot data - ordering by scHelper cell type and then by hclust
    plot_data <- PrepPeakModuleHeatmap(seacell_filtered_normalised_matrix, seacell_filtered_metadata, 
                                      col_order = c('scHelper_cell_type_by_proportion'), custom_order_column = "scHelper_cell_type_by_proportion", custom_order = order, 
                                      hclust_SEACells = TRUE, hclust_SEACells_within_groups = TRUE,
                                      peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE,
                                      log_path = NULL, scale_data = TRUE)
    
    # Plot heatmap
    plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                    show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                    row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                    bottom_annotation = CreateCellTypeAnnotation(plot_data, scHelper_cell_type_cols),
                    top_annotation = CreateStageAnnotation(plot_data, stage_colours[i]), 
                    col = PurpleAndYellow())
    
    png(paste0(temp_plot_path, 'All_peak_modules_mapped_NC_SEACells_ordered_by_cell_type.png'), width = 10, height = 20, res = 400, units = 'cm')
    print(plot)
    graphics.off()

  }
  
  ## Finished plots
  print(paste0(temp_stage, " peak modules plotted!"))
  
}

print("Stage peak modules plotted!")