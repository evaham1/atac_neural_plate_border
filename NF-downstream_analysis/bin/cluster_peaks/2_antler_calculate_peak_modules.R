#!/usr/bin/env Rscript

print("Calculate Peak modules using Antler")

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
    # for peak matrix - normalised and raw
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/"
    # for combined SEACell metadata
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/0_combining_outputs/csv_files/"
    # output paths:
    rds_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/plots/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

########################       CELL STATE COLOURS    ########################################
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak', 
                              'PGC', 'BI', 'meso', 'endo', 'Unmapped')

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

## Function to generate plot data for Complex Heatmap with ordering of data by cell metadata and peak modules
PrepPeakModuleHeatmap <- function (peak_normalised_matrix, metadata, col_order,
                                    custom_order_column = metadata[1], custom_order = NULL,
                                    peak_modules, peak_row_annotation = TRUE,
                                    scale_data = TRUE) 
{
  
  ### Cell-level ordering and annotations ###
  
  # Initiated column anndata
  col_ann <- metadata %>% mutate_if(is.character, as.factor)
  
  # If 'custom_order' is set use this to reorder cells
  if (!is.null(custom_order)) {
    if (!setequal(custom_order, unique(col_ann[[custom_order_column]]))) {
      stop("custom_order factors missing from custom_order_column \n\n")
    }
    levels(col_ann[[custom_order_column]]) <- custom_order
    col_ann <- col_ann[order(col_ann[[custom_order_column]]),
                       , drop = FALSE]
  }
  
  # Automatically order cells by the columns specified in 'col_order'
  col_ann <- col_ann[do.call("order", c(col_ann[col_order],
                                        list(decreasing = FALSE))), , drop = FALSE]
  
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

########################################################################################################
#                                 Read in data and clean up                               #
########################################################################################################

########## RAW COUNTS MATRIX -> FOR ANTLER #############

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

########## NORMALISED COUNTS MATRIX -> FOR PLOTTING #############

# read in SEACells data
SEACells_normalised_summarised <- fread(paste0(data_path, "Filtered_normalised_summarised_counts.csv"), header = TRUE)
print("Normalised data read in!")

# Extract SEACell IDs from first column
SEACells_IDs <- SEACells_normalised_summarised$V1
print(head(SEACells_IDs))
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

########## COMBINED SEACELL METADATA #############

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

print(head(metadata))


########################################################################################################
#                                 Generate Antler Inputs                               #
########################################################################################################

# using raw count matrix for this as Antler has its own normalisation step

# generate fake metadata needed for Antler
pheno_data <- data.frame(row.names = rownames(SEACells_summarised),
                         "timepoint" = rep(1, nrow(SEACells_summarised)),
                         "treatment" = rep("null", nrow(SEACells_summarised)),
                         "replicate_id" = rep(1, nrow(SEACells_summarised))
)

# create antler folder
antler_path = "./antler/"
dir.create(antler_path)

# save pheno data
write.table(pheno_data, file = paste0(antler_path, "phenoData.csv"), row.names = T, sep = "\t", col.names = T)

# save count data
write.table(t(SEACells_summarised), file = paste0(antler_path, "assayData.csv"), row.names = T, sep = "\t", col.names = T, quote = F)

########################################################################################################
#                            Load Antler data and generate correlation matrix                          #
########################################################################################################

# Create Antler object
antler_data <- Antler$new(output_folder = plot_path, num_cores = ncores)
antler_data$load_dataset(folder_path = antler_path)

# Normalize data
antler_data$normalize(method = 'CPM')

###########################
# Calculate GMs unbiasedly
antler_data$gene_modules$identify(
  name                  = "unbiasedPMs",
  corr_t                = 0.3,  # the Spearman correlation threshold
  corr_min              = 3,    # min. number of genes a gene must correlate with
  mod_min_cell          = 10,   # min. number of cells expressing the module
  mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
  process_plots         = TRUE)

# rename peak modules
names(antler_data$gene_modules$lists$unbiasedPMs$content) <- paste0("PM", 1:length(antler_data$gene_modules$lists$unbiasedPMs$content))

# how many peak modules were generated
print(paste0("Number of peak modules made: ", length(antler_data$gene_modules$lists$unbiasedPMs$content)))

########## DE PMs ############## <- NEED TO WORK ON THIS ONCE HAVE MY METADATA
# # Subset peak modules that have at least 50% of peaks DE > 0.25 logFC & FDR < 0.001 between scHelperCellTypes
# gms <- DEGeneModules(seurat_data, antler_data$gene_modules$get("unbiasedGMs"), logfc = 0.25, pval = 0.001, selected_gene_proportion = 0.5, active_ident = meta_col,
#                      ident_1 = c('aPPR', 'pPPR', 'PPR'), ident_2 = c('NC', 'dNC'))
# 
# # save unbiasedGMs_DE in antler object
# antler_data$gene_modules$set(name= "unbiasedGMs_DE", content = gms)

########## Write GMs ##############
export_antler_modules <- function(antler_object, publish_dir, names_list){
  for(gm_list in names_list){
    mods = antler_data$gene_modules$lists[[gm_list]]$content
    for (i in seq(length(mods))) {
      modname = base::names(mods)[i]
      if (is.null(modname)) {
        modname = paste0("PM: ", i)
      }
      write(paste0(modname, "; ", paste0(mods[[i]], collapse = ", ")), file = paste0(publish_dir, '/', gm_list, '.txt'), append = TRUE)
    }
  }
}

export_antler_modules(antler_data, publish_dir = rds_path, names_list = 'unbiasedGMs')

########## Save Antler object ##############
saveRDS(antler_data, paste0(rds_path, 'antler_out.RDS'))


########################################################################################################
#                                        Visualisations                                                #
########################################################################################################

# read in antler RDS file if working interactively:
# antler_data <- readRDS(paste0(rds_path, 'antler_out.RDS'))

## subset matrix to only include peaks that are in PMs
peaks <- unlist(antler_data$gene_modules$lists$unbiasedPMs$content)
length(peaks)
filtered_normalised_matrix <- SEACells_normalised_summarised[, peaks]

# prepare scHelper_cell_type order and colors so by subsetting based on what is in the matrix
order <- scHelper_cell_type_order[scHelper_cell_type_order %in% metadata$scHelper_cell_type]
cols <- scHelper_cell_type_colours[order]

## prepare plot data - ordering by stage and then within that by scHelper_cell_type with custom order
plot_data <- PrepPeakModuleHeatmap(filtered_normalised_matrix, metadata, col_order = c('stage', 'scHelper_cell_type'),
                                    custom_order_column = "scHelper_cell_type", custom_order = order,
                                    peak_modules = antler_data$gene_modules$lists$unbiasedPMs$content, peak_row_annotation = TRUE)

## Create bottom annotation - scHelper_cell_type labels
bottom_annotation = HeatmapAnnotation(scHelper_cell_type = anno_simple(x = as.character(plot_data$col_ann$scHelper_cell_type),
                                                                       col = cols, height = unit(0.5, "cm")), show_annotation_name = FALSE,
                                      labels = anno_mark(at = cumsum(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths) - floor(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths/2),
                                                         labels = rle(as.character(plot_data$col_ann$scHelper_cell_type))$values,
                                                         which = "column", side = 'bottom',
                                                         labels_gp = gpar(fontsize = 10), lines_gp = gpar(lwd=2)))

## Create top annotation - stage
top_annotation <- HeatmapAnnotation(stage = anno_block(gp = gpar(fill = stage_colours),
                                                       labels = levels(plot_data$col_ann$stage),
                                                       labels_gp = gpar(col = "white", fontsize = 20, fontface='bold')),
                                    simple_anno_size = unit(1, "cm"),
                                    annotation_label = "stage", gp = gpar(fontsize = 20))

###### Complex Heatmap
plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = FALSE,
                show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                row_split = plot_data$row_ann$`Peak Modules`, column_split = plot_data$col_ann$stage,
                bottom_annotation = bottom_annotation,
                top_annotation = top_annotation)
png('temp.png', width = 60, height = 40, res = 400, units = 'cm')
plot
graphics.off()


## still need to add:
# column annotations ie colours and labels for stages on top
# row annotations ie colours and labels for PMs
# option to order cells (columns) by scHelper cell type within stage
# option to order PMs by accessibility









