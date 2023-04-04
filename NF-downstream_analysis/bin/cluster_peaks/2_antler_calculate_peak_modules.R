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
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/"
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/0_combining_outputs/csv_files/"
    #data_path = "./output/NF-downstream_analysis/Processing/FullData/Peak_call/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/rds_files/"
    rds_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/plots/"
    
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

########################################################################################################
#                                 Read in data and clean up                               #
########################################################################################################

# read in seurat data - which seurat data?? would need to be full data object of SEACells which I havent made yet...
# seurat_data <- readRDS(data_path, full.names = TRUE))
# seurat_data <- readRDS(list.files(data_path, full.names = TRUE))

# read in SEACells data
SEACells_summarised <- fread(paste0(data_path, "Filtered_summarised_counts.csv"), header = TRUE)
print("Data read in!")

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
rownames(SEACells_summarised_numeric) <- gsub('#', '', rownames(SEACells_summarised_numeric))

# Check resulting matrix
print(dim(SEACells_summarised_numeric))
print("Preview of summarised count df:")
print(SEACells_summarised_numeric[1:4, 1:4])

# Overwrite cleaned data
SEACells_summarised <- SEACells_summarised_numeric

print("Data read in!")


######## read in metadata

metadata <- read.csv(paste0(data_path, "Combined_SEACell_integrated_metadata.csv"), row.names = 'ATAC')


########################################################################################################
#                                 Generate Antler Inputs                               #
########################################################################################################

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
names(antler_data$gene_modules$lists$unbiasedPMs$content) <- paste0("GM", 1:length(antler_data$gene_modules$lists$unbiasedPMs$content))

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
        modname = paste0("GM: ", i)
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
# antler <- readRDS(paste0(rds_path, 'antler_out.RDS'))

antler_data$gene_modules$lists$unbiasedPMs$content

peaks <- unlist(antler_data$gene_modules$lists$unbiasedPMs$content)
length(peaks)

## subset matrix to only include peaks that are in PMs
SEACells_summarised[1:2, 1:2]
dim(SEACells_summarised)

filtered_matrix <- SEACells_summarised[, peaks]

dim(filtered_matrix)
filtered_matrix[1:2, 1:2]

## scale matrix
#plot_data <- round(plot_data, digits = 3)
# plot_data <- t(scale(filtered_matrix))

plot_data[1:2, 1:2]


## plot data
pheatmap(plot_data, show_colnames = FALSE, show_rownames = FALSE, border_color = NA)


# 
# GeneModuleOrder <- function (seurat_obj, gene_modules, metadata_1 = NULL, order_1 = NULL, 
#                              metadata_2 = NULL, order_2 = NULL, rename_modules = NULL, 
#                              plot_path = "scHelper_log/GM_classification/") 
# {
#   classified_gms_1 <- GeneModuleClassify(seurat_obj, gene_modules, 
#                                          metadata = metadata_1, plot_path = plot_path)
#   classified_gms_1 <- classified_gms_1 %>% arrange(match(!!sym(metadata_1), 
#                                                          order_1)) %>% group_by(!!sym(metadata_1)) %>% mutate(pos = 1:n()) %>% 
#     mutate(new_name = paste(!!sym(metadata_1), pos, sep = "-")) %>% 
#     dplyr::select(gene_module, !!sym(metadata_1), new_name)
#   ordered_gms <- gene_modules[order(match(names(gene_modules), 
#                                           classified_gms_1$gene_module))] # ordered on metadata_1 but not renamed
#   if (is.null(metadata_2) || is.na(metadata_2) || is.nan(metadata_2)) {
#     print(paste0("Gene modules ordered only on ", metadata_1))
#   }
#   else {
#     print(paste0("Gene modules ordered on ", metadata_1, 
#                  " AND ", metadata_2))
#     temp_seurat <- SplitObject(seurat_data, split.by = metadata_1)
#     classified_gms_2 <- c()
#     for (i in order_1) {
#       print(i)
#       subset_gms <- gene_modules[classified_gms_1 %>% filter(!!sym(metadata_1) == 
#                                                                i) %>% dplyr::pull(gene_module)]
#       if (length(subset_gms) != 0) {
#         temp <- GeneModuleClassify(seurat_data, subset_gms, 
#                                    metadata = metadata_2, plot_path = paste0(plot_path, 
#                                                                              i, "/"))
#         classified_gms_2 <- rbind(classified_gms_2, temp)
#       }
#     }
#     
#     classified_gms_2 <- classified_gms_2 %>% add_column(!!sym(metadata_1) := classified_gms_1[[metadata_1]])
#     
#     classified_gms_2 <- classified_gms_2 %>% arrange(!!sym(metadata_1), (match(!!sym(metadata_2), order_2))) %>% 
#       mutate(pos = 1:n()) %>% mutate(new_name = paste(!!sym(metadata_2), pos, sep = "-")) %>% 
#       dplyr::select(gene_module, !!sym(metadata_2), !!sym(metadata_1), new_name)
#     
#     ordered_gms <- gene_modules[order(match(names(gene_modules), classified_gms_2$gene_module))]
#     
#     if (!is.null(rename_modules) && rename_modules == metadata_2) {
#       names(ordered_gms) <- classified_gms_2$new_name
#     }
#   }
#   if (!is.null(rename_modules) && rename_modules == metadata_1) {
#     names(ordered_gms) <- classified_gms_1$new_name
#   }
#   return(ordered_gms)
# }

### TEMP - overwrite pheatmap function to return dfs to pass to complex heatmap
GeneModulePheatmap <- function (peak_normalised_matrix, metadata, custom_order = NULL, 
                                custom_order_column = metadata[1],
                                gene_modules, gm_row_annotation = TRUE,
                                annotation_names_row = FALSE, scale_data = TRUE) 
{

  metadata <- metadata %>% mutate_if(is.character, as.factor)
  
  col_ann <- col_ann[do.call("order", c(col_ann[custom_order_column],
                                        list(decreasing = FALSE))), , drop = FALSE]
  if (!is.null(custom_order)) {
    if (!setequal(custom_order, unique(col_ann[[custom_order_column]]))) {
      stop("custom_order factors missing from custom_order_column \n\n")
    }
    levels(col_ann[[custom_order_column]]) <- custom_order
    col_ann <- col_ann[order(col_ann[[custom_order_column]]),
                       , drop = FALSE]
  }
  
  if (gm_row_annotation == TRUE) {
    row_ann <- stack(selected_GM) %>% rename(`Gene Modules` = ind) %>%
      column_to_rownames("values")
  }
  else {
    row_ann <- NA
  }
  
  plot_data <- t(peak_normalised_matrix)
  
  if (scale_data) {
    cat("Scaling data \n")
    plot_data <- t(scale(t(plot_data)))
  }
  
  plot_data <- replace(plot_data, plot_data >= 2, 2)
  plot_data <- replace(plot_data, plot_data <= -2, -2)
  
  output <- list(plot_data = plot_data,
                row_ann = row_ann,
                col_ann = col_ann,
                ann_colours = ann_colours)
  return(output)
}







########################################################################################################################################################
##################################################   Setting metadata and colours for heatmaps:   ################################################

# # Set metadata and order cells in heatmap
# metadata <- if(length(unique(seurat_data@meta.data$stage)) == 1){meta_col
# }else{c("stage", meta_col)}
# metadata <- if(length(unique(seurat_data@meta.data$run)) == 1){metadata
# }else{c(metadata, "run")}
# 
# # Allow manual setting of metadata using --force_order command line arg
# if(!is.null(opt$force_order)){metadata <- unlist(str_split(opt$force_order, pattern = ','))}

# 
# ########################################################################################################################################################
# ##################################################   Plotting heatmaps with ordered GMs:   #############################################################
# 
# # Extract ordering of gms from metadata
# labels <- c("stage", "scHelper_cell_type", "seurat_clusters")
# stage_order <- levels(seurat_data@meta.data$stage)
# scHelper_cell_type_order <- levels(seurat_data@meta.data$scHelper_cell_type)
# 
# metadata_1 <- NULL
# order_1 <- NULL
# metadata_2 <- NULL
# order_2 <- NULL
# 
# if(sum(labels %in% metadata) !=0){
#   metadata_1 <- labels[labels %in% metadata][1]
#   order_1 <- get(paste0(metadata_1, '_order'))
#   
#   if(sum(labels %in% metadata) == 2){
#     metadata_2 <- labels[labels %in% metadata][2]
#     order_2 <- get(paste0(metadata_2, '_order'))
#   }
# } else {
#   print(paste(c(labels, 'not found in metadata. GMs will not be ordered'), collapse = ' '))
# }
# 
# ##########################################################################################
# # plot all gene modules (unbiasedGMs)
# 
# # Order gms
# if (!is.null(metadata_1)){
#   antler_data$gene_modules$lists$unbiasedGMs$content <- GeneModuleOrder(seurat_obj = seurat_data, gene_modules = antler_data$gene_modules$lists$unbiasedGMs$content,
#                                                                         metadata_1 = metadata_1, order_1 = order_1,
#                                                                         metadata_2 = metadata_2, order_2 = order_2,
#                                                                         plot_path = "scHelper_log/GM_classification/unbiasedGMs/")
# }
# 
# # if stage not present in metadata, add it here so all heatmaps have the stage annotation
# if (!("stage" %in% metadata)){
#   metadata <- c("stage", metadata)}
# 
# # generate plot data
# plot_data <- GeneModulePheatmap(seurat_obj = seurat_data,  metadata = metadata, gene_modules = antler_data$gene_modules$lists$unbiasedGMs$content,
#                                 show_rownames = FALSE, col_order = metadata, col_ann_order = metadata, gaps_col = ifelse('stage' %in% metadata, 'stage', meta_col), fontsize = 15, fontsize_row = 10,
#                                 return = "plot_data")
# plot_data$ann_colours$scHelper_cell_type <- scHelper_cell_type_colours[names(plot_data$ann_colours$scHelper_cell_type)]
# plot_data$ann_colours$stage <- stage_colours[names(plot_data$ann_colours$stage)]
# 
# # Set annotations for heatmap
# if (!is.null(plot_data$ann_colours$run)){
#   if(length(plot_data$ann_colours$stage) > 1){
#     top_annotation <- HeatmapAnnotation(stage = anno_block(gp = gpar(fill = plot_data$ann_colours$stage),
#                                                            labels = levels(plot_data$col_ann$stage),
#                                                            labels_gp = gpar(col = "white", fontsize = 50, fontface='bold')),
#                                         run = anno_simple(x = as.character(plot_data$col_ann$run),
#                                                           col = plot_data$ann_colours$run, height = unit(1, "cm")),
#                                         simple_anno_size = unit(1, "cm"),
#                                         annotation_label = "Run", gp = gpar(fontsize = 35))
#   } else {
#     top_annotation <- HeatmapAnnotation(run = anno_simple(x = as.character(plot_data$col_ann$run),
#                                                           col = plot_data$ann_colours$run, height = unit(1, "cm")),
#                                         simple_anno_size = unit(1, "cm"),
#                                         annotation_label = "Run", gp = gpar(fontsize = 35))
#   }
# } else {
#   top_annotation = NULL
# }
# 
# png(paste0(plot_path, 'unbiasedGMs.png'), height = 150, width = 100, units = 'cm', res = 400)
# Heatmap(t(plot_data$plot_data), col = PurpleAndYellow(), cluster_columns = FALSE, cluster_rows = FALSE,
#         show_column_names = FALSE, column_title = NULL, show_row_names = FALSE, row_title_gp = gpar(fontsize = 45), row_title_rot = 90,
#         row_split = plot_data$row_ann$`Gene Modules`, column_split = if(length(plot_data$ann_colours$stage) > 1){
#           plot_data$col_ann$stage
#         }else{
#           plot_data$col_ann$scHelper_cell_type
#         },
#         heatmap_legend_param = list(
#           title = "Scaled expression", at = c(-2, 0, 2), 
#           labels = c(-2, 0, 2),
#           legend_height = unit(11, "cm"),
#           grid_width = unit(1.5, "cm"),
#           title_gp = gpar(fontsize = 35, fontface = 'bold'),
#           labels_gp = gpar(fontsize = 30, fontface = 'bold'),
#           title_position = 'lefttop-rot',
#           x = unit(13, "npc")
#         ),
#         top_annotation = top_annotation,
#         bottom_annotation = HeatmapAnnotation(scHelper_cell_type = anno_simple(x = as.character(plot_data$col_ann$scHelper_cell_type),
#                                                                                col = plot_data$ann_colours$scHelper_cell_type, height = unit(1, "cm")), show_annotation_name = FALSE,
#                                               labels = anno_mark(at = cumsum(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths) - floor(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths/2),
#                                                                  labels = rle(as.character(plot_data$col_ann$scHelper_cell_type))$values,
#                                                                  which = "column", side = 'bottom',
#                                                                  labels_gp = gpar(fontsize = 40), lines_gp = gpar(lwd=8))),
#         raster_quality = 8
# )
# graphics.off()
# 
# 
# 
