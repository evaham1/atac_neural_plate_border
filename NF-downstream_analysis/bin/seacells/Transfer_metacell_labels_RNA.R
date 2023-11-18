##### Transferring metacell labels onto single RNA cells in seurat object

# load libraries
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
  make_option(c("-k", "--k"), action = "store", type = "integer", help = "How many clusters to split metacells into"),
  make_option(c("-t", "--categories"), action = "store", type = "character", help = "Which categories to use to check for purity", default = "stage,clusters,scHelper_cell_type,scHelper_cell_type_broad"),
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
    addArchRThreads(threads = 1)
    
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


############################## Functions #######################################

### need to update schelper with this function
SEACells_MetacellFrequencies <- function(input_data, input_data_type, metacell_slot = "SEACell_ID", category = "clusters", calc_proportions = FALSE){
  
  # check if input data is in ArchR or Seurat format
  ifelse(input_data_type %in% c("seurat", "ArchR"), print(""), stop("Input data type should be seurat or ArchR!"))
  
  # if data input type is seurat...
  if (input_data_type == "seurat"){
    print("input data: seurat")
    df <- data.frame(FetchData(object = seurat, vars = c(metacell_slot, category)))
    colnames(df) <- c("Metacell", "Category")
  }
  
  # if data input type is ArchR...
  if (input_data_type == "ArchR"){
    print("input data: ArchR")
    df <- data.frame(getCellColData(ArchR, select = c(metacell_slot, category)))
    colnames(df) <- c("Metacell", "Category")
  }
  
  df <- df %>%
    group_by(Metacell, Category) %>% 
    dplyr::summarize(count = n()) %>%
    mutate(Metacell = str_split(Metacell, "-", simplify = TRUE)[ , 2]) %>%
    mutate(Metacell = as.numeric(Metacell)) %>% 
    arrange(Metacell) %>% 
    mutate(count = as.numeric(count))
  
  print(paste0("Number of metacells: ", length(unique(df$Metacell))))
  print(paste0("Number of categories: ", length(unique(df$Category))))
  
  if (calc_proportions){
    
    # calculate total cell counts per metacell
    totals_df <- aggregate(count ~ Metacell, data=df, FUN=sum)
    totals <- totals_df$count
    print(length(totals))
    
    # calculate proportions per metacell
    prop_table <- df %>%
      mutate(prop = count/totals[Metacell+1])
    
    df <- prop_table
  }
  
  return(df)
  
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

stage_order <- c("HH4", "HH5", "HH6", "HH7", "ss4", "ss8")
############################################################################################

############################## Set colours #######################################

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

###### schelper cell type colours ~ 29 cell states
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

############################## Read in seurat object + Metacell assignments #######################################

# read in full seurat object
seurat <- readRDS(paste0(data_path, "seurat.RDS"))
seurat
head(seurat@meta.data)

# read in metacell seurat object
seurat_metacells <- readRDS(paste0(data_path, "/rds_files/", "Classified_metacells.RDS"))
seurat_metacells
head(seurat_metacells@meta.data)

############################## Add Metacell schelper cell type labels to full seurat object #######################################

# extract single cell to metacell dictionary
metacell_dictionary <- dplyr::select(seurat@meta.data, c("SEACell"))
metacell_dictionary <- rownames_to_column(metacell_dictionary, var = "Cell_ID")

# extract metacell to schelper cell type dictionary
metacell_idents <- dplyr::select(seurat_metacells@meta.data, c("scHelper_cell_type"))
metacell_idents <- rownames_to_column(metacell_idents, var = "SEACell")

# project the de novo metacell annotations onto single cells
metacell_identity_dictionary <- merge(metacell_dictionary, metacell_idents, by = "SEACell")
metacell_identity_dictionary <- metacell_identity_dictionary[match(rownames(seurat@meta.data), metacell_identity_dictionary$Cell_ID),]
rownames(metacell_identity_dictionary) <- NULL
metacell_identity_dictionary <- column_to_rownames(metacell_identity_dictionary, var = "Cell_ID")

# add this back to seurat object
seurat <- AddMetaData(seurat, metacell_identity_dictionary$scHelper_cell_type, col.name = "SEACell_scHelper_cell_type")

############################## Check single cell identity across 3 annotations #######################################

# count matching labels across single cells
all_labels <- dplyr::select(seurat@meta.data, c("scHelper_cell_type", "SEACell_identity", "SEACell_scHelper_cell_type"))

assesment_proportions_rows <- all_labels %>%
  filter(scHelper_cell_type == SEACell_identity) %>%
  nrow()
assesment_denovo_rows <- all_labels %>%
  filter(scHelper_cell_type == SEACell_scHelper_cell_type) %>%
  nrow()

# print out final matching labels at a single cell level
# count how many mixed single cells
df <- data.frame(c("Total number of single cells", 
                   "Number of single cells with correct label - proportion based",
                   "Number of single cells with correct label - de novo based"),
                 c(nrow(all_labels), assesment_proportions_rows, assesment_denovo_rows)
)
colnames(df) <- NULL

png(paste0(plot_path, 'singlecell_label_correspondance.png'), height = 10, width = 30, units = 'cm', res = 400)
grid.arrange(tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# UMAP for singlecell cell states
scHelper_cols <- scHelper_cell_type_colours[levels(droplevels(unique(seurat@meta.data$scHelper_cell_type)))]

png(paste0(plot_path, "scHelper_original_umap.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'scHelper_cell_type', label = TRUE, 
        label.size = ifelse(length(unique(seurat$stage)) == 1, 9, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(seurat$stage)) == 1, 6, 6), 
        cols = scHelper_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

# UMAP for proportion based identity
scHelper_cols <- scHelper_cell_type_colours[unique(seurat@meta.data$SEACell_identity)]

png(paste0(plot_path, "metacell_scHelper_proportion_umap_metacell.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'SEACell_identity', label = TRUE, 
        label.size = ifelse(length(unique(seurat$stage)) == 1, 9, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(seurat$stage)) == 1, 6, 6), 
        cols = scHelper_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

# UMAP for metacell cell states
scHelper_cols <- scHelper_cell_type_colours[levels(droplevels(unique(seurat@meta.data$SEACell_scHelper_cell_type)))]

png(paste0(plot_path, "metacell_scHelper_denovo_umap_metacell.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat, group.by = 'SEACell_scHelper_cell_type', label = TRUE, 
        label.size = ifelse(length(unique(seurat$stage)) == 1, 9, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(seurat$stage)) == 1, 6, 6), 
        cols = scHelper_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()
