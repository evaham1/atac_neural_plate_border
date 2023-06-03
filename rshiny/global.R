library(shiny)
library(bs4Dash)
library(data.table)
library(tidyverse)
library(viridis)
library(mgcv)
library(patchwork)
library(Seurat)
library(ComplexHeatmap)
library(scHelper)
#library(shinycssloaders)

source('./custom_functions.R')

# Don't know what this does
options(scipen = 1)
options(digits = 2)

##########################################################################
###################       Read in data      ##############################

# Read in normalised matrix of SEACells x all peaks (temp with only 100 peaks)
SEACells_peak_matrix <- fread('../output/Rshiny_input_SEACells_matrix.csv')
SEACells_IDs <- SEACells_peak_matrix$V1
table(duplicated(SEACells_IDs))
SEACells_peak_matrix <- SEACells_peak_matrix[,-1]
SEACells_peak_matrix <- as.matrix(sapply(SEACells_peak_matrix, as.numeric))
rownames(SEACells_peak_matrix) <- SEACells_IDs
SEACells_peak_matrix[1:2, 1:2]
dim(SEACells_peak_matrix)

# Read in metadata for all SEACells
SEACells_metadata <- fread('../output/Rshiny_input_metadata.csv')
SEACells_metadata <- SEACells_metadata %>% mutate(stage = sub(".*-", "", ATAC))
head(SEACells_metadata)

# Read in test SEACells seurat object
SEACells_seurat <- readRDS('../output/NF-downstream_analysis/Processing/ss8/SEACELLS_ATAC_WF/5_Classify_metacells/rds_files/Classified_metacells.RDS')

# Read in peak modules
paths <- c("./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/PMs/FullData/unbiasedPMs.txt",
           "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/PMs/HH5/unbiasedPMs.txt",
           "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/PMs/HH6/unbiasedPMs.txt",
           "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/PMs/HH7/unbiasedPMs.txt",
           "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/PMs/ss4/unbiasedPMs.txt",
           "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/PMs/ss8/unbiasedPMs.txt")

dfs <- list()

for (path in paths){
  
  # Read the file line by line
  lines <- readLines(path)
  
  # Split each line into separate values
  data <- lapply(strsplit(lines, ","), function(x) as.data.frame(t(x), stringsAsFactors = FALSE))
  
  # Determine the maximum number of columns
  max_cols <- max(sapply(data, ncol))
  
  # Pad the data frames with empty values to ensure consistent number of columns
  data <- lapply(data, function(df) {
    if (ncol(df) < max_cols) {
      extra_cols <- max_cols - ncol(df)
      df <- cbind(df, matrix("", nrow = nrow(df), ncol = extra_cols))
    }
    colnames(df) <- paste0("V", seq_len(ncol(df)))  # Assign unique column names
    df
  })
  
  # Combine the list of data frames into a single data frame
  df <- do.call(rbind, data)
  
  df <- df %>% separate(V1, c("PM", "1"), sep = ";")
  df <- column_to_rownames(df, "PM")
  
  # Add to list of dfs
  dfs[[substr(path, 91, 93)]] <- df
  
}

PM_df <- bind_rows(dfs)

# Read in peaks found from HiChip
PPR_hichip_peaks <- read.csv(("../output/Rshiny_PPR_input_peaks.csv"))
PPR_hichip_peaks <- PPR_hichip_peaks$x
NC_hichip_peaks <- read.csv(("../output/Rshiny_NC_input_peaks.csv"))
NC_hichip_peaks <- NC_hichip_peaks$x

############################################################################
##################       Read in aesthetic params      #####################

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
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
names(stage_colours) <- stage_order

my_theme <- theme(axis.text=element_text(size=14),
                  axis.title=element_text(size=16))

############################################################################
####################       Read in UI options      ##########################

# Potential ways to subset the SEACells matrix to make heatmaps
data_subsets = c("Full Data", "HH5", "HH6", "HH7", "ss4", "ss8")

# Potential ways to group UMAP
dimplot_groupby_options = c("seurat_clusters", "scHelper_cell_type")

