#!/usr/bin/env Rscript

print("Calculate average peak module accessibility across metacells and plots featureplots on UMAPs and co-accessibility plots")

############################## Load libraries #######################################
library(optparse)
library(future)
library(pheatmap)
library(tidyverse)
library(Antler)
library(RColorBrewer)
library(scHelper)
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

########## PEAK MODULE AVERAGE SCORES DF #############
PM_avg_scores <- read_csv(paste0(data_path, "rds_files/PM_avg_scores.csv"))
PM_avg_scores <- column_to_rownames(PM_avg_scores, "SEACell_ID")
print("PM average scores: ")
print(head(PM_avg_scores))

######### SEACell METADATA #############
metadata <- read_csv(paste0(data_path, "rds_files/Combined_SEACell_integrated_metadata_latent_time.csv"))
print("SEACell metadata:")
print(head(metadata))

########## ATAC METACELL SEURAT OBJECT ############# 
label <- setdiff(sub('_.*', '', list.files(data_path)), "seacells_seurat_integrated.RDS")
print(label)
seurat <- readRDS(paste0(data_path, label, "_seacells_seurat_integrated.RDS"))

print("Data read in!")

########################################################################################################
#                                 Add PMs as features to seurat object                               #
########################################################################################################

## split up the PM df by stage
split_df <- split(PM_avg_scores, substring(rownames(PM_avg_scores), nchar(rownames(PM_avg_scores)) - 3))
names(split_df) <- substr(names(split_df), 2, 4)

## detect which stage is the seurat object and extract the correct df
df <- split_df[[unique(seurat_data@meta.data$stage)]]
head(df)

## reorder df so the metacells are in the same order as they appear in the seurat metadata
df <- rownames_to_column(df, "SEACell_ID")
df <- df %>%
  mutate(SEACell_ID = gsub("_", "-", SEACell_ID)) %>%
  mutate(SEACell_ID = substr(SEACell_ID, 1, nchar(SEACell_ID) - 4))
df <- df %>% arrange(SEACell_ID)
head(df)

## check that the seacells match and are in the same order
if (sum(rownames(seurat_data@meta.data) == df$SEACell_ID) == nrow(seurat_data@meta.data)){
  "SEACell IDs match!"
} else {stop("ERROR! SEACell IDs dont match!")}

## add each average PM score as metadata
df <- df[,-1]
for (module in colnames(df)){
  seurat_data@meta.data[[module]] <- df[[module]]
  
  # feature plot of that PM
  png(paste0(plot_path, module, '_feature_plot.png'), width = 15, height = 15, units='cm', res=200)
  print(
    FeaturePlot(seurat_data, features = module, pt.size = 1.4) +
      theme_void() +
      theme(plot.title = element_blank(),
            legend.text = element_text(size=16),
            legend.key.size = unit(1, 'cm'))
  ) 
  graphics.off()
}



########################################################################################################
#                                 Plot co-accessibility of PMs                       #
########################################################################################################






########################################################################################################
#                                 Plot latent time on SEACells UMAPs                       #
########################################################################################################









