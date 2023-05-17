library(shiny)
library(bs4Dash)
library(data.table)
library(tidyverse)
library(viridis)
library(mgcv)
library(patchwork)
#library(shinycssloaders)

#source('./custom_functions.R')

# Don't know what this does
options(scipen = 1)
options(digits = 2)

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

############################################################################################
# Read in normalised matrix of SEACells x all peaks (temp with only 100 peaks)
SEACells_peak_matrix <- fread('../output/Rshiny_input_SEACells_matrix.csv')
SEACells_peak_matrix <- column_to_rownames(SEACells_peak_matrix, var = "V1")
SEACells_peak_matrix[1:2, 1:2]

# Read in metadata for all SEACells
SEACells_metadata <- fread('../output/Rshiny_input_metadata.csv')
SEACells_metadata <- SEACells_metadata %>% mutate(stage = sub(".*-", "", ATAC))
head(SEACells_metadata)
############################################################################################

# Potential ways to subset the SEACells matrix to make heatmaps
data_subsets = c("Full Data", "HH5", "HH6", "HH7", "ss4", "ss8")

my_theme <- theme(axis.text=element_text(size=14),
                  axis.title=element_text(size=16))