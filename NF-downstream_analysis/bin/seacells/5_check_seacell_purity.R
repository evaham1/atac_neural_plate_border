##### Putting the SEACell assignments onto ArchR and checking their purity
# also adding metadata to the SEACells themselves and exporting this for visualisations of peak modules


## basically it looks like the metacells are pretty pure in terms of stage but mixed in terms of cell type
# is this reflected in the purity plots made by seacells in python?
# can we explore this further, eg if combine PPR cell states are metacells pretty pure?
# does this mean we shouldnt use them or need to be careful?? need to look into this


# load libraries
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(hexbin)
library(gridExtra)
library(grid)
library(parallel)
library(data.table)
library(doSNOW)

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
    addArchRThreads(threads = 1) 
    
    # for transfer labels archr object
    data_path = "./output/NF-downstream_analysis/Processing/TransferLabels/3_peak_call/rds_files/"
    # for metacell assignments
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/SEACells/4_exported_SEACells_data/rds_files/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


############################## Functions #######################################

# function to make heatmap showing contribution of cell groups to other cell groups
cell_counts_heatmap <- function(ArchR = ArchR, group1 = "scHelper_cell_type_new", group2 = "clusters") {
  group1_data <- getCellColData(ArchR, select = group1)[,1]
  group2_data <- getCellColData(ArchR, select = group2)[,1]
  cM <- confusionMatrix(paste0(group2_data), paste0(group1_data))
  cM <- cM / Matrix::rowSums(cM)
  
  p <- pheatmap::pheatmap(
    mat = cM,
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
  )
}

############################## Set colours #######################################

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

###### schelper cell type colours
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

############################## Read in ArchR TL project + Metacell assignments #######################################

ArchR <- loadArchRProject(path = paste0(data_path, "TransferLabels_Save-ArchR"), force = FALSE, showLogo = TRUE)

SEACells_metadata <- read_csv(paste0(data_path, "AnnData_metacells_assigned_cell_metadata.csv"))

############################## Read in SEACells assignments and add to ArchR object #######################################

# extract just SEACell assignments for each cell id
SEACells_cell_assignments <- SEACells_metadata %>% select(index, SEACell)

# check these match
dim(SEACells_cell_assignments)
dim(getCellColData(ArchR))

# add SEACell assignments to ArchR object - these are a cell id
ArchR <- addCellColData(ArchRProj = ArchR, data = SEACells_cell_assignments$SEACell,
                        cells = SEACells_cell_assignments$index, name = "SEACell")

# convert the seacell assignment into a unique number for ease
seacell_dictionary <- data.frame(unique(SEACells_cell_assignments$SEACell), 1:length(unique(SEACells_cell_assignments$SEACell)))
colnames(seacell_dictionary) <- c("SEACell", "SEACell_ID")
SEACells_cell_assignments_and_ids <- merge(SEACells_cell_assignments, seacell_dictionary, by = "SEACell")
ArchR <- addCellColData(ArchRProj = ArchR, data = SEACells_cell_assignments_and_ids$SEACell_ID,
                        cells = SEACells_cell_assignments$index, name = "SEACell_ID")

# add scHelpercelltype col that is same as stage_scHelper_cell_type but without stage label
ArchR <- addCellColData(ArchRProj = ArchR, data = sub(".*_", "", getCellColData(ArchR)$stage_scHelper_cell_type_old),
                        cells = rownames(getCellColData(ArchR)), name = "scHelper_cell_type")

getCellColData(ArchR)

############################## Visualise how metacells spread across data #######################################

###### plot seacells on UMAP next to scHelper_cell_type and stage
p1 <- plotEmbedding(ArchR, 
                    name = "SEACell",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    randomize = TRUE)
p2 <- plotEmbedding(ArchR, 
                    name = "stage",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0,
                    randomize = TRUE, pal = stage_colours)
p3 <- plotEmbedding(ArchR, 
                    name = "scHelper_cell_type",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0,
                    randomize = TRUE, pal = scHelper_cell_type_colours)

png(paste0(plot_path_temp, "UMAPs.png"), width=100, height=40, units = 'cm', res = 200)
ggAlignPlots(p1, p2, p3, type = "h")
graphics.off()

############################## Visualise seacell purity - confusion matrices #######################################

png(paste0(plot_path, "SEACell_by_stage_distribution.png"), width=25, height=20, units = 'cm', res = 200)
cell_counts_heatmap(ArchR = ArchR, group1 = "SEACell", group2 = "stage")
graphics.off()

png(paste0(plot_path, "SEACell_by_cell_type_distribution.png"), width=25, height=20, units = 'cm', res = 200)
cell_counts_heatmap(ArchR = ArchR, group1 = "SEACell", group2 = "scHelper_cell_type")
graphics.off()


############################## Explore seacell purity - stages #######################################

SEACells_stages <- data.frame(getCellColData(ArchR)$SEACell_ID, getCellColData(ArchR)$stage)
colnames(SEACells_stages) <- c("SEACell_ID", "stage")
dim(SEACells_stages)
length(unique(SEACells_stages$SEACell_ID))
length(unique(SEACells_stages$stage))

stages_freq <- SEACells_stages %>%
  group_by(SEACell_ID, stage) %>% 
  dplyr::summarize(count = n())
dim(stages_freq)
length(unique(stages_freq$SEACell_ID))
length(unique(stages_freq$stage))

# these are the metacells that include cells of more than one stage - theres only 4!
mixed_metacells <- stages_freq$SEACell_ID[duplicated(stages_freq$SEACell_ID)]

# and these are their relative distributions - they are generally of neibouring stages and mainly one stage
subset(stages_freq, SEACell_ID %in% mixed_metacells)

############################## Explore seacell purity - scHelper_cell_type #######################################

SEACells_cell_type <- data.frame(getCellColData(ArchR)$SEACell_ID, getCellColData(ArchR)$scHelper_cell_type)
colnames(SEACells_cell_type) <- c("SEACell_ID", "cell_type")
dim(SEACells_cell_type)
length(unique(SEACells_cell_type$SEACell_ID))
length(unique(SEACells_cell_type$cell_type))

cell_type_freq <- SEACells_cell_type %>%
  group_by(SEACell_ID, cell_type) %>% 
  dplyr::summarize(count = n())
dim(cell_type_freq)
length(unique(cell_type_freq$SEACell_ID))
length(unique(cell_type_freq$cell_type))

# these are the metacells that include cells of more than one stage - theres only 4!
mixed_metacells <- cell_type_freq$SEACell_ID[duplicated(cell_type_freq$SEACell_ID)]
length(unique(mixed_metacells))

# and these are their relative distributions - they are generally of neibouring stages and mainly one stage
subset(cell_type_freq, SEACell_ID %in% mixed_metacells)

# they are pretty mixed by cell type!!


