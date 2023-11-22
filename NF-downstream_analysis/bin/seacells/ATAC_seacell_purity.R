##### Putting the SEACell assignments onto ArchR and checking their purity
# also adding metadata to the SEACells themselves and exporting this for visualisations of peak modules

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
library(readr)
library(scHelper)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-t", "--categories"), action = "store", type = "character", help = "Which categories to use to check for purity", default = "clusters,scHelper_cell_type_old"),
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
    
    # for ArchR object
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Integrated_SEACells_label_transfer/rds_files/"
    # for metacell IDs csv
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Integrated_SEACells_label_transfer/rds_files/ss8_ATAC_singlecell_integration_map.csv"

    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    label = "AnnData_summarised_by_metacells_peak_counts.csv"
    
    addArchRThreads(threads = ncores)
    
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
                              'PGC', 'BI', 'meso', 'endo', 'MIXED', 'Unmapped',
                              'Neural', 'Placodal', 'Non-neural', 'Contam')
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

############################## Read in ArchR TL project + Metacell assignments #######################################

# If files are not in rds_files subdirectory look in input dir 
label <- unique(sub('_.*', '', list.files(data_path)))
print(label) 

if (length(label) == 0){
  data_path = "./input/"
  label <- sub('_.*', '', list.files(data_path))
  print(label)
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
} else {
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
}

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
head(getCellColData(ArchR))
colnames(getCellColData(ArchR))
head(getPeakSet(ArchR))
getAvailableMatrices(ArchR)

# read in metacell assignments
SEACells_metadata <- as.data.frame(read_delim(paste0(data_path, label, "_ATAC_singlecell_integration_map.csv"), delim = ","))
print(head(SEACells_metadata))

############################## Add Metacell data to seurat object #######################################

# check these match
dim(SEACells_metadata)
dim(getCellColData(ArchR)) 

# add ATAC SEACell assignments to ArchR object
ArchR <- addCellColData(ArchRProj = ArchR, data = paste0(SEACells_metadata$SEACell, "-", label),
                        cells = SEACells_metadata$index, name = "SEACell_ID", force = TRUE)

# add closest RNA SEACell ID to ArchR object
ArchR <- addCellColData(ArchRProj = ArchR, data = paste0(SEACells_metadata$Integrated_RNA_SEACell_ID, "-", label),
                        cells = SEACells_metadata$index, name = "RNA_SEACell_ID", force = TRUE)

# add SEACell integrated labels to ArchR object
ArchR <- addCellColData(ArchRProj = ArchR, data = SEACells_metadata$scHelper_cell_type_by_proportion,
                        cells = SEACells_metadata$index, name = "SEACell_scHelper_cell_type", force = TRUE)

# add broad SEACell integrated labels to ArchR object
ArchR <- addCellColData(ArchRProj = ArchR, data = SEACells_metadata$scHelper_cell_type_broad_by_proportion,
                        cells = SEACells_metadata$index, name = "SEACell_scHelper_cell_type_broad", force = TRUE)

print("Cell metadata added to ArchR object!")
print(getCellColData(ArchR))

############################## Visualise SEACell labels on single cells #######################################

p1 <- plotEmbedding(ArchR, 
                    name = "SEACell_ID",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    randomize = TRUE)
png(paste0(plot_path, "SEACell_IDs_UMAPs.png"), width=60, height=100, units = 'cm', res = 200)
print(p1)
graphics.off()


p1 <- plotEmbedding(ArchR, 
                    name = "SEACell_scHelper_cell_type",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    pal = scHelper_cell_type_colours, randomize = TRUE)
p2 <- plotEmbedding(ArchR, 
                    name = "transferred_scHelper_cell_type",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0,
                    pal = scHelper_cell_type_colours, randomize = TRUE)

png(paste0(plot_path, "scHelper_cell_state_UMAPs.png"), width=60, height=40, units = 'cm', res = 200)
ggAlignPlots(p1, p2, type = "h")
graphics.off()

p1 <- plotEmbedding(ArchR, 
                    name = "SEACell_scHelper_cell_type_broad",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    pal = scHelper_cell_type_colours, randomize = TRUE)
p2 <- plotEmbedding(ArchR, 
                    name = "transferred_scHelper_cell_type_broad",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0,
                    pal = scHelper_cell_type_colours, randomize = TRUE)

png(paste0(plot_path, "scHelper_cell_state_broad_UMAPs.png"), width=60, height=40, units = 'cm', res = 200)
ggAlignPlots(p1, p2, type = "h")
graphics.off()

############################## Explore SEACell purity #######################################

# categories <- strsplit(opt$categories, ",")[[1]]

# ## loop through categories and check each of them for purity in metacells
# for (cat in categories) {
  
#   print(paste0("Checking purity of ", cat))
#   if (!(cat %in% colnames(getCellColData(ArchR)))){
#     stop("Category not in ArchR cell col data!")}
  
#   plot_path_temp = paste0(plot_path, cat, "/")
#   dir.create(plot_path_temp, recursive = T)
  
#   # calculate proportions
#   prop_table <- SEACells_MetacellFrequencies(input_data = ArchR, input_data_type = "ArchR", metacell_slot = "SEACell_ID", category = cat, calc_proportions = TRUE)
  
#   # plot the relative proportions of labels in each metacell
#   png(paste0(plot_path_temp, "Hist_all_proportions.png"), width=25, height=20, units = 'cm', res = 200)
#   hist(prop_table$prop)
#   graphics.off()
  
#   # plot the max relative proportions of labels in each metacell
#   max_prop_table <- prop_table %>% group_by(Metacell) %>% dplyr::summarise(prop = max(prop))
#   png(paste0(plot_path_temp, "Hist_max_proportions_per_metacell.png"), width=25, height=20, units = 'cm', res = 200)
#   hist(max_prop_table$prop)
#   graphics.off()
  
#   ## how many metacells have > 50% of their cells from same label
#   png(paste0(plot_path_temp, "Pie_prop_over_0.5.png"), width=25, height=20, units = 'cm', res = 200)
#   SEACells_PiechartProportionThreshold(prop_table, threshold = 0.5)
#   graphics.off()
  
#   ## how many metacells have > 75% of their cells from same label
#   png(paste0(plot_path_temp, "Pie_prop_over_0.75.png"), width=25, height=20, units = 'cm', res = 200)
#   SEACells_PiechartProportionThreshold(prop_table, threshold = 0.75)
#   graphics.off()
  
#   ## how many metacells have > 90% of their cells from same label
#   png(paste0(plot_path_temp, "Pie_prop_over_0.9.png"), width=25, height=20, units = 'cm', res = 200)
#   SEACells_PiechartProportionThreshold(prop_table, threshold = 0.9)
#   graphics.off()
  
# }


############################## Save ArchR with new obs #######################################

paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)
print("ArchR object saved")