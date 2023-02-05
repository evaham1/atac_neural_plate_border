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
    # output from SEACells - summarised by metacells
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/SEACells/4_exported_SEACells_data/rds_files/"
    label = "AnnData_summarised_by_metacells_peak_counts.csv"
    
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

# function to take ArchR object and a category and make a freq table of how frequently metacells found in each category
calculate_seacells_frequencies <- function(ArchR, metacell_slot = "SEACell_ID", category = "stage"){
  
  df <- data.frame(getCellColData(ArchR, select = metacell_slot), getCellColData(ArchR, select = category))
  colnames(df) <- c("Metacell", "Category")
  
  freq <- df %>%
    group_by(Metacell, Category) %>% 
    dplyr::summarize(count = n())
  
  print(paste0("Number of metacells: ", length(unique(freq$Metacell))))
  print(paste0("Number of categories: ", length(unique(freq$Category))))
  
  return(freq)
}

# function to plot a piechart showing how many cells are mixed identity and how many are not
mixed_metacells_piechart <- function(mixed_metacell_number, total_metacell_number = 1000, plot_path = plot_path){
  slices <- c(mixed_metacell_number, total_metacell_number-mixed_metacell_number)
  lbls <- c(paste0("Mixed: ", slices[1]), paste0("Not mixed: ", slices[2]))
  png(paste0(plot_path, "Pie_number_of_mixed_metacells.png"), width=25, height=20, units = 'cm', res = 200)
  pie(slices, labels = lbls, main="How many metacells with mixed identities")
  graphics.off()
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

# add more general scHelpercelltype labels (ie group together all PPRs, neurals, etc)
scHelper_cell_types <- data.frame(getCellColData(ArchR, select = "scHelper_cell_type"))
broad <- scHelper_cell_types %>%
  mutate(broad = mapvalues(scHelper_cell_type, 
                                               from=c("NP", "aNP", "iNP", "pNP", "eN", "vFB", "FB", "MB", "HB",
                                                      'PPR', 'aPPR', 'pPPR',
                                                      'eNPB', 'NPB', 'aNPB', 'pNPB',
                                                      'NC', 'dNC'),
                                               to=c(rep("Neural", 9), rep("Placodal", 3), rep("NPB", 4), rep("NC", 2))))


ArchR <- addCellColData(ArchRProj = ArchR, data = broad$broad, cells = rownames(getCellColData(ArchR)), name = "scHelper_cell_type_broad")


print("Cell metadata added to ArchR object!")
print(getCellColData(ArchR))

############################## Save ArchR with new obs #######################################
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, "TransferLabels_Save-ArchR"), load = FALSE)
print("ArchR object saved")

############################## Read and write summarised counts by metacells #######################################

SEACells_summarised <- as.matrix(fread(paste0(data_path, label), header = TRUE), rownames = 1)
write.csv(SEACells_summarised, paste0(rds_path, "SEACells_summarised.csv"))


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


############################## Explore seacell purity - stages #######################################
plot_path = "./plots/stages/"
dir.create(plot_path)

# calculate frequencies in which metacells are in each stage
stage_freq <- calculate_seacells_frequencies(ArchR, metacell_slot = "SEACell_ID", category = "stage")

# these are the metacells that include cells of more than one stage - theres only 4!
mixed_metacells_stages <- stage_freq$Metacell[duplicated(stage_freq$Metacell)]
print(paste0("Number of metacells that have mixed stage identity: ", length(unique(mixed_metacells_stages))))

mixed_metacells_piechart(mixed_metacell_number = length(unique(mixed_metacells_stages)), plot_path = plot_path)

# and these are their relative distributions - they are generally of neighbouring stages and mainly one stage
print(subset(stage_freq, Metacell %in% mixed_metacells_stages))
mixed_cell_table_stages <- subset(stage_freq, Metacell %in% mixed_metacells_stages)
how_many_cell_types_stages <- table(mixed_cell_table_stages$Metacell)
print(how_many_cell_types_stages)

png(paste0(plot_path, "Hist_mixed_stages.png"), width=25, height=20, units = 'cm', res = 200)
hist(how_many_cell_types_stages)
graphics.off()

############################## Explore seacell purity - scHelper_cell_type #######################################

plot_path = "./plots/scHelper_cell_type/"
dir.create(plot_path)

# calculate frequencies in which metacells are in each stage
celltype_freq <- calculate_seacells_frequencies(ArchR, metacell_slot = "SEACell_ID", category = "scHelper_cell_type")

# these are the metacells that include cells of more than one cell type - its all 1000 of them!
mixed_metacells_celltype <- celltype_freq$Metacell[duplicated(celltype_freq$Metacell)]
print(paste0("Number of metacells that have mixed scHelper_cell_type identity: ", length(unique(mixed_metacells_celltype))))

mixed_metacells_piechart(mixed_metacell_number = length(unique(mixed_metacells_celltype)), plot_path = plot_path)

# and these are their relative distributions - they are generally of neibouring stages and mainly one stage
mixed_cell_table_celltype <- subset(celltype_freq, Metacell %in% mixed_metacells_celltype)
how_many_cell_types_celltype <- table(mixed_cell_table_celltype$Metacell)

png(paste0(plot_path, "Hist_mixed_cell_types.png"), width=25, height=20, units = 'cm', res = 200)
hist(how_many_cell_types_celltype)
graphics.off()

# they are pretty mixed by cell type!!

############################## Explore seacell purity - scHelper_cell_type #######################################

plot_path = "./plots/scHelper_cell_type_broad/"
dir.create(plot_path)

# calculate frequencies in which metacells are in each stage
celltype_broad_freq <- calculate_seacells_frequencies(ArchR, metacell_slot = "SEACell_ID", category = "scHelper_cell_type_broad")

# these are the metacells that include cells of more than one cell type - its all 1000 of them!
mixed_metacells_celltype_broad <- celltype_broad_freq$Metacell[duplicated(celltype_broad_freq$Metacell)]
print(paste0("Number of metacells that have mixed scHelper_cell_type_broad identity: ", length(unique(mixed_metacells_celltype_broad))))

mixed_metacells_piechart(mixed_metacell_number = length(unique(mixed_metacells_celltype_broad)), plot_path = plot_path)

# and these are their relative distributions - they are generally of neibouring stages and mainly one stage
mixed_cell_broad_table_celltype <- subset(celltype_broad_freq, Metacell %in% mixed_metacells_celltype_broad)
how_many_cell_types_celltype <- table(mixed_cell_broad_table_celltype$Metacell)

png(paste0(plot_path, "Hist_mixed_cell_types.png"), width=25, height=20, units = 'cm', res = 200)
hist(how_many_cell_types_celltype)
graphics.off()

# they are pretty mixed by cell type!!


