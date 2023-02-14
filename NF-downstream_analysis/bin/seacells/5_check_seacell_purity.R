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
    
    # for transfer labels archr object
    data_path = "./output/NF-downstream_analysis/Processing/TransferLabels/3_peak_call/rds_files/"
    # for metacell assignments
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/SEACells/4_exported_SEACells_data/rds_files/"
    # output from SEACells - summarised by metacells
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/SEACells/4_exported_SEACells_data/rds_files/"
    label = "summarised_counts_1000.csv"
    
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

# function to take ArchR object and a category and make a freq table of how frequently metacells found in each category
calculate_metacell_frequencies <- function(ArchR, metacell_slot = "SEACell_ID", category = "stage"){
  
  df <- data.frame(getCellColData(ArchR, select = metacell_slot), getCellColData(ArchR, select = category))
  colnames(df) <- c("Metacell", "Category")
  
  if (metacell_slot == "SEACell_ID"){
    freq <- df %>%
      group_by(Metacell, Category) %>% 
      dplyr::summarize(count = n()) %>%
      mutate(Metacell = as.numeric(Metacell)) %>% 
      arrange(Metacell) %>% 
      mutate(count = as.numeric(count))
  } else {
    if (metacell_slot == "SEACell_cluster"){
      freq <- df %>%
        group_by(Metacell, Category) %>% 
        dplyr::summarize(count = n()) %>%
        mutate(Metacell = str_split(Metacell, "_", simplify = TRUE)[ , 3]) %>%
        mutate(Metacell = as.numeric(Metacell)) %>% 
        arrange(Metacell) %>% 
        mutate(count = as.numeric(count))
    } else {
      stop("metacell_slot must be either SEACell_ID or SEACell_cluster!")
    }
  }
  
  print(paste0("Number of metacells: ", length(unique(freq$Metacell))))
  print(paste0("Number of categories: ", length(unique(freq$Category))))
  
  return(freq)
}

## Function to take freq table and turn it into proportions per metacell
calculate_metacell_proportions <- function(freq_table){
  
  # calculate total cell counts per metacell
  totals_df <- aggregate(count ~ Metacell, data=freq_table, FUN=sum)
  totals <- totals_df$count
  print(length(totals))
  
  # calculate proportions per metacell
  prop_table <- freq_table %>%
    mutate(prop = count/totals[Metacell])
  
  return(prop_table)
}

## Function to plot piechart of how many metacells pass threshold for proportion of cells coming from one label
piechart_proportion_threshold <- function(prop_table, threshold = 0.5){
  
  # filter cells to only include those that pass threshold
  passed_cells <- prop_table %>% filter(prop > threshold)
  # number of cells that pass threshold:
  passed_cells <- length(unique(passed_cells$Metacell))
  # number of cells that didn't pass threshold:
  failed_cells <- length(unique(prop_table$Metacell)) - passed_cells
  # plot piechart
  slices <- c(passed_cells, failed_cells)
  lbls <- c(paste0("Passed: ", passed_cells), paste0("Didn't pass: ", failed_cells))
  
  return(pie(slices, labels = lbls, main = paste0("Number of cells that passed threshold (", threshold, ") of label proportions")))
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



############################## Explore individual seacell purity #######################################

categories <- strsplit(opt$categories, ",")[[1]]

## loop through categories and check each of them for purity in metacells
for (cat in categories) {
  
  print(paste0("Checking purity of ", cat))
  if (!(cat %in% colnames(getCellColData(ArchR)))){
    stop("Category not in ArchR cell col data!")}
  
  plot_path_temp = paste0(plot_path, cat, "/")
  dir.create(plot_path_temp, recursive = T)
  
  # calculate frequencies
  freq_table <- calculate_metacell_frequencies(ArchR, metacell_slot = "SEACell_ID", category = cat)
  
  # calculate proportions of labels in metacells
  prop_table <- calculate_metacell_proportions(freq_table)
  
  # plot the relative proportions of labels in each metacell
  png(paste0(plot_path_temp, "Hist_all_proportions.png"), width=25, height=20, units = 'cm', res = 200)
  hist(prop_table$prop)
  graphics.off()
  
  # plot the max relative proportions of labels in each metacell
  max_prop_table <- prop_table %>% group_by(Metacell) %>% dplyr::summarise(prop = max(prop))
  png(paste0(plot_path_temp, "Hist_max_proportions_per_metacell.png"), width=25, height=20, units = 'cm', res = 200)
  hist(max_prop_table$prop)
  graphics.off()
  
  ## how many metacells have > 50% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.5.png"), width=25, height=20, units = 'cm', res = 200)
  piechart_proportion_threshold(prop_table, threshold = 0.5)
  graphics.off()
  
  ## how many metacells have > 75% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.75.png"), width=25, height=20, units = 'cm', res = 200)
  piechart_proportion_threshold(prop_table, threshold = 0.75)
  graphics.off()
  
  ## how many metacells have > 90% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.9.png"), width=25, height=20, units = 'cm', res = 200)
  piechart_proportion_threshold(prop_table, threshold = 0.9)
  graphics.off()
  
}

############################## Read summarised counts by metacells #######################################

SEACells_summarised <- as.matrix(fread(paste0(data_path, label), header = TRUE), rownames = 1)
print("Summarised count data read in!")

####################################  Cluster metacells #######################################

plot_path = "plots/seacell_clusters/"
dir.create(plot_path, recursive = T)

## normalise each metacell by the total number of cut sites
normalised_counts <- t(apply(SEACells_summarised, 1, function(x) x/sum(x))) * 1000

## calculate new variance and plot 
variance <- apply(SEACells_summarised, 2, var)
print("Before normalising:")
print(summary(variance))

png(paste0(plot_path, "hist_variance_before_normalising.png"), width=60, height=40, units = 'cm', res = 200)
hist(variance, breaks = 1000)
graphics.off()

variance <- apply(normalised_counts, 2, var)
print("After normalising:")
print(summary(variance))

png(paste0(plot_path, "hist_variance_after_normalising.png"), width=60, height=40, units = 'cm', res = 200)
hist(variance, breaks = 1000)
graphics.off()

## cluster metacells by hclust
corr_mat <- cor(t(normalised_counts), method = "spearman")
dim(corr_mat)

diss_matrix <- as.dist(1 - corr_mat)
tree <- hclust(diss_matrix, method="complete")

png(paste0(plot_path, "hclust_tree.png"), width=100, height=40, units = 'cm', res = 200)
plot(tree)
graphics.off()

# split tree into clusters based on k
metacell_clusters <- as.data.frame(cutree(tree, k = opt$k))
metacell_clusters <- rownames_to_column(metacell_clusters)
colnames(metacell_clusters) <- c("SEACell", "SEACell_cluster")
metacell_clusters <- metacell_clusters %>%
  mutate(SEACell_cluster = paste0("seacell_cluster_", SEACell_cluster))
head(metacell_clusters)
dim(metacell_clusters)

# add new metacell cluster IDs to ArchR object
metacell_clusters <- merge(metacell_clusters, SEACells_cell_assignments)
ArchR <- addCellColData(ArchRProj = ArchR, data = metacell_clusters$SEACell_cluster,
                        cells = SEACells_cell_assignments$index, name = "SEACell_cluster", force = TRUE)

print("SEACells clustered!")
print(paste0("Number of clusters made: ", length(unique(metacell_clusters$SEACell_cluster))))

############################## Plot seacell clusters on UMAP #######################################

p1 <- plotEmbedding(ArchR, 
                    name = "SEACell_cluster",
                    plotAs = "points", size = 1,
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    randomize = TRUE)

png(paste0(plot_path, "SEACell_cluster_UMAP.png"), width=60, height=40, units = 'cm', res = 200)
print(p1)
graphics.off()

############################## Plot seacell clusters size distribution #######################################

seacell_cluster_sizes <- as.data.frame(table(getCellColData(ArchR, select = "SEACell_cluster"))) %>%
  mutate(Freq = as.numeric(Freq))

png(paste0(plot_path, "SEACell_cluster_sizes_hist.png"), width=60, height=40, units = 'cm', res = 200)
hist(seacell_cluster_sizes$Freq, breaks = opt$k+10)
graphics.off()

############################## Explore seacell cluster purity #######################################

## loop through categories and check each of them for purity in metacells
for (cat in categories) {
  
  print(paste0("Checking purity of ", cat))
  if (!(cat %in% colnames(getCellColData(ArchR)))){
    stop("Category not in ArchR cell col data!")}
  
  plot_path_temp = paste0(plot_path, cat, "/")
  dir.create(plot_path_temp, recursive = T)
  
  # calculate frequencies
  freq_table <- calculate_metacell_frequencies(ArchR, category = cat, metacell_slot = "SEACell_cluster")
  
  # calculate proportions of labels in metacells
  prop_table <- calculate_metacell_proportions(freq_table)
  
  # plot the relative proportions of labels in each metacell
  png(paste0(plot_path_temp, "Hist_all_proportions.png"), width=25, height=20, units = 'cm', res = 200)
  hist(prop_table$prop)
  graphics.off()
  
  # plot the max relative proportions of labels in each metacell
  max_prop_table <- prop_table %>% group_by(Metacell) %>% dplyr::summarise(prop = max(prop))
  png(paste0(plot_path_temp, "Hist_max_proportions_per_metacell.png"), width=25, height=20, units = 'cm', res = 200)
  hist(max_prop_table$prop)
  graphics.off()
  
  ## how many metacells have > 50% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.5.png"), width=25, height=20, units = 'cm', res = 200)
  piechart_proportion_threshold(prop_table, threshold = 0.5)
  graphics.off()
  
  ## how many metacells have > 75% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.75.png"), width=25, height=20, units = 'cm', res = 200)
  piechart_proportion_threshold(prop_table, threshold = 0.75)
  graphics.off()
  
  ## how many metacells have > 90% of their cells from same label
  png(paste0(plot_path_temp, "Pie_prop_over_0.9.png"), width=25, height=20, units = 'cm', res = 200)
  piechart_proportion_threshold(prop_table, threshold = 0.9)
  graphics.off()
  
}

############################## Save ArchR with new obs #######################################

paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, "TransferLabels_Save-ArchR"), load = FALSE)
print("ArchR object saved")

############################## Write summarised counts by metacells (unchanged) #######################################

write.csv(SEACells_summarised, paste0(rds_path, "SEACells_summarised.csv"))
print("Summarised count data written out!")