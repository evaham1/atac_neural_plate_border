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
  
  freq <- df %>%
    group_by(Metacell, Category) %>% 
    dplyr::summarize(count = n()) %>%
    mutate(Metacell = as.numeric(Metacell)) %>% 
    arrange(Metacell) %>% 
    mutate(count = as.numeric(count))
  
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

SEACells_metadata <- read_csv(paste0(data_path, "rds_files/AnnData_metacells_assigned_cell_metadata.csv"))

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

############################## Read summarised counts by metacells #######################################

SEACells_summarised <- as.matrix(fread(paste0(data_path, "rds_files/", label), header = TRUE), rownames = 1)
print("Summarised count data read in!")

############################## Explore individual seacell purity #######################################

########### Stages

print("Stages purity...")

# calculate frequencies in which metacells are in each stage
stage_freq_table <- calculate_metacell_frequencies(ArchR, category = "stage")
head(stage_freq_table)

# calculate proportions of labels in metacells
stage_prop_table <- calculate_metacell_proportions(stage_freq_table)
head(stage_prop_table)

# plot the relative proportions of labels in each metacell
png(paste0(plot_path, "Hist_mixed_stages_proportions.png"), width=25, height=20, units = 'cm', res = 200)
hist(stage_prop_table$prop)
graphics.off()

########### scHelper_cell_type

print("scHelper_cell_type purity...")

# calculate frequencies in which metacells are in each stage
celltype_freq_table <- calculate_metacell_frequencies(ArchR, category = "scHelper_cell_type")
head(celltype_freq_table)

# calculate proportions of labels in metacells
celltype_prop_table <- calculate_metacell_proportions(celltype_freq_table)
head(celltype_prop_table)

# plot the relative proportions of labels in each metacell
png(paste0(plot_path, "Hist_mixed_celltype_proportions.png"), width=25, height=20, units = 'cm', res = 200)
hist(celltype_prop_table$prop)
graphics.off()

## how many metacells have > 50% of their cells from same label
high_proportion_cells <- celltype_prop_table %>% filter(prop > 0.5)
head(high_proportion_cells)
print(paste0("Number of metacells with more than 50% of their cells from same scHelper_cell_type: ", length(unique(high_proportion_cells$Metacell))))

########### scHelper_cell_type_broad

print("scHelper_cell_type_broad purity...")

# calculate frequencies in which metacells are in each stage
broad_celltype_freq_table <- calculate_metacell_frequencies(ArchR, category = "scHelper_cell_type_broad")
head(broad_celltype_freq_table)

# calculate proportions of labels in metacells
broad_celltype_prop_table <- calculate_metacell_proportions(broad_celltype_freq_table)
head(broad_celltype_prop_table)

# plot the relative proportions of labels in each metacell
png(paste0(plot_path, "Hist_mixed_broad_celltype_proportions.png"), width=25, height=20, units = 'cm', res = 200)
hist(broad_celltype_prop_table$prop)
graphics.off()

## how many metacells have > 50% of their cells from same label
high_proportion_cells <- broad_celltype_prop_table %>% filter(prop > 0.5)
head(high_proportion_cells)
print(paste0("Number of metacells with more than 50% of their cells from same broad_scHelper_cell_type: ", length(unique(high_proportion_cells$Metacell))))

########### clusters

# calculate frequencies in which metacells are in each stage
clusters_freq_table <- calculate_metacell_frequencies(ArchR, category = "clusters")
head(clusters_freq_table)

# calculate proportions of labels in metacells
clusters_prop_table <- calculate_metacell_proportions(clusters_freq_table)
head(clusters_prop_table)

# plot the relative proportions of labels in each metacell
png(paste0(plot_path, "Hist_mixed_clusters_proportions.png"), width=25, height=20, units = 'cm', res = 200)
hist(clusters_prop_table$prop)
graphics.off()

## how many metacells have > 50% of their cells from same label
high_proportion_cells <- clusters_prop_table %>%
  filter(prop > 0.5)
print(paste0("Number of metacells with more than 50% of their cells from same cluster: ", length(unique(high_proportion_cells$Metacell)))

print("SEACell purity plots made!")

####################################  Cluster metacells #######################################

## normalise each metacell by the total number of cut sites
normalised_counts <- t(apply(SEACells_summarised, 1, function(x) x/sum(x))) * 1000
normalised_counts[1:2, 1:2]
dim(normalised_counts)

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
plot(tree)

# split tree into clusters based on k
metacell_clusters <- as.data.frame(cutree(tree, k = 29))
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

############################## Plot seacell clusters on UMAP #######################################

plot_path = "plot_path/seacell_clusters/"
dir.create(plot_path, recursive = T)

p1 <- plotEmbedding(ArchR, 
                    name = "SEACell_cluster",
                    plotAs = "points", size = 1,
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    randomize = TRUE)

png(paste0(plot_path_temp, "SEACell_cluster_UMAP.png"), width=60, height=40, units = 'cm', res = 200)
print(p1)
graphics.off()

############################## Plot seacell clusters size distribution #######################################

seacell_cluster_sizes <- as.data.frame(table(getCellColData(ArchR, select = "SEACell_cluster"))) %>%
  mutate(Freq = as.numeric(Freq))

png(paste0(plot_path_temp, "SEACell_cluster_sizes_hist.png"), width=60, height=40, units = 'cm', res = 200)
hist(seacell_cluster_sizes$Freq, breaks = 30)
graphics.off()

############################## Explore seacell cluster purity #######################################

########### Stages

# calculate frequencies in which metacells are in each stage
stage_freq_table <- calculate_metacell_frequencies(ArchR, category = "stage", metacell_slot = "SEACell_cluster")
head(stage_freq_table)

# calculate proportions of labels in metacells
stage_prop_table <- calculate_metacell_proportions(stage_freq_table)
head(stage_prop_table)

# plot the relative proportions of labels in each metacell
png(paste0(plot_path, "Hist_mixed_stages_proportions.png"), width=25, height=20, units = 'cm', res = 200)
hist(stage_prop_table$prop)
graphics.off()

########### scHelper_cell_type

# calculate frequencies in which metacells are in each stage
celltype_freq_table <- calculate_metacell_frequencies(ArchR, category = "scHelper_cell_type")
head(celltype_freq_table)

# calculate proportions of labels in metacells
celltype_prop_table <- calculate_metacell_proportions(celltype_freq_table)
head(celltype_prop_table)

# plot the relative proportions of labels in each metacell
png(paste0(plot_path, "Hist_mixed_celltype_proportions.png"), width=25, height=20, units = 'cm', res = 200)
hist(celltype_prop_table$prop)
graphics.off()

## how many metacells have > 50% of their cells from same label
high_proportion_cells <- celltype_prop_table %>%
  filter(prop > 0.5)
print(paste0("Number of metacells with more than 50% of their cells from same scHelper_cell_type: ", length(unique(high_proportion_cells$Metacell)))

########### scHelper_cell_type_broad

# calculate frequencies in which metacells are in each stage
broad_celltype_freq_table <- calculate_metacell_frequencies(ArchR, category = "scHelper_cell_type_broad")
head(broad_celltype_freq_table)

# calculate proportions of labels in metacells
broad_celltype_prop_table <- calculate_metacell_proportions(broad_celltype_freq_table)
head(broad_celltype_prop_table)

# plot the relative proportions of labels in each metacell
png(paste0(plot_path, "Hist_mixed_broad_celltype_proportions.png"), width=25, height=20, units = 'cm', res = 200)
hist(broad_celltype_prop_table$prop)
graphics.off()

## how many metacells have > 50% of their cells from same label
high_proportion_cells <- broad_celltype_prop_table %>%
  filter(prop > 0.5)
print(paste0("Number of metacells with more than 50% of their cells from same broad_scHelper_cell_type: ", length(unique(high_proportion_cells$Metacell)))


########### clusters

# calculate frequencies in which metacells are in each stage
clusters_freq_table <- calculate_metacell_frequencies(ArchR, category = "clusters")
head(clusters_freq_table)

# calculate proportions of labels in metacells
clusters_prop_table <- calculate_metacell_proportions(clusters_freq_table)
head(clusters_prop_table)

# plot the relative proportions of labels in each metacell
png(paste0(plot_path, "Hist_mixed_clusters_proportions.png"), width=25, height=20, units = 'cm', res = 200)
hist(clusters_prop_table$prop)
graphics.off()

## how many metacells have > 50% of their cells from same label
high_proportion_cells <- clusters_prop_table %>%
  filter(prop > 0.5)
print(paste0("Number of metacells with more than 50% of their cells from same cluster: ", length(unique(high_proportion_cells$Metacell)))

print("SEACells clusters purity plots made!")

############################## Save ArchR with new obs #######################################

paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, "TransferLabels_Save-ArchR"), load = FALSE)
print("ArchR object saved")

############################## Write summarised counts by metacells (unchanged) #######################################

write.csv(SEACells_summarised, paste0(rds_path, "SEACells_summarised.csv"))
print("Summarised count data written out!")

