#!/usr/bin/env Rscript

print("late differences across all timepoints")
# calculates differential peaks between clusters, scHelper_cell_types_old and stage -> plots on heatmaps

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)
library(clustree)
library(plyr)
library(ComplexHeatmap)
library(scHelper)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-g", "--group_by"), action = "store", type = "character", help = "How to group cells to call peaks", default = "clusters",),
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
    
    # input folder - none of these work
    data_path = "./output/NF-downstream_analysis/Processing/ss8/ARCHR_INTEGRATING_WF/Single_cell_integration/rds_files/"
    data_path = "./output/NF-downstream_analysis/Processing/HH6/ARCHR_INTEGRATING_WF/Single_cell_integration/rds_files/"
    data_path = "./output/NF-downstream_analysis/Processing/HH7/ARCHR_INTEGRATING_WF/Single_cell_integration/rds_files/"
    data_path = "./output/NF-downstream_analysis/Processing/ss4/ARCHR_INTEGRATING_WF/Single_cell_integration/rds_files/"
    data_path = "./output/NF-downstream_analysis/Processing/ss4/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/"
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Peak_call/rds_files/"
    
    # works but no peaks have to recalculate interactively
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Clustering/rds_files/"
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
}

############################### FUNCTIONS ####################################

### Function to Log2normalise matrix (must be run before subsetting!)
Log2norm <- function(mat, scaleTo = 10^4) {
  mat <- log2(t(t(mat)/colSums(mat)) * scaleTo + 1) # normalising means for depth of cluster
  return(mat)
}

### Function to subset normalised matrix using IDs
subset_matrix <- function(mat, ids) {
  subsetted_matrix <- mat[ids, ]
  return(subsetted_matrix)
}

###########################################################################################
############################## Read in ArchR project #####################################

# If files are not in rds_files subdirectory look in input dir
label <- unique(sub('_.*', '', list.files(data_path)))
print(label)

if (length(label) == 0){
  data_path = "./input/"
  label <- unique(sub('_.*', '', list.files(data_path)))
  print(label)
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
} else {
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
}

getAvailableMatrices(ArchR)
ArchR@peakSet

###########################################################################################
########################## Calculate se across all cell groups ###############################

print("Calculating se across all cell groups...")
print(paste0("Cells grouped by: ", opt$group_by))

se <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = opt$group_by)
se <- scHelper::ArchRAddUniqueIdsToSe(se, ArchR, matrix_type = "PeakMatrix")

#saveRDS(se, file = paste0(rds_path, label, "_SE.RDS"))

print("se RDS saved")

seMarker <- se

###########################################################################################
############################## Read in data #####################################

# # Retrieve object label
# label <- unique(sub('_.*', '', list.files(data_path)))
# print(label)

# # load ArchR object using its retrieved name
# ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
# paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
# print('data read in')
# print(ArchR)

# # check that gene score matrix and gene integration matrix are available
# getAvailableMatrices(ArchRProj = ArchR)

# # read in se object
# seMarker <- readRDS(paste0(data_path, label, "_SE.RDS"))

# print("data read in!")

############################################################################################
############################## COLOURS #######################################

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
# set colour palettes for UMAPs
atac_scHelper_old_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_old)]

# set pal for heatmaps
pal = paletteContinuous(set = "solarExtra", n = 100)


##############################################################################
############################# Peak Heatmaps ##################################

plot_path_temp <- paste0(plot_path, "diff_peaks_heatmaps/")
dir.create(plot_path_temp, recursive = T)
  
# prepare for plotting
normalised_matrix <- ArchR_ExtractMeansFromSe(seMarker, Log2norm = TRUE, scaleTo = 10^4) # extract means df from se object and log2norm all features in each cell group
print("matrix for plotting made")
  
# Heatmap of positive markers which pass cutoff thresholds
ids <- ArchR_ExtractIds(seMarker, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE) # extract ids
if (length(ids) < 2){
  print(paste0(length(ids), " features passed cutoff - not enough to make heatmap"))
} else {
  print(paste0(length(ids), " features passed cutoff - now plotting heatmap"))
  subsetted_matrix <- normalised_matrix[ids, ]
    
  png(paste0(plot_path_temp, 'diff_cutoff_heatmap.png'), height = 40, width = 20, units = 'cm', res = 400)
  print(ArchR_PlotMarkerHeatmap(subsetted_matrix, pal = pal, fontSizeCols = 20))
  graphics.off()
}
  
# Heatmap of positive markers top 10 per cell group
ids <- ArchR_ExtractIds(seMarker, cutOff = "FDR <= 0.05 & Log2FC >= 0", top_n = TRUE, n = 10) # extract ids
print(paste0(length(ids), " features passed cutoff for top 10 heatmap"))
subsetted_matrix <- normalised_matrix[ids, ]
  
png(paste0(plot_path_temp, 'diff_top10_heatmap.png'), height = 70, width = 20, units = 'cm', res = 400)
ArchR_PlotMarkerHeatmap(subsetted_matrix, labelRows = TRUE, pal = pal, cluster_columns = FALSE, cluster_rows = FALSE, fontSizeCols = 20)
graphics.off()

print("heatmaps made")

############################# Promoter Peaks #######################################

# extract diff peaks that are promoters
markersPeaks_promoter <- subset(seMarker, rowData(seMarker)$peakType == "Promoter")
marker_tables_promoter = markersPeaks_promoter %>% getMarkers(cutOff = "FDR <= 0.01 & Log2FC >= 1")

# make df of them and write
marker_tables_promoter_tmp = marker_tables_promoter %>%  as.data.frame()
write.csv(marker_tables_promoter_tmp, paste0(plot_path, "differentially_accessible_peaks_promoter.csv"), row.names = FALSE)

# plot them in heatmap
ids <- marker_tables_promoter_tmp$unique_id
if (length(ids) < 2){
  print(paste0(length(ids), " features passed cutoff - not enough to make heatmap"))
} else {
  print(paste0(length(ids), " features passed cutoff - now plotting heatmap"))
  subsetted_matrix <- normalised_matrix[ids, ]
    
  png(paste0(plot_path_temp, 'diff_cutoff_promoter_peaks_heatmap.png'), height = 40, width = 20, units = 'cm', res = 400)
  print(ArchR_PlotMarkerHeatmap(subsetted_matrix, pal = pal, fontSizeCols = 20))
  graphics.off()
}

############################# Distal Peaks #######################################

# extract diff peaks that are promoters
markersPeaks_distal <- subset(seMarker, rowData(seMarker)$peakType == "Distal")
marker_tables_distal = markersPeaks_distal %>% getMarkers(cutOff = "FDR <= 0.01 & Log2FC >= 1")

# make df of them and write
marker_tables_distal_tmp = marker_tables_distal %>%  as.data.frame()
write.csv(marker_tables_distal_tmp, paste0(plot_path, "differentially_accessible_peaks_distal.csv"), row.names = FALSE)

# plot them in heatmap
ids <- marker_tables_distal_tmp$unique_id
if (length(ids) < 2){
  print(paste0(length(ids), " features passed cutoff - not enough to make heatmap"))
} else {
  print(paste0(length(ids), " features passed cutoff - now plotting heatmap"))
  subsetted_matrix <- normalised_matrix[ids, ]
    
  png(paste0(plot_path_temp, 'diff_cutoff_distal_peaks_heatmap.png'), height = 40, width = 20, units = 'cm', res = 400)
  print(ArchR_PlotMarkerHeatmap(subsetted_matrix, pal = pal, fontSizeCols = 20))
  graphics.off()
}


###################### Boxplots showing distribution of FDR and Logf2c values #############################

png(paste0(plot_path, 'Log2FC_boxplot.png'), height = 20, width = 20, units = 'cm', res = 400)
boxplot(assays(seMarker)$Log2FC)
graphics.off()

png(paste0(plot_path, 'FDR_boxplot.png'), height = 20, width = 20, units = 'cm', res = 400)
boxplot(assays(seMarker)$FDR)
graphics.off()


###################### Extract Log2FC and p-vals of significant peaks from each cell group #############################

group_df <- getCellColData(ArchR, select = opt$group_by)
groups <- unique(group_df[,1])

print("Looping over cell groups:")
print(groups)

final_df <- data.frame()

for (group in groups){
  
  print(group)
  
  # extract peaks from this group using pseudoreplicate names
  length(grep(group, rowData(seMarker)$GroupReplicate))
  length(rowData(seMarker)$GroupReplicate)
  
  subset_se <- seMarker[grep(group, rowData(seMarker)$GroupReplicate), ]
  
  # extract logFC and p-vals for these peaks and put in df
  df <- data.frame(LogFC = c(t(assays(subset_se)$Log2FC)), FDR = c(t(assays(subset_se)$FDR)),
                   LogFDR = log10(c(t(assays(subset_se)$FDR))), stringsAsFactors=FALSE)
  df <- df %>% mutate(Passed = as.factor(ifelse(FDR < 0.01 & LogFC > 1 | FDR < 0.01 & LogFC < -1,"passed", "failed")))
  # add group name to this df
  df <- df %>% mutate(Cell_group = group)
  
  # combine to final df
  final_df <- rbind(final_df, df)
}

###################### Volcano plots #############################

set.seed(42)
rows <- sample(nrow(final_df))
final_df <- final_df[rows, ]

png(paste0(plot_path, 'FDR_Log2FC_scatterplot_sig_colour.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(final_df, aes(x = -LogFDR, y = LogFC, color = Passed, shape = Passed)) +
  geom_point() +
  scale_color_manual(values=c("black", "4f7942")) +
  scale_shape_manual(values=c(16, 17)) +
  xlim(0, 50) +
  ylim(-15, 15) +
  theme_minimal(text = element_text(size = 20))
graphics.off()

png(paste0(plot_path, 'FDR_Log2FC_scatterplot_cell_group_col.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(final_df, aes(x = -LogFDR, y = LogFC, color = Cell_group, shape = Passed)) +
  geom_point() +
  # scale_color_manual(values=c("black", "4f7942")) +
  scale_shape_manual(values=c(16, 17)) +
  xlim(0, 50) +
  ylim(-15, 15) +
  theme_minimal(text = element_text(size = 20))
graphics.off()

for (group in groups){
  print(group)
  
  # add temp column to df to only highlight these peaks
  group_df <- final_df %>% mutate(Highlight = as.factor(ifelse(Cell_group == group,"passed", "failed")))
  
  # plot with highlight on this cell group
  png(paste0(plot_path, 'FDR_Log2FC_scatterplot_cell_group_col_', group, '.png'), height = 23, width = 20, units = 'cm', res = 400)
  print(
    ggplot(group_df, aes(x = -LogFDR, y = LogFC, color = Highlight, shape = Passed)) +
    geom_point() +
    scale_color_manual(values=c("grey", "red")) +
    scale_shape_manual(values=c(16, 17)) +
    xlim(0, 50) +
    ylim(-15, 15) +
    theme_minimal(text = element_text(size = 20))
  )
  graphics.off()
}

###################### Barcharts of nPeaks per cluster and how many are significant #############################

freq_table <- table(final_df$Cell_group, final_df$Passed)

# Convert the data frame to a tidy format
freq_table <- as.data.frame(freq_table)
colnames(freq_table) <- c("Group_name", "Status", "Freq")

# Stacked
png(paste0(plot_path, 'Freq_peaks_passed_by_cell_group_stacked_boxplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(freq_table, aes(fill=Status, y=Freq, x=Group_name)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("black", "4f7942")) +
  ylim(0, 100000) +
  theme_minimal()
graphics.off()