#!/usr/bin/env Rscript

print("transfer peak set and labels from one ArchR object onto another ArchR object")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
# library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)
library(ComplexHeatmap)
# library(BSgenome.Ggallus.UCSC.galGal6)
library(scHelper)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-l", "--labels"), action = "store", type = "character", help = "Which labels to transfer", default = "clusters,scHelper_cell_type"),
  make_option(c("-g", "--group_by"), action = "store", type = "character", help = "How to group cells to call peaks", default = "clusters"),
  make_option(c("-v", "--verbose"), action = "store", type = "logical", help = "Verbose", default = FALSE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    
    # donor object - full data that has undergone peak calling and integration
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Single_cell_integration/rds_files/"
    
    # recipient object - stage data
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Clustering/rds_files/"
    
    # all the labels I want to transfer:
    opt$labels <- "clusters,scHelper_cell_type,scHelper_cell_type_broad,predictedCell_Un,predictedScore_Un"
    
    addArchRThreads(threads = 1)
    
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

############################## Read in ArchR project #######################################

# Extract stage name
label <- unique(sub('_.*', '', list.files(data_path)))
print(paste0("Label: ", label))

# Read in target ArchR object (in rds_files because it came from previous process)
print("reading in target ArchR object in rds_files folder")
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getPeakSet(ArchR)
getAvailableMatrices(ArchR)

# Read in donor ArchR object (in input folder)
print("reading in donor ArchR object in input folder")
ArchR_donor <- loadArchRProject(path = "./input/FullData_Save-ArchR", force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# see what is in the ArchR donor object already
print("ArchR donor object info: ")
print(ArchR_donor)
getPeakSet(ArchR_donor)
getAvailableMatrices(ArchR_donor)

# how cells are grouped for pseudobulk replicates + peak calling
print(paste0("Cells grouped by: ", opt$group_by))

##########################################################################################
############################## FUNCTIONS ###############################

ArchR_TransferLabels <- function(donor, recipient, transfer_label){
  cell_ids <- donor$cellNames
  transfer_ids <- getCellColData(donor, select = transfer_label)
  transfer_ids <- as.character(transfer_ids[,1])
  recipient <- addCellColData(
    recipient, force = TRUE,
    data = transfer_ids,
    cells = cell_ids, 
    name = paste0("transferred_", transfer_label)
  )
  return(recipient)
}

##########################################################################################
############################## Transfer integrated labels ###############################

print("Transferring labels...")

## first subset donor object to only include cells in recipient object
idxSample <- BiocGenerics::which(ArchR_donor$stage %in% label)
cellsSample <- ArchR_donor$cellNames[idxSample]
ArchR_donor_subset <- ArchR_donor[cellsSample, ]

## then transfer each label in turn
print(opt$labels)
labels <- unlist(strsplit(opt$labels, ",", fixed = TRUE))

for (label in labels){
  print(paste0("Transferring: ", label))
  
  ArchR <- ArchR_TransferLabels(donor = ArchR_donor_subset,
                                recipient = ArchR, 
                                transfer_label = label)
}

##########################################################################
############################## Plot labels ###############################

print("Label plots...")

###### schelper cell type colours
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                "#A5718D", "#3F918C", "#ed5e5f", "#9792A3")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                       'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo',
                                       'Neural', 'Placodal', 'Non-neural', 'Contam')
cols <- scHelper_cell_type_colours[unique(ArchR$transferred_scHelper_cell_type)]
cols_broad <- scHelper_cell_type_colours[unique(ArchR$transferred_scHelper_cell_type_broad)]

### Plot transferred clusters
png(paste0(plot_path, "transferred_clusterse.png"), width=30, height=40, units = 'cm', res = 200)
plotEmbedding(ArchR,
              name = "transferred_clusters",
              plotAs = "points", size = 1.8,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

### Plot scHelper_cell_type
png(paste0(plot_path, "transferred_scHelper_cell_type_cell_type.png"), width=30, height=40, units = 'cm', res = 200)
plotEmbedding(ArchR,
              name = "transferred_scHelper_cell_type", pal = cols,
              plotAs = "points", size = 1.8,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

### Plot scHelper_cell_type broad
png(paste0(plot_path, "transferred_scHelper_cell_type_cell_type_broad.png"), width=30, height=40, units = 'cm', res = 200)
plotEmbedding(ArchR,
              name = "transferred_scHelper_cell_type_broad", pal = cols_broad,
              plotAs = "points", size = 1.8,
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)
graphics.off()

### Plot integration scores
png(paste0(plot_path, 'Integration_Scores_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "transferred_predictedScore_Un", plotAs = "points", size = 1.8, baseSize = 0, 
              legendSize = 10)
graphics.off()

##################### Distribution of labels across clusters ##################################

plot_path = "./plots/label_by_cluster_distribution/"
dir.create(plot_path, recursive = T)

# visualise distribution across clusters: table of cell counts
png(paste0(plot_path, 'label_by_cluster_cell_number_table.png'), height = 25, width = 40, units = 'cm', res = 400)
scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "transferred_scHelper_cell_type", group2 = "clusters", print_table = TRUE, scHelper_cell_type_order = scHelper_cell_type_order)
graphics.off()

# visualise distribution across clusters: confusion matrix
png(paste0(plot_path, "label_by_cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
scHelper::ArchRCellCountsHeatmap(ArchR = ArchR, group1 = "transferred_scHelper_cell_type", group2 = "clusters")
graphics.off()

# visualise distribution across clusters: table of cell percentages
cell_counts <- scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "transferred_scHelper_cell_type", group2 = "clusters", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
percentage_counts <- as.data.frame(lapply(cell_counts, function(x) (x / sum(x))*100))
rownames(percentage_counts) <- rownames(cell_counts)

png(paste0(plot_path, 'label_by_cluster_cell_percentage_table.png'), height = 25, width = 40, units = 'cm', res = 400)
grid.arrange(tableGrob(round(percentage_counts, 2), theme = ttheme_minimal()))
graphics.off()

# visualise distribution across clusters: piecharts
counts <- scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
png(paste0(plot_path, "label_by_cluster_piecharts.png"), width=50, height=40, units = 'cm', res = 200)
scHelper::CellLabelPieCharts(counts, col = scHelper_cell_type_colours)
graphics.off()


################################################################################
############################## Transfer peak set ###############################

print("Transferring peak set...")

# extract peakset granges from donor ArchR object
peakset_granges <- getPeakSet(ArchR_donor)

# add the peak set to the target ArchR object and calculate the matrix
ArchR <- addPeakSet(ArchR, peakset_granges, force = TRUE)
ArchR <- addPeakMatrix(ArchR, force = TRUE)

#################################################################################
############################## Save ArchR project ###############################

print("Saving data...")

output_directory <- paste0(rds_path, label, "_Save-ArchR")
print(paste0("Output directory of object with labels and peak set transferred onto it: ", output_directory))
saveArchRProject(ArchRProj = ArchR, outputDirectory = output_directory, load = FALSE)

########################################################################
##################### How many cut sites per peak ######################

plot_path_temp = "./plots/cut_sites_per_peak/"
dir.create(plot_path, recursive = T)

peak_data <- getMatrixFromProject(ArchR_peaks, useMatrix = "PeakMatrix")
peak_matrix <- t(assays(peak_data)[[1]])
colnames(peak_matrix) <- rowData(peak_data)$name

png(paste0(plot_path_temp, 'cutsites_per_peak_histogram.png'), height = 20, width = 40, units = 'cm', res = 400)
hist(colSums(peak_matrix), breaks = 1000, main = "Histogram of cut sites per peak", 
     xlab = "Number of cut sites", ylab = "Frequency")
graphics.off()

png(paste0(plot_path_temp, 'cutsites_per_peak_boxplot.png'), height = 30, width = 15, units = 'cm', res = 400)
boxplot(colSums(peak_matrix), breaks = 1000, main = "Boxplot of cut sites per peak", 
        xlab = "Number of cut sites", ylab = "Frequency")
graphics.off()

cutsites <- summary(colSums(peak_matrix))
table <- data.frame(Stats = names(cutsites), Value = as.vector(cutsites))
png(paste0(plot_path_temp, 'cutsites_per_peak_table.png'), height = 30, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Cut sites per peak", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(table, rows=NULL, theme = ttheme_minimal()))
graphics.off()

print("cut sites per peak calculated")