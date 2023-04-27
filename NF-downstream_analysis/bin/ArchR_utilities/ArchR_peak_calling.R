#!/usr/bin/env Rscript

print("peak_calling_ArchR ")
# creates pseudorepicates and calls peaks on user-defined groups of cells

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)
library(clustree)
library(ComplexHeatmap)
library(BSgenome.Ggallus.UCSC.galGal6)
library(scHelper)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-g", "--group_by"), action = "store", type = "character", help = "How to group cells to call peaks", default = "clusters",),
  make_option(c("", "--heatmaps_stage"), action = "store", type = "logical", help = "Whether to plot heatmap of diff peaks for stage data", default = FALSE,),
  make_option(c("", "--heatmaps_full"), action = "store", type = "logical", help = "Whether to plot heatmap of diff peaks for full data", default = FALSE,),
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
    
    # peaks already called on ss8
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Peak_call/rds_files/"
    
    # peaks already called on FullData
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Peak_call/rds_files/"
    
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

# If files are not in rds_files subdirectory look in input dir
label <- sub('_.*', '', list.files(data_path))
print(paste0("Label: ", label))

if (length(label) == 0){
  print("ArchR object not in rds_files folder, checking input folder...")
  data_path = "./input/"
  label <- sub('_.*', '', list.files(data_path))
  print(paste0("Label: ", label))
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
} else {
  print("ArchR object in rds_files folder")
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
}

print(paste0("Cells grouped by: ", opt$group_by))

##############################################################################################
############################## Generating pseudo-replicates ##################################

plot_path_temp <- paste0(plot_path, "pseudoreplicates/")
dir.create(plot_path_temp, recursive = T)

# If there is more than one sample, plot number of cells in each group that come from each sample
if (length(unique(ArchR$Sample)) > 1){
  png(paste0(plot_path_temp, 'cell_counts_by_sample_table.png'), height = 25, width = 30, units = 'cm', res = 400)
  ArchR_cell_counting(ArchR = ArchR, group1 = opt$group_by, group2 = "Sample")
  graphics.off()
}

# Make pseudo replicates
pseudo_replicates <- addGroupCoverages(ArchR, groupBy = opt$group_by, returnGroups = TRUE, force = TRUE)

# Plot table to see which samples and groups the pseudo replicate cells come from 
png(paste0(plot_path_temp, 'pseudoreplicate_cell_counts_per_sample_table.png'), height = 40, width = 30, units = 'cm', res = 400)
ArchR_pseudoreplicate_counts(ArchR, pseudo_replicates, group_by = "Sample")
graphics.off()

png(paste0(plot_path_temp, 'pseudoreplicate_cell_counts_per_group_table.png'), height = 40, width = 70, units = 'cm', res = 400)
ArchR_pseudoreplicate_counts(ArchR, pseudo_replicates, group_by = opt$group_by)
graphics.off()

#####  Make actual pseudo-replicates for peak calling:
ArchR <- addGroupCoverages(ArchR, groupBy = opt$group_by, returnGroups = FALSE, force = TRUE)
print("pseudo replicates created")

##############################################################################################
############################## Call peaks on pseudo-replicates ###############################

ArchR <- addReproduciblePeakSet(
  ArchRProj = ArchR,
  groupBy = opt$group_by,
  pathToMacs2 = "/opt/conda/bin/macs2",
  force = TRUE,
  genomeSize = 1230258557, # copied from Grace's paper, need to check this
)
print("peaks called using Macs2")

getPeakSet(ArchR)
getAvailableMatrices(ArchR)

## add unique names to peaks
ps_df <- data.frame(ArchR@peakSet)
ps <- ArchR@peakSet
ps$name <- paste0(ps_df$seqnames, "-", ps_df$start, "-", ps_df$end)
ArchR@peakSet <- ps
ArchR_peaks <- addPeakMatrix(ArchR, force = TRUE)

getPeakSet(ArchR_peaks)
getAvailableMatrices(ArchR_peaks)

#################################################################################
############################## Save ArchR project ###############################
saveArchRProject(ArchRProj = ArchR_peaks, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)
print("ArchR project saved")

##############################################################################
##################### How many peaks per cell group ###########################

## Create df of peaks
peaks_granges <- getPeakSet(ArchR_peaks)
ids <- names(peaks_granges@ranges)
names(peaks_granges@ranges) <- c(1:length(peaks_granges))

peaks_df <- as.data.frame(peaks_granges)
peaks_df <- peaks_df %>% mutate(ID = ids)

counts <- as.data.frame(table(peaks_df$ID))
colnames(counts) <- c("ID", "nPeaks")

# order rows by cluster number or scHelper cell state or stage_clusters???
if (opt$group_by == "clusters") {
  counts <- counts %>%
    mutate(ID = substr(counts$ID, 2, nchar(as.character(counts$ID)))) %>%
    mutate(ID = as.numeric(as.character(ID))) %>%
    arrange(ID) %>%
    mutate(ID = as.character(ID))
}
if (opt$group_by == "scHelper_cell_type_old") {
  order <- intersect(scHelper_cell_type_order, counts$ID)
  counts <- counts[match(order, counts$ID),]
}

# Add total peak counts
counts <- counts %>% tibble::add_row(ID = "Total", nPeaks = sum(counts$nPeaks))
print(counts)

## Plot how many peaks found per cluster
png(paste0(plot_path, 'peak_counts_per_group.png'), height = 45, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Peak Counts per group", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

print("peaks per group calculated")

########################################################################
##################### How many cut sites per peak ######################

plot_path_temp <- paste0(plot_path, "cut_sites_per_peak/")
dir.create(plot_path_temp, recursive = T)

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


##############################################################################
############################# Peak Annotations ###############################

plot_path_temp <- paste0(plot_path, "annotations/")
dir.create(plot_path_temp, recursive = T)

## What are these peaks annotated too
counts <- as.data.frame(table(peaks_df$peakType))
colnames(counts) <- c("Peak Annotation", "Number of peaks")
counts <- counts %>% mutate('Percentage' = round(100*(`Number of peaks`/sum(counts$`Number of peaks`))))

png(paste0(plot_path_temp, 'peak_counts_per_type.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Peak Counts per type", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path_temp, 'peak_counts_per_type_piechart.png'), height = 10, width = 20, units = 'cm', res = 400)
ggplot(counts, aes(x="", y=`Number of peaks`, fill=`Peak Annotation`)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  geom_text(aes(label = paste0(Percentage, "%")),
            position = position_stack(vjust = 0.5))
graphics.off()


## What is the distance of these peaks to nearest TSS
plot_path_temp <- paste0(plot_path, "annotations/distance_to_nearest_TSS/")
dir.create(plot_path_temp, recursive = T)

# ggplot(peaks_df, aes(x = distToTSS, fill = peakType)) + 
#   geom_histogram(binwidth=1000, alpha=0.5)

png(paste0(plot_path_temp, 'dist_to_TSS.png'), height = 40, width = 30, units = 'cm', res = 400)
ggplot(peaks_df, aes(x = distToTSS)) + 
  geom_histogram(binwidth=1000) + 
  facet_grid(peakType ~ .) +
  geom_vline(aes(xintercept = 500), color = "black", linetype = "dashed", size = 1)
graphics.off()

peaks_df <- mutate(peaks_df, log_distToTSS = log(peaks_df$distToTSS))
png(paste0(plot_path_temp, 'log_dist_to_TSS.png'), height = 40, width = 30, units = 'cm', res = 400)
ggplot(peaks_df, aes(x = log_distToTSS)) + 
  geom_histogram(binwidth = 0.05) + 
  facet_grid(peakType ~ .) +
  geom_vline(aes(xintercept = log(500)), color = "black", linetype = "dashed", size = 1)
graphics.off()

## What is the distance of these peaks to nearest gene
plot_path_temp <- paste0(plot_path, "annotations/distance_to_nearest_gene/")
dir.create(plot_path_temp, recursive = T)

png(paste0(plot_path_temp, 'dist_to_gene_start.png'), height = 40, width = 30, units = 'cm', res = 400)
ggplot(peaks_df, aes(x = distToGeneStart)) + 
  geom_histogram(binwidth=1000) + 
  facet_grid(peakType ~ .) +
  geom_vline(aes(xintercept = 500), color = "black", linetype = "dashed", size = 1)
graphics.off()

peaks_df <- mutate(peaks_df, log_distToGeneStart = log(peaks_df$distToGeneStart))
png(paste0(plot_path_temp, 'log_dist_to_gene_start.png'), height = 40, width = 30, units = 'cm', res = 400)
ggplot(peaks_df, aes(x = log_distToGeneStart)) + 
  geom_histogram(binwidth = 0.05) + 
  facet_grid(peakType ~ .) +
  geom_vline(aes(xintercept = log(500)), color = "black", linetype = "dashed", size = 1)
graphics.off()

## What proportion of these peaks from each annotation are less or more than 500bp from a TSS
plot_path_temp <- paste0(plot_path, "annotations/how_many_peaks_within_ranges_of_nearest_TSS/")
dir.create(plot_path_temp, recursive = T)

#Distal
distal_peaks <- peaks_df[which(peaks_df$peakType == "Distal"), ]
distal_peaks <- mutate(distal_peaks, near_TSS = ifelse(distal_peaks$distToTSS < 500, "Less", "More"))
counts <- table(distal_peaks$near_TSS)
png(paste0(plot_path_temp, 'Distal_500bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Distal peaks more or less than 500bp from TSS")
graphics.off

distal_peaks <- mutate(distal_peaks, near_TSS = ifelse(distal_peaks$distToTSS < 1000, "Less", "More"))
counts <- table(distal_peaks$near_TSS)
png(paste0(plot_path_temp, 'Distal_1000bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Distal peaks more or less than 1000bp from TSS")
graphics.off

#Promoter
promoter_peaks <- peaks_df[which(peaks_df$peakType == "Promoter"), ]
promoter_peaks <- mutate(promoter_peaks, near_TSS = ifelse(promoter_peaks$distToTSS < 500, "Less", "More"))
counts <- table(promoter_peaks$near_TSS)
png(paste0(plot_path_temp, 'Promoter_500bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Promoter peaks more or less than 500bp from TSS")
graphics.off

promoter_peaks <- mutate(promoter_peaks, near_TSS = ifelse(promoter_peaks$distToTSS < 1000, "Less", "More"))
counts <- table(promoter_peaks$near_TSS)
png(paste0(plot_path_temp, 'Promoter_1000bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Promoter peaks more or less than 1000bp from TSS")
graphics.off

#Intronic
intronic_peaks <- peaks_df[which(peaks_df$peakType == "Intronic"), ]
intronic_peaks <- mutate(intronic_peaks, near_TSS = ifelse(intronic_peaks$distToTSS < 500, "Less", "More"))
counts <- table(intronic_peaks$near_TSS)
png(paste0(plot_path_temp, 'Intronic_500bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Intronic peaks more or less than 500bp from TSS")
graphics.off

intronic_peaks <- mutate(intronic_peaks, near_TSS = ifelse(intronic_peaks$distToTSS < 1000, "Less", "More"))
counts <- table(intronic_peaks$near_TSS)
png(paste0(plot_path_temp, 'Intronic_1000bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Intronic peaks more or less than 1000bp from TSS")
graphics.off

#Exonic
exonic_peaks <- peaks_df[which(peaks_df$peakType == "Exonic"), ]
exonic_peaks <- mutate(exonic_peaks, near_TSS = ifelse(exonic_peaks$distToTSS < 500, "Less", "More"))
counts <- table(exonic_peaks$near_TSS)
png(paste0(plot_path_temp, 'Exonic_500bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Exonic peaks more or less than 500bp from TSS")
graphics.off

exonic_peaks <- mutate(exonic_peaks, near_TSS = ifelse(exonic_peaks$distToTSS < 1000, "Less", "More"))
counts <- table(exonic_peaks$near_TSS)
png(paste0(plot_path_temp, 'Exonic_1000bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Exonic peaks more or less than 1000bp from TSS")
graphics.off

print("annotation plots made")

##############################################################################
############################# Peak Heatmaps ###############################

run_heatmaps <- ifelse(length(unique(ArchR_peaks$stage)) == 1 & isTRUE(opt$heatmaps_stage) | length(unique(ArchR_peaks$stage)) > 1 & isTRUE(opt$heatmaps_full),
                       TRUE, FALSE)

if (isTRUE(run_heatmaps)) {

  print("Running peak heatmaps...")

  plot_path_temp <- paste0(plot_path, "diff_peaks_heatmaps/")
  dir.create(plot_path_temp, recursive = T)
  
  seMarker <- getMarkerFeatures(
    ArchRProj = ArchR_peaks, 
    useMatrix = "PeakMatrix", 
    groupBy = opt$group_by)
  seMarker <- ArchR_add_unique_ids_to_se(seMarker, ArchR, matrix_type = "PeakMatrix")
  print("seMarker object made")
  
  # prepare for plotting
  normalised_matrix <- ArchR_extract_means_from_se(seMarker, Log2norm = TRUE, scaleTo = 10^4) # extract means df from se object and log2norm all features in each cell group
  print("matrix for plotting made")
  
  # heatmap palette 
  pal <- paletteContinuous(set = "solarExtra", n = 100)
  
  # Heatmap of positive markers which pass cutoff thresholds
  ids <- ArchR_extract_ids(seMarker, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE) # extract ids
  if (length(ids) < 2){
    print(paste0(length(ids), " features passed cutoff - not enough to make heatmap"))
  } else {
    print(paste0(length(ids), " features passed cutoff - now plotting heatmap"))
    subsetted_matrix <- normalised_matrix[ids, ]
    
    png(paste0(plot_path_temp, 'diff_cutoff_heatmap.png'), height = 40, width = 20, units = 'cm', res = 400)
    print(ArchR_marker_heatmap(subsetted_matrix, pal = pal))
    graphics.off()
  }
  
  # Heatmap of positive markers top 10 per cell group
  ids <- ArchR_extract_ids(seMarker, cutOff = "FDR <= 0.05 & Log2FC >= 0", top_n = TRUE, n = 10) # extract ids
  print(paste0(length(ids), " features passed cutoff for top 10 heatmap"))
  subsetted_matrix <- normalised_matrix[ids, ]
  
  png(paste0(plot_path_temp, 'diff_top10_heatmap.png'), height = 70, width = 20, units = 'cm', res = 400)
  print(ArchR_marker_heatmap(subsetted_matrix, labelRows = TRUE, pal = pal, cluster_columns = FALSE, cluster_rows = FALSE))
  graphics.off()
  
}

print("heatmaps made")