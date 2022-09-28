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

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-g", "--group_by"), action = "store", type = "character", help = "How to group cells to call peaks", default = "clusters",),
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
    
    #plot_path = "./output/NF-downstream_analysis/ArchR_peak_calling/FullData/plots/"
    #rds_path = "./output/NF-downstream_analysis/ArchR_peak_calling/FullData/rds_files/"
    rds_path = "./output/NF-downstream_analysis/ArchR_peak_calling/ss8/rds_files/"
    plot_path = "./output/NF-downstream_analysis/ArchR_peak_calling/ss8/plots/"
    
    # peaks already called
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/FILTERING/ss8/postfiltering/peak_calling/rds_files/"
    #data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/7_ArchR_clustering_postfiltering_twice/rds_files/"
    
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

############################### FUNCTIONS ####################################
cell_counts <- function(ArchR = ArchR, group1 = "clusters", group2 = "Sample") {
  group1_data <- getCellColData(ArchR, select = group1)[,1]
  group1_cell_counts <- as.data.frame(table(group1_data))
  colnames(group1_cell_counts) <- c("ID", "Total_count")
  
  group2_cell_counts <- data.frame()
  group2_data <- getCellColData(ArchR, select = group2)[,1]
  for (i in unique(group1_data)) {
    data_group1 <- getCellColData(ArchR, select = group1)[,1]
    cells <- ArchR$cellNames[BiocGenerics::which(data_group1 == i)]
    ArchR_subset <- ArchR[cells, ]
    data_group2 <- getCellColData(ArchR_subset, select = group2)[,1]
    group2_cell_counts_i <- as.data.frame(table(data_group2)) %>%
      pivot_wider(names_from = data_group2, values_from = Freq) %>% 
      add_column(ID = !!i)
    group2_cell_counts <- rbind.fill(group2_cell_counts, group2_cell_counts_i)
  }
  
  cell_counts <- merge(group1_cell_counts, group2_cell_counts)
  cell_counts[is.na(cell_counts)] <- 0
  
  ## Ordering rows and columns to better visualise
  if (group1 == "clusters"){
    cell_counts <- cell_counts %>% 
      mutate(ID = as.numeric(gsub('^.', '', ID))) %>%
      arrange(ID)
    }
  
  if (group2 == "clusters"){
    new_names <- as.numeric(gsub('^.', '', colnames(cell_counts)[3:length(colnames(cell_counts))]))
    colnames(cell_counts)[3:length(colnames(cell_counts))] <- new_names
    cell_counts <- cell_counts[, c("ID", "Total_count", 1:max(new_names))]
  }
  
  grid.arrange(tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
}


pseudoreplicate_counts <- function(ArchR = ArchR, pseudo_replicates) {
  unlisted <- unlist(pseudo_replicates, recursive=FALSE)
  
  pseudoreplicate_cell_counts_samples <- data.frame()
  for (i in c(1:length(unlisted))){
    group_name <- names(unlisted[i])
    cell_IDs <- unlisted[i][[1]]
    ArchR_pseudo_replicate <- ArchR[cell_IDs, ]
    sample_cell_counts <- as.data.frame(table(ArchR_pseudo_replicate$Sample)) %>%
      pivot_wider(names_from = Var1, values_from = Freq) %>% 
      add_column(pseudo_replicate_ID = group_name)
    pseudoreplicate_cell_counts_samples <- rbind.fill(pseudoreplicate_cell_counts_samples, sample_cell_counts)
  }

  pseudoreplicate_cell_counts_samples <- pseudoreplicate_cell_counts_samples %>%
    mutate(Cluster_ID = gsub("\\..*","", pseudo_replicate_ID)) %>%
    mutate(Cluster_ID = as.numeric(parse_number(Cluster_ID))) %>%
    arrange(Cluster_ID) %>% 
    mutate(sum_contributing_samples = rowSums(is.na(pseudoreplicate_cell_counts_samples[-which(names(pseudoreplicate_cell_counts_samples) %in% c("Cluster_ID", "pseudo_replicate_ID"))]) == FALSE))
  
  pseudoreplicate_cell_counts_samples <- pseudoreplicate_cell_counts_samples[,order(colnames(pseudoreplicate_cell_counts_samples))]
  pseudoreplicate_cell_counts_samples[is.na(pseudoreplicate_cell_counts_samples)] <- 0
  
  grid.arrange(tableGrob(pseudoreplicate_cell_counts_samples, rows=NULL, theme = ttheme_minimal()))
}

scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak', 
                              'PGC', 'BI', 'meso', 'endo')

############################## Read in ArchR project #######################################

# If files are not in rds_files subdirectory look in input dir
label <- sub('_.*', '', list.files(data_path))
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

print(paste0("Cells grouped by: ", opt$group_by))

##############################################################################################
############################## Generating pseudo-replicates ##################################

plot_path_1 <- paste0(plot_path, "pseudoreplicates/")
dir.create(plot_path_1, recursive = T)

# Plot number of cells in each group that come from each sample
png(paste0(plot_path_1, 'cell_counts_by_sample_table.png'), height = 25, width = 30, units = 'cm', res = 400)
cell_counts(ArchR = ArchR, group1 = opt$group_by, group2 = "Sample")
graphics.off()

# Make pseudo replicates and see which samples these cells come from
pseudo_replicates <- addGroupCoverages(ArchR, groupBy = opt$group_by, returnGroups = TRUE)

png(paste0(plot_path_1, 'cell_counts_by_pseudoreplicate_table.png'), height = 80, width = 30, units = 'cm', res = 400)
pseudoreplicate_counts(ArchR, pseudo_replicates)
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
colnames(counts) <- c("ID", "Number of peaks")

# order rows by cluster number or scHelper cell state
if (opt$group_by == "clusters") {
  counts <- counts %>%
    mutate(ID = substr(counts$ID, 2, nchar(as.character(counts$ID)))) %>%
    mutate(ID = as.numeric(as.character(ID))) %>%
    arrange(ID)
}
if (opt$group_by == "scHelper_cell_type_old") {
  order <- intersect(scHelper_cell_type_order, counts$ID)
  counts <- counts[match(order, counts$ID),]
}

## Plot how many peaks found per cluster
png(paste0(plot_path, 'peak_counts_per_group.png'), height = 45, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Peak Counts per group", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

########################################################################
##################### How many cut sites per peak ###########################

peak_data <- getMatrixFromProject(ArchR_peaks, useMatrix = "PeakMatrix")
peak_matrix <- t(assays(peak_data)[[1]])
colnames(peak_matrix) <- rowData(peak_data)$name

png(paste0(plot_path, 'cutsites_per_peak.png'), height = 45, width = 10, units = 'cm', res = 400)
hist(colSums(peak_matrix), breaks = 1000, main = "Histogram of cut sites per peak", 
     xlab = "Number of cut sites", ylab = "Frequency")
graphics.off()

table <- data.frame(Stats = names(cutsites), Value = as.vector(cutsites))
grid.arrange(top=textGrob("Cut sites per peak", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(table, rows=NULL, theme = ttheme_minimal()))


##############################################################################
############################# Peak Annotations ###############################

plot_path <- paste0(plot_path, "annotations/")

## What are these peaks annotated too
counts <- as.data.frame(table(peaks_df$peakType))
colnames(counts) <- c("Peak Annotation", "Number of peaks")
counts <- counts %>% mutate('Percentage' = round(100*(`Number of peaks`/sum(counts$`Number of peaks`))))

png(paste0(plot_path, 'peak_counts_per_type.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Peak Counts per type", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path, 'peak_counts_per_type_piechart.png'), height = 10, width = 20, units = 'cm', res = 400)
ggplot(counts, aes(x="", y=`Number of peaks`, fill=`Peak Annotation`)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  geom_text(aes(label = paste0(Percentage, "%")),
            position = position_stack(vjust = 0.5))
graphics.off()


## What is the distance of these peaks to nearest TSS
ggplot(peaks_df, aes(x = distToTSS, fill = peakType)) + 
  geom_histogram(binwidth=1000, alpha=0.5)

ggplot(peaks_df, aes(x = log(distToTSS), fill = peakType)) + 
  geom_histogram(binwidth=1000, alpha=0.5)

png(paste0(plot_path, 'dist_to_TSS.png'), height = 40, width = 30, units = 'cm', res = 400)
ggplot(peaks_df, aes(x = distToTSS)) + 
  geom_histogram(binwidth=1000) + 
  facet_grid(peakType ~ .) +
  geom_vline(aes(xintercept = 500), color = "black", linetype = "dashed", size = 1)
graphics.off()

peaks_df <- mutate(peaks_df, log_distToTSS = log(peaks_df$distToTSS))
png(paste0(plot_path, 'log_dist_to_TSS.png'), height = 40, width = 30, units = 'cm', res = 400)
ggplot(peaks_df, aes(x = log_distToTSS)) + 
  geom_histogram(binwidth = 0.05) + 
  facet_grid(peakType ~ .) +
  geom_vline(aes(xintercept = log(500)), color = "black", linetype = "dashed", size = 1)
graphics.off()

## What is the distance of these peaks to nearest gene
png(paste0(plot_path, 'dist_to_gene_start.png'), height = 40, width = 30, units = 'cm', res = 400)
ggplot(peaks_df, aes(x = distToGeneStart)) + 
  geom_histogram(binwidth=1000) + 
  facet_grid(peakType ~ .) +
  geom_vline(aes(xintercept = 500), color = "black", linetype = "dashed", size = 1)
graphics.off()

peaks_df <- mutate(peaks_df, log_distToGeneStart = log(peaks_df$distToGeneStart))
png(paste0(plot_path, 'log_dist_to_gene_start.png'), height = 40, width = 30, units = 'cm', res = 400)
ggplot(peaks_df, aes(x = log_distToGeneStart)) + 
  geom_histogram(binwidth = 0.05) + 
  facet_grid(peakType ~ .) +
  geom_vline(aes(xintercept = log(500)), color = "black", linetype = "dashed", size = 1)
graphics.off()

## What proportion of these peaks from each annotation are less or more than 500bp from a TSS
distal_peaks <- peaks_df[which(peaks_df$peakType == "Distal"), ]
distal_peaks <- mutate(distal_peaks, near_TSS = ifelse(distal_peaks$distToTSS < 500, "Less", "More"))
counts <- table(distal_peaks$near_TSS)
png(paste0(plot_path, 'Distal_500bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Distal peaks more or less than 500bp from TSS")
graphics.off

promoter_peaks <- peaks_df[which(peaks_df$peakType == "Promoter"), ]
promoter_peaks <- mutate(promoter_peaks, near_TSS = ifelse(promoter_peaks$distToTSS < 500, "Less", "More"))
counts <- table(promoter_peaks$near_TSS)
png(paste0(plot_path, 'Promoter_500bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Promoter peaks more or less than 500bp from TSS")
graphics.off

intronic_peaks <- peaks_df[which(peaks_df$peakType == "Intronic"), ]
intronic_peaks <- mutate(intronic_peaks, near_TSS = ifelse(intronic_peaks$distToTSS < 500, "Less", "More"))
counts <- table(intronic_peaks$near_TSS)
png(paste0(plot_path, 'Intronic_500bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Intronic peaks more or less than 500bp from TSS")
graphics.off

exonic_peaks <- peaks_df[which(peaks_df$peakType == "Exonic"), ]
exonic_peaks <- mutate(exonic_peaks, near_TSS = ifelse(exonic_peaks$distToTSS < 500, "Less", "More"))
counts <- table(exonic_peaks$near_TSS)
png(paste0(plot_path, 'Exonic_500bp_TSS_threshold.png'), height = 20, width = 25, units = 'cm', res = 400)
pie(counts, main = "Exonic peaks more or less than 500bp from TSS")
graphics.off