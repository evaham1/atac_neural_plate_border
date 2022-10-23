#!/usr/bin/env Rscript

print("compare stages se objects")
# compares stages using their se objects

#### TO DO: need to sort out barcharts to just plot nPeaks on Y axis and colours peaks based on if they have a FDR<0.01 for any of the differential tests
### can also keep the volcano plots but clean them up and plot them together so have same y axis

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

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-m", "--matrix"), action = "store", type = "character", help = "Matrix to use", default = "PeakMatrix",),
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
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
}

####################### Read in RDS objects ##########################

HH5_se <- readRDS("./input/HH5_SE.RDS")
HH6_se <- readRDS("./input/HH6_SE.RDS")
HH7_se <- readRDS("./input/HH7_SE.RDS")
ss4_se <- readRDS("./input/ss4_SE.RDS")
ss8_se <- readRDS("./input/ss8_SE.RDS")

dim(rowData(HH5_se)) # 92,628 peaks
dim(rowData(HH6_se)) # 80,775 peaks
dim(rowData(HH7_se)) # 119,929 peaks
dim(rowData(ss4_se)) # 116,420 peaks
dim(rowData(ss8_se)) # 146,059 peaks

#### TO DO: need to figure out why there aren't all the peaks in these se objects
## Expected total peak counts:
# HH5: 117,799
# HH6: 100,322
# HH7: 131,613
# ss4: 147,124
# ss8: 176,488

########################################################################################################
########################################################################################################
###### Barchart of total number of PEAKS and PEAKS that have a FDR < 0.01 for at least one cluster

dir.create(paste0(plot_path, "Barcharts/"), recursive = T)

HH5_FDRs <- as.data.frame(assays(HH5_se)$FDR)
HH5 <- table(apply(HH5_FDRs, 1, function(x) sum(x < 0.01) >= 1))
HH6_FDRs <- as.data.frame(assays(HH6_se)$FDR)
HH6 <- table(apply(HH6_FDRs, 1, function(x) sum(x < 0.01) >= 1))
HH7_FDRs <- as.data.frame(assays(HH7_se)$FDR)
HH7 <- table(apply(HH7_FDRs, 1, function(x) sum(x < 0.01) >= 1))
ss4_FDRs <- as.data.frame(assays(ss4_se)$FDR)
ss4 <- table(apply(ss4_FDRs, 1, function(x) sum(x < 0.01) >= 1))
ss8_FDRs <- as.data.frame(assays(ss8_se)$FDR)
ss8 <- table(apply(ss8_FDRs, 1, function(x) sum(x < 0.01) >= 1))

all_peak_counts <- rbind(HH5, HH6, HH7, ss4, ss8)
plot_data <- rownames_to_column(as.data.frame(all_peak_counts), var = "stage")
colnames(plot_data) <- c("Stage", "Failed", "Passed")
plot_data <- plot_data %>% gather(FDR_less_than_0.01, Count, Failed:Passed)

png(paste0(plot_path, 'Barcharts/All_peaks_pass_or_fail.png'), height = 26, width = 30, units = 'cm', res = 400)
ggplot(plot_data, aes(x = Stage, y = Count, fill = FDR_less_than_0.01)) +
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label=Count), vjust=1.6, color="white", size=3.5)+
  theme_minimal()
graphics.off()


########################################################################################################
########################################################################################################
###### Barchart of total number of TESTS and TESTS that have a FDR < 0.01 for at least one cluster

HH5_df <- data.frame(LogFC = c(t(assays(HH5_se)$Log2FC)), FDR = c(t(assays(HH5_se)$FDR)), stage = "HH5", stringsAsFactors=FALSE) %>% mutate(Passed = ifelse(FDR < 0.01, "Passed", "Failed"))
HH5 <- table(HH5_df$Passed)
HH6_df <- data.frame(LogFC = c(t(assays(HH6_se)$Log2FC)), FDR = c(t(assays(HH6_se)$FDR)), stage = "HH6", stringsAsFactors=FALSE) %>% mutate(Passed = ifelse(FDR < 0.01, "Passed", "Failed"))
HH6 <- table(HH6_df$Passed)
HH7_df <- data.frame(LogFC = c(t(assays(HH7_se)$Log2FC)), FDR = c(t(assays(HH7_se)$FDR)), stage = "HH7", stringsAsFactors=FALSE) %>% mutate(Passed = ifelse(FDR < 0.01, "Passed", "Failed"))
HH7 <- table(HH7_df$Passed)
ss4_df <- data.frame(LogFC = c(t(assays(ss4_se)$Log2FC)), FDR = c(t(assays(ss4_se)$FDR)), stage = "ss4", stringsAsFactors=FALSE) %>% mutate(Passed = ifelse(FDR < 0.01, "Passed", "Failed"))
ss4 <- table(ss4_df$Passed)
ss8_df <- data.frame(LogFC = c(t(assays(ss8_se)$Log2FC)), FDR = c(t(assays(ss8_se)$FDR)), stage = "ss8", stringsAsFactors=FALSE) %>% mutate(Passed = ifelse(FDR < 0.01, "Passed", "Failed"))
ss8 <- table(ss8_df$Passed)

all_counts <- rbind(HH5, HH6, HH7, ss4, ss8)
plot_data <- rownames_to_column(as.data.frame(all_counts), var = "stage")
plot_data <- plot_data %>% gather(FDR_less_than_0.01, Count, Failed:Passed)

png(paste0(plot_path, 'Barcharts/All_diff_accessibility_tests_pass_or_fail.png'), height = 26, width = 30, units = 'cm', res = 400)
ggplot(plot_data, aes(x = stage, y = Count, fill = FDR_less_than_0.01)) +
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label=Count), vjust=1.6, color="white", size=3.5)+
  theme_minimal()
graphics.off()


########################################################################################################
########################################################################################################
###### Volcano plots of TESTS LogFC and FDRs and TESTS that have a FDR < 0.01

dir.create(paste0(plot_path, "Volcano_plots/"), recursive = T)

#### adjust these plots so colour points based on pass/fail threshold
### also maybe plot them together so they have the same y axes??

png(paste0(plot_path, 'Volcano_plots/HH5_FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(HH5_df, aes(x = -log(FDR), y = LogFC, colour = Passed)) + 
  geom_point(alpha = 0.5) +
  theme_minimal()
graphics.off()

png(paste0(plot_path, 'Volcano_plots/HH6_FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(HH6_df, aes(x = -log(FDR), y = LogFC, colour = Passed)) + 
  geom_point(alpha = 0.5) + 
  theme_minimal()
graphics.off()

png(paste0(plot_path, 'Volcano_plots/HH7_FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(HH7_df, aes(x = -log(FDR), y = LogFC, colour = Passed)) + 
  geom_point(alpha = 0.5) + 
  theme_minimal()
graphics.off()

png(paste0(plot_path, 'Volcano_plots/ss4_FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(ss4_df, aes(x = -log(FDR), y = LogFC, colour = Passed)) + 
  geom_point(alpha = 0.5) + 
  theme_minimal()
graphics.off()

png(paste0(plot_path, 'Volcano_plots/ss8_FDR_Log2FC_scatterplot.png'), height = 23, width = 20, units = 'cm', res = 400)
ggplot(ss8_df, aes(x = -log(FDR), y = LogFC, colour = Passed)) + 
  geom_point(alpha = 0.5) + 
  theme_minimal()
graphics.off()

## all plots at once
all_data <- rbind(HH5_df, HH6_df, HH7_df, ss4_df, ss8_df)

png(paste0(plot_path, 'Volcano_plots/all_FDR_Log2FC_scatterplot.png'), height = 25, width = 100, units = 'cm', res = 400)
ggplot(all_data, aes(x = -log(FDR), y = LogFC, colour = Passed)) + 
  geom_point(alpha = 0.5) + 
  theme_minimal() +
  facet_grid(. ~ stage) + 
  theme(text = element_text(size = 20)) 
graphics.off()

# ###################### Plots showing how many features pass different thresholds #############################
# 
# sequence <- seq(from = 1, to = 0.01, by = -0.01)
# 
# ############. LogFC >= 0
# number_of_features <- c()
# for (i in sequence){
#   cutOff <- paste0("FDR <= ", i, " & Log2FC >= 0")
#   ids <- extract_ids(ss8_se, cutOff = cutOff, top_n = FALSE) # extract ids
#   print(length(ids))
#   number_of_features <- c(number_of_features, length(ids))
# }
# ss8_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "ss8")
# 
# number_of_features <- c()
# for (i in sequence){
#   cutOff <- paste0("FDR <= ", i, " & Log2FC >= 0")
#   ids <- extract_ids(ss4_se, cutOff = cutOff, top_n = FALSE) # extract ids
#   print(length(ids))
#   number_of_features <- c(number_of_features, length(ids))
# }
# ss4_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "ss4")
# 
# number_of_features <- c()
# for (i in sequence){
#   cutOff <- paste0("FDR <= ", i, " & Log2FC >= 0")
#   ids <- extract_ids(HH7_se, cutOff = cutOff, top_n = FALSE) # extract ids
#   print(length(ids))
#   number_of_features <- c(number_of_features, length(ids))
# }
# HH7_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH7")
# 
# number_of_features <- c()
# for (i in sequence){
#   cutOff <- paste0("FDR <= ", i, " & Log2FC >= 0")
#   ids <- extract_ids(HH6_se, cutOff = cutOff, top_n = FALSE) # extract ids
#   print(length(ids))
#   number_of_features <- c(number_of_features, length(ids))
# }
# HH6_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH6")
# 
# number_of_features <- c()
# for (i in sequence){
#   cutOff <- paste0("FDR <= ", i, " & Log2FC >= 0")
#   ids <- extract_ids(HH5_se, cutOff = cutOff, top_n = FALSE) # extract ids
#   print(length(ids))
#   number_of_features <- c(number_of_features, length(ids))
# }
# HH5_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH5")
# 
# df <- do.call("rbind", list(HH5_df, HH6_df, HH7_df, ss4_df, ss8_df))
# 
# png(paste0(plot_path, 'changing_FDR_cutoffs_logFC_0_linegraph.png'), height = 20, width = 20, units = 'cm', res = 400)
# ggplot(df, aes(x=FDR_cutoff, y=number_of_features, group=stage, colour=stage)) +
#   geom_line() + scale_color_manual(values=stage_colours) + theme_minimal()
# graphics.off()
# 
# df_cut <- df %>% group_by(stage) %>% filter(FDR_cutoff < 0.1)
# png(paste0(plot_path, 'changing_FDR_cutoffs_zoom_logFC_0_linegraph.png'), height = 20, width = 20, units = 'cm', res = 400)
# ggplot(df_cut, aes(x=FDR_cutoff, y=number_of_features, group=stage, colour=stage)) +
#   geom_line() + scale_color_manual(values=stage_colours) + theme_minimal()
# graphics.off()
# 
# 
# ############. LogFC >= 1
# number_of_features <- c()
# for (i in sequence){
#   cutOff <- paste0("FDR <= ", i, " & Log2FC >= 1")
#   ids <- extract_ids(ss8_se, cutOff = cutOff, top_n = FALSE) # extract ids
#   print(length(ids))
#   number_of_features <- c(number_of_features, length(ids))
# }
# ss8_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "ss8")
# 
# number_of_features <- c()
# for (i in sequence){
#   cutOff <- paste0("FDR <= ", i, " & Log2FC >= 1")
#   ids <- extract_ids(ss4_se, cutOff = cutOff, top_n = FALSE) # extract ids
#   print(length(ids))
#   number_of_features <- c(number_of_features, length(ids))
# }
# 
# ss4_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "ss4")
# 
# number_of_features <- c()
# for (i in sequence){
#   cutOff <- paste0("FDR <= ", i, " & Log2FC >= 1")
#   ids <- extract_ids(HH7_se, cutOff = cutOff, top_n = FALSE) # extract ids
#   print(length(ids))
#   number_of_features <- c(number_of_features, length(ids))
# }
# HH7_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH7")
# 
# number_of_features <- c()
# for (i in sequence){
#   cutOff <- paste0("FDR <= ", i, " & Log2FC >= 1")
#   ids <- extract_ids(HH6_se, cutOff = cutOff, top_n = FALSE) # extract ids
#   print(length(ids))
#   number_of_features <- c(number_of_features, length(ids))
# }
# HH6_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH6")
# 
# number_of_features <- c()
# for (i in sequence){
#   cutOff <- paste0("FDR <= ", i, " & Log2FC >= 1")
#   ids <- extract_ids(HH5_se, cutOff = cutOff, top_n = FALSE) # extract ids
#   print(length(ids))
#   number_of_features <- c(number_of_features, length(ids))
# }
# HH5_df <- data.frame(FDR_cutoff = sequence, number_of_features = number_of_features, stage = "HH5")
# 
# df <- do.call("rbind", list(HH5_df, HH6_df, HH7_df, ss4_df, ss8_df))
# 
# png(paste0(plot_path, 'changing_FDR_cutoffs_logFC_1_linegraph.png'), height = 20, width = 20, units = 'cm', res = 400)
# ggplot(df, aes(x=FDR_cutoff, y=number_of_features, group=stage, colour=stage)) +
#   geom_line() + scale_color_manual(values=stage_colours) + theme_minimal()
# graphics.off()
# 
# df_cut <- df %>% group_by(stage) %>% filter(FDR_cutoff < 0.1)
# png(paste0(plot_path, 'changing_FDR_cutoffs_zoom_logFC_1_linegraph.png'), height = 20, width = 20, units = 'cm', res = 400)
# ggplot(df_cut, aes(x=FDR_cutoff, y=number_of_features, group=stage, colour=stage)) +
#   geom_line() + scale_color_manual(values=stage_colours) + theme_minimal()
# graphics.off()