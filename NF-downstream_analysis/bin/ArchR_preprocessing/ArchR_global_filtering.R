#!/usr/bin/env Rscript

### Script to preprocess in ArchR
print("global filtering")
# generates TSS enrichment, nucleosome score and number of fragment plots for each sample in dataset
# optionally filters whole ArchR object on nFrags using user-defined parameters

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ArchR)
library(GenomicFeatures)
library(parallel)
library(gridExtra)
library(grid)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-f", "--filter"), action = "store", type = "logical", help = "whether to filter data", default = FALSE),
  make_option(c("-m", "--factor"), action = "store", type = "double", help = "how many times to multiply SD by to create upper limit", default = 1),
  make_option(c("-v", "--verbose"), action = "store", type = "logical", help = "Verbose", default = TRUE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    
    data_path = "./output/NF-downstream_analysis/Upstream_processing/PREPROCESSING/preprocess/rds_files/"
    
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

set.seed(42)

############################## Read in ArchR project #######################################
label <- sub('_.*', '', list.files(data_path))
print(label)

ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

############################## QC plots across samples #######################################
##############################################################################################

############################## Plot TSS Enrichment #######################################
p2 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  baseSize = 25,
  pal = stage_colours
)
png(paste0(plot_path, 'TSS_enrichment_vln.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p2)
graphics.off()

p3 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges",
  baseSize = 25,
  pal = stage_colours
)
png(paste0(plot_path, 'TSS_enrichment_ridge.png'), height = 15, width = 21, units = 'cm', res = 400)
print(p3)
graphics.off()

p2 <- plotTSSEnrichment(ArchRProj = ArchR, groupBy = "stage", pal = stage_colours)
png(paste0(plot_path, 'TSS_enrichment_plot.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p2)
graphics.off()

print(paste0("Minimum TSS Enrichment score:", min(ArchR$TSSEnrichment)))
print(paste0("Maximum TSS Enrichment score:", max(ArchR$TSSEnrichment)))

############################## Plot Unique Fragments - log 10 #######################################
p4 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  baseSize = 25,
  pal = stage_colours
)
png(paste0(plot_path, 'fragment_log10_count_vln.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p4)
graphics.off()

p3 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges",
  baseSize = 25,
  pal = stage_colours
)
png(paste0(plot_path, 'fragment_log10_count_ridge.png'), height = 15, width = 21, units = 'cm', res = 400)
print(p3)
graphics.off()

############################## Plot Unique Fragments #######################################
p4 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "nFrags",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  baseSize = 25,
  pal = stage_colours
)
png(paste0(plot_path, 'fragment_count_vln.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p4)
graphics.off()

p <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "nFrags",
  plotAs = "ridges",
  baseSize = 25,
  pal = stage_colours)
png(paste0(plot_path, 'fragment_count_ridge.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p)
graphics.off()

print(paste0("Minimum number of fragments", min(ArchR$nFrags)))
print(paste0("Maximum number of fragments:", max(ArchR$nFrags)))

############################## Plot nucleosome banding #######################################

p2 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "NucleosomeRatio",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE,
  baseSize = 25,
  pal = stage_colours
)
png(paste0(plot_path, 'Nucleosome_ratio_vln.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p2)
graphics.off()

p3 <- plotGroups(
  ArchRProj = ArchR, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "NucleosomeRatio",
  plotAs = "ridges",
  baseSize = 25,
  pal = stage_colours
)
png(paste0(plot_path, 'Nucleosome_ratio_ridge.png'), height = 15, width = 21, units = 'cm', res = 400)
print(p3)
graphics.off()

p1 <- plotFragmentSizes(ArchRProj = ArchR, groupBy = "stage", pal = stage_colours,
                        threads = 1)
png(paste0(plot_path, 'nucleosome_banding_plot.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p1)
graphics.off()

############# Plot log10(Unique Fragments) vs TSS enrichment score #######################
df <- getCellColData(ArchR, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed") +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=25), legend.text=element_text(size=12),
        legend.key.size = unit(1.4, 'cm'), legend.title=element_text(size=18))

png(paste0(plot_path, 'fragments_vs_TSS.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p)
graphics.off()

############################## FILTERING #######################################
################################################################################
### we want the metrics to be normally distributed within each sample! ###

plot_path <- paste0(plot_path, "filtering/")
dir.create(plot_path, recursive = T)

### calculate max and min thresholds based on standard deviations per group
df <- tibble(cell_ids = getCellNames(ArchR), stage = ArchR$stage, nFrags = as.numeric(ArchR$nFrags))
df %>% group_by(stage)
medians <- aggregate(df$nFrags, list(df$stage), median)
sds <- aggregate(df$nFrags, list(df$stage), sd)
limits_df <- tibble(stage = medians$Group.1, median = medians$x,
                    upper_limit = medians$x + (sds$x * opt$factor))

if (opt$filter == FALSE) {
  
  print("data not filtered")
  ArchR_filtered <- ArchR
  
} else {
  
  ## plot thresholds
  p <- plotGroups(
    ArchRProj = ArchR, groupBy = "stage", colorBy = "cellColData", alpha = 0.4,
    name = "nFrags", plotAs = "violin", baseSize = 25, pal = stage_colours)
  p1 <- p + 
    geom_segment(aes(x = 0.5, xend = 1.5, y = limits_df$upper_limit[1], , yend = limits_df$upper_limit[1]), linetype = "dashed", colour = "black") +
    geom_segment(aes(x = 1.5, xend = 2.5, y = limits_df$upper_limit[2], , yend = limits_df$upper_limit[2]), linetype = "dashed", colour = "black") +
    geom_segment(aes(x = 2.5, xend = 3.5, y = limits_df$upper_limit[3], , yend = limits_df$upper_limit[3]), linetype = "dashed", colour = "black") +
    geom_segment(aes(x = 3.5, xend = 4.5, y = limits_df$upper_limit[4], , yend = limits_df$upper_limit[4]), linetype = "dashed", colour = "black") +
    geom_segment(aes(x = 4.5, xend = 5.5, y = limits_df$upper_limit[5], , yend = limits_df$upper_limit[5]), linetype = "dashed", colour = "black")
  png(paste0(plot_path, 'fragment_count_violin_thresholds.png'), height = 25, width = 25, units = 'cm', res = 400)
  print(p1)
  graphics.off()
  
  # filter cells using stage-specific upper threshold
  cell_ids <- c()
  for (i in c(1:5)) {
    print(i)
    stage <- limits_df$stage[i]
    print(stage)
    df_stage <- df[which(df$stage == stage), ]
    print(length(df_stage$cell_ids))
    df_stage <- df_stage[df_stage$nFrags < limits_df$upper_limit[i], ]
    print(length(df_stage$cell_ids))
    cell_ids <- c(cell_ids, df_stage$cell_ids)
  }
  ArchR_filtered <- ArchR[cell_ids, ]
  
}

# save ArchR project
saveArchRProject(ArchRProj = ArchR_filtered, outputDirectory = paste0(rds_path, "FullData_Save-ArchR"), load = FALSE)

############################ POST-FILTERING ####################################
################################################################################

p <- plotGroups(
  ArchRProj = ArchR_filtered, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "nFrags",
  plotAs = "violin",
  baseSize = 25,
  pal = stage_colours,
  alpha = 0.4)
png(paste0(plot_path, 'fragment_count_violin_filtered.png'), height = 25, width = 25, units = 'cm', res = 400)
print(p)
graphics.off()

p3 <- plotGroups(
  ArchRProj = ArchR_filtered, 
  groupBy = "stage", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges",
  baseSize = 25,
  pal = stage_colours)
png(paste0(plot_path, 'fragment_log10_count_violin_filtered.png'), height = 15, width = 21, units = 'cm', res = 400)
print(p3)
graphics.off()

unfiltered <- table(ArchR$stage)
filtered <- table(ArchR_filtered$stage)
cell_counts <- as_tibble(rbind(unfiltered, filtered))
cell_counts <- cbind(cell_counts, Total = rowSums(cell_counts))

png(paste0(plot_path, 'cell_counts_table.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()