#!/usr/bin/env Rscript

### Script to create seurat object, filter data, remove poor quality clusters and integrate across batches

############################## Load libraries #######################################
library(getopt)
library(Signac)
library(Seurat)
library(future)
library(tidyverse)
library(grid)
library(gridExtra)
library(clustree)
library(ggplot2)
library(dplyr)
library(scHelper)


############################## Set up script options #######################################
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)
test = TRUE

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    setwd("~/NF-downstream_analysis")
    ncores = 8
    ref_path = "../output/NF-luslab_sc_multiomic/reference/"
    
    if(test == TRUE){
      plot_path = "../output/NF-downstream_analysis/filtering/TEST/plots/"
      rds_path = "../output/NF-downstream_analysis/filtering/TEST/rds_files/"
      data_path = "../output/NF-downstream_analysis/preprocessing/rds_files/"
      }else{
      plot_path = "../output/NF-downstream_analysis/filtering/plots/"
      rds_path = "../output/NF-downstream_analysis/filtering/rds_files/"
      data_path = "../output/NF-downstream_analysis/preprocessing/rds_files/"}
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 75* 1024^3)
    plan()
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Read merged Seurat RDS object and make plots #######################################

seurat_all <- readRDS(paste0(data_path, "seurat_all.RDS"))

print(seurat_all)


############################## Set filtering thresholds #######################################

# These have been adjusted based on simulated plots below
filter_thresholds <- data.frame(pct_reads_in_peaks = c(0, 40, 50, 60), 
                                TSS.enrichment = c(0, 0.4, 2, 3), 
                                nucleosome_signal = c(Inf, 1.25, 1, 0.7),
                                peak_region_fragments_min = c(0, 210, 500, 1000),
                                peak_region_fragments_max = c(Inf, 1250, 800, 500),
                                row.names = c("unfilt", "low", "med", "high"))



png(paste0(plot_path, 'filter_thresholds.png'), height = 6, width = 30, units = 'cm', res = 400)
grid.arrange(top=textGrob("Selected Filter Thresholds", gp=gpar(fontsize=8, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(filter_thresholds, rows=NULL, theme = ttheme_minimal()))
graphics.off()


############################## Plot simulated thresholds to QC #######################################

simulated_plot_path = paste0(plot_path, "simulated_thresholds/")
dir.create(simulated_plot_path, recursive = T)

####    Number of fragments in peaks (peak_region_fragments) - Minimum cut off  ####

# Max per sample
maximums <- 
  seurat_all@meta.data %>%
  group_by(orig.ident) %>%
  summarise(max = max(peak_region_fragments, na.rm = TRUE))

# png(paste0(simulated_plot_path, 'peak_region_fragments_maximums.png'), height = 10, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Maximum Number of Fragments in Peaks", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(maximums, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# Simulation
filter_qc <- lapply(seq(from = 0, to = min(maximums$max), by = 10), function(cutoff){
  seurat_all@meta.data %>%
    filter(peak_region_fragments > cutoff) %>%
    group_by(orig.ident) %>%
    summarise(median = median(peak_region_fragments, na.rm = TRUE)) %>%
    mutate(median = as.integer(median)) %>%
    dplyr::rename(!! paste(cutoff) := median)
})
filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc) %>% reshape2::melt() %>% mutate(variable = as.integer(variable)*10)

png(paste0(simulated_plot_path, 'peak_region_fragments_simulation_min.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc, aes(x=variable, y=value, group=orig.ident)) +
  geom_line(aes(colour = orig.ident)) +
  xlab("Minimum cut off") +
  ylab("Median number of fragments in peaks") +
  ggtitle("Median number of fragments in peaks at simulated minimum filter thresholds") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x = filter_thresholds$peak_region_fragments_min[1])) +
  geom_line(aes(x = filter_thresholds$peak_region_fragments_min[2])) +
  geom_line(aes(x = filter_thresholds$peak_region_fragments_min[3])) +
  geom_line(aes(x = filter_thresholds$peak_region_fragments_min[4]))
graphics.off()

####    Number of fragments in peaks (peak_region_fragments) - Maximum cut off  ####

# Mins per sample
minimums <- 
  seurat_all@meta.data %>%
  group_by(orig.ident) %>%
  summarise(min = min(peak_region_fragments, na.rm = TRUE))

# png(paste0(simulated_plot_path, 'peak_region_fragments_maximums.png'), height = 10, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Maximum Number of Fragments in Peaks", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(maximums, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# Simulation
filter_qc <- lapply(seq(from = min(maximums$max), to = max(minimums$min), by = -10), function(cutoff){
  seurat_all@meta.data %>%
    filter(peak_region_fragments < cutoff) %>%
    group_by(orig.ident) %>%
    summarise(median = median(peak_region_fragments, na.rm = TRUE)) %>%
    mutate(median = as.integer(median)) %>%
    dplyr::rename(!! paste(cutoff) := median)
})
filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc) %>% reshape2::melt() %>% mutate(variable = as.integer(variable)*10)

png(paste0(simulated_plot_path, 'peak_region_fragments_simulation_max.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc, aes(x=variable, y=value, group=orig.ident)) +
  geom_line(aes(colour = orig.ident)) +
  xlab("Maximum cut off") +
  ylab("Median number of fragments in peaks") +
  ggtitle("Median number of fragments in peaks at simulated maximum filter thresholds") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x = filter_thresholds$peak_region_fragments_max[1])) +
  geom_line(aes(x = filter_thresholds$peak_region_fragments_max[2])) +
  geom_line(aes(x = filter_thresholds$peak_region_fragments_max[3])) +
  geom_line(aes(x = filter_thresholds$peak_region_fragments_max[4]))
graphics.off()

####    Percentage of reads in peaks (pct_reads_in_peaks) - Minimum cut off ####

# Max per sample
maximums <- 
  seurat_all@meta.data %>%
  group_by(orig.ident) %>%
  summarise(max = max(pct_reads_in_peaks, na.rm = TRUE))

# png(paste0(simulated_plot_path, 'pct_reads_in_peaks_maximums.png'), height = 10, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Maximum Percentage of Fragments in Peaks", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(maximums, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# Simulation
filter_qc <- lapply(seq(from = 0, to = min(maximums$max), by = 1), function(cutoff){
  seurat_all@meta.data %>%
    filter(pct_reads_in_peaks > cutoff) %>%
    group_by(orig.ident) %>%
    summarise(median = median(pct_reads_in_peaks, na.rm = TRUE)) %>%
    mutate(median = as.integer(median)) %>%
    dplyr::rename(!! paste(cutoff) := median)
})
filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc) %>% reshape2::melt() %>% mutate(variable = as.integer(variable))

png(paste0(simulated_plot_path, 'pct_reads_in_peaks_simulation.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc, aes(x=variable, y=value, group=orig.ident)) +
  geom_line(aes(colour = orig.ident)) +
  xlab("Minimum cut off") +
  ylab("Median % fragments in peaks") +
  ggtitle("Meadian percentage of fragments in peaks at simulated minimum filter thresholds") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x = filter_thresholds$pct_reads_in_peaks[1])) +
  geom_line(aes(x = filter_thresholds$pct_reads_in_peaks[2])) +
  geom_line(aes(x = filter_thresholds$pct_reads_in_peaks[3])) +
  geom_line(aes(x = filter_thresholds$pct_reads_in_peaks[4]))
graphics.off()

####    TSS enrichment (TSS.enrichment)  - Minimum cut off ####

# Max per sample
maximums <- 
  seurat_all@meta.data %>%
  group_by(orig.ident) %>%
  summarise(max = max(TSS.enrichment, na.rm = TRUE))

# png(paste0(simulated_plot_path, 'TSS_enrichment_maximums.png'), height = 10, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Maximum TSS Enrichment Score", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(maximums, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# Simulation
filter_qc <- lapply(seq(from = 0, to = min(maximums$max), by = 1), function(cutoff){
  seurat_all@meta.data %>%
    filter(TSS.enrichment > cutoff) %>%
    group_by(orig.ident) %>%
    summarise(median = median(TSS.enrichment, na.rm = TRUE)) %>%
    mutate(median = as.integer(median)) %>%
    dplyr::rename(!! paste(cutoff) := median)
})
filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc) %>% reshape2::melt() %>% mutate(variable = as.integer(variable)/10)

png(paste0(simulated_plot_path, 'TSS_enrichment_simulation.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc, aes(x=variable, y=value, group=orig.ident)) +
  geom_line(aes(colour = orig.ident)) +
  xlab("Minimum cut off") +
  ylab("Median TSS enrichment score") +
  ggtitle("TSS enrichment scores at simulated minimum filter thresholds") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x = filter_thresholds$TSS.enrichment[1])) +
  geom_line(aes(x = filter_thresholds$TSS.enrichment[2])) +
  geom_line(aes(x = filter_thresholds$TSS.enrichment[3])) +
  geom_line(aes(x = filter_thresholds$TSS.enrichment[4]))
graphics.off()

####    Nucleosome signal (nucleosome_signal) - Maximum cut off ####

# Max per sample
maximums <- 
  seurat_all@meta.data %>%
  group_by(orig.ident) %>%
  summarise(max = max(nucleosome_signal, na.rm = TRUE))

# png(paste0(simulated_plot_path, 'nucleosome_signal_maximums.png'), height = 10, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Maximum Nucleosome Signal Score", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(maximums, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# Min per sample
minimums <- 
  seurat_all@meta.data %>%
  group_by(orig.ident) %>%
  summarise(min = min(nucleosome_signal, na.rm = TRUE))

# png(paste0(simulated_plot_path, 'nucleosome_signal_minimums.png'), height = 10, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Minimum Nucleosome Signal Score", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(minimums, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# Simulation
filter_qc <- lapply(seq(from = min(maximums$max), to = max(minimums$min), by = -0.01), function(cutoff){
  seurat_all@meta.data %>%
    filter(nucleosome_signal < cutoff) %>%
    group_by(orig.ident) %>%
    summarise(median = median(nucleosome_signal, na.rm = TRUE)) %>%
    dplyr::rename(!! paste(cutoff) := median)
})
filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc) %>% reshape2::melt() %>% mutate(variable = as.integer(variable)/100)

png(paste0(simulated_plot_path, 'nucleosome_signal_simulation.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc, aes(x=variable, y=value, group=orig.ident)) +
  geom_line(aes(colour = orig.ident)) +
  xlab("Minimum cut off") +
  ylab("Median nucleosome signal score") +
  ggtitle("Nucleosome signal at simulated minimum filter thresholds") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x = filter_thresholds$nucleosome_signal[1])) +
  geom_line(aes(x = filter_thresholds$nucleosome_signal[2])) +
  geom_line(aes(x = filter_thresholds$nucleosome_signal[3])) +
  geom_line(aes(x = filter_thresholds$nucleosome_signal[4]))
graphics.off()

############################## Plot Vln plots with different thresholds #######################################

# Plot violins for nCount, nFeature and percent.mt at different filtering thresholds
filter_qc <- lapply(rownames(filter_thresholds), function(condition){
  seurat_all@meta.data %>%
    filter(peak_region_fragments > filter_thresholds[condition,'peak_region_fragments_min']) %>%
    filter(peak_region_fragments < filter_thresholds[condition,'peak_region_fragments_max']) %>%
    filter(pct_reads_in_peaks > filter_thresholds[condition,'pct_reads_in_peaks']) %>%
    filter(TSS.enrichment > filter_thresholds[condition,'TSS.enrichment']) %>%
    filter(nucleosome_signal < filter_thresholds[condition,'nucleosome_signal']) %>%
    dplyr::select(orig.ident, peak_region_fragments, pct_reads_in_peaks, TSS.enrichment, nucleosome_signal) %>%
    mutate(filter_condition = !!condition)
})

filter_qc <- do.call(rbind, filter_qc) %>%
  mutate(filter_condition = factor(filter_condition, rownames(filter_thresholds))) %>%
  reshape2::melt()

png(paste0(plot_path, 'violins_filter_thresholds.png'), height = 18, width = 40, units = 'cm', res = 400)
ggplot(filter_qc, aes(x = filter_condition, y = value, fill = orig.ident)) +
  geom_violin() +
  facet_wrap(~ variable, nrow = 4, strip.position = "left", scales = "free_y") +
  theme_minimal() +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(strip.placement = "outside")
graphics.off()

############################## Try low/med/high filtering thresholds to QC #######################################

# Calculate remaining cells following different filter thresholds
filter_qc <- lapply(rownames(filter_thresholds), function(condition){
  seurat_all@meta.data %>%
    filter(pct_reads_in_peaks > filter_thresholds[condition,'pct_reads_in_peaks']) %>%
    filter(TSS.enrichment > filter_thresholds[condition,'TSS.enrichment']) %>%
    filter(nucleosome_signal < filter_thresholds[condition,'nucleosome_signal']) %>%
    group_by(orig.ident) %>%
    tally() %>%
    dplyr::rename(!!condition := n)
})

filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc)

# Plot remaining cell counts
filter_qc <-  filter_qc %>% column_to_rownames('orig.ident')
filter_qc <- rbind(filter_qc, Total = colSums(filter_qc)) %>% rownames_to_column("orig.ident")

png(paste0(plot_path, 'remaining_cell_table.png'), height = 10, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(filter_qc, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path, 'remaining_cell_bar.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc[filter_qc$orig.ident != "Total",] %>% reshape2::melt(), aes(x=variable, y=value, fill=orig.ident)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  xlab("Filter Condition") +
  ylab("Cell Count") +
  ggtitle("Cell count after filtering") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()

# Calculate median fragment count per cell following different filter thresholds
filter_qc <- lapply(rownames(filter_thresholds), function(condition){
  seurat_all@meta.data %>%
    filter(pct_reads_in_peaks > filter_thresholds[condition,'pct_reads_in_peaks']) %>%
    filter(TSS.enrichment > filter_thresholds[condition,'TSS.enrichment']) %>%
    filter(nucleosome_signal < filter_thresholds[condition,'nucleosome_signal']) %>%
    group_by(orig.ident) %>%
    summarise(median = median(passed_filters, na.rm = TRUE)) %>%
    mutate(median = as.integer(median)) %>%
    dplyr::rename(!!condition := median)
})

# Plot median fragment count per cell
filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc)

png(paste0(plot_path, 'median_fragment_count_table.png'), height = 10, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Median Fragment Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(filter_qc, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path, 'median_fragment_count_line.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc %>% reshape2::melt(), aes(x=variable, y=value, group=orig.ident)) +
  geom_line(aes(colour = orig.ident)) +
  xlab("Filter Condition") +
  ylab("Median Fragment Count") +
  ggtitle("Median fragment count per cell after filtering") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()

############################## Filter data #######################################

seurat_all <- subset(seurat_all, subset = pct_reads_in_peaks > filter_thresholds['med','pct_reads_in_peaks'] &
                                                 peak_region_fragments > filter_thresholds['med','peak_region_fragments_min'] &
                                                 peak_region_fragments < filter_thresholds['med','peak_region_fragments_max'] &
                                                 TSS.enrichment > filter_thresholds['med','TSS.enrichment'] &
                                                 nucleosome_signal < filter_thresholds['med','nucleosome_signal'])

print(seurat_all)

print("seurat object filtered based on QC thresholds")

############################## Normalization and linear dimensional reduction #######################################

normalising_plot_path = paste0(plot_path, "normalising/")
dir.create(normalising_plot_path, recursive = T)

seurat_all <- RunTFIDF(seurat_all)
seurat_all <- FindTopFeatures(seurat_all, min.cutoff = 'q0')
seurat_all <- RunSVD(seurat_all)

png(paste0(normalising_plot_path, 'DepthCor.png'), height = 15, width = 21, units = 'cm', res = 400)
DepthCor(seurat_all)
graphics.off()
# just remove first LSI component

############################## Non-linear dimension reduction and clustering #######################################

clustering_plot_path = paste0(plot_path, "clustering/")
dir.create(clustering_plot_path, recursive = T)

seurat_all <- RunUMAP(object = seurat_all, reduction = 'lsi', dims = 2:30)
seurat_all <- FindNeighbors(object = seurat_all, reduction = 'lsi', dims = 2:30)
seurat_all <- FindClusters(object = seurat_all, verbose = FALSE, algorithm = 3)
DimPlot(object = seurat_all, label = TRUE) + NoLegend()

# Find optimal cluster resolution
png(paste0(clustering_plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = seurat_all, by = 0.2, prefix = "peaks_snn_res.") ## might need to change algorithm??
graphics.off()


############################## Identify poor quality clusters #######################################

# Use higher cluster resolution for filtering poor clusters
seurat_all <- FindClusters(object = seurat_all, verbose = FALSE, algorithm = 3, resolution = 2)

# Plot UMAP for clusters and developmental stage
png(paste0(clustering_plot_path, "UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(seurat_all)
graphics.off()

# Plot QC for each cluster
png(paste0(clustering_plot_path, "QCPlot.png"), width=32, height=28, units = 'cm', res = 200)
QCPlot(seurat_all, quantiles = c(0.25, 0.75), y_elements = c("pct_reads_in_peaks", "peak_region_fragments", 
                                                             "TSS.enrichment", "nucleosome_signal"),
       x_lab = c("% fragments in peaks", "Number of fragments in peaks", "TSS enrichment score", "Nucleosome signal score"))
graphics.off()

# Automatically find poor quality clusters
poor_clusters <- IdentifyOutliers(seurat_all, metrics = c("pct_reads_in_peaks", "peak_region_fragments", 
                                                                  "TSS.enrichment", "nucleosome_signal"), quantiles = c(0.25, 0.75))

# Plot UMAP for poor quality clusters
png(paste0(clustering_plot_path, "PoorClusters.png"), width=60, height=20, units = 'cm', res = 200)
ClusterDimplot(seurat_all, clusters = poor_clusters, plot_title = 'poor quality clusters')
graphics.off()

# Filter poor quality clusters
seurat_all_filtered <- subset(seurat_all, cells = rownames(filter(seurat_all@meta.data, seurat_clusters %in% poor_clusters)), invert = T)

print("seurat object filtered based on poor quality clusters")

# Plot table with remaining cell counts after full filtering
cell_counts <- data.frame(unfilt = summary(seurat_all@meta.data$orig.ident),
                          filtered = summary(seurat_all_filtered@meta.data$orig.ident))

cell_counts <- rbind(cell_counts, Total = colSums(cell_counts)) %>% rownames_to_column("orig.ident")

png(paste0(clustering_plot_path, 'final_remaining_cell_table.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# Save RDS output
saveRDS(seurat_all_filtered, paste0(rds_path, "seurat_all_filtered.RDS"), compress = FALSE)

