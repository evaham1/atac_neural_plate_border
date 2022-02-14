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
test = FALSE

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
      data_path = "../output/NF-downstream_analysis/test_input/"
    }else{
      plot_path = "../output/NF-downstream_analysis/filtering/plots/"
      rds_path = "../output/NF-downstream_analysis/filtering/rds_files/"
      data_path = "../output/NF-downstream_analysis/test_input/"}
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 155* 1024^3)
    plan()
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## FUNCTIONS #######################################

QC_metric_hist <- function(seurat_obj, QC_metric, identity = "stage", bin_width = 10, ident_cols = NULL, 
                           title = "QC Metric", xmin = NULL, xmax = NULL, intercept = NULL){
  df <- data.frame(ident = FetchData(object = seurat_obj, vars = c(identity)),
                   value = FetchData(object = seurat_obj, vars = c(QC_metric)))
  colnames(df) <- c("ident", "value")
  if(is.null(xmin) == TRUE){
    xmin <- min(df$value)
  }
  if(is.null(xmax) == TRUE){
    xmax <- max(df$value)
  }
  if(is.null(ident_cols) == TRUE){
    ident_cols <- palette(rainbow(length(unique(df$ident))))
  }
  ggplot(df, aes(x=value, fill=ident, color=ident)) + 
    geom_histogram(binwidth = bin_width) +
    facet_wrap(~ ident) +
    xlab("QC Metric") +
    ylab("Frequency") +
    ggtitle(title) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(strip.text = element_blank()) +
    scale_fill_manual(values = ident_cols) +
    scale_color_manual(values = ident_cols) +
    xlim(xmin, xmax) +
    geom_vline(xintercept = intercept, linetype="dashed", size=1) 
}


############################## Read in Seurat RDS object and fragment files #######################################

seurat_all <- readRDS(paste0(data_path, "rds_files/seurat_all.RDS"))
print(seurat_all)

# create a vector with all the new paths, in the correct order for your list of fragment objects
paths <- list.dirs(paste0(data_path, "cellranger_atac_output/"), recursive = FALSE, full.names = TRUE)
input <- data.frame(sample = sub('.*/', '', paths), 
                    matrix_path = paste0(paths, "/outs/filtered_peak_bc_matrix.h5"),
                    metadata_path = paste0(paths, "/outs/singlecell.csv"),
                    fragments_path = paste0(paths, "/outs/fragments.tsv.gz"))
new.paths <- as.list(input$fragments_path)
frags <- Fragments(seurat_all)  # get list of fragment objects
Fragments(seurat_all) <- NULL  # remove fragment information from assay

for (i in seq_along(frags)) {
  frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
}

Fragments(seurat_all) <- frags # assign updated list back to the object

############################## Set stage as metadata and colour them #######################################
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

seurat_all <- AddMetaData(seurat_all, substr(seurat_all@meta.data$orig.ident, 1, 3), col.name = "stage")
seurat_all@meta.data$stage <- factor(seurat_all@meta.data$stage, levels = stage_order)

stage_cols <- stage_colours[levels(droplevels(seurat_all@meta.data$stage))]

############################## TSS Enrichment Plot before filtering ######################################

before_plot_path = paste0(plot_path, "before_filtering/")
dir.create(before_plot_path, recursive = T)

png(paste0(before_plot_path, 'QC_TSS.png'), height = 15, width = 21, units = 'cm', res = 400)
TSSPlot(seurat_all, group.by = 'stage') + NoLegend()
graphics.off()

############################## Nucleosome Banding Plot before filtering #######################################

# just using first 100bps of chromosome one as takes a long time - still takes too long!
#png(paste0(before_plot_path, 'QC_Nucleosome_banding.png'), height = 15, width = 21, units = 'cm', res = 400)
#FragmentHistogram(object = seurat_all, region = "chr1-1-100", group.by = 'stage')
#graphics.off()

############################## Histograms of QC metrics before filtering #######################################

### Percentage of reads in peaks
png(paste0(before_plot_path, 'pct_reads_in_peaks_hist.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "pct_reads_in_peaks", bin_width = 0.1,
               ident_cols = stage_cols, title = "Percentage of fragments in peaks")
graphics.off()

png(paste0(before_plot_path, 'pct_reads_in_peaks_hist_zoom.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "pct_reads_in_peaks", bin_width = 0.1,
               ident_cols = stage_cols, title = "Percentage of fragments in peaks - minimum threshold 50%", xmin = 0, xmax = 60, intercept = 50)
graphics.off()

### TSS enrichment score
png(paste0(before_plot_path, 'TSS.enrichment_hist.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "TSS.enrichment", bin_width = 0.1,
               ident_cols = stage_cols, title = "TSS enrichment score")
graphics.off()

png(paste0(before_plot_path, 'TSS.enrichment_hist_zoom.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "TSS.enrichment", bin_width = 0.1,
               ident_cols = stage_cols, title = "TSS enrichment score - minimum threshold 2", xmin = 0, xmax = 4, intercept = 2)
graphics.off()

### Nucleosome signal score
png(paste0(before_plot_path, 'nucleosome_signal_hist.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "nucleosome_signal", bin_width = 0.1,
               ident_cols = stage_cols, title = "Nucleosome signal score")
graphics.off()

png(paste0(before_plot_path, 'nucleosome_signal_hist_zoom.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "nucleosome_signal", bin_width = 0.1,
               ident_cols = stage_cols, title = "Nucleosome signal score - maximum threshold 1.5", xmin = 0, xmax = 3, intercept = 1.5)
graphics.off()

### Number of reads in peaks - NB: thresholds vary depending on whether data is test data or not!
png(paste0(before_plot_path, 'peak_region_fragments_hist.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "peak_region_fragments", bin_width = 100,
               ident_cols = stage_cols, title = "Number of fragments in peaks")
graphics.off()

if (length(unique(seurat_all@meta.data$stage)) > 2){
  png(paste0(before_plot_path, 'peak_region_fragments_hist_zoom_min.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all, QC_metric = "peak_region_fragments", bin_width = 100,
                 ident_cols = stage_cols, title = "Number of fragments in peaks - minimum threshold 2000", xmin = 0, xmax = 5000, intercept = 2000)
  )
  graphics.off()
  
  png(paste0(before_plot_path, 'peak_region_fragments_hist_zoom_max.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all, QC_metric = "peak_region_fragments", bin_width = 100,
                 ident_cols = stage_cols, title = "Number of fragments in peaks - maximum threshold 50000", xmin = 20000, xmax = 100000, intercept = 50000)
  )
  graphics.off()
} else {
  png(paste0(before_plot_path, 'peak_region_fragments_hist_zoom_min.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all, QC_metric = "peak_region_fragments", bin_width = 100,
                 ident_cols = stage_cols, title = "Number of fragments in peaks - minimum threshold 200", xmin = 0, xmax = 2000, intercept = 200)
  )
  graphics.off()
  
  png(paste0(before_plot_path, 'peak_region_fragments_hist_zoom_max.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all, QC_metric = "peak_region_fragments", bin_width = 100,
                 ident_cols = stage_cols, title = "Number of fragments in peaks - maximum threshold 1800", xmin = 1000, xmax = 2000, intercept = 1800)
  )
  graphics.off()
}

### All QC metrics
png(paste0(before_plot_path, 'All_QC_VlnPlot.png'), height = 15, width = 28, units = 'cm', res = 400)
VlnPlot(
  object = seurat_all, group.by = "mitochondrial",
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0, ncol = 4)
graphics.off()

############################## ########### #######################################
############################## Filter data #######################################

after_plot_path = paste0(plot_path, "after_filtering/")
dir.create(after_plot_path, recursive = T)

if (length(unique(seurat_all@meta.data$stage)) > 2){
  seurat_all_filtered <- subset(seurat_all, subset = 
                       pct_reads_in_peaks > 50 &
                       TSS.enrichment > 2 &
                       nucleosome_signal < 1.5 &
                       peak_region_fragments < 50000 &
                       peak_region_fragments > 2000)
} else {
  seurat_all_filtered <- subset(seurat_all, subset = 
                                  pct_reads_in_peaks > 50 &
                                  TSS.enrichment > 2 &
                                  nucleosome_signal < 1.5 &
                                  peak_region_fragments < 1800 &
                                  peak_region_fragments > 200)
}

print(seurat_all_filtered)

print("seurat object filtered based on QC thresholds")

# Plot table with remaining cell counts after full filtering
cell_counts <- data.frame(unfilt = summary(seurat_all@meta.data$stage),
                          filtered = summary(seurat_all_filtered@meta.data$stage))

cell_counts <- rbind(cell_counts, Total = colSums(cell_counts)) %>% rownames_to_column("stage")

png(paste0(after_plot_path, 'remaining_cell_table.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# set as seurat_all for further filtering
seurat_all <- seurat_all_filtered

############################## TSS Enrichment Plot after filtering #######################################

png(paste0(after_plot_path, 'QC_TSS.png'), height = 15, width = 21, units = 'cm', res = 400)
TSSPlot(seurat_all, group.by = 'stage') + NoLegend()
graphics.off()

############################## Histograms of QC metrics after filtering #######################################

### Percentage of reads in peaks
png(paste0(after_plot_path, 'pct_reads_in_peaks_hist.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "pct_reads_in_peaks", bin_width = 0.1,
               ident_cols = stage_cols, title = "Percentage of fragments in peaks")
graphics.off()

png(paste0(after_plot_path, 'pct_reads_in_peaks_hist_zoom.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "pct_reads_in_peaks", bin_width = 0.1,
               ident_cols = stage_cols, title = "Percentage of fragments in peaks - minimum threshold 50%", xmin = 0, xmax = 60, intercept = 50)
graphics.off()

### TSS enrichment score
png(paste0(after_plot_path, 'TSS.enrichment_hist.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "TSS.enrichment", bin_width = 0.1,
               ident_cols = stage_cols, title = "TSS enrichment score")
graphics.off()

png(paste0(after_plot_path, 'TSS.enrichment_hist_zoom.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "TSS.enrichment", bin_width = 0.1,
               ident_cols = stage_cols, title = "TSS enrichment score - minimum threshold 2", xmin = 0, xmax = 4, intercept = 2)
graphics.off()

### Nucleosome signal score
png(paste0(after_plot_path, 'nucleosome_signal_hist.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "nucleosome_signal", bin_width = 0.1,
               ident_cols = stage_cols, title = "Nucleosome signal score")
graphics.off()

png(paste0(after_plot_path, 'nucleosome_signal_hist_zoom.png'), height = 15, width = 21, units = 'cm', res = 400)
QC_metric_hist(seurat_obj = seurat_all, QC_metric = "nucleosome_signal", bin_width = 0.1,
               ident_cols = stage_cols, title = "Nucleosome signal score - maximum threshold 1.5", xmin = 0, xmax = 3, intercept = 1.5)
graphics.off()

### Number of reads in peaks
if (length(unique(seurat_all@meta.data$stage)) > 2){
  png(paste0(after_plot_path, 'peak_region_fragments_hist_zoom_min.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all, QC_metric = "peak_region_fragments", bin_width = 100,
                   ident_cols = stage_cols, title = "Number of fragments in peaks - minimum threshold 2000", xmin = 0, xmax = 5000, intercept = 2000)
  )
  graphics.off()
  
  png(paste0(after_plot_path, 'peak_region_fragments_hist_zoom_max.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all, QC_metric = "peak_region_fragments", bin_width = 100,
                   ident_cols = stage_cols, title = "Number of fragments in peaks - maximum threshold 50000", xmin = 20000, xmax = 100000, intercept = 50000)
  )
  graphics.off()
} else {
  png(paste0(after_plot_path, 'peak_region_fragments_hist_zoom_min.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all, QC_metric = "peak_region_fragments", bin_width = 100,
                   ident_cols = stage_cols, title = "Number of fragments in peaks - minimum threshold 200", xmin = 0, xmax = 2000, intercept = 200)
  )
  graphics.off()
  
  png(paste0(after_plot_path, 'peak_region_fragments_hist_zoom_max.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all, QC_metric = "peak_region_fragments", bin_width = 100,
                   ident_cols = stage_cols, title = "Number of fragments in peaks - maximum threshold 1800", xmin = 1000, xmax = 2000, intercept = 1800)
  )
  graphics.off()
}

### All QC metrics
png(paste0(after_plot_path, 'All_QC_VlnPlot.png'), height = 15, width = 28, units = 'cm', res = 400)
VlnPlot(
  object = seurat_all, group.by = "mitochondrial",
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0, ncol = 4)
graphics.off()

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

# Find optimal cluster resolution -- need to change algorithm to match the one above (SLM)??
png(paste0(clustering_plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = seurat_all, by = 0.2, prefix = "peaks_snn_res.")
graphics.off()

############################## UMAP Visualtions before cluster filtering #######################################

png(paste0(clustering_plot_path, "UMAP.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(object = seurat_all, label = TRUE) + NoLegend()
graphics.off()

# UMAP for clusters and developmental stage
png(paste0(clustering_plot_path, "ClustStagePlot_UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(seurat_all)
graphics.off()

png(paste0(clustering_plot_path, "stage_umap.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(seurat_all, group.by = 'stage', label = TRUE, label.size = 12,
        label.box = TRUE, repel = TRUE,
        pt.size = 0.9, cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none",
                 plot.title = element_blank())
graphics.off()

############################## Identify poor quality clusters #######################################

# Don't use pct reads in peaks as it may be biological

# nucleosome_signal
png(paste0(clustering_plot_path, "QCPlot_nucleosome_signal.png"), width=20, height=40, units = 'cm', res = 200)
QCPlot(seurat_all, stage = "stage", quantiles = c(0, 0.85), y_elements = c("nucleosome_signal"),
       x_lab = c("Cluster"))
graphics.off()

poor_clusters_nucleosome_signal <- IdentifyOutliers(seurat_all, metrics = c("nucleosome_signal"), quantiles = c(0, 0.85))

png(paste0(clustering_plot_path, "PoorClusters_nucleosome_signal.png"), width=60, height=20, units = 'cm', res = 200)
ClusterDimplot(seurat_all, clusters = poor_clusters_nucleosome_signal, plot_title = 'poor quality clusters')
graphics.off()

# TSS.enrichment
png(paste0(clustering_plot_path, "QCPlot_TSS.enrichment.png"), width=20, height=40, units = 'cm', res = 200)
QCPlot(seurat_all, stage = "stage", quantiles = c(0.15, 1), y_elements = c("TSS.enrichment"),
       x_lab = c("Cluster"))
graphics.off()

poor_clusters_TSS_enrichment <- IdentifyOutliers(seurat_all, metrics = c("TSS.enrichment"), quantiles = c(0.15, 1))

png(paste0(clustering_plot_path, "PoorClusters_TSS.enrichment.png"), width=60, height=20, units = 'cm', res = 200)
ClusterDimplot(seurat_all, clusters = poor_clusters_TSS_enrichment, plot_title = 'poor quality clusters')
graphics.off()

# peak_region_fragments
# QCPlot(seurat_all, stage = "stage", quantiles = c(0.2, 0.8), y_elements = c("peak_region_fragments"),
#        x_lab = c("Cluster"))

# poor_clusters_peak_region_fragments <- IdentifyOutliers(seurat_all, metrics = c("peak_region_fragments"), quantiles = c(0.2, 0.8))

# png(paste0(clustering_plot_path, "PoorClusters_peak_region_fragments.png"), width=60, height=20, units = 'cm', res = 200)
# ClusterDimplot(seurat_all, clusters = poor_clusters, plot_title = 'poor quality clusters')
# graphics.off()

############################## Filter out poor quality clusters and count remaining cells #######################################

seurat_all_filtered <- subset(seurat_all, cells = rownames(filter(seurat_all@meta.data, seurat_clusters %in% poor_clusters_nucleosome_signal)), invert = T)
print("seurat object filtered based on poor quality clusters")

# Plot table with remaining cell counts after full filtering
cell_counts <- data.frame(unfilt = summary(seurat_all@meta.data$orig.ident),
                          filtered = summary(seurat_all_filtered@meta.data$orig.ident))

cell_counts <- rbind(cell_counts, Total = colSums(cell_counts)) %>% rownames_to_column("orig.ident")

png(paste0(clustering_plot_path, 'final_remaining_cell_table.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

############################## UMAP Visulations after cluster filtering #######################################

clustering_plot_path_filtered = paste0(plot_path, "clustering_filtered/")
dir.create(clustering_plot_path_filtered, recursive = T)

# Recluster after removing cells
seurat_all_filtered <- RunUMAP(object = seurat_all_filtered, reduction = 'lsi', dims = 2:30)
seurat_all_filtered <- FindNeighbors(object = seurat_all_filtered, reduction = 'lsi', dims = 2:30)
seurat_all_filtered <- FindClusters(object = seurat_all_filtered, verbose = FALSE, algorithm = 3)

png(paste0(clustering_plot_path_filtered, "UMAP.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(object = seurat_all_filtered, label = TRUE) + NoLegend()
graphics.off()

# UMAP for clusters and developmental stage
png(paste0(clustering_plot_path_filtered, "ClustStagePlot_UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(seurat_all_filtered)
graphics.off()

png(paste0(clustering_plot_path_filtered, "stage_umap.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(seurat_all_filtered, group.by = 'stage', label = TRUE, label.size = 12,
        label.box = TRUE, repel = TRUE,
        pt.size = 0.9, cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none",
                 plot.title = element_blank())
graphics.off()

############################## Save output #######################################

# Save RDS output
saveRDS(seurat_all_filtered, paste0(rds_path, "seurat_all_filtered.RDS"), compress = FALSE)


############################## ########### ####################################### ############ ############ ######
############################## Filtering data without using fragment counts #######################################

after_plot_path = paste0(plot_path, "after_filtering_no_frag_count/")
dir.create(after_plot_path, recursive = T)

seurat_all_filtered_pct <- subset(seurat_all, subset = 
                       pct_reads_in_peaks > 50 )
seurat_all_filtered_pct_TSS <- subset(seurat_all, subset = 
                       pct_reads_in_peaks > 50 &
                       TSS.enrichment > 2 )
seurat_all_filtered_pct_TSS_nucleosome <- subset(seurat_all, subset = 
                       pct_reads_in_peaks > 50 &
                       TSS.enrichment > 2 &
                       nucleosome_signal < 1.5 )
                       
# Plot table with remaining cell counts after filtering
cell_counts <- data.frame(unfilt = summary(seurat_all@meta.data$stage),
                          pct = summary(seurat_all_filtered_pct@meta.data$stage),
                          pct_and_TSS = summary(seurat_all_filtered_pct_TSS@meta.data$stage),
                          pct_and_TSS_and_nucleosome = summary(seurat_all_filtered_pct_TSS_nucleosome@meta.data$stage))

cell_counts <- rbind(cell_counts, Total = colSums(cell_counts)) %>% rownames_to_column("stage")

png(paste0(after_plot_path, 'remaining_cell_table_no_fragment_count.png'), height = 20, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

print("seurat object filtered with metrics other than fragment count")

### Number of reads in peaks
if (length(unique(seurat_all_filtered_pct_TSS_nucleosome@meta.data$stage)) > 2){
  png(paste0(after_plot_path, 'peak_region_fragments_hist_zoom_min_no_fragment_count.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all_filtered_pct_TSS_nucleosome, QC_metric = "peak_region_fragments", bin_width = 100,
                   ident_cols = stage_cols, title = "Number of fragments in peaks - minimum threshold 2000", xmin = 0, xmax = 5000, intercept = 2000)
  )
  graphics.off()
  
  png(paste0(after_plot_path, 'peak_region_fragments_hist_zoom_max_no_fragment_count.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all_filtered_pct_TSS_nucleosome, QC_metric = "peak_region_fragments", bin_width = 100,
                   ident_cols = stage_cols, title = "Number of fragments in peaks - maximum threshold 50000", xmin = 20000, xmax = 100000, intercept = 50000)
  )
  graphics.off()
} else {
  png(paste0(after_plot_path, 'peak_region_fragments_hist_zoom_min_no_fragment_count.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all_filtered_pct_TSS_nucleosome, QC_metric = "peak_region_fragments", bin_width = 100,
                   ident_cols = stage_cols, title = "Number of fragments in peaks - minimum threshold 200", xmin = 0, xmax = 2000, intercept = 200)
  )
  graphics.off()
  
  png(paste0(after_plot_path, 'peak_region_fragments_hist_zoom_max_no_fragment_count.png'), height = 15, width = 21, units = 'cm', res = 400)
  print(
    QC_metric_hist(seurat_obj = seurat_all_filtered_pct_TSS_nucleosome, QC_metric = "peak_region_fragments", bin_width = 100,
                   ident_cols = stage_cols, title = "Number of fragments in peaks - maximum threshold 1800", xmin = 1000, xmax = 2000, intercept = 1800)
  )
  graphics.off()
}

### All QC metrics
png(paste0(after_plot_path, 'All_QC_VlnPlot_no_fragment_count.png'), height = 15, width = 28, units = 'cm', res = 400)
VlnPlot(
  object = seurat_all_filtered_pct_TSS_nucleosome, group.by = "mitochondrial",
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0, ncol = 4)
graphics.off()


############################## Pre-processing #######################################

clustering_plot_path = paste0(plot_path, "clustering_no_fragment_count/")
dir.create(clustering_plot_path, recursive = T)

seurat_all <- RunTFIDF(seurat_all)
seurat_all <- FindTopFeatures(seurat_all, min.cutoff = 'q0')
seurat_all <- RunSVD(seurat_all)

png(paste0(clustering_plot_path, 'DepthCor.png'), height = 15, width = 21, units = 'cm', res = 400)
DepthCor(seurat_all)
graphics.off()

seurat_all <- RunUMAP(object = seurat_all, reduction = 'lsi', dims = 2:30)
seurat_all <- FindNeighbors(object = seurat_all, reduction = 'lsi', dims = 2:30)
seurat_all <- FindClusters(object = seurat_all, verbose = FALSE, algorithm = 3)

############################## UMAP Visualtions before cluster filtering #######################################

png(paste0(clustering_plot_path, "UMAP.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(object = seurat_all, label = TRUE) + NoLegend()
graphics.off()

# UMAP for clusters and developmental stage
png(paste0(clustering_plot_path, "ClustStagePlot_UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(seurat_all)
graphics.off()

png(paste0(clustering_plot_path, "stage_umap.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(seurat_all, group.by = 'stage', label = TRUE, label.size = 12,
        label.box = TRUE, repel = TRUE,
        pt.size = 0.9, cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none",
                 plot.title = element_blank())
graphics.off()

############################## Identify poor quality clusters #######################################

# Don't use pct reads in peaks as it may be biological

# nucleosome_signal
png(paste0(clustering_plot_path, "QCPlot_nucleosome_signal.png"), width=20, height=40, units = 'cm', res = 200)
QCPlot(seurat_all, stage = "stage", quantiles = c(0, 0.85), y_elements = c("nucleosome_signal"),
       x_lab = c("Cluster"))
graphics.off()

poor_clusters_nucleosome_signal <- IdentifyOutliers(seurat_all, metrics = c("nucleosome_signal"), quantiles = c(0, 0.85))

png(paste0(clustering_plot_path, "PoorClusters_nucleosome_signal.png"), width=60, height=20, units = 'cm', res = 200)
ClusterDimplot(seurat_all, clusters = poor_clusters_nucleosome_signal, plot_title = 'poor quality clusters')
graphics.off()

# TSS.enrichment
png(paste0(clustering_plot_path, "QCPlot_TSS.enrichment.png"), width=20, height=40, units = 'cm', res = 200)
QCPlot(seurat_all, stage = "stage", quantiles = c(0.15, 1), y_elements = c("TSS.enrichment"),
       x_lab = c("Cluster"))
graphics.off()

poor_clusters_TSS_enrichment <- IdentifyOutliers(seurat_all, metrics = c("TSS.enrichment"), quantiles = c(0.15, 1))

png(paste0(clustering_plot_path, "PoorClusters_TSS.enrichment.png"), width=60, height=20, units = 'cm', res = 200)
ClusterDimplot(seurat_all, clusters = poor_clusters_TSS_enrichment, plot_title = 'poor quality clusters')
graphics.off()

# peak_region_fragments
# QCPlot(seurat_all, stage = "stage", quantiles = c(0.2, 0.8), y_elements = c("peak_region_fragments"),
#        x_lab = c("Cluster"))

# poor_clusters_peak_region_fragments <- IdentifyOutliers(seurat_all, metrics = c("peak_region_fragments"), quantiles = c(0.2, 0.8))

# png(paste0(clustering_plot_path, "PoorClusters_peak_region_fragments.png"), width=60, height=20, units = 'cm', res = 200)
# ClusterDimplot(seurat_all, clusters = poor_clusters, plot_title = 'poor quality clusters')
# graphics.off()

############################## Filter out poor quality clusters and count remaining cells #######################################

seurat_all_filtered <- subset(seurat_all, cells = rownames(filter(seurat_all@meta.data, seurat_clusters %in% poor_clusters_nucleosome_signal)), invert = T)
print("seurat object filtered based on poor quality clusters")

# Plot table with remaining cell counts after full filtering
cell_counts <- data.frame(unfilt = summary(seurat_all@meta.data$orig.ident),
                          filtered = summary(seurat_all_filtered@meta.data$orig.ident))

cell_counts <- rbind(cell_counts, Total = colSums(cell_counts)) %>% rownames_to_column("orig.ident")

png(paste0(clustering_plot_path, 'final_remaining_cell_table.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

############################## UMAP Visulations after cluster filtering #######################################

clustering_plot_path_filtered = paste0(plot_path, "clustering_filtered_no_frag_count/")
dir.create(clustering_plot_path_filtered, recursive = T)

# Recluster after removing cells
seurat_all_filtered <- RunUMAP(object = seurat_all_filtered, reduction = 'lsi', dims = 2:30)
seurat_all_filtered <- FindNeighbors(object = seurat_all_filtered, reduction = 'lsi', dims = 2:30)
seurat_all_filtered <- FindClusters(object = seurat_all_filtered, verbose = FALSE, algorithm = 3)

png(paste0(clustering_plot_path_filtered, "UMAP.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(object = seurat_all_filtered, label = TRUE) + NoLegend()
graphics.off()

# UMAP for clusters and developmental stage
png(paste0(clustering_plot_path_filtered, "ClustStagePlot_UMAP.png"), width=40, height=20, units = 'cm', res = 200)
ClustStagePlot(seurat_all_filtered)
graphics.off()

png(paste0(clustering_plot_path_filtered, "stage_umap.png"), width=20, height=20, units = 'cm', res = 200)
DimPlot(seurat_all_filtered, group.by = 'stage', label = TRUE, label.size = 12,
        label.box = TRUE, repel = TRUE,
        pt.size = 0.9, cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none",
                 plot.title = element_blank())
graphics.off()




############################## **** ARCHIVED **** #######################################
#########################################################################################
#########################################################################################

############################## Set filtering thresholds #######################################

# These have been adjusted based on simulated plots below
# filter_thresholds <- data.frame(pct_reads_in_peaks = c(0, 50, 57, 62), 
#                                 TSS.enrichment = c(0, 2.3, 2.8, 3.2), 
#                                 nucleosome_signal = c(Inf, 2, 1.8, 1.6),
#                                 peak_region_fragments_min = c(0, 250, 500, 1000),
#                                 peak_region_fragments_max = c(Inf, 40000, 35000, 30000),
#                                 row.names = c("unfilt", "low", "med", "high"))
# 
# 
# 
# png(paste0(plot_path, 'filter_thresholds.png'), height = 6, width = 30, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Selected Filter Thresholds", gp=gpar(fontsize=8, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(filter_thresholds, rows=NULL, theme = ttheme_minimal()))
# graphics.off()


############################## Plot simulated thresholds to QC #######################################

# simulated_plot_path = paste0(plot_path, "simulated_thresholds/")
# dir.create(simulated_plot_path, recursive = T)
# 
# ####    Number of fragments in peaks (peak_region_fragments) - Minimum cut off  ####
# 
# ######    NEED TO SORT OUT x/y AXES AND ALSO WHY MAX CUT OFF GOES DOWNWARDS
# 
# # # Max per sample
# maximums <-
#   seurat_all@meta.data %>%
#   group_by(orig.ident) %>%
#   summarise(max = max(peak_region_fragments, na.rm = TRUE))
# 
# # Median simulation
# filter_qc <- lapply(seq(from = 0, to = min(maximums$max), by = 10), function(cutoff){
#   seurat_all@meta.data %>%
#     filter(peak_region_fragments > cutoff) %>%
#     group_by(orig.ident) %>%
#     summarise(median = median(peak_region_fragments, na.rm = TRUE)) %>%
#     mutate(median = as.integer(median)) %>%
#     dplyr::rename(!! paste(cutoff) := median)
# })
# filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc) %>% reshape2::melt() %>% mutate(variable = as.integer(variable))
# 
# png(paste0(simulated_plot_path, 'peak_region_fragments_simulation_min.png'), height = 15, width = 21, units = 'cm', res = 400)
# ggplot(filter_qc, aes(x=variable*10, y=value, group=orig.ident)) +
#   geom_line(aes(colour = orig.ident)) +
#   xlab("Minimum cut off") +
#   ylab("Median number of fragments in peaks") +
#   ggtitle("Median number of fragments in peaks at simulated minimum filter thresholds") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5))
# graphics.off()
# 
# ####    Number of fragments in peaks (peak_region_fragments) - Maximum cut off  ####
# 
# # Mins per sample
# minimums <- 
#   seurat_all@meta.data %>%
#   group_by(orig.ident) %>%
#   summarise(min = min(peak_region_fragments, na.rm = TRUE))
# 
# # Median simulation
# filter_qc <- lapply(seq(from = min(maximums$max), to = max(minimums$min), by = -10), function(cutoff){
#   seurat_all@meta.data %>%
#     filter(peak_region_fragments < cutoff) %>%
#     group_by(stage) %>%
#     summarise(median = median(peak_region_fragments, na.rm = TRUE)) %>%
#     mutate(median = as.integer(median)) %>%
#     dplyr::rename(!! paste(cutoff) := median)
# })
# filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc) %>% reshape2::melt() %>% mutate(variable = as.integer(variable))
# 
# png(paste0(simulated_plot_path, 'peak_region_fragments_simulation_max.png'), height = 15, width = 21, units = 'cm', res = 400)
# ggplot(filter_qc, aes(x=variable*10, y=value, group=orig.ident)) +
#   geom_line(aes(colour = orig.ident)) +
#   xlab("Maximum cut off") +
#   ylab("Median number of fragments in peaks") +
#   ggtitle("Median number of fragments in peaks at simulated maximum filter thresholds") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5))
# graphics.off()

####    Percentage of reads in peaks (pct_reads_in_peaks) - Minimum cut off ####

# # Max per sample
# maximums <- 
#   seurat_all@meta.data %>%
#   group_by(stage) %>%
#   summarise(max = max(pct_reads_in_peaks, na.rm = TRUE))
# 
# # Median simulation
# png(paste0(simulated_plot_path, 'pct_reads_in_peaks_medians_simulation.png'), height = 15, width = 21, units = 'cm', res = 400)
# simulation_medians(seurat_all, QC_metric = "pct_reads_in_peaks", idents = "stage", from = 0, to = min(maximums$max), by = 1, 
#                    ident_cols = stage_colours)
# graphics.off()
# 
# # Cell count simulation
# png(paste0(simulated_plot_path, 'pct_reads_in_peaks_cell_count_simulation.png'), height = 15, width = 21, units = 'cm', res = 400)
# simulation_cell_count(seurat_all, QC_metric = "pct_reads_in_peaks", idents = "stage", from = 0, to = min(maximums$max), by = 1, 
#                       ident_cols = stage_colours)
# graphics.off()
# 
# 
# ####    TSS enrichment (TSS.enrichment)  - Minimum cut off ####
# 
# # Max per sample
# maximums <- 
#   seurat_all@meta.data %>%
#   group_by(stage) %>%
#   summarise(max = max(TSS.enrichment, na.rm = TRUE))
# 
# # Median simulation
# png(paste0(simulated_plot_path, 'TSS.enrichment_medians_simulation.png'), height = 15, width = 21, units = 'cm', res = 400)
# simulation_medians(seurat_all, QC_metric = "TSS.enrichment", idents = "stage", from = 0, to = min(maximums$max), by = 0.1, 
#                    ident_cols = stage_colours)
# graphics.off()
# 
# # Cell count simulation
# png(paste0(simulated_plot_path, 'TSS.enrichment_cell_count_simulation.png'), height = 15, width = 21, units = 'cm', res = 400)
# simulation_cell_count(seurat_all, QC_metric = "TSS.enrichment", idents = "stage", from = 0, to = min(maximums$max), by = 0.1, 
#                       ident_cols = stage_colours)
# graphics.off()
# 
# 
# ####    Nucleosome signal (nucleosome_signal) - Maximum cut off ####
# 
# # Max per sample
# maximums <- 
#   seurat_all@meta.data %>%
#   group_by(orig.ident) %>%
#   summarise(max = max(nucleosome_signal, na.rm = TRUE))
# 
# # Min per sample
# minimums <- 
#   seurat_all@meta.data %>%
#   group_by(orig.ident) %>%
#   summarise(min = min(nucleosome_signal, na.rm = TRUE))
# 
# # Median simulation
# png(paste0(simulated_plot_path, 'nucleosome_signal_medians_simulation.png'), height = 15, width = 21, units = 'cm', res = 400)
# simulation_medians(seurat_all, QC_metric = "nucleosome_signal", idents = "stage", 
#                    from = min(maximums$max), to = max(minimums$min), by = -0.01, 
#                    ident_cols = stage_colours, cutoff_type = "maximum")
# graphics.off()
# 
# # Cell count simulation
# png(paste0(simulated_plot_path, 'nucleosome_signal_cell_count_simulation.png'), height = 15, width = 21, units = 'cm', res = 400)
# simulation_cell_count(seurat_all, QC_metric = "nucleosome_signal", idents = "stage", 
#                       from = min(maximums$max), to = max(minimums$min), by = -0.01, 
#                       ident_cols = stage_colours, cutoff_type = "maximum")
# graphics.off()

############################## Try low/med/high filtering thresholds to QC #######################################

# # Plot violins for nCount, nFeature and percent.mt at different filtering thresholds
# filter_qc <- lapply(rownames(filter_thresholds), function(condition){
#   seurat_all@meta.data %>%
#     filter(peak_region_fragments > filter_thresholds[condition,'peak_region_fragments_min']) %>%
#     filter(peak_region_fragments < filter_thresholds[condition,'peak_region_fragments_max']) %>%
#     filter(pct_reads_in_peaks > filter_thresholds[condition,'pct_reads_in_peaks']) %>%
#     filter(TSS.enrichment > filter_thresholds[condition,'TSS.enrichment']) %>%
#     filter(nucleosome_signal < filter_thresholds[condition,'nucleosome_signal']) %>%
#     dplyr::select(orig.ident, peak_region_fragments, pct_reads_in_peaks, TSS.enrichment, nucleosome_signal) %>%
#     mutate(filter_condition = !!condition)
# })
# filter_qc <- do.call(rbind, filter_qc) %>%
#   mutate(filter_condition = factor(filter_condition, rownames(filter_thresholds))) %>%
#   reshape2::melt()
# 
# png(paste0(plot_path, 'violins_filter_thresholds.png'), height = 18, width = 40, units = 'cm', res = 400)
# ggplot(filter_qc, aes(x = filter_condition, y = value, fill = orig.ident)) +
#   geom_violin() +
#   facet_wrap(~ variable, nrow = 4, strip.position = "left", scales = "free_y") +
#   theme_minimal() +
#   theme(axis.title.x=element_blank()) +
#   theme(axis.title.y=element_blank()) +
#   theme(strip.placement = "outside")
# graphics.off()
# 
# #### TO ADD: NEW TSS ENRICHMENT PLOT AND NUCLEOSOME BANDING PLOT WITH DIFFERENT THRESHOLDS  #####
# 
# # Calculate remaining cells following different filter thresholds
# filter_qc <- lapply(rownames(filter_thresholds), function(condition){
#   seurat_all@meta.data %>%
#     filter(peak_region_fragments > filter_thresholds[condition,'peak_region_fragments_min']) %>%
#     filter(peak_region_fragments < filter_thresholds[condition,'peak_region_fragments_max']) %>%
#     filter(pct_reads_in_peaks > filter_thresholds[condition,'pct_reads_in_peaks']) %>%
#     filter(TSS.enrichment > filter_thresholds[condition,'TSS.enrichment']) %>%
#     filter(nucleosome_signal < filter_thresholds[condition,'nucleosome_signal']) %>%
#     group_by(orig.ident) %>%
#     tally() %>%
#     dplyr::rename(!!condition := n)
# })
# filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc)
# 
# # Plot remaining cell counts
# filter_qc <-  filter_qc %>% column_to_rownames('orig.ident')
# filter_qc <- rbind(filter_qc, Total = colSums(filter_qc)) %>% rownames_to_column("orig.ident")
# 
# png(paste0(plot_path, 'remaining_cell_table.png'), height = 10, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(filter_qc, rows=NULL, theme = ttheme_minimal()))
# graphics.off()
# 
# png(paste0(plot_path, 'remaining_cell_bar.png'), height = 15, width = 21, units = 'cm', res = 400)
# ggplot(filter_qc[filter_qc$orig.ident != "Total",] %>% reshape2::melt(), aes(x=variable, y=value, fill=orig.ident)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   xlab("Filter Condition") +
#   ylab("Cell Count") +
#   ggtitle("Cell count after filtering") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5))
# graphics.off()
# 
# # Calculate median fragment count per cell following different filter thresholds
# filter_qc <- lapply(rownames(filter_thresholds), function(condition){
#   seurat_all@meta.data %>%
#     filter(pct_reads_in_peaks > filter_thresholds[condition,'pct_reads_in_peaks']) %>%
#     filter(TSS.enrichment > filter_thresholds[condition,'TSS.enrichment']) %>%
#     filter(nucleosome_signal < filter_thresholds[condition,'nucleosome_signal']) %>%
#     filter(peak_region_fragments > filter_thresholds[condition,'peak_region_fragments_min']) %>%
#     filter(peak_region_fragments < filter_thresholds[condition,'peak_region_fragments_max']) %>%
#     group_by(orig.ident) %>%
#     summarise(median = median(passed_filters, na.rm = TRUE)) %>%
#     mutate(median = as.integer(median)) %>%
#     dplyr::rename(!!condition := median)
# })
# 
# # Plot median fragment count per cell
# filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc)
# 
# png(paste0(plot_path, 'median_fragment_count_table.png'), height = 10, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Median Fragment Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(filter_qc, rows=NULL, theme = ttheme_minimal()))
# graphics.off()
# 
# png(paste0(plot_path, 'median_fragment_count_line.png'), height = 15, width = 21, units = 'cm', res = 400)
# ggplot(filter_qc %>% reshape2::melt(), aes(x=variable, y=value, group=orig.ident)) +
#   geom_line(aes(colour = orig.ident)) +
#   xlab("Filter Condition") +
#   ylab("Median Fragment Count") +
#   ggtitle("Median fragment count per cell after filtering") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5))
# graphics.off()

#### Functions
# x axes wrong
# simulation_cell_count <- function(seurat_obj, QC_metric, idents, 
#                                   from, to, by, cutoff_type = c("minimum", "maximum"), 
#                                   ident_cols = NULL){
#   
#   if (cutoff_type == "minimum"){
#     xlab <- "Minimum cut off"
#     filter_cells <- lapply(seq(from = from, to, by = by), function(cutoff){
#       seurat_obj@meta.data %>%
#         filter(!!rlang::sym(QC_metric) > cutoff) %>%
#         group_by(!!rlang::sym(idents)) %>%
#         tally() %>%
#         dplyr::rename(!! paste(cutoff) := n)
#     })} else {
#       xlab <- "Maximum cut off"
#       filter_cells <- lapply(seq(from = from, to, by = by), function(cutoff){
#         seurat_obj@meta.data %>%
#           filter(!!rlang::sym(QC_metric) < cutoff) %>%
#           group_by(!!rlang::sym(idents)) %>%
#           tally() %>%
#           dplyr::rename(!! paste(cutoff) := n)
#       })}
#   filter_cells <- Reduce(function(x, y) merge(x, y), filter_cells) %>% reshape2::melt() %>% mutate(variable = as.integer(variable))
#   
#   if(is.null(ident_cols) == TRUE){
#     ident_cols <- palette(rainbow(length(unique(filter_cells[,1]))))
#   }
#   
#   ggplot(filter_cells, aes(x=variable, y=value, group=!!rlang::sym(idents))) +
#     geom_line(aes(colour = !!rlang::sym(idents))) +
#     xlab(xlab) +
#     ylab("Number of cells") +
#     theme_classic() +
#     scale_color_manual(values = ident_cols)
# }
# 
# ## not working properly
# simulation_medians <- function(seurat_obj, QC_metric, idents, 
#                                from, to, by, cutoff_type = c("minimum", "maximum"), 
#                                ident_cols = NULL){
#   
#   if (cutoff_type == "minimum"){
#     xlab <- "Minimum cut off"
#     filter_cells <- lapply(seq(from = from, to, by = by), function(cutoff){
#       seurat_obj@meta.data %>%
#         filter(!!rlang::sym(QC_metric) > cutoff) %>%
#         group_by(!!rlang::sym(idents)) %>%
#         summarise(median = median(!!rlang::sym(QC_metric), na.rm = TRUE)) %>%
#         mutate(median = as.integer(median)) %>%
#         dplyr::rename(!! paste(cutoff) := median)
#     })} else {
#       xlab <- "Maximum cut off"
#       filter_cells <- lapply(seq(from = from, to, by = by), function(cutoff){
#         seurat_obj@meta.data %>%
#           filter(!!rlang::sym(QC_metric) < cutoff) %>%
#           group_by(!!rlang::sym(idents)) %>%
#           summarise(median = median(!!rlang::sym(QC_metric), na.rm = TRUE)) %>%
#           mutate(median = as.integer(median)) %>%
#           dplyr::rename(!! paste(cutoff) := median)
#       })}
#   filter_cells <- Reduce(function(x, y) merge(x, y), filter_cells) %>% reshape2::melt() %>% mutate(variable = as.integer(variable))
#   
#   if(is.null(ident_cols) == TRUE){
#     ident_cols <- palette(rainbow(length(unique(filter_cells[,1]))))
#   }
#   
#   ggplot(filter_cells, aes(x=variable, y=value, group=!!rlang::sym(idents))) +
#     geom_line(aes(colour = !!rlang::sym(idents))) +
#     xlab(xlab) +
#     ylab("Median QC Metric") +
#     theme_classic() +
#     scale_color_manual(values = ident_cols)
# }

# seurat_all <- subset(seurat_all, subset = pct_reads_in_peaks > filter_thresholds['med','pct_reads_in_peaks'] &
#                        peak_region_fragments > filter_thresholds['med','peak_region_fragments_min'] &
#                        peak_region_fragments < filter_thresholds['med','peak_region_fragments_max'] &
#                        TSS.enrichment > filter_thresholds['med','TSS.enrichment'] &
#                        nucleosome_signal < filter_thresholds['med','nucleosome_signal'])

