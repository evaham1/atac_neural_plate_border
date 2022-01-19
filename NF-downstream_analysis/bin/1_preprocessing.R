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
library(GenomeInfoDb)
library(ggplot2)
library(dplyr)
library(rtracklayer)

# add this to dockerfile
#devtools::install_github('alexthiery/scHelper@v0.2.4', dependencies = TRUE, force = TRUE)
#library(scHelper)


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
      plot_path = "../output/NF-downstream_analysis/TEST/plots/"
      rds_path = "../output/NF-downstream_analysis/TEST/rds_files/"
      data_path = "../output/NF-luslab_sc_multiomic/test/cellranger_atac_output/"
      }else{
      plot_path = "../output/NF-downstream_analysis/1_preprocessing/plots/"
      rds_path = "../output/NF-downstream_analysis/1_preprocessing/rds_files/"
      data_path = "../output/NF-luslab_sc_multiomic/full/cellranger_atac_output/"}
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/cellranger_atac_output/"
    ref_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 16gb
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Read in data and set up Signac object #######################################

# Make dataframe with stage and replicate info extracted from path
paths <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)

input <- data.frame(sample = sub('.*/', '', paths), 
                   matrix_path = paste0(paths, "/outs/filtered_peak_bc_matrix.h5"),
                   metadata_path = paste0(paths, "/outs/singlecell.csv"),
                   fragments_path = paste0(paths, "/outs/fragments.tsv.gz"))

# Read in the 3 files needed in list format
counts_list <- apply(input, 1, function(x) Read10X_h5(filename = x[["matrix_path"]]))
metadata_list <- apply(input, 1, function(x) read.csv(file = x[["metadata_path"]], header = TRUE, row.names = 1))
fragments_list <- as.list(input$fragments_path)

# Build list of assays using these files
chrom_assays <- lapply(1:nrow(input), function(x) CreateChromatinAssay(
                                                    counts = counts_list[[x]],
                                                    sep = c(":", "-"),
                                                    fragments = fragments_list[[x]],
                                                    min.cells = 10,
                                                    min.features = 200))
signac_datas <- lapply(1:nrow(input), function(x) CreateSeuratObject(
  counts = chrom_assays[[x]],
  project = input$sample[x],
  assay = "peaks",
  meta.data = metadata_list[[x]]))

# add annotations using chick gtf
gtf <- rtracklayer::import(paste0(ref_path, "genes.gtf.gz"))
gene.coords <- gtf[gtf$type == 'gene']

signac_list <- lapply(signac_datas, function(x) SetAssayData(x, slot = "annotation", new.data = gene.coords))

# Init list of signac objects then merge
names(signac_list) <- input$sample
seurat_all <- merge(x = signac_list[[1]], y=signac_list[-1], add.cell.ids = names(signac_list), project = "chick.10x.atac")

# Add metadata col for stage and flow cell
seurat_all@meta.data[["stage"]] <- substr(seurat_all@meta.data$orig.ident, 1, 3)
seurat_all@meta.data[["flow_cell"]] <- substr(seurat_all@meta.data$orig.ident, 5, 5)

# Convert metadata character cols to factors
seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)] <- lapply(seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)], as.factor)

# to test: save RDS
saveRDS(seurat_all, paste0(rds_path, "seurat_all.RDS"), compress = FALSE)

############################## Try low/med/high filtering thresholds to QC #######################################
seurat_all <- NucleosomeSignal(object = seurat_all)
seurat_all <- TSSEnrichment(object = seurat_all, fast = FALSE)
seurat_all$pct_reads_in_peaks <- seurat_all$peak_region_fragments / seurat_all$passed_filters * 100

# make dataframe with different filtering parameters which can be put into a loop for carrying out downstream analysis
# this doesnt work if strictest threshold has no cells in it? how do I choose these numbers?
filter_thresholds <- data.frame(pct_reads_in_peaks = c(0, 40, 50, 60), TSS.enrichment = c(0, 0.3, 2.5, 3), nucleosome_signal = c(Inf, 1.5, 0.8, 0.6), row.names = c("unfilt", "low", "med", "high"))

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


##!!! Median gene counts --> adapt to median peak counts??
# Calculate median peak count per cell following different filter thresholds
# filter_qc <- lapply(rownames(filter_thresholds), function(condition){
#   seurat_all@meta.data %>%
#     filter(nFeature_peaks > filter_thresholds[condition,'peaks_min']) %>%
#     filter(nFeature_peaks < filter_thresholds[condition,'peaks_max']) %>%
#     filter(percent.mt < filter_thresholds[condition,'MT_max']) %>%
#     group_by(orig.ident) %>%
#     summarise(median = median(nFeature_RNA, na.rm = TRUE)) %>%
#     mutate(median = as.integer(median)) %>%
#     rename(!!condition := median)
# })
# 
# # Plot median gene count per cell
# filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc)
# 
# png(paste0(plot_path, 'median_gene_count_table.png'), height = 10, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Median Gene Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(filter_qc, rows=NULL, theme = ttheme_minimal()))
# graphics.off()
# 
# png(paste0(plot_path, 'median_gene_count_bar.png'), height = 15, width = 21, units = 'cm', res = 400)
# ggplot(filter_qc %>% reshape2::melt(), aes(x=variable, y=value, group=orig.ident)) +
#   geom_line(aes(colour = orig.ident)) +
#   xlab("Filter Condition") +
#   ylab("Median Gene Count") +
#   ggtitle("Median gene count per cell after filtering") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5))
# graphics.off()
# 
# 
# # Plot median gene count on simulated minimum filter thresholds
# tests <- data.frame(gene_cutoff = seq(from = 0, to = 3000, by = 10))
# 
# filter_qc <- lapply(seq(from = 0, to = 3000, by = 10), function(cutoff){
#   seurat_all@meta.data %>%
#     filter(nFeature_RNA > cutoff) %>%
#     group_by(orig.ident) %>%
#     summarise(median = median(nFeature_RNA, na.rm = TRUE)) %>%
#     mutate(median = as.integer(median)) %>%
#     rename(!! paste(cutoff) := median)
# })
# 
# filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc) %>% reshape2::melt() %>% mutate(variable = as.integer(variable)*10)
# 
# png(paste0(plot_path, 'median_gene_count_simulation.png'), height = 15, width = 21, units = 'cm', res = 400)
# ggplot(filter_qc, aes(x=variable, y=value, group=orig.ident)) +
#   geom_line(aes(colour = orig.ident)) +
#   xlab("Lower Gene Threshold") +
#   ylab("Median Gene Count") +
#   ggtitle("Median gene counts at simulated minimum filter thresholds") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5))
# graphics.off()


# Plot violins for nCount, nFeature and percent.mt at different filtering thresholds
filter_qc <- lapply(rownames(filter_thresholds), function(condition){
  seurat_all@meta.data %>%
    filter(nFeature_RNA > filter_thresholds[condition,'gene_min']) %>%
    filter(nFeature_RNA < filter_thresholds[condition,'gene_max']) %>%
    filter(percent.mt < filter_thresholds[condition,'MT_max']) %>%
    dplyr::select(orig.ident, nCount_RNA, nFeature_RNA, percent.mt) %>%
    mutate(filter_condition = !!condition)
})

filter_qc <- lapply(rownames(filter_thresholds), function(condition){
  seurat_all@meta.data %>%
    filter(pct_reads_in_peaks > filter_thresholds[condition,'pct_reads_in_peaks']) %>%
    filter(TSS.enrichment > filter_thresholds[condition,'TSS.enrichment']) %>%
    filter(nucleosome_signal < filter_thresholds[condition,'nucleosome_signal']) %>%
    dplyr::select(orig.ident, pct_reads_in_peaks, TSS.enrichment, nucleosome_signal) %>%
    mutate(filter_condition = !!condition)
})

filter_qc <- do.call(rbind, filter_qc) %>%
  mutate(filter_condition = factor(filter_condition, rownames(filter_thresholds))) %>%
  #mutate(nCount_RNA = ifelse(nCount_RNA >= 100000, 100000, nCount_RNA)) %>% # limit max RNA to 100k for plotting
  reshape2::melt()

png(paste0(plot_path, 'violins_filter_thresholds.png'), height = 18, width = 30, units = 'cm', res = 400)
ggplot(filter_qc, aes(x = filter_condition, y = value, fill = orig.ident)) +
  geom_violin() +
  facet_wrap(~ variable, nrow = 3, strip.position = "left", scales = "free_y") +
  theme_minimal() +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(strip.placement = "outside")
graphics.off()

##########################################################################
############# Remove data which do not pass filter threshold #############
seurat_split <- subset(seurat_all, subset = pct_reads_in_peaks > filter_thresholds['med','pct_reads_in_peaks'] &
                                                TSS.enrichment > filter_thresholds['med','TSS.enrichment'] &
                                                nucleosome_signal < filter_thresholds['med','nucleosome_signal'])
# Normalise
seurat_split <- RunTFIDF(seurat_split)
seurat_split <- FindTopFeatures(seurat_split, min.cutoff = 'q0')
seurat_split <- RunSVD(seurat_split)

# Plot LSI components
png(paste0(plot_path, 'LSI_components.png'), height = 18, width = 30, units = 'cm', res = 400)
DepthCor(seurat_split)
graphics.off()

# Dimensionality reduction and clustering
seurat_split <- RunUMAP(object = seurat_split, reduction = 'lsi', dims = 2:30)
seurat_split <- FindNeighbors(object = seurat_split, reduction = 'lsi', dims = 2:30)
seurat_split <- FindClusters(object = seurat_split, verbose = FALSE, algorithm = 3)

# Plot UMAP
png(paste0(plot_path, 'UMAP.png'), height = 18, width = 18, units = 'cm', res = 400)
DimPlot(object = seurat_split, label = TRUE) + NoLegend()
graphics.off()

####!!! Do I need to pick optimal LSI components like PCA components?

####!!!! Find optimal cluster resolution - need scHelper for this
#png(paste0(plot_path, "clustree_run_", run, ".png"), width=70, height=35, units = 'cm', res = 200)
#ClustRes(seurat_object = seurat_split[[run]], by = 0.2)
#graphics.off()

############################## Identify poor quality clusters #######################################

# Use higher cluster resolution for filtering poor clusters
#seurat_split <- lapply(seurat_split, FindClusters, resolution = 2, verbose = FALSE)

####!!! Plot UMAP for clusters and developmental stage - need scHelper
#png(paste0(plot_path, "UMAP_run_", run, ".png"), width=40, height=20, units = 'cm', res = 200)
#ClustStagePlot(seurat_split)
#graphics.off()

####!!! Need scHelper- Plot QC for each cluster
# png(paste0(plot_path, "QCPlot_run_", run, ".png"), width=32, height=28, units = 'cm', res = 200)
# QCPlot(seurat_split[[run]], quantiles = c(0.25, 0.75))
# graphics.off()

# # Automatically find poor quality clusters
# poor_clusters <- lapply(seurat_split, IdentifyOutliers, metrics = c('nCount_RNA', 'nFeature_RNA'), quantiles = c(0.25, 0.75))
# 
# # Plot UMAP for poor quality clusters
# for(run in names(seurat_split)){
#   png(paste0(plot_path, "PoorClusters_run_", run, ".png"), width=60, height=20, units = 'cm', res = 200)
#   ClusterDimplot(seurat_split[[run]], clusters = poor_clusters[[run]], plot_title = 'poor quality clusters')
#   graphics.off()
# }
# 
# # Filter poor quality clusters
# preprocessing_data <- lapply(names(seurat_split), function(x) {
#   subset(seurat_split[[x]], cells = rownames(filter(seurat_split[[x]]@meta.data, seurat_clusters %in% poor_clusters[[x]])), invert = T)
# })
# names(preprocessing_data) <- names(seurat_split)
# 
# 
# # Plot table with remaining cell counts after full filtering
# cell_counts <- data.frame(unfilt = summary(seurat_all@meta.data$orig.ident),
#                           filtered = lapply(preprocessing_data, function(x) summary(x@meta.data$orig.ident)) %>% do.call(cbind.data.frame, .) %>% rowSums())
# 
# cell_counts <- rbind(cell_counts, Total = colSums(cell_counts)) %>% rownames_to_column("orig.ident")
# 
# png(paste0(plot_path, 'final_remaining_cell_table.png'), height = 10, width = 10, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
# graphics.off()
# 
# # Save RDS output
# saveRDS(preprocessing_data, paste0(rds_path, "preprocessing_data.RDS"), compress = FALSE)