#!/usr/bin/env Rscript

print("ArchR_gene_scores")

############################## Load libraries #######################################
library(getopt)
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
library(presto)

############################## Set up script options #######################################
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    
    plot_path = "./output/NF-downstream_analysis/ArchR_preprocessing/QC_MED/ss8/postfiltering/gene_scores/plots/"
    rds_path = "./output/NF-downstream_analysis/ArchR_preprocessing/QC_MED/ss8/postfiltering/gene_scores/rds_files/"
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/QC_MED/ss8/postfiltering/clustering/rds_files/"
    
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/QC_MED/HH6/postfiltering/clustering/rds_files/"
    
    #data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ss8/1_ArchR_clustering_prefiltering/rds_files/"
    #plot_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ss8/1.5_ArchR_gene_scores_unfiltered/plots/"

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

# Feature plot function to create grid of feature plots
feature_plot_grid <- function(ArchRProj = ArchR, matrix = "GeneScoreMatrix", gene_list) {
  p <- plotEmbedding(ArchRProj, colorBy = matrix, name = gene_list, 
                    plotAs = "points", size = 1.8, baseSize = 0, labelSize = 8, legendSize = 10)
  p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
      )
  })
  do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
}

## Bubble plot function adapted from https://github.com/NoemieL/bubble-plot-ArchR/blob/main/Script
bubble_plot <- function(ArchRProj = ArchR, matrix = "GeneScoreMatrix", gene_list) {
  
  extracted_matrix <- getMatrixFromProject(ArchRProj, useMatrix = matrix)
  data <- as.data.frame(as.matrix(t(assay(extracted_matrix)))) # extract expression matrix
  colnames(data) <- rowData(extracted_matrix)[,5] # add gene names
  for(i in rownames(data)){ data[i,"clusters"] = ArchRProj$clusters[which(ArchRProj$cellNames==i)] } # add cluster IDs
  data$clusters <- as.numeric(gsub('^.', '', data$clusters)) # remove the Cs so clustered are ordered
  data$clusters  = factor(data$clusters, levels = c(1 : max(data$clusters))) # order clusters

  bubble_plot_info = data.frame()

  for(i in gene_list) {
  
    for(k in 1:length(unique(data$clusters))) {
      a = nrow(bubble_plot_info)
      l = unique(data$clusters)[k]
      bubble_plot_info[a+1,"gene_name"] = i
      bubble_plot_info[a+1,"clusters"] = l
      eval(parse(text=(paste("bubble_plot_info[",a,"+1,'pct_exp'] = (length(data[(data$",i,">0 & data$clusters == '",l,"'),'",i,"'])/nrow(data[data$clusters=='",l,"',]))*100", sep=""))))
      eval(parse(text=(paste("bubble_plot_info[",a,"+1,'avg_exp'] = mean(data[data$",i,">0 & data$clusters == '",l,"','",i,"'])", sep=""))))
    }
  }

  ggplot(data = bubble_plot_info, mapping = aes_string(x = 'gene_name', y = 'clusters')) +
    geom_point(mapping = aes_string(size = 'pct_exp', color = "avg_exp")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = 'Percent Expressed')) +
    labs(
      x = 'gene_name',
      y = 'Clusters'
    )+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 0))
}

### function to extract top n logFC for each group from summarized experiment object
# returns a summarised experiment object with only the top n features per cell group
extract_top_features <- function(markers, n = 10) {
  markerList <- getMarkers(markers, cutOff = "FDR <= 1")
  df <- data.frame()
  for (i in 1:length(names(markerList))) {
    print(i)
    df_i <- as.data.frame(markerList[i])
    df <- rbind(df, df_i)
  }
  df <- df %>%
    group_by(group_name) %>%
    top_n(n, Log2FC) %>%
    dplyr::arrange(Log2FC, .by_group = TRUE)
  top_markers <- df$idx # top 15 markers by logF2C between cell groups
  top_markers <- unique(top_markers) # some of these markers are shared??
  coords <-  rownames(markers)[rownames(markers) %in% top_markers]
  top_markers_se <- markers[coords, ]
  return(top_markers_se)
}

############################## Read in ArchR project #######################################
# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

############################## Read in marker genes #################################

# Contaminating markers
contaminating_markers <- c(
  'DAZL', #PGC
  'CDH5', 'TAL1', 'HBZ', # Blood island
  'CDX2', 'GATA6', 'ALX1', 'PITX2', 'TWIST1', 'TBXT', 'MESP1', #mesoderm
  'SOX17', 'CXCR4', 'FOXA2', 'NKX2-2', 'GATA6' #endoderm
)

# Late marker genes
late_markers <- c(
  "GATA3", "DLX5", "SIX1", "EYA2", #PPR
  "MSX1", "TFAP2A", "TFAP2B", #mix
  "PAX7", "CSRNP1", "SNAI2", "SOX10", #NC
  "SOX2", "SOX21" # neural
)

# look for ap marker genes
ap_markers <- c(
  "PAX2", "WNT4", "SIX3", "SHH" # no GBX2 in matrix
)

# look for early markers
early_markers <- c(
  "EPAS1", "BMP4", "YEATS4", "SOX3", "HOXB1", "ADMP", "EOMES"
)

dotplot_1_genes <- c("EPAS1", "GATA3", "SIX1", "EYA2",
                  "DLX5", "BMP4", "MSX1", "TFAP2A", "TFAP2B",
                  "PAX3", "PAX7", "SOX2", "OTX2", "YEATS4",
                  "SOX11", "SOX3", "SOX21", "HOXB1", "GBX2",
                  "SIX3", "ADMP", "EOMES")

dotplot_2_genes <- c("GATA3", "DLX5", "SIX1", "EYA2",
                  "MSX1", "TFAP2A", "TFAP2B", "PAX3",
                  "PAX7", "CSRNP1", "SNAI2", "SOX10",
                  "SOX2", "SOX21", "GBX2", "PAX2",
                  "WNT4", "SIX3", "SHH")

feature_plot_genes <- c("SIX1", "PAX7", "DLX5", "CSRNP1", "SOX10",
           "SOX21", "SOX2", "BMP4", "HOXB1")

all_genes <- unique(c(contaminating_markers, late_markers, early_markers, ap_markers, dotplot_1_genes, dotplot_2_genes))

############################## Dot Plots #################################

# addArchRThreads(threads = 1) 

# png(paste0(plot_path, 'Contaminating_markers_DotPlots.png'), height = 15, width = 15, units = 'cm', res = 400)
# bubble_plot(ArchR, gene_list = contaminating_markers)
# graphics.off()

# png(paste0(plot_path, 'Late_markers_DotPlots.png'), height = 15, width = 15, units = 'cm', res = 400)
# bubble_plot(ArchR, gene_list = late_markers)
# graphics.off()

# png(paste0(plot_path, 'AP_markers_DotPlots.png'), height = 15, width = 15, units = 'cm', res = 400)
# bubble_plot(ArchR, gene_list = ap_markers)
# graphics.off()

# png(paste0(plot_path, 'Early_markers_DotPlots.png'), height = 15, width = 15, units = 'cm', res = 400)
# bubble_plot(ArchR, gene_list = early_markers)
# graphics.off()

# png(paste0(plot_path, 'DotPlot_1.png'), height = 15, width = 15, units = 'cm', res = 400)
# bubble_plot(ArchR, gene_list = dotplot_1_genes)
# graphics.off()

# png(paste0(plot_path, 'DotPlot_2.png'), height = 15, width = 15, units = 'cm', res = 400)
# bubble_plot(ArchR, gene_list = dotplot_2_genes)
# graphics.off()

# print("Dotplots done")

############################## Feature Plots #################################

# impute weights using MAGIC to plot better feature plots
ArchR <- addImputeWeights(ArchR)

png(paste0(plot_path, 'Contaminating_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = contaminating_markers)
graphics.off()

png(paste0(plot_path, 'Late_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = late_markers)
graphics.off()

png(paste0(plot_path, 'AP_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = ap_markers)
graphics.off()

png(paste0(plot_path, 'Early_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = early_markers)
graphics.off()

png(paste0(plot_path, 'Useful_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = feature_plot_genes)
graphics.off()

print("Feature plots done")

############################## Calculate top gene markers and plot heatmap #################################

markers <- getMarkerFeatures(
  ArchRProj = ArchR,
  useMatrix = "GeneScoreMatrix",
  groupBy = "clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
print("marker genes calculated")

## heatmap of all calculated markers that pass thresholds
heatmap <- markerHeatmap(
  seMarker = markers,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  nLabel = 3
)
png(paste0(plot_path, 'cutoff_heatmap.png'), height = 30, width = 40, units = 'cm', res = 400)
draw(heatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
graphics.off()

## heatmap of top 50 differentially expressed genes per cluster
# adapted from: https://github.com/maehrlab/pharyngeal_endoderm_development/blob/main/scATAC_qc_analysis/scATAC_analysis/scATAC_7_gene_score_and_peaks.ipynb

marker_genes <- markers
handpicked_genes <- all_genes

add_name = function(X, c) {
  if(nrow(X)==0) return(NULL)
  X$clusters = c
  X
}
marker_tables = marker_genes %>% getMarkers(cutOff = "FDR <= 1 & Log2FC >= 0")
marker_tables = mapply(add_name, marker_tables, names(marker_tables), SIMPLIFY = F) 
marker_tables_clean <- c()
for (i in 1:length(marker_tables)){
  print(i)
  if(!is.null(marker_tables[[i]])){
    add <- marker_tables[[i]]
    marker_tables_clean <- c(marker_tables_clean, add)
  } else { print("this cluster has no markers!") }
}
marker_tables <- marker_tables_clean %>% Reduce(f = rbind)

mixedrank = function(x) order(gtools::mixedorder(x))
markers_top_table_S2 = marker_tables %>%  
  as.data.frame() %>%
  dplyr::group_by(clusters) %>%
  dplyr::mutate(rank = rank(-Log2FC, ties = "first")) %>%
  dplyr::top_n(Log2FC, n = 100) %>%
  dplyr::arrange(mixedrank(clusters), desc(Log2FC))

# subsetting markers to only include top 50 per cluster
markers_include = markers_top_table_S2 %>% subset(rank<=50, select = "name", drop = T)
marker_subset = marker_genes[rowData(marker_genes)$name %in% markers_include,]

# genes to label
markers_label = intersect(handpicked_genes, markers_include)
if (length(markers_label) == 0){markers_label <- NULL}

# plot heatmap
heatmapGS <- plotMarkerHeatmap(
  seMarker = marker_subset, 
  cutOff = "FDR <= 1 & Log2FC >= 0",
  labelMarkers = markers_label,
  clusterCols = TRUE,
  pal = viridis::magma(30),
  nLabel = 3
)
png(paste0(plot_path, 'top50_heatmap.png'), height = 30, width = 30, units = 'cm', res = 400)
draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
graphics.off()
