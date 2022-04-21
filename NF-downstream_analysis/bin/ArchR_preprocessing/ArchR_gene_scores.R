#!/usr/bin/env Rscript

print("5_ArchR_gene_scores")

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
    
    #plot_path = "./output/NF-downstream_analysis/8_ArchR_gene_scores/plots/"
    #rds_path = "./output/NF-downstream_analysis/8_ArchR_gene_scores/rds_files/"
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ss8/ArchR_clustering/rds_files/"

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
# add a function here to extract top differentially expressed genes per cluster

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
  data <- as.data.frame(t(as.data.frame(assay(extracted_matrix)))) # extract expression matrix
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

############################## Read in ArchR project #######################################
# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

############################## Calculate top gene markers and plot heatmap #################################

markers <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
print("marker genes calculated")

markerList <- getMarkers(markers) # could make more stringent in future
top_markers <- tibble()
print(top_markers)
for (i in 1:length(markerList)){
  table <- as.tibble(markerList[[i]]) 
  print(table)
  table <- table %>% top_n(5, Log2FC) %>% mutate(cluster = i)
  top_markers <- rbind(top_markers, table)
}
if(nrow(top_markers) != 0){
  print("significant markers found")
  
  png(paste0(plot_path, 'top_genes.png'), height = 100, width = 30, units = 'cm', res = 400)
  grid.arrange(tableGrob(top_markers))
  dev.off()
  
  markerGenes <- top_markers$name
  heatmap <- markerHeatmap(
    seMarker = markers, 
    cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
    transpose = TRUE
  )
  png(paste0(plot_path, 'heatmap.png'), height = 30, width = 40, units = 'cm', res = 400)
  ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  graphics.off()
  
} else { print("No markers found that passed thresholds")}


############################## Dot plots and Feature plots of marker genes #################################

# impute weights using MAGIC to plot better feature plots
ArchR <- addImputeWeights(ArchR)

# Contaminating markers
contaminating_markers <- c(
  'DAZL', #PGC
  'CDH5', 'TAL1', 'HBZ', # Blood island
  'CDX2', 'GATA6', 'ALX1', 'PITX2', 'TWIST1', 'TBXT', 'MESP1', #mesoderm
  'SOX17', 'CXCR4', 'FOXA2', 'NKX2-2', 'GATA6' #endoderm
)

png(paste0(plot_path, 'Contaminating_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = contaminating_markers)
graphics.off()

png(paste0(plot_path, 'Contaminating_markers_DotPlots.png'), height = 15, width = 15, units = 'cm', res = 400)
bubble_plot(ArchR, gene_list = contaminating_markers)
graphics.off()

# Late marker genes
late_markers <- c(
  "GATA3", "DLX5", "SIX1", "EYA2", #PPR
  "MSX1", "TFAP2A", "TFAP2B", #mix
  "PAX7", "CSRNP1", "SNAI2", "SOX10", #NC
  "SOX2", "SOX21" # neural
  )

png(paste0(plot_path, 'Late_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = late_markers)
graphics.off()

png(paste0(plot_path, 'Late_markers_DotPlots.png'), height = 15, width = 15, units = 'cm', res = 400)
bubble_plot(ArchR, gene_list = late_markers)
graphics.off()

# look for ap marker genes
ap_markers <- c(
  "PAX2", "WNT4", "SIX3", "SHH" # no GBX2 in matrix
)

png(paste0(plot_path, 'AP_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = ap_markers)
graphics.off()

png(paste0(plot_path, 'AP_markers_DotPlots.png'), height = 15, width = 15, units = 'cm', res = 400)
bubble_plot(ArchR, gene_list = ap_markers)
graphics.off()

# look for early markers
early_markers <- c(
  "EPAS1", "BMP4", "YEATS4", "SOX3", "HOXB1", "ADMP", "EOMES"
)

png(paste0(plot_path, 'Early_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = early_markers)
graphics.off()

png(paste0(plot_path, 'Early_markers_DotPlots.png'), height = 15, width = 15, units = 'cm', res = 400)
bubble_plot(ArchR, gene_list = early_markers)
graphics.off()

############################## Dot plots from RNAseq paper #################################

plot_1_genes <- c("EPAS1", "GATA3", "SIX1", "EYA2",
                  "DLX5", "BMP4", "MSX1", "TFAP2A", "TFAP2B",
                  "PAX3", "PAX7", "SOX2", "OTX2", "YEATS4",
                  "SOX11", "SOX3", "SOX21", "HOXB1", "GBX2",
                  "SIX3", "ADMP", "EOMES")
png(paste0(plot_path, 'DotPlot_1.png'), height = 15, width = 15, units = 'cm', res = 400)
bubble_plot(ArchR, gene_list = plot_1_genes)
graphics.off()

plot_2_genes <- c("GATA3", "DLX5", "SIX1", "EYA2",
                  "MSX1", "TFAP2A", "TFAP2B", "PAX3",
                  "PAX7", "CSRNP1", "SNAI2", "SOX10",
                  "SOX2", "SOX21", "GBX2", "PAX2",
                  "WNT4", "SIX3", "SHH")
png(paste0(plot_path, 'DotPlot_2.png'), height = 15, width = 15, units = 'cm', res = 400)
bubble_plot(ArchR, gene_list = plot_2_genes)
graphics.off()

############################## More informative Feature Plots #################################

genes <- c("SIX1", "PAX7", "DLX5", "CSRNP1", "SOX10",
           "SOX21", "SOX2", "BMP4", "HOXB1")

png(paste0(plot_path, 'Useful_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = genes)
graphics.off()