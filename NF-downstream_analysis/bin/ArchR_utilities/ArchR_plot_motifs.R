#!/usr/bin/env Rscript

print("plot motifs and TF footprinting - need to split cells by a grouping")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(parallel)
# library(scHelper)
library(BSgenome.Ggallus.UCSC.galGal6)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-g", "--group_by"), action = "store", type = "character", help = "How to group cells", default = "clusters",),
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
    
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Metacell_to_singlecell/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/ss8/motif_analysis/plots/"
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    #addArchRThreads(threads = ncores)
    addArchRThreads(threads = 1) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

set.seed(42)

############################## Read in ArchR project #######################################

# Extract stage name by removing anything with the word 'peak' in it
labels <- unique(sub('_.*', '', list.files(data_path)))
print(labels)
label <- labels[labels != "Peak"]
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

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
# getPeakSet(ArchR)
getAvailableMatrices(ArchR)


#############################   MOTIF ANALYSIS PLOTS    #####################################

############################## Motifs in differentially accessible peaks 

plot_path = "./plots/diff_accessible_peaks_motifs/"
dir.create(plot_path, recursive = T)

# calculate differentially accessible peaks between clusters
se <- getMarkerFeatures(
  ArchRProj = ArchR,
  useMatrix = "PeakMatrix",
  groupBy = opt$group_by)
# se <- scHelper::ArchRAddUniqueIdsToSe(se, ArchR, matrix_type = "PeakMatrix")

motif_marker_peaks <- peakAnnoEnrichment(
  seMarker = se,
  ArchRProj = ArchR,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
motif_marker_peaks

# make dataframe of results
df <- data.frame(TF = rownames(motif_marker_peaks), mlog10Padj = assay(motif_marker_peaks)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

head(df)

 # plot results
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

png(paste0(plot_path, 'diff_peaks_motif_enrichment_plot.png'), height = 10, width = 10, units = 'cm', res = 400)
print(ggUp)
graphics.off()

png(paste0(plot_path, 'diff_peaks_motif_enrichment_plot_bigger.png'), height = 30, width = 30, units = 'cm', res = 400)
print(ggUp)
graphics.off()

# plot as heatmap
if (label == "ss8"){topn = 10} # 7 was ok so trying higher
if (label == "ss4"){topn = 10} # 7 was ok so trying higher
if (label == "HH7"){topn = 10} # 7 was ok so trying higher
if (label == "HH6"){topn = 10} # 7 was ok so trying higher
if (label == "HH5"){topn = 1} # 1 was too much!

if (topn > 1){
  heatmapEM <- plotEnrichHeatmap(motif_marker_peaks, n = topn, transpose = TRUE)
  png(paste0(plot_path, 'diff_peaks_motif_enrichment_heatmap.png'), height = 10, width = 18, units = 'cm', res = 400)
  ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  graphics.off()
}

############################## Extract TFs of interest #######################################

# set genes of interest
TFs <- c("SIX1", "IRF6", "DLX5", "DLX6", "GATA2", "GATA3", "TFAP2A", "TFAP2B", "TFAP2C", "PITX1", "PITX2",
           "PAX7", "MSX1", "ETS1", "SOX9", "SOX8", "SOX10", "SOX5", "SOX21", "NKX6-2")
# CTNRP and LMX1B and ZEB2 not found

ArchR <- addImputeWeights(ArchR)

# Plot ridge plot of each TF deviation
for (TF in TFs){
  print(TF)
  markerMotif <- getFeatures(ArchR, select = TF, useMatrix = "MotifMatrix")
  if(length(markerMotif) == 0){stop("Motif of that TF not found!")}
  
  p <- plotGroups(ArchR, 
                  groupBy = opt$group_by, 
                  colorBy = "MotifMatrix", 
                  name = markerMotif,
                  imputeWeights = getImputeWeights(ArchR))

  # plot distribution of chromvar deviation score for each cluster
  png(paste0(plot_path, TF, '_chromvar_ridge_plot.png'), height = 12, width = 10, units = 'cm', res = 400)
  print(p)
  graphics.off()
  
  # Plot chromvar scores on UMAP
  p <- plotEmbedding(ArchR, colorBy = "MotifMatrix", name = markerMotif, embedding = "UMAP", 
                     imputeWeights = getImputeWeights(ArchR), plotAs = "points", size = 1.8,)
  png(paste0(plot_path, TF, '_chromvar_UMAP.png'), height = 12, width = 10, units = 'cm', res = 400)
  print(p)
  graphics.off()
  
}

print("Chromvar plots made!")

#######################################################################################
#############################   TF FOOTPRINTING    #####################################
#######################################################################################

plot_path = "./plots/Footprinting/"
dir.create(plot_path, recursive = T)

# slightly smaller TFs list:
TFs <- c("SIX1", "IRF6", "DLX5", "DLX6", "GATA2", "GATA3", 
        "TFAP2A", "TFAP2B", "TFAP2C", "PITX1", "PITX2", 
        "PAX7", "MSX1", "ETS1", "SOX9", "SOX8", "SOX10", "SOX21", "NKX6-2")

# extract all motif positions
motifPositions <- getPositions(ArchR)

# compute pseudoreplicates
ArchR <- addGroupCoverages(ArchR, groupBy = opt$group_by)

# loop through each TF to do footprinting
for (TF in TFs){
  print(paste0("Calculating and plotting footprint for: ", TF))
  
  # extract motif positions for that TF
  positions <- motifPositions[TF]
  print(paste0(length(positions[[1]]), " motif positions found"))
  
  # compute footprints for TFs of interest
  seFoot <- getFootprints(ArchR, positions = positions, groupBy = opt$group_by)
  seFoot
  
  ####### Plot TF footprints
  p <- plotFootprints(seFoot, names = TF, normMethod = "Subtract", plotName = "Footprints-Subtract-Bias",
                      smoothWindow = 10, baseSize = 16, plot = FALSE)
  
  png(paste0(plot_path, TF, '_TF_footprint.png'), height = 20, width = 20, units = 'cm', res = 400)
  grid::grid.newpage()
  grid::grid.draw(p[[1]])
  graphics.off()
}