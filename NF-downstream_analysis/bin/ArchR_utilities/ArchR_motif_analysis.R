#!/usr/bin/env Rscript

print("calculate motifs info and save archr object with this - cell grouping agnostic")

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
library(TFBSTools)
library(JASPAR2020)
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
    
    data_path = "./output/NF-downstream_analysis/Processing/ss8/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/"
    rds_path = "./output/NF-downstream_analysis/Processing/ss8/motif_analysis/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/ss8/motif_analysis/plots/"
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    csv_path = "./csv_files/"
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
  dir.create(csv_path, recursive = T)
}

set.seed(42)


# TFs <- c("SIX1", "IRF6", "DLX5", "DLX6", "GATA2", "GATA3", "TFAP2A", "TFAP2B", "TFAP2C", "PITX1", "PITX2",
#            "PAX7", "MSX1", "ETS1", "SOX9", "SOX8", "SOX10", "SOX5", "SOX21", "NKX6-2",
#            "EPAS1", "RARB", "HRKB1", "IRX1",
#            "TEAD3", "TEAD4", "TEAD2", "TEAD1",
#             "SNAI2", "SNAI1", "SNAI3", "TCF12",
#             "TCF3", "HMBOX1", "FOXK1", "ETV1","REL", "KLF6", "THRB")

TFs <- c("SIX1", "IRF6", "DLX5", "DLX6", "GATA2", "GATA3", "PITX1", "PITX2",
           "TFAP2A", "TFAP2B", "TFAP2C", "TFAP2E",
           "PAX7", "MSX1", "ETS1", "SOX2", "SOX9", "SOX8", "SOX10", "SOX5", "SOX21", "NKX6-2", "NKX2-3",
           "TEAD1", "TEAD2", "TEAD3", "TEAD4",
           "TCF3", "TCF4", "TCF12", "LEF1",
           "FOXK1", "FOXK2", "ZEB1",
           "ETV1","REL", "KLF6", "THRB", "TBXT", "EOMES")

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
getPeakSet(ArchR)
getAvailableMatrices(ArchR)

#############################   MOTIF ANALYSIS    #####################################

############################## Add motif annotations to ArchR #######################################
# creates a binary matrix where the prescence of a motif in each peak is indicated numerically

# download motif database
motifList <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))

# rename each motif to have TF name
name_vector <- c()
for (i in 1:length(motifList)){
  name <- name(motifList[[i]])
  name_vector <- c(name_vector, name)
}
names(motifList) <- name_vector

# annotate peaks in ArchR object with these motifs
ArchR <- addMotifAnnotations(ArchR, name = "Motif", motifPWMs = motifList, cutOff = 1e-05, force = T)
print("Motifs matrix added to ArchR object!")

############################## Run ChromVar #######################################

# plot_path = "./plots/ChromVar/"
# dir.create(plot_path, recursive = T)

# # background peaks used to compute motif deviations
# ArchR <- addBgdPeaks(ArchR)

# # compute per-cell deviations across all motif annotations
# ArchR <- addDeviationsMatrix(ArchR, peakAnnotation = "Motif", force = TRUE)

# # save results as df
# df <- getVarDeviations(ArchR, name = "MotifMatrix", plot = FALSE)
# write.csv(df, file = paste0(csv_path, "Chromvar_df.csv"), row.names = FALSE)

# # plot results
# plotVarDev <- getVarDeviations(ArchR, name = "MotifMatrix", plot = TRUE)
# png(paste0(plot_path, 'ChromVar_results.png'), height = 15, width = 15, units = 'cm', res = 400)
# plotVarDev
# graphics.off()

# print("Chromvar scores calculated!")

# ############################## Save object #######################################
# paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
# saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)
# print("ArchR object saved")

## this object doesnt seem like it can be used for footprinting...

############################## Extract target peaks of key TFs #######################################

# # can see where the motif x peak matrix is stored here
# anno <- getPeakAnnotation(ArchR, name = "Motif")
# anno$Matches

# # need to read in this rds file to extract peak annotations with motifs - dont know if this will work in nextflow...
# #motifs_peaks <- readRDS("~/output/NF-downstream_analysis/Processing/ss8/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss8_Save-ArchR/Annotations/Motif-Matches-In-Peaks.rds")
# motifs_peaks <- readRDS(anno$Matches)
# dim(motifs_peaks) # 268007 peaks x 746 motifs

# # then need to extract the sparse matrix from the ranged experiment object
# motif_matrix <- assays(motifs_peaks)$matches
# rownames(motif_matrix) <- rowData(motifs_peaks)$name
# colnames(motif_matrix) <- colnames(motifs_peaks)

# # then can subset the matrix to only include TFs and only peaks with a motif in at least one of these
# subsetted_motif_matrix <- motif_matrix[ , TFs]
# subsetted_motif_matrix <- subsetted_motif_matrix[rowSums(subsetted_motif_matrix[])>0,]

# # now want to extract the list of target peak IDs for each TF
# TF_targets <- list()
# for (TF in colnames(subsetted_motif_matrix)){
#   print(TF)
#   peak_ids <- rownames(subsetted_motif_matrix)[subsetted_motif_matrix[ , TF]]
#   TF_targets[[TF]] <- peak_ids
# }
# length(TF_targets)

# # now write out the targets into a txt file
# for (TF in names(TF_targets)) {
#   print(TF)
#   write(paste0(TF, "; ", paste0(TF_targets[[TF]], collapse = ", ")), 
#           file = paste0(csv_path, "/", "TF_targets.txt"), 
#           append = TRUE)
# }

# print("TF peak targets written to file!")

# #############################   MOTIF ANALYSIS PLOTS    #####################################

# ############################## Motifs in differentially accessible peaks 

# plot_path = "./plots/diff_accessible_peaks_motifs/"
# dir.create(plot_path, recursive = T)

# # calculate differentially accessible peaks between clusters
# se <- getMarkerFeatures(
#   ArchRProj = ArchR,
#   useMatrix = "PeakMatrix",
#   groupBy = opt$group_by)
# # se <- scHelper::ArchRAddUniqueIdsToSe(se, ArchR, matrix_type = "PeakMatrix")

# motif_marker_peaks <- peakAnnoEnrichment(
#   seMarker = se,
#   ArchRProj = ArchR,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
# )
# motif_marker_peaks

# # make dataframe of results
# df <- data.frame(TF = rownames(motif_marker_peaks), mlog10Padj = assay(motif_marker_peaks)[,1])
# df <- df[order(df$mlog10Padj, decreasing = TRUE),]
# df$rank <- seq_len(nrow(df))

# head(df)

#  # plot results
# ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
#   geom_point(size = 1) +
#   ggrepel::geom_label_repel(
#     data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
#     size = 1.5,
#     nudge_x = 2,
#     color = "black"
#   ) + theme_ArchR() + 
#   ylab("-log10(P-adj) Motif Enrichment") + 
#   xlab("Rank Sorted TFs Enriched") +
#   scale_color_gradientn(colors = paletteContinuous(set = "comet"))

# png(paste0(plot_path, 'diff_peaks_motif_enrichment_plot.png'), height = 10, width = 10, units = 'cm', res = 400)
# print(ggUp)
# graphics.off()

# png(paste0(plot_path, 'diff_peaks_motif_enrichment_plot_bigger.png'), height = 30, width = 30, units = 'cm', res = 400)
# print(ggUp)
# graphics.off()

# # plot as heatmap
# if (label == "ss8"){topn = 10} # 7 was ok so trying higher
# if (label == "ss4"){topn = 10} # 7 was ok so trying higher
# if (label == "HH7"){topn = 10} # 7 was ok so trying higher
# if (label == "HH6"){topn = 10} # 7 was ok so trying higher
# if (label == "HH5"){topn = 1} # 1 was too much!

# if (topn > 1){
#   heatmapEM <- plotEnrichHeatmap(motif_marker_peaks, n = topn, transpose = TRUE)
#   png(paste0(plot_path, 'diff_peaks_motif_enrichment_heatmap.png'), height = 10, width = 18, units = 'cm', res = 400)
#   ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#   graphics.off()
# }

############################## Plot Gene Integration values of TFs of interest #######################################

print("plotting gene integration counts...")

png(paste0(plot_path, 'TEST_gene_integration_UMAP.png'), height = 12, width = 10, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "SIX1",
                plotAs = "points", size = 1.8,
                colorBy = "GeneIntegrationMatrix", continuousSet = "blueYellow")
graphics.off()

plot_path = "./plots/GeneIntegrationCounts/"
dir.create(plot_path, recursive = T)

for (TF in TFs){
  print(TF)
  
  # Plot GeneIntegration values on UMAP
  png(paste0(plot_path, TF, '_gene_integration_UMAP.png'), height = 12, width = 10, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = TF,
                plotAs = "points", size = 1.8,
                colorBy = "GeneIntegrationMatrix", continuousSet = "blueYellow") + 
    theme_ArchR(legendTextSize = 17, baseSize = 17, plotMarginCm = 0.5))
  graphics.off()
  
}


############################## Extract TFs of interest #######################################

# features <- getFeatures(ArchR, useMatrix = "MotifMatrix", select = NULL, ignoreCase = TRUE)
# print(features)

# ArchR <- addImputeWeights(ArchR)

# # # Plot ridge plot of each TF deviation
# for (TF in TFs){
#   print(TF)
#   # markerMotif <- getFeatures(ArchR, select = TF, useMatrix = "MotifMatrix")
#   # if(length(markerMotif) == 0){stop("Motif of that TF not found!")}

#   markerMotif <- paste0("z:", TF)

# #   p <- plotGroups(ArchR, 
# #                   groupBy = opt$group_by, 
# #                   colorBy = "MotifMatrix", 
# #                   name = markerMotif,
# #                   imputeWeights = getImputeWeights(ArchR))

# #   # plot distribution of chromvar deviation score for each cluster
# #   png(paste0(plot_path, TF, '_chromvar_ridge_plot.png'), height = 12, width = 10, units = 'cm', res = 400)
# #   print(p)
# #   graphics.off()
  
#   # Plot chromvar scores on UMAP
#   p <- plotEmbedding(ArchR, colorBy = "MotifMatrix", name = markerMotif, embedding = "UMAP", 
#                      imputeWeights = getImputeWeights(ArchR), plotAs = "points", size = 1.8,)
#   png(paste0(plot_path, TF, '_chromvar_UMAP.png'), height = 12, width = 10, units = 'cm', res = 400)
#   print(p)
#   graphics.off()
  
# }

# print("Chromvar plots made!")

#######################################################################################
#############################   TF FOOTPRINTING    #####################################
#######################################################################################

print("Running footprinting...")

TFs <- c("SIX1", "DLX5", "DLX6", "GATA2", "GATA3",
           "TFAP2A", "TFAP2B", "TFAP2C", "TFAP2E",
           "PAX7", "MSX1", "SOX2", "SOX10",
           "TEAD3",
           "TCF3",
           "FOXK2", "ZEB1")

plot_path = "./plots/Footprinting/"
dir.create(plot_path, recursive = T)

# # set up for footprinting
# motifPositions <- getPositions(ArchR)
# ArchR <- addGroupCoverages(ArchRProj = ArchR, groupBy = opt$group_by)

# # loop through each TF to do footprinting
# for (TF in TFs){
#   print(paste0("Calculating and plotting footprint for: ", TF))
  
#   # compute footprints for TFs of interest
#   seFoot <- getFootprints(ArchR, positions = motifPositions[TF], groupBy = opt$group_by)
#   seFoot
  
#   ####### Plot TF footprints
#   p <- plotFootprints(seFoot, names = TF, normMethod = "Subtract", plotName = "Footprints-Subtract-Bias",
#                       smoothWindow = 10, baseSize = 16, plot = FALSE)
  
#   png(paste0(plot_path, TF, '_TF_footprint.png'), height = 20, width = 20, units = 'cm', res = 400)
#   grid::grid.newpage()
#   grid::grid.draw(p[[1]])
#   graphics.off()
# }

# set up for footprinting
motifPositions <- getPositions(ArchR)
ArchR <- addGroupCoverages(ArchRProj = ArchR, groupBy = "stage")

# loop through each TF to do footprinting
for (TF in TFs){
  print(paste0("Calculating and plotting footprint for: ", TF))
  
  # compute footprints for TFs of interest
  seFoot <- getFootprints(ArchR, positions = motifPositions[TF], groupBy = "stage")
  seFoot
  
  ####### Plot TF footprints
  p <- plotFootprints(seFoot, names = TF, normMethod = "Subtract", plotName = "Footprints-Subtract-Bias",
                      smoothWindow = 10, baseSize = 16, plot = FALSE)
  
  png(paste0(plot_path, TF, '_TF_footprint_stage.png'), height = 20, width = 20, units = 'cm', res = 400)
  grid::grid.newpage()
  grid::grid.draw(p[[1]])
  graphics.off()
}

# set up for footprinting
motifPositions <- getPositions(ArchR)
ArchR <- addGroupCoverages(ArchRProj = ArchR, groupBy = "SEACell_scHelper_cell_type_broad")

# loop through each TF to do footprinting
for (TF in TFs){
  print(paste0("Calculating and plotting footprint for: ", TF))
  
  # compute footprints for TFs of interest
  seFoot <- getFootprints(ArchR, positions = motifPositions[TF], groupBy = "SEACell_scHelper_cell_type_broad")
  seFoot
  
  ####### Plot TF footprints
  p <- plotFootprints(seFoot, names = TF, normMethod = "Subtract", plotName = "Footprints-Subtract-Bias",
                      smoothWindow = 10, baseSize = 16, plot = FALSE)
  
  png(paste0(plot_path, TF, '_TF_footprint_SEACell_scHelper_cell_type_broad.png'), height = 20, width = 20, units = 'cm', res = 400)
  grid::grid.newpage()
  grid::grid.draw(p[[1]])
  graphics.off()
}