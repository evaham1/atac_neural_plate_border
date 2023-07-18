#!/usr/bin/env Rscript

print("coaccessibility ArchR")

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
library(scHelper)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
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


############################## FUNCTIONS #######################################



############################## Read in ArchR project #######################################

# If files are not in rds_files subdirectory look in input dir 
label <- unique(sub('_.*', '', list.files(data_path)))
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



############################## Add motif annotations to ArchR #######################################
# creates a binary matrix where the prescence of a motif in each peak is indicated numerically

# BiocManager::install("JASPAR2020")
# ArchR <- addMotifAnnotations(ArchRProj = ArchR, motifSet = "JASPAR2020", name = "Motif")
    # Error in motif_mats[[x]] : subscript out of bounds
    # In addition: Warning message:
    #   In .get_IDlist_by_query(x, opts) :
    #   Warning: Zero matrices returned with current critera

ArchR <- addMotifAnnotations(ArchR, motifSet = "cisbptest", annoName = "test", species = "Gallusgallus", force = TRUE)

# alex used: JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
library(TFBSTools)
library(JASPAR2020)

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
ArchR <- addMotifAnnotations(ArchR, name = "Motif", motifPWMs = motifList, force = T)


############################## Motifs in differentially accessible peaks #######################################

# se <- getMarkerFeatures(
#   ArchRProj = ArchR, 
#   useMatrix = "PeakMatrix", 
#   groupBy = "clusters")
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

ggUp

# plot as heatmap
heatmapEM <- plotEnrichHeatmap(motif_marker_peaks, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")



############################## Run ChromVar #######################################

# background peaks used to compute motif deviations
ArchR <- addBgdPeaks(ArchR)

# compute per-cell deviations across all motif annotations
ArchR <- addDeviationsMatrix(ArchR, peakAnnotation = "Motif", force = TRUE)

# extract results
df <- getVarDeviations(ArchR, name = "MotifMatrix", plot = FALSE)

# plot results
plotVarDev <- getVarDeviations(ArchR, name = "MotifMatrix", plot = TRUE)
png(paste0(plot_path, 'ChromVar_results.png'), height = 40, width = 40, units = 'cm', res = 400)
plotVarDev
graphics.off()

# save ArchR object
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)
print("ArchR object saved")

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
  
  p <- plotGroups(ArchR, groupBy = "clusters", 
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
  
  # Plot chromvar scores on UMAP
  p <- plotEmbedding(ArchR, colorBy = "GeneIntegrationMatrix", name = TF, embedding = "UMAP", continuousSet = "blueYellow", 
                     imputeWeights = getImputeWeights(ArchR), plotAs = "points", size = 1.8,)
  png(paste0(plot_path, TF, '_gene_integration_UMAP.png'), height = 12, width = 10, units = 'cm', res = 400)
  print(p)
  graphics.off()
  
}


p <- plotGroups(ArchR, groupBy = "clusters", 
                colorBy = "GeneIntegrationMatrix", 
                name = TF)
p


############################## Extract target peaks of key TFs #######################################

# can see where the motif x peak matrix is stored here
anno <- getPeakAnnotation(ArchR, name = "Motif")
anno$Matches

# need to read in this rds file to extract peak annotations with motifs
motifs_peaks <- readRDS("~/output/NF-downstream_analysis/Processing/ss8/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss8_Save-ArchR/Annotations/Motif-Matches-In-Peaks.rds")
dim(motifs_peaks) # 268007 peaks x 746 motifs

# then need to extract the sparse matrix from the ranged experiment object
motif_matrix <- assays(motifs_peaks)$matches
rownames(motif_matrix) <- rowData(motifs_peaks)$name
colnames(motif_matrix) <- colnames(motifs_peaks)

# then can subset the matrix to only include TFs and only peaks with a motif in at least one of these
subsetted_motif_matrix <- motif_matrix[ , c("SIX1", "IRF6", "DLX5", "DLX6", "GATA2", "GATA3", 
                  "TFAP2A", "TFAP2B", "TFAP2C", "PITX1", "PITX2", 
                  "PAX7", "MSX1", "ETS1", "SOX9", "SOX8", "SOX10", "SOX21", "NKX6-2")]
subsetted_motif_matrix <- subsetted_motif_matrix[rowSums(subsetted_motif_matrix[])>0,]

# now want to extract the list of target peak IDs for each TF
TF_targets <- list()
for (TF in colnames(subsetted_motif_matrix)){
  print(TF)
  peak_ids <- rownames(subsetted_motif_matrix)[subsetted_motif_matrix[ , TF]]
  TF_targets[[TF]] <- peak_ids
}

length(TF_targets)

length(TF_targets[["SIX1"]])
tail(TF_targets[["SIX1"]])

for (TF in names(TF_targets)) {
  print(TF)
  write(paste0(TF, "; ", paste0(TF_targets[[TF]], collapse = ", ")), 
          file = paste0(rds_path, "/", "TF_targets.txt"), 
          append = TRUE)
  }