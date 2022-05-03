#!/usr/bin/env Rscript

print("differential peaks ArchR")
# calculates differential peaks between clusters, scHelper_cell_types_old and stage -> plots on heatmaps

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
    
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/QC_HIGH/ss8/prefiltering/peak_calling/rds_files/"
    #plot_path = "./output/NF-downstream_analysis/ArchR_peak_calling/ss8/2_differential_peaks/plots/"
    
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

# If files are not in rds_files subdirectory look in input dir
label <- sub('_.*', '', list.files(data_path))
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

getPeakSet(ArchR)
getAvailableMatrices(ArchR)

# ############################## Subset to speed up interactive work #######################################
# ArchR_big <- ArchR
# 
# idxSample <- BiocGenerics::which(ArchR_big$scHelper_cell_type_old %in% c("NC", "dNC", "aPPR", "pPPR"))
# cellsSample <- ArchR_big$cellNames[idxSample]
# 
# ArchR <- subsetArchRProject(ArchR_big, cellsSample, outputDirectory = "ArchRSubset",
#   dropCells = TRUE, logFile = NULL, force = TRUE)
# 
# ArchR <- loadArchRProject(path = "./ArchRSubset", force = TRUE, showLogo = TRUE)


############################# CLUSTER Marker Peaks #######################################

markersPeaks <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "clusters")
markersPeaks

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.3 & Log2FC >= 0.5",
  nLabel = 3)

png(paste0(plot_path, 'clusters_diff_peak_cutoff_heatmap.png'), height = 50, width = 40, units = 'cm', res = 400)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
graphics.off()

# https://github.com/maehrlab/pharyngeal_endoderm_development/blob/main/scATAC_qc_analysis/scATAC_analysis/scATAC_10_Fig3a_3b.ipynb
add_name = function(X, c) {
  if(nrow(X)==0) return(NULL)
  X$cluster = c
  X
}
## All peaks:
marker_tables = markersPeaks %>% getMarkers(cutOff = "FDR <= 1 & Log2FC >= 0.1")
marker_tables_S = mapply(add_name, marker_tables, names(marker_tables), SIMPLIFY = F) %>% Reduce(f = rbind)
mixedrank = function(x) order(gtools::mixedorder(x))
# markers_top_S = marker_tables_S %>%  
#   as.data.frame() %>%
#   #dplyr::mutate( highest_logfc_by_gene = ave(Log2FC, name, FUN = max)) %>%
#   #dplyr::mutate( is_highest = Log2FC == highest_logfc_by_gene) %>%
#   #subset(is_highest) %>%
#   dplyr::group_by(cluster) %>%
#   dplyr::mutate(rank = rank(-Log2FC, ties = "first")) %>%
#   # dplyr::top_n(Log2FC, n = 1000000) %>%
#   dplyr::arrange(mixedrank(cluster), desc(Log2FC))
# #write.csv(markers_top_S,"all_differentially_expressed_peaks.csv", row.names = FALSE)

# Plot top peaks:
markers_top = marker_tables_S %>%  
  as.data.frame() %>%
  dplyr::mutate( highest_logfc_by_gene = ave(Log2FC, name, FUN = max)) %>%
  dplyr::mutate( is_highest = Log2FC == highest_logfc_by_gene) %>%
  subset(is_highest) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(rank = rank(-Log2FC, ties = "first")) %>%
  dplyr::top_n(Log2FC, n = 50) %>%
  dplyr::arrange(mixedrank(cluster), desc(Log2FC))
markers_label   = markers_top %>% subset(rank<=2,  select = "name", drop = T) 
markers_include = markers_top %>% subset(rank<=4000000, select = "unique_id", drop = T) 
#length(markers_include)
marker_subset = markersPeaks[rowData(markersPeaks)$unique_id %in% markers_include,]

png(paste0(plot_path, 'clusters_diff_peak_top50_heatmap.png'), height = 50, width = 40, units = 'cm', res = 400)
heatmapGS <- plotMarkerHeatmap(
  seMarker = marker_subset, 
  cutOff = "FDR <= 1 & Log2FC >= 0.1",
  labelMarkers = labelMarkers,
  pal = viridisLite::mako(256)
)
graphics.off()

# ######################
# # subset for proximal/distal/etc
# 
# tmp_peaks = data.frame(ArchR@peakSet)
# tmp_diff_peaks = data.frame(rowData(markersPeaks))
# 
# diff_peaks_join_peakset = left_join(tmp_diff_peaks, tmp_peaks, 
#                                     by = c("seqnames" = "seqnames", "start" = "start", "end" = "end"))
# diff_peaks_join_peakset$name = paste(diff_peaks_join_peakset$nearestGene, diff_peaks_join_peakset$distToTSS,sep="_")
# diff_peaks_join_peakset$unique_id = paste(diff_peaks_join_peakset$seqnames, diff_peaks_join_peakset$start, diff_peaks_join_peakset$end, sep=":")
# rowData(markersPeaks) = diff_peaks_join_peakset
# 
# # extract labels of peaks closest to handpicked genes:
# handpicked_genes  <- c(
#   "GATA3", "DLX5", "SIX1", "EYA2", #PPR
#   "MSX1", "TFAP2A", "TFAP2B", #mix
#   "PAX7", "CSRNP1", "SNAI2", "SOX10", #NC
#   "SOX2", "SOX21" # neural
# )
# markers_label_handpicked = diff_peaks_join_peakset[as.vector(diff_peaks_join_peakset$nearestGene %in% handpicked_genes),]$name
# ## Promoter peaks::
# markersPeaks_promoter <- subset(markersPeaks, rowData(markersPeaks)$peakType == "Promoter")
# marker_tables_promoter = markersPeaks_promoter %>% getMarkers(cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
# #marker_tables_promoter_tmp = marker_tables_promoter %>%  as.data.frame() 
# #write.csv(marker_tables_promoter_tmp,"all_differentially_expressed_peaks_promoter.csv", row.names = FALSE)
# 
# ## Distal peaks:
# markersPeaks_distal <- subset(markersPeaks, rowData(markersPeaks)$peakType == "Distal")
# marker_tables_distal = markersPeaks_distal %>% getMarkers(cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
# #marker_tables_distal_tmp = marker_tables_distal %>%  as.data.frame() 
# #write.csv(marker_tables_distal_tmp,"all_differentially_expressed_peaks_distal.csv", row.names = FALSE)

# ## Plot top distal peaks:
# marker_tables_distal = markersPeaks_distal %>% getMarkers(cutOff = "FDR <= 0.01 & Log2FC >= 0.5")
# marker_tables_distal = mapply(add_name, marker_tables_distal, names(marker_tables_distal), SIMPLIFY = F) 
# marker_tables_clean <- c()
# for (i in 1:length(marker_tables_distal)){
#   print(i)
#   if(!is.null(marker_tables_distal[[i]])){
#     add <- marker_tables_distal[[i]]
#     marker_tables_clean <- c(marker_tables_clean, add)
#   } else { print("this cluster has no markers!") }
# }
# marker_tables_distal <- marker_tables_clean %>% Reduce(f = rbind)
# mixedrank = function(x) order(gtools::mixedorder(x))
# markers_top_distal = marker_tables_distal %>%  
#   as.data.frame() %>%
#   dplyr::mutate( highest_logfc_by_gene = ave(Log2FC, name, FUN = max)) %>%
#   dplyr::mutate( is_highest = Log2FC == highest_logfc_by_gene) %>%
#   subset(is_highest) %>%
#   dplyr::group_by(cluster) %>%
#   dplyr::mutate(rank = rank(-Log2FC, ties = "first")) %>%
#   #  dplyr::top_n(Log2FC, n = 100) %>%
#   dplyr::arrange(mixedrank(cluster), desc(Log2FC))
# #dim(marker_tables_distal)
# markers_label   = markers_top_distal %>% subset(rank<=2,  select = "name", drop = T) 
# markers_include = markers_top_distal %>% subset(rank<=4000000, select = "unique_id", drop = T) 
# #length(markers_include)
# marker_subset = markersPeaks_distal[rowData(markersPeaks_distal)$unique_id %in% markers_include,]
# #dim(marker_subset)
# 
# labelMarkers = markers_label_handpicked %>% union(markers_label)
# #length(labelMarkers)
# 
# heatmapGS <- plotMarkerHeatmap(
#   seMarker = marker_subset, 
#   cutOff = "FDR <= 1 & Log2FC >= 0.1",
#   #labelMarkers = labelMarkers,
#   pal = viridisLite::mako(256),
#   returnMatrix = TRUE
# )
# 
# heatmap(heatmapGS)
# pheatmap(heatmapGS, show_rownames = FALSE)

############################# SC_HELPER_CELL_TYPE Marker Peaks #######################################

if (is.null(ArchR$scHelper_cell_type_old) == FALSE) {
  ArchR$scHelper_cell_type_old <- as.character(ArchR$scHelper_cell_type_old)
  
  markersPeaks <- getMarkerFeatures(
    ArchRProj = ArchR, 
    useMatrix = "PeakMatrix", 
    groupBy = "scHelper_cell_type_old")
  markersPeaks
  
  heatmapPeaks <- plotMarkerHeatmap(
    seMarker = markersPeaks, 
    cutOff = "FDR <= 0.3 & Log2FC >= 0.5",
    nLabel = 3)
  
  png(paste0(plot_path, 'scHelper_cell_type_diff_peak_cutoff_heatmap.png'), height = 50, width = 40, units = 'cm', res = 400)
  print(draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot"))
  graphics.off()
}

############################# STAGE Marker Peaks #######################################

if (length(unique(ArchR$stage)) > 1) {
  
  markersPeaks <- getMarkerFeatures(
    ArchRProj = ArchR, 
    useMatrix = "PeakMatrix", 
    groupBy = "stage")
  markersPeaks
  
  heatmapPeaks <- plotMarkerHeatmap(
    seMarker = markersPeaks, 
    cutOff = "FDR <= 0.3 & Log2FC >= 0.5",
    nLabel = 3)
  
  png(paste0(plot_path, 'stage_diff_peak_cutoff_heatmap.png'), height = 50, width = 40, units = 'cm', res = 400)
  print(draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot"))
  graphics.off()
  
}

# ############################## Individual group plots #######################################

# ### make for loop so make MA/volcano plot for each group in group_by
# for (i in )
# pma <- plotMarkers(seMarker = markersPeaks, name = "aPPR", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
# pma

# pv <- plotMarkers(seMarker = markersPeaks, name = "C6", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
# pv


# ############################## Browser tracks #######################################

# p <- plotBrowserTrack(
#   ArchRProj = ArchR, 
#   groupBy = "clusters", 
#   geneSymbol = c("SIX1"),
#   features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["C6"],
#   upstream = 50000,
#   downstream = 50000
# )

# grid::grid.newpage()
# grid::grid.draw(p$SIX1)