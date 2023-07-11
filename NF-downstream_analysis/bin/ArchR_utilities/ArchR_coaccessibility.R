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
    
    # already clustered
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Clustering/rds_files/"
    # smaller object
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Peak_call//rds_files/"
    
    rds_path = "./output/NF-downstream_analysis/Processing/ss8/Coaccessibility/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/ss8/Coaccessibility/plots/"
    
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


############################## FUNCTIONS #######################################

# extracts the chromosome and TSS coordinate of gene of interest from P2G linkage data attached to ArchR object
ArchR_ExtractTss <- function(ArchR_obj, gene){
  gene <- toupper(gene)
  gene_metadata <- metadata(getPeak2GeneLinks(ArchR_obj, corCutOff = corCutOff, returnLoops = FALSE))$geneSet
  gene_coord <- as.numeric(start(ranges(gene_metadata[which(gene_metadata$name %in% gene)])))
  chr <- as.character(seqnames(gene_metadata[which(gene_metadata$name %in% gene)]))
  return(c(chr, gene_coord))
}


# This function creates a GRanges object that can be used to plot genome browser plots 
# (i.e. one end of the range is Anchor I and the other end of the range is Anchor J)
# It creates this object from a target gene and the peaks it interacts with
# It find the centre coordinate of the TSS and all the peaks, make the ranges and then uses these
# to filter the original P2G linkage data associated with the ArchR object

ArchR_ExtractLoopsToPlot <- function(ArchR_obj, gene, interacting_peaks, corCutOff = 0.5){
  
  # extract TSS coordinate of gene of interest
  chr <- ArchR_ExtractTss(ArchR_obj, gene)[1]
  gene_coord <- ArchR_ExtractTss(ArchR_obj, gene)[2]
  
  # extract interacting peak coordinates from unique peak IDs
  split <- unlist(lapply(interacting_peaks, FUN = function(x) strsplit(x, split = "-")))
  peak_coord <- as.numeric(split[1 : (length(split) / 3) * 3]) - 250
  
  # combine TSS and peak coords into one df
  df <- data.frame(temp_start=gene_coord, temp_end=peak_coord)
  df <- mutate_all(df, function(x) as.numeric(as.character(x)))
  df_ordered <- df %>% 
    mutate(start = apply(df, 1, function(x) min(x))) %>%
    mutate(end = apply(df, 1, function(x) max(x))) %>%
    mutate(chr = chr) %>% 
    dplyr::select(chr, start, end)
  
  # extract p2gl data with user-defined corCutOff - return loops
  all_interactions <- getPeak2GeneLinks(ArchR_obj, corCutOff = corCutOff, returnLoops = TRUE)[[1]]
  
  # filter all interactions from P2GL for these ones
  filtered_granges <- subset(all_interactions, start %in% df_ordered$start & end %in% df_ordered$end)
  
  # check have extracted all interactions
  if (length(filtered_granges) == nrow(df_ordered)){"Check passed"}else{stop("Not all interactions have been extracted!")}
  
  # return interactions
  return(filtered_granges)
}


# This function takes gene and interactions_granges, as well as ArchR object and makes genome browser plot showing interactions
# Can select how to group cells for browser using group_by
# will automatically zoom out so all interactions are seen and centre plot on the gene, can extend further using extend_by

ArchR_PlotInteractions <- function(ArchR_obj, gene, interactions_granges, extend_by = 500, max_dist = Inf, group_by = "clusters", highlight_granges = NULL, return_plot = TRUE){
  
  # extract gene TSS coord
  gene_coord <- as.numeric(ArchR_ExtractTss(ArchR_obj, gene)[2])
  
  # identify the range you need to plot to see all these interactions
  max_coordinate <- as.numeric(max(end(interactions_granges)))
  min_coordinate <- as.numeric(min(start(interactions_granges)))
  
  # set plotting distance limits depending on the size of the loops
  distance <- max(width(interactions_granges)) + extend_by
  if (distance > max_dist){
    print("Distance above max, setting distance from gene as max_dist...")
    distance <- max_dist}
  print(paste0("Distance plotting: ", distance))
  
  # make plot of all the interactions pertaining to one gene around that gene
  if(is.null(highlight_granges)){
    p <- plotBrowserTrack(
      ArchRProj = ArchR_obj,
      groupBy = group_by,
      geneSymbol = gene, 
      upstream = distance,
      downstream = distance,
      loops = interactions_granges,
      title = paste0(gene, "locus - ", length(interactions_granges), " interactions found - distance around: ", distance, "bp")
    )
  } else {
    p <- plotBrowserTrack(
      ArchRProj = ArchR_obj,
      groupBy = group_by,
      geneSymbol = gene, 
      upstream = distance,
      downstream = distance,
      loops = interactions_granges,
      highlight = highlight_granges,
      title = paste0(gene, "locus - ", length(interactions_granges), " interactions found - distance around: ", distance, "bp")
    )
  }
  
  

  
  # optionally return plot
  if (return_plot){
    return(p)
  } else {
    grid::grid.newpage()
    grid::grid.draw(p[[1]])
  }

}


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


###############################################################################################
############################## CO-ACCESSIBILITY BETWEEN PEAKS #################################

# calculate co-accessibility between all peaks
ArchR <- addCoAccessibility(ArchR)

# extract interactions - returns indexes of queryHits and subjectHits
cA <- getCoAccessibility(ArchR, corCutOff = 0.5, returnLoops = FALSE)
cA
  # DataFrame with 120270 rows and 11 columns
  # queryHits subjectHits seqnames correlation Variability1 Variability2     TStat        Pval         FDR VarQuantile1 VarQuantile2
  # <integer>   <integer>    <Rle>   <numeric>    <numeric>    <numeric> <numeric>   <numeric>   <numeric>    <numeric>    <numeric>
  #   1              3           4     chr1    0.548725   0.00437754   0.00683964   14.5441 4.15759e-40 4.52151e-38     0.911185     0.965430
  # 2              4           3     chr1    0.548725   0.00683964   0.00437754   14.5441 4.15759e-40 4.52151e-38     0.965430     0.911185
  # 3              4           5     chr1    0.517190   0.00683964   0.00356568   13.3901 4.49249e-35 3.64027e-33     0.965430     0.870967
  # 4              5           4     chr1    0.517190   0.00356568   0.00683964   13.3901 4.49249e-35 3.64027e-33     0.870967     0.965430
  # 5             27          40     chr1    0.761607   0.01690577   0.00855042   26.0418 1.47916e-94 2.12498e-91     0.995825     0.978303
coacessibility_df <- as.data.frame(cA)

# Need to use indices from df to extract granges and therefore informative peak IDs
coacessibility_df <- coacessibility_df %>% 
  mutate(query_PeakID = paste0(seqnames(metadata(cA)[[1]][queryHits]), "-", start(metadata(cA)[[1]][queryHits]), "-", end(metadata(cA)[[1]][queryHits]))) %>%
  mutate(subject_PeakID = paste0(seqnames(metadata(cA)[[1]][subjectHits]), "-", start(metadata(cA)[[1]][subjectHits]), "-", end(metadata(cA)[[1]][subjectHits])))

head(coacessibility_df)

# sanity check that all interaction Peak IDs are in the ArchR peakset
table(coacessibility_df$subject_PeakID %in% getPeakSet(ArchR)$name)

# save df
write.csv(coacessibility_df, file = paste0(rds_path, "Peak_coaccessibility_df.csv"), row.names = FALSE)

# #### Browser tracks
# p <- plotBrowserTrack(
#   ArchRProj = ArchR,
#   groupBy = "clusters", 
#   geneSymbol = "SIX1", 
#   upstream = 50000,
#   downstream = 50000,
#   loops = getCoAccessibility(ArchR)
# )
# grid::grid.newpage()
# grid::grid.draw(p[[1]])

#########################################################################################################
############################## CO-ACCESSIBILITY BETWEEN PEAKS AND GENES #################################

# calculate gene-to-peak co-accessibility using GeneIntegrationMatrix
#ArchR <- addPeak2GeneLinks(ArchR)
ArchR <- addPeak2GeneLinks(ArchR, maxDist = 90000000000)
# biggest chrom chrom 1 size: 197608386

# extract resulting interactions
p2g <- getPeak2GeneLinks(ArchR, corCutOff = 0.5, returnLoops = FALSE)
p2g_df <- as.data.frame(p2g)
head(p2g)

# need to add correct Peak IDs and gene names to df
p2g_df <- p2g_df %>% 
  mutate(PeakID = paste0(seqnames(metadata(p2g)$peakSet[idxATAC]), "-", start(metadata(p2g)$peakSet[idxATAC]), "-", end(metadata(p2g)$peakSet[idxATAC]))) %>%
  mutate(gene_name = metadata(p2g)$geneSet[idxRNA]$name)

head(p2g_df)

# sanity check that all interaction Peak IDs are in the ArchR peakset
table(p2g_df$PeakID %in% getPeakSet(ArchR)$name)

# save df
write.csv(p2g_df, file = paste0(rds_path, "Peak_to_gene_linkage_df.csv"), row.names = FALSE)

################################################################################
############################## HEATMAPS OF P2L #################################

## Heatmap of linkage across clusters
p <- plotPeak2GeneHeatmap(ArchRProj = ArchR, groupBy = "clusters")
png(paste0(plot_path, 'Peak_to_gene_linkage_clusters_heatmap.png'), height = 80, width = 60, units = 'cm', res = 400)
print(p)
graphics.off()

p <- plotPeak2GeneHeatmap(ArchRProj = ArchR, groupBy = "stage")
png(paste0(plot_path, 'Peak_to_gene_linkage_stage_heatmap.png'), height = 80, width = 60, units = 'cm', res = 400)
print(p)
graphics.off()


###########################################################################################
############################## BROWSER TRACKS P2G LINKAGE #################################

# read in P2G dataframe
p2g_df <- read.csv(paste0(rds_path, "Peak_to_gene_linkage_df.csv"))
head(p2g_df)

# set genes of interest
genes <- c("SIX1", "IRF6", "DLX5", "DLX6", "GATA2", "GATA3", "TFAP2A", "TFAP2B", "TFAP2C", "PITX1", "PITX2",
           "PAX7", "MSX1", "CSRNP1", "ETS1", "SOX9", "SOX8", "SOX10", "SOX5",
           "LMX1B", "ZEB2", "SOX21", "NKX6-2")

for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest
  interactions <- p2g_df %>% filter(gene_name %in% gene)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene
  interacting_peaks <- unique(interactions$PeakID)
  
  # extract loops between the gene and these peaks
  extracted_loops <- ArchR_ExtractLoopsToPlot(ArchR, gene = gene, interacting_peaks = interacting_peaks)
  
  # make plot of these interactions
  p <- ArchR_PlotInteractions(ArchR, gene = gene, interactions_granges = extracted_loops, return_plot = TRUE,
                              extend_by = 500, max_dist = Inf)
  grid::grid.newpage()
  
  # plot
  png(paste0(plot_path, gene, '_interactions_browser_plot.png'), height = 15, width = 18, units = 'cm', res = 400)
  grid::grid.draw(p[[1]])
  graphics.off()
}

# ### Overlay known enhancers
# SIX1_enhancers_df <- data.frame(
#   chr = c("chr5", "chr5", "chr5", "chr5"),
#   start = c(54587474, 54589267, 54590423, 54596822),
#   end = c(54587537, 54589478, 54590562, 54597223)
#   )
# SIX1_enhancers_granges <- makeGRangesFromDataFrame(SIX1_enhancers_df)

# SOX2_enhancers_df <- data.frame(
#   chr = c("chr9", "chr9", "chr9", "chr9", "chr9"),
#   start = c(17013525, 17029120, 17041213, 17003823, 17018130),
#   end = c(17013822, 17029653, 17041793, 17004302, 17018498)
# )
# SOX2_enhancers_granges <- makeGRangesFromDataFrame(SOX2_enhancers_df)

# SOX10_enhancers_df <- data.frame(
#   chr = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
#   start = c(51032535, 51036237, 51046546, 51048966, 51051673, 51065028),
#   end = c(51035291, 51038859, 51049401, 51050564, 51054928, 51068292)
# )
# SOX10_enhancers_granges <- makeGRangesFromDataFrame(SOX10_enhancers_df)

# ETS1_enhancers_df <- data.frame(
#   chr = c("chr24", "chr24", "chr24", "chr24", "chr24", "chr24"),
#   start = c(267288, 385057, 397366, 495221, 1213638, 1717736),
#   end = c(269227, 386795, 399331, 497158, 1215537, 1719643)
# )
# ETS1_enhancers_granges <- makeGRangesFromDataFrame(ETS1_enhancers_df)


# enhancers = ETS1_enhancers_granges
# gene = "ETS1"

# ArchR_ExtractTss(ArchR, gene)

# # plot
# grid::grid.newpage()
# p <- ArchR_PlotInteractions(ArchR, gene = gene, interactions_granges = extracted_loops, return_plot = TRUE,
#                             extend_by = 500, max_dist = Inf, highlight_granges = enhancers)
# png(paste0(plot_path, gene, '_interactions_browser_plot_enhancers.png'), height = 15, width = 18, units = 'cm', res = 400)
# grid::grid.draw(p[[1]])
# graphics.off()
