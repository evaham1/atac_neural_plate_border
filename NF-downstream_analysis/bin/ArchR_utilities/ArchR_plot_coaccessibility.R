#!/usr/bin/env Rscript

print("plot coaccessibility ArchR - needs cell grouping")

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
    
    # ArchR object with no coaccessibility data
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Transfer_peaks/"
    # co-accessibility csv and RDS files
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Single_cell_integration/csv_files/"
    # # ArchR object with coaccessibility data
    # data_path = "./output/NF-downstream_analysis/Processing/FullData/Single_cell_integration/"
    
    rds_path = "./output/NF-downstream_analysis/Processing/ss8/Coaccessibility/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/ss8/Coaccessibility/plots/"
    
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

# extracts the chromosome and TSS coordinate of gene of interest from P2G linkage data attached to ArchR object
ArchR_ExtractTss <- function(gene_metadata, gene){
  gene <- toupper(gene)
  gene_coord <- as.numeric(start(ranges(gene_metadata[which(gene_metadata$name %in% gene)])))
  chr <- as.character(seqnames(gene_metadata[which(gene_metadata$name %in% gene)]))
  return(c(chr, gene_coord))
}


# This function creates a GRanges object that can be used to plot genome browser plots 
# (i.e. one end of the range is Anchor I and the other end of the range is Anchor J)
# It creates this object from a target gene and the peaks it interacts with
# It find the centre coordinate of the TSS and all the peaks, make the ranges and then uses these
# to filter the original P2G linkage data associated with the ArchR object

ArchR_ExtractLoopsToPlot <- function(ArchR_obj, gene, gene_locations, interacting_peaks, interactions_granges, corCutOff = 0.5){
  
  # extract TSS coordinate of gene of interest
  chr <- ArchR_ExtractTss(gene_locations, gene)[1]
  gene_coord <- ArchR_ExtractTss(gene_locations, gene)[2]
  
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
  
  # subset granges of interactions by corCutOff
  all_interactions <- interactions_granges[(elementMetadata(interactions_granges)[,1] > corCutOff)]
  
  # filter all interactions from P2G for these ones
  filtered_granges <- subset(all_interactions, start %in% df_ordered$start & end %in% df_ordered$end)
  
  # check have extracted all interactions
  if (length(filtered_granges) == nrow(df_ordered)){"Check passed"}else{stop("Not all interactions have been extracted!")}
  
  # return interactions
  return(filtered_granges)
}


# This function takes gene and interactions_granges, as well as ArchR object and makes genome browser plot showing interactions
# Can select how to group cells for browser using group_by
# will automatically zoom out so all interactions are seen and centre plot on the gene, can extend further using extend_by

ArchR_PlotInteractions <- function(ArchR_obj, gene, gene_locations, interactions_granges, extend_by = 500, max_dist = Inf, group_by = "clusters", highlight_granges = NULL, return_plot = TRUE){
  
  # extract gene TSS coord
  gene_coord <- as.numeric(ArchR_ExtractTss(gene_locations, gene)[2])
  
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
label <- unique(sub('_.*', '', list.files(paste0(data_path, "rds_files/"))))
print(label) 

ArchR <- loadArchRProject(path = paste0(data_path, "rds_files/", label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getPeakSet(ArchR)
getAvailableMatrices(ArchR)

############################## Read in interactions data #######################################

# read in P2G dataframe
p2g_df <- read.csv(paste0(data_path, "Peak_to_gene_linkage_df_250000_distance.csv"))
head(p2g_df)

# read in P2G granges
p2g_granges <- readRDS(paste0(data_path, "Peak_to_gene_linkage_df_250000_distance.RDS"))
head(p2g_granges)

# read in gene locations
gene_locations <- readRDS(paste0(data_path, "Gene_locations.RDS"))
head(gene_locations)

###########################################################################################
############################## BROWSER TRACKS P2G LINKAGE #################################

# set genes of interest
genes <- c("SIX1", "IRF6", "DLX5", "DLX6", "GATA2", "GATA3", "TFAP2A", "TFAP2B", "TFAP2C", "PITX1", "PITX2",
           "PAX7", "MSX1", "CSRNP1", "ETS1", "SOX9", "SOX8", "SOX10", "SOX5",
           "LMX1B", "ZEB2", "SOX21", "NKX6-2")


# make list of known enhancers organised by genes
SIX1_enhancers_df <- data.frame(
  chr = c("chr5", "chr5", "chr5", "chr5"),
  start = c(54587474, 54589267, 54590423, 54596822),
  end = c(54587537, 54589478, 54590562, 54597223)
)
SOX2_enhancers_df <- data.frame(
  chr = c("chr9", "chr9", "chr9", "chr9", "chr9", "chr9"),
  start = c(17013525, 17029120, 17041213, 17003823, 17018130, 16883843),
  end = c(17013822, 17029653, 17041793, 17004302, 17018498, 16885306)
)
SOX10_enhancers_df <- data.frame(
  chr = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1"),
  start = c(51032535, 51036237, 51046546, 51048966, 51051673, 51065028),
  end = c(51035291, 51038859, 51049401, 51050564, 51054928, 51068292)
)
ETS1_enhancers_df <- data.frame(
  chr = c("chr24", "chr24", "chr24", "chr24", "chr24", "chr24"),
  start = c(267288, 385057, 397366, 495221, 1213638, 1717736),
  end = c(269227, 386795, 399331, 497158, 1215537, 1719643)
)
FOXI3_enhancers_df <- data.frame(
  chr = c("chr4", "chr4", "chr4", "chr4", "chr4"),
  start = c(85595131, 85594801, 85595131, 85593696, 86017750),
  end = c(85596131, 85596939, 85595778, 85597534, 86019361)
)
FOXD3_enhancers_df <- data.frame(
  chr = c("chr8", "chr8", "chr8", "chr8", "chr8", "chr8", "chr8", "chr8"),
  start = c(27926327, 27934827, 27979648, 28017350, 28087080, 28106380, 28139635, 28399561),
  end = c(27928266, 27936786, 27981550, 28019278, 28088418, 28108358, 28141548, 28401557)
)
LMX1A_enhancers_df <- data.frame(
  chr = c("chr8", "chr8"),
  start = c(5172341, 5570901),
  end = c(5173737, 5572297)
)
LMX1B_enhancers_df <- data.frame(
  chr = c("chr17"),
  start = c(10470326),
  end = c(10472156)
)
MSX1_enhancers_df <- data.frame(
  chr = c("chr4", "chr4", "chr4", "chr4", "chr4", "chr4"),
  start = c(78295854, 78663629, 78711612, 78728706, 78895539, 79799413),
  end = c(78297772, 78665953, 78713525, 78730951, 78897524, 79801348)
)
PAX2_enhancers_df <- data.frame(
  chr = c("chr6", "chr6", "chr6", "chr6"),
  start = c(17836767, 17850104, 18094913, 18079378),
  end = c(17837601, 17850840, 18095104, 18079989)
)
PAX7_enhancers_df <- data.frame(
  chr = c("chr21", "chr21", "chr21", "chr21"),
  start = c(3826990, 4238685, 4523315, 4767796),
  end = c(3828981, 4240684, 4525246, 4769743)
)
TFAP2A_enhancers_df <- data.frame(
  chr = c("chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2"),
  start = c(63311515, 63057851, 63310862, 63365989, 63389615, 63453876, 63522766, 63714903, 63772784, 63870453, 63878567, 64233006, 64421924),
  end = c(63312508, 63059773, 63312573, 63367907, 63391605, 63455810, 63524480, 63716810, 63774778, 63871022, 63880501, 64234932, 64423822)
)
TFAP2B_enhancers_df <- data.frame(
  chr = c("chr3", "chr3", "chr3", "chr3", "chr3", "chr3", "chr3", "chr3", "chr3", "chr3"),
  start = c(107523019, 107735860, 107849673, 107872557, 107894338, 107927002, 108021169, 108046754, 108072186, 108223899),
  end = c(107524948, 107737597, 107851588, 107874552, 107896273, 107928971, 108021490, 108048735, 108074077, 108225829)
)

enhancers_df_list <- list(
  SIX1_enhancers_df, SOX2_enhancers_df, SOX10_enhancers_df, ETS1_enhancers_df, FOXI3_enhancers_df, FOXD3_enhancers_df,
  LMX1A_enhancers_df, LMX1B_enhancers_df, MSX1_enhancers_df, PAX2_enhancers_df, PAX7_enhancers_df, TFAP2A_enhancers_df, TFAP2B_enhancers_df
)
names(enhancers_df_list) <- c("SIX1", "SOX2", "SOX10", "ETS1", "FOXI3", "FOXD3", "LMX1A", "LMX1B", "MSX1", "PAX2", "PAX7", "TFAP2A", "TFAP2B")


# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest - from the csv file
  interactions <- p2g_df %>% filter(gene_name %in% gene)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # extract loops between the gene and these peaks - now from archr but make it from csv file?
  extracted_loops <- ArchR_ExtractLoopsToPlot(ArchR, 
                                              gene = gene, gene_locations = gene_locations,
                                              interacting_peaks = interacting_peaks, 
                                              interactions_granges = p2g_granges,
                                              corCutOff = 0.5)
  
  # make plot of these interactions
  p <- ArchR_PlotInteractions(ArchR, gene = gene, gene_locations = gene_locations,
                              interactions_granges = extracted_loops, return_plot = TRUE,
                              group_by = opt$group_by,
                              extend_by = 500, max_dist = Inf)
  
  
  grid::grid.newpage()
  
  # plot
  png(paste0(plot_path, gene, '_interactions_browser_plot.png'), height = 15, width = 18, units = 'cm', res = 400)
  grid::grid.draw(p[[1]])
  graphics.off()
  
  # overlay known enhancers
  if (gene %in% names(enhancers_df_list)){
    enhancers_granges <- makeGRangesFromDataFrame(enhancers_df_list[[gene]])
    grid::grid.newpage()
    p <- ArchR_PlotInteractions(ArchR, gene = gene, gene_locations = gene_locations,
                                interactions_granges = extracted_loops, 
                                return_plot = TRUE,
                                extend_by = 500, max_dist = Inf, 
                                highlight_granges = enhancers_granges,
                                group_by = opt$group_by)
    png(paste0(plot_path, gene, '_interactions_browser_plot_enhancers.png'), height = 15, width = 18, units = 'cm', res = 400)
    grid::grid.draw(p[[1]])
    graphics.off()
  }
  
  
}