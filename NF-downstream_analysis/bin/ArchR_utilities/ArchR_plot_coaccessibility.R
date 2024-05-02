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
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Metacell_ID_purity/rds_files/"
    # co-accessibility csv and RDS files
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Single_cell_integration/"
    # # ArchR object with coaccessibility data
    # data_path = "./output/NF-downstream_analysis/Processing/FullData/Single_cell_integration/"
    
    # rds_path = "./output/NF-downstream_analysis/Processing/ss8/Coaccessibility/rds_files/"
    # plot_path = "./output/NF-downstream_analysis/Processing/ss8/Coaccessibility/plots/"
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    txt_path = "./text_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    #addArchRThreads(threads = ncores)
    addArchRThreads(threads = 1) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  dir.create(txt_path, recursive = T)
  
  
}

set.seed(42)

########################       CELL STATE COLOURS    ########################################
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak', 
                              'PGC', 'BI', 'meso', 'endo', 'MIXED', 'Unmapped',
                              'Neural', 'Placodal', 'Non-neural', 'Contam')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                "#A5718D", "#3F918C", "#ed5e5f", "#9792A3",
                                "#7C8483", "#EAEAEA")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                       'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo',
                                       'Neural', 'Placodal', 'Non-neural', 'Contam',
                                       'MIXED', 'Unmapped')

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

ArchR_ExtractLoopsToPlot <- function(gene, gene_locations, interacting_peaks, interactions_granges){
  
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

  # filter all interactions from P2G for these ones
  filtered_granges <- subset(interactions_granges, start %in% df_ordered$start & end %in% df_ordered$end)
  
  # return interactions
  return(filtered_granges)
}


# This function takes gene and interactions_granges, as well as ArchR object and makes genome browser plot showing interactions
# Can select how to group cells for browser using group_by
# will automatically zoom out so all interactions are seen and centre plot on the gene, can extend further using extend_by

ArchR_PlotInteractions <- function(ArchR_obj, gene, gene_locations, interactions_granges, extend_by = 500, max_dist = Inf, group_by = "clusters", highlight_granges = NULL, return_plot = TRUE, pal = pal){
  
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
    p <- ArchR::plotBrowserTrack(
      ArchRProj = ArchR_obj,
      groupBy = group_by,
      geneSymbol = gene, 
      upstream = distance,
      downstream = distance,
      loops = interactions_granges,
      pal = pal,
      title = paste0(gene, "locus - ", length(interactions_granges), " interactions found - distance around: ", distance, "bp")
    )
  } else {
    p <- ArchR::plotBrowserTrack(
      ArchRProj = ArchR_obj,
      groupBy = group_by,
      geneSymbol = gene, 
      upstream = distance,
      downstream = distance,
      loops = interactions_granges,
      highlight = highlight_granges,
      pal = pal,
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

print("reading in interactions data...")

# read in P2G dataframe
p2g_df <- read.csv(paste0(data_path, "csv_files/", "Peak_to_gene_linkage_df_250000_distance.csv"))
head(p2g_df)

# read in P2G granges
p2g_granges <- readRDS(paste0(data_path, "csv_files/", "Peak_to_gene_linkage_250000_distance.RDS"))
p2g_granges

# read in gene locations
gene_locations <- readRDS(paste0(data_path, "csv_files/", "Gene_locations.RDS"))
gene_locations

###########################################################################################
############################## BROWSER TRACKS P2G LINKAGE #################################

plot_path = "./plots/key_genes/"
dir.create(plot_path, recursive = T)

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

# init list of interacting peaks
peaks <- list()

# filter p2g dataframe in the same way did for GRNi
P2G_filt <- p2g_df %>% dplyr::filter(FDR < 0.01)

# set colour scheme

pal <- scHelper_cell_type_colours[unique(getCellColData(ArchR, select = opt$group_by)[,1])]

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- P2G_filt %>% 
    filter(gene_name %in% gene)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
    # make plot of these interactions
    p <- ArchR_PlotInteractions(ArchR, gene = gene, gene_locations = gene_locations,
                                interactions_granges = extracted_loops, return_plot = TRUE,
                                group_by = opt$group_by,
                                extend_by = 500, max_dist = Inf,
                                pal = pal)
    
    
    grid::grid.newpage()
    
    # plot
    png(paste0(plot_path, gene, '_interactions_browser_plot.png'), height = 15, width = 18, units = 'cm', res = 400)
    grid::grid.draw(p[[1]])
    graphics.off()
    
    #overlay known enhancers - not working atm
    #Error in ArchR::plotBrowserTrack(ArchRProj = ArchR_obj, groupBy = group_by,  :
    #unused argument (highlight = highlight_granges)
    #Calls: ArchR_PlotInteractions
    #Execution halted
    if (gene %in% names(enhancers_df_list)){
      enhancers_granges <- makeGRangesFromDataFrame(enhancers_df_list[[gene]])
      grid::grid.newpage()
      p <- ArchR_PlotInteractions(ArchR, gene = gene, gene_locations = gene_locations,
                                  interactions_granges = extracted_loops,
                                  return_plot = TRUE,
                                  extend_by = 500, max_dist = Inf,
                                  highlight_granges = enhancers_granges,
                                  group_by = opt$group_by, pal = pal)
      png(paste0(plot_path, gene, '_interactions_browser_plot_enhancers.png'), height = 15, width = 18, units = 'cm', res = 400)
      grid::grid.draw(p[[1]])
      graphics.off()
    }
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "key_genes_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)

###########################################################################################
############################## ALEX'S GMs P2G LINKAGE #################################

# from Alex's GMs analysis, GMs made on ss8 data, these are the ones used for coexpression and temporal analysis

# PPR specific GMs
# GM12; AKR1D1, ATP1A1, ATP1B1, ATP2B1, B3GNT7, CCDC25, CGNL1, DAG1, EMILIN2, ENSGALG00000004814, EPCAM, FAM184B, FAM89A, GATA2, GATA3, IRAK2, IRF6, MAP7, METRNL, MPZL3, NET1, PLEKHA5, POGLUT2, PPFIBP1, RGN, SLC16A10, SLC25A4, TSPAN13, TTC39A, UNC5B, Z-TJP2
# GM14; BRINP1, CITED4, DLX5, DLX6, ENSGALG00000023936, ENSGALG00000040010, FN1, HESX1, HIF1A, KCNAB1, LAMB1, LRP11, NFKB1, PAX6, PITX1, PITX2, SEM1, SFRP1, SH3D19, SHISA2, SIX3, SPON1, SST, WDR1, Z-HAPLN1
# GM13; AKAP12, ASS1, BASP1, CD99, ENSGALG00000011296, ENSGALG00000041054, ENSGALG00000042443, EYA2, FERMT2, LGMN, METTL24, NR2F2, NUCKS1, SIX1


# PPR_genes <- c("AKR1D1", "ATP1A1", "ATP1B1", "ATP2B1", "B3GNT7", "CCDC25", "CGNL1", "DAG1", "EMILIN2", "ENSGALG00000004814", "EPCAM", "FAM184B", "FAM89A", "GATA2", "GATA3", "IRAK2", "IRF6", "MAP7", "METRNL", "MPZL3", "NET1", "PLEKHA5", "POGLUT2", "PPFIBP1", "RGN", "SLC16A10", "SLC25A4", "TSPAN13", "TTC39A", "UNC5B", "Z-TJP2",
#                "BRINP1", "CITED4", "DLX5", "DLX6", "ENSGALG00000023936", "ENSGALG00000040010", "FN1", "HESX1", "HIF1A", "KCNAB1", "LAMB1", "LRP11", "NFKB1", "PAX6", "PITX1", "PITX2", "SEM1", "SFRP1", "SH3D19", "SHISA2", "SIX3", "SPON1", "SST", "WDR1", "Z-HAPLN1",
#                "AKAP12", "ASS1", "BASP1", "CD99", "ENSGALG00000011296", "ENSGALG00000041054", "ENSGALG00000042443", "EYA2", "FERMT2", "LGMN", "METTL24", "NR2F2", "NUCKS1", "SIX1")
# 
# # NC specific GMs
# # GM40; AGTRAP, BRINP2, CDH11, CMTM8, ENSGALG00000001136, FRZB, GADD45A, HUNK, LARP7, LMX1B, MRAS, MSX1, SFRP2, SOX11, SPSB4, Z-FST, ZEB2, ZIC1
# # GM42; BMP5, CDH6, CSRNP1, DRAXIN, EN1, ENSGALG00000013505, FOXD3, NKD1, NRP2, OLFM1, PAX7, SNAI2, SOX9, TFAP2B, TMEM132C, TSPAN18, WLS, WNT6, Z-ENC1, ZFHX4, ZNF423
# # GM43; COL9A3, ENSGALG00000037717, ENSGALG00000053185, ERMN, ETS1, LMO4, PPP1R1C, RASL11B, RFTN2, SOX10, SOX8, TNC, Z-MEF2C, Z-PLK2
# # GM44; CXCR4, ENSGALG00000030512, ENSGALG00000031427, FABP7, FKBP11, GLIPR2, ID1, ID2, MEOX1, MYL4, OLFML3, SOX5, WNT1
# 
# NC_genes <- c("AGTRAP", "BRINP2", "CDH11", "CMTM8", "ENSGALG00000001136", "FRZB", "GADD45A", "HUNK", "LARP7", "LMX1B", "MRAS", "MSX1", "SFRP2", "SOX11", "SPSB4", "Z-FST", "ZEB2", "ZIC1",
#               "BMP5", "CDH6", "CSRNP1", "DRAXIN", "EN1", "ENSGALG00000013505", "FOXD3", "NKD1", "NRP2", "OLFM1", "PAX7", "SNAI2", "SOX9", "TFAP2B", "TMEM132C", "TSPAN18", "WLS", "WNT6", "Z-ENC1", "ZFHX4", "ZNF423",
#               "CDON", "COTL1", "ENSGALG00000048488", "MFAP2", "PRTG", "TUBB3",
#               "COL9A3", "ENSGALG00000037717", "ENSGALG00000053185", "ERMN", "ETS1", "LMO4", "PPP1R1C", "RASL11B", "RFTN2", "SOX10", "SOX8", "TNC", "Z-MEF2C", "Z-PLK2",
#               "CXCR4", "ENSGALG00000030512", "ENSGALG00000031427", "FABP7", "FKBP11", "GLIPR2", "ID1", "ID2", "MEOX1", "MYL4", "OLFML3", "SOX5", "WNT1")


#####################################################################
############################## ss8 PPR GM12 #################################

plot_path = "./plots/ss8_PPR_GM12/"
dir.create(plot_path, recursive = T)

# set genes of interest
genes <- c("AKR1D1", "ATP1A1", "ATP1B1", "ATP2B1", "B3GNT7", "CCDC25", "CGNL1", "DAG1", "EMILIN2", "ENSGALG00000004814", "EPCAM", "FAM184B", "FAM89A", "GATA2", "GATA3", "IRAK2", "IRF6", "MAP7", "METRNL", "MPZL3", "NET1", "PLEKHA5", "POGLUT2", "PPFIBP1", "RGN", "SLC16A10", "SLC25A4", "TSPAN13", "TTC39A", "UNC5B", "Z-TJP2")

# init list of interacting peaks
peaks <- list()

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- p2g_df %>% 
    filter(gene_name %in% gene) %>%
    filter(Correlation > 0.5)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
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
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "ss8_PPR_GM12_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)

#####################################################################
############################## ss8 PPR GM14 #################################

plot_path = "./plots/ss8_PPR_GM14/"
dir.create(plot_path, recursive = T)

# set genes of interest
genes <- c("BRINP1", "CITED4", "DLX5", "DLX6", "ENSGALG00000023936", "ENSGALG00000040010", "FN1", "HESX1", "HIF1A", "KCNAB1", "LAMB1", "LRP11", "NFKB1", "PAX6", "PITX1", "PITX2", "SEM1", "SFRP1", "SH3D19", "SHISA2", "SIX3", "SPON1", "SST", "WDR1", "Z-HAPLN1")

# init list of interacting peaks
peaks <- list()

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- p2g_df %>% 
    filter(gene_name %in% gene) %>%
    filter(Correlation > 0.5)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
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
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "ss8_PPR_GM14_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)

#####################################################################
############################## ss8 PPR GM13 #################################

plot_path = "./plots/ss8_PPR_GM13/"
dir.create(plot_path, recursive = T)

# set genes of interest
genes <- c("AKAP12", "ASS1", "BASP1", "CD99", "ENSGALG00000011296", "ENSGALG00000041054", "ENSGALG00000042443", "EYA2", "FERMT2", "LGMN", "METTL24", "NR2F2", "NUCKS1", "SIX1")


# init list of interacting peaks
peaks <- list()

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- p2g_df %>% 
    filter(gene_name %in% gene) %>%
    filter(Correlation > 0.5)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
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
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "ss8_PPR_GM13_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)

# NC_genes <- c("AGTRAP", "BRINP2", "CDH11", "CMTM8", "ENSGALG00000001136", "FRZB", "GADD45A", "HUNK", "LARP7", "LMX1B", "MRAS", "MSX1", "SFRP2", "SOX11", "SPSB4", "Z-FST", "ZEB2", "ZIC1",
#               "BMP5", "CDH6", "CSRNP1", "DRAXIN", "EN1", "ENSGALG00000013505", "FOXD3", "NKD1", "NRP2", "OLFM1", "PAX7", "SNAI2", "SOX9", "TFAP2B", "TMEM132C", "TSPAN18", "WLS", "WNT6", "Z-ENC1", "ZFHX4", "ZNF423",
#               "CDON", "COTL1", "ENSGALG00000048488", "MFAP2", "PRTG", "TUBB3",
#               "COL9A3", "ENSGALG00000037717", "ENSGALG00000053185", "ERMN", "ETS1", "LMO4", "PPP1R1C", "RASL11B", "RFTN2", "SOX10", "SOX8", "TNC", "Z-MEF2C", "Z-PLK2",
#               "CXCR4", "ENSGALG00000030512", "ENSGALG00000031427", "FABP7", "FKBP11", "GLIPR2", "ID1", "ID2", "MEOX1", "MYL4", "OLFML3", "SOX5", "WNT1")

#####################################################################
############################## ss8_NC_GM40 #################################

plot_path = "./plots/ss8_NC_GM40/"
dir.create(plot_path, recursive = T)

# set genes of interest
genes <- c("AGTRAP", "BRINP2", "CDH11", "CMTM8", "ENSGALG00000001136", "FRZB", "GADD45A", "HUNK", "LARP7", "LMX1B", "MRAS", "MSX1", "SFRP2", "SOX11", "SPSB4", "Z-FST", "ZEB2", "ZIC1")

# init list of interacting peaks
peaks <- list()

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- p2g_df %>% 
    filter(gene_name %in% gene) %>%
    filter(Correlation > 0.5)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
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
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "ss8_NC_GM40_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)

#####################################################################
############################## ss8_NC_GM42 #################################

plot_path = "./plots/ss8_NC_GM42/"
dir.create(plot_path, recursive = T)

# set genes of interest
genes <- c("BMP5", "CDH6", "CSRNP1", "DRAXIN", "EN1", "ENSGALG00000013505", "FOXD3", "NKD1", "NRP2", "OLFM1", "PAX7", "SNAI2", "SOX9", "TFAP2B", "TMEM132C", "TSPAN18", "WLS", "WNT6", "Z-ENC1", "ZFHX4", "ZNF423")

# init list of interacting peaks
peaks <- list()

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- p2g_df %>% 
    filter(gene_name %in% gene) %>%
    filter(Correlation > 0.5)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
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
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "ss8_NC_GM42_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)

#####################################################################
############################## ss8_NC_GM43 #################################

plot_path = "./plots/ss8_NC_GM43/"
dir.create(plot_path, recursive = T)

# set genes of interest
genes <- c("COL9A3", "ENSGALG00000037717", "ENSGALG00000053185", "ERMN", "ETS1", "LMO4", "PPP1R1C", "RASL11B", "RFTN2", "SOX10", "SOX8", "TNC", "Z-MEF2C", "Z-PLK2")

# init list of interacting peaks
peaks <- list()

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- p2g_df %>% 
    filter(gene_name %in% gene) %>%
    filter(Correlation > 0.5)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
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
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "ss8_NC_GM43_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)

#####################################################################
############################## ss8_NC_GM42 #################################

plot_path = "./plots/ss8_NC_GM42/"
dir.create(plot_path, recursive = T)

# set genes of interest
genes <- c("CXCR4", "ENSGALG00000030512", "ENSGALG00000031427", "FABP7", "FKBP11", "GLIPR2", "ID1", "ID2", "MEOX1", "MYL4", "OLFML3", "SOX5", "WNT1")

# init list of interacting peaks
peaks <- list()

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- p2g_df %>% 
    filter(gene_name %in% gene) %>%
    filter(Correlation > 0.5)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
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
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "ss8_NC_GM42_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)



# Full data GMs!!!!
# PPR:
#   GM23 <- c("ASS1", "BMP4", "BMP6", "CLDN3", "CSRP2", "DLX5", "DLX6", "EMILIN2", "ENSGALG00000001885", "ENSGALG00000023936", "ENSGALG00000040010", "ENSGALG00000042443", "ENSGALG00000050334", "ENSGALG00000052786", "EYA2", "FABP3", "FAM184B", "FAM89A", "FN1", "GATA2", "GATA3", "HAS2", "KCNAB1", "KRT18", "KRT19", "KRT7", "LAMB1", "NEDD9", "NET1", "PITX1", "SIX1", "SPON1", "TFAP2A", "TUBAL3", "UNC5B", "Z-HAPLN1")
# Neural:
#   GM9 <- c("CDH2", "CNTNAP1", "ENSGALG00000054487", "FEZF2", "LMO1", "LRRN1", "MSN", "MYLK", "NUAK1", "SOX21", "TOX", "Z-GRP", "Z-RAX")
# Epiblast:
#   GM21 <- c("ACSL1", "ACTN1", "AKR1D1", "ATP1A1", "ATP1B1", "B3GNT7", "BAMBI", "CCDC25", "CDK6", "CITED4", "CTSB", "CXCL12", "DAG1", "ENSGALG00000002988", "ENSGALG00000004814", "ENSGALG00000005572", "ENSGALG00000008518", "ENSGALG00000016570", "ENSGALG00000027805", "ENSGALG00000051984", "FIBIN", "HOMER2", "ID3", "IRF6", "KLF3", "LY6E", "MAP7", "METRNL", "MYL3", "MYL9", "NDRG1", "NT5DC2", "PDZK1IP1", "RGN", "S100A11", "SMIM4", "TFAP2C", "TTC39A")
#   GM24 <- c("ADD3", "AIG1", "AMY1A", "AP1S3", "ARFGAP1", "ARL8BL", "ATP2B1", "CACNG3", "CD24", "CLDN1", "ENSGALG00000026754", "ENSGALG00000040886", "ENSGALG00000043421", "EPCAM", "G3BP1", "HSPA9", "IVNS1ABP", "MPZL3", "PDGFA", "SALL4", "SLC25A4", "SUCLG1", "TSPAN13", "Z-CTSV", "Z-PLPP1")

#####################################################################
############################## FullData_PPR_GM23 #################################

plot_path = "./plots/FullData_PPR_GM23/"
dir.create(plot_path, recursive = T)

# set genes of interest
genes <- c("ASS1", "BMP4", "BMP6", "CLDN3", "CSRP2", "DLX5", "DLX6", "EMILIN2", "ENSGALG00000001885", "ENSGALG00000023936", "ENSGALG00000040010", "ENSGALG00000042443", "ENSGALG00000050334", "ENSGALG00000052786", "EYA2", "FABP3", "FAM184B", "FAM89A", "FN1", "GATA2", "GATA3", "HAS2", "KCNAB1", "KRT18", "KRT19", "KRT7", "LAMB1", "NEDD9", "NET1", "PITX1", "SIX1", "SPON1", "TFAP2A", "TUBAL3", "UNC5B", "Z-HAPLN1")

# init list of interacting peaks
peaks <- list()

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- p2g_df %>% 
    filter(gene_name %in% gene) %>%
    filter(Correlation > 0.5)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
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
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "FullData_PPR_GM23_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)

#####################################################################
############################## FullData_Neural_GM9 #################################

plot_path = "./plots/FullData_Neural_GM9/"
dir.create(plot_path, recursive = T)

# set genes of interest
genes <- c("CDH2", "CNTNAP1", "ENSGALG00000054487", "FEZF2", "LMO1", "LRRN1", "MSN", "MYLK", "NUAK1", "SOX21", "TOX", "Z-GRP", "Z-RAX")

# init list of interacting peaks
peaks <- list()

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- p2g_df %>% 
    filter(gene_name %in% gene) %>%
    filter(Correlation > 0.5)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
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
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "FullData_Neural_GM9_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)

#####################################################################
############################## FullData_Epiblast_GM21 #################################

plot_path = "./plots/FullData_Epiblast_GM21/"
dir.create(plot_path, recursive = T)

# set genes of interest
genes <- c("ACSL1", "ACTN1", "AKR1D1", "ATP1A1", "ATP1B1", "B3GNT7", "BAMBI", "CCDC25", "CDK6", "CITED4", "CTSB", "CXCL12", "DAG1", "ENSGALG00000002988", "ENSGALG00000004814", "ENSGALG00000005572", "ENSGALG00000008518", "ENSGALG00000016570", "ENSGALG00000027805", "ENSGALG00000051984", "FIBIN", "HOMER2", "ID3", "IRF6", "KLF3", "LY6E", "MAP7", "METRNL", "MYL3", "MYL9", "NDRG1", "NT5DC2", "PDZK1IP1", "RGN", "S100A11", "SMIM4", "TFAP2C", "TTC39A")

# init list of interacting peaks
peaks <- list()

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- p2g_df %>% 
    filter(gene_name %in% gene) %>%
    filter(Correlation > 0.5)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
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
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "FullData_Epiblast_GM21_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)

#####################################################################
############################## FullData_Epiblast_GM24 #################################

plot_path = "./plots/FullData_Epiblast_GM24/"
dir.create(plot_path, recursive = T)

# set genes of interest
genes <- c("ADD3", "AIG1", "AMY1A", "AP1S3", "ARFGAP1", "ARL8BL", "ATP2B1", "CACNG3", "CD24", "CLDN1", "ENSGALG00000026754", "ENSGALG00000040886", "ENSGALG00000043421", "EPCAM", "G3BP1", "HSPA9", "IVNS1ABP", "MPZL3", "PDGFA", "SALL4", "SLC25A4", "SUCLG1", "TSPAN13", "Z-CTSV", "Z-PLPP1")

# init list of interacting peaks
peaks <- list()

# loop through genes and make plots (+ with enhancers highlighted if have that data)
for (gene in genes){
  
  print(gene)
  
  # extract interactions to gene of interest above a correlation cut off
  interactions <- p2g_df %>% 
    filter(gene_name %in% gene) %>%
    filter(Correlation > 0.5)
  print(paste0(nrow(interactions), " interactions found"))
  
  # extract the peak IDs that interact with that gene - from the csv file
  interacting_peaks <- unique(interactions$PeakID)
  
  # write the interacting peaks to list
  peaks[[gene]] <- interacting_peaks
  
  # only run this bit if there are interacting peaks
  if (length(interacting_peaks) > 0){
    # extract the loops between the gene and these peaks for plotting
    extracted_loops <- ArchR_ExtractLoopsToPlot(gene = gene, gene_locations = gene_locations,
                                                interacting_peaks = interacting_peaks, 
                                                interactions_granges = p2g_granges)
    extracted_loops
    
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
    
  }
  
}

# Open the file for writing
file_conn <- file(paste0(txt_path, "FullData_Epiblast_GM24_interacting_peaks.txt"), "w")

# Write the list to the file with the desired format
for (key in names(peaks)) {
  cat(paste(key, ";", paste0("\"", peaks[[key]], "\"", collapse = ", "), "\n"), file = file_conn)
}

# Close the file connection
close(file_conn)