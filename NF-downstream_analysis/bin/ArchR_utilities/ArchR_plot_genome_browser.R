#!/usr/bin/env Rscript

print("plot genome browser plots of known enhancers - needs cell grouping")

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
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Metacell_to_singlecell/"
    
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

###### STAGE COLOURS ####################
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order


############################## FUNCTIONS #######################################

## function which takes a table of regions (chr | start | end) and plots them individually on genome browser plots
# the given regions are highlighted the the plot is extended by argument extend_by
# the accessibility is shown over cells grouped by group_by argument 

ArchR_genome_browser_plot <- function(ArchR, regions_df, group_by = "clusters", extend_by = 500, pal, show_peaks = TRUE){
  
  GRanges_highlight <- makeGRangesFromDataFrame(regions_df)
  df_extended <- regions_df %>% 
    mutate(start = start - extend_by) %>%
    mutate(end = end + extend_by)
  GRanges_plot <- makeGRangesFromDataFrame(df_extended)
  
  if (show_peaks){
    p <- plotBrowserTrack(
      ArchRProj = ArchR, 
      groupBy = group_by, 
      region = GRanges_plot,
      highlight = GRanges_highlight,
      features = getPeakSet(ArchR),
      plotSummary = c("bulkTrack", "geneTrack", "featureTrack"),
      facetbaseSize = 12,
      baseSize = 12,
      title = "",
      pal = pal)
  } else {
    p <- plotBrowserTrack(
      ArchRProj = ArchR, 
      groupBy = group_by, 
      region = GRanges_plot,
      highlight = GRanges_highlight,
      plotSummary = c("bulkTrack", "geneTrack"),
      facetbaseSize = 12,
      baseSize = 12,
      title = "",
      pal = pal)
  }
  return(p)
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
getAvailableMatrices(ArchR)
colnames(getCellColData(ArchR))

###########################################################################################
############################## Set up enhancers df ########################################

# make list of known enhancers organised by genes
AXUD1_enhancers_df <- data.frame(
  chr = c("chr2"),
  start = c(5263169),
  end = c(5263696)
)
CHAMP1_enhancers_df <- data.frame(
  chr = c("chr1"),
  start = c(137543885),
  end = c(137545827)
)
DACT2_enhancers_df <- data.frame(
  chr = c("chr3"),
  start = c(42085878),
  end = c(42087663)
)
DRAXIN_enhancers_df <- data.frame(
  chr = c("chr21"),
  start = c(5765452),
  end = c(5767088)
)
EPHA4_enhancers_df <- data.frame(
  chr = c("chr9"),
  start = c(7800500),
  end = c(7802007)
)
ETS1_enhancers_df <- data.frame(
  chr = c("chr24", "chr24", "chr24", "chr24"),
  start = c(385057, 397366, 1213638, 1717736),
  end = c(386795, 399331, 1215537, 1719643)
)
FOXD3_enhancers_df <- data.frame(
  chr = c("chr8", "chr8", "chr8", "chr8", "chr8", "chr8", "chr8", "chr8"),
  start = c(27926327, 27934827, 27979648, 28017350, 28087080, 28106380, 28139635, 28399561),
  end = c(27928266, 27936786, 27981550, 28019278, 28088418, 28108358, 28141548, 28401557)
)
FOXI1_enhancers_df <- data.frame(
  chr = c("chr13", "chr13", "chr13"),
  start = c(5077007, 5378232, 5531347),
  end = c(5078794, 5380149, 5533263)
)
FOXI3_enhancers_df <- data.frame(
  chr = c("chr4"),
  start = c(86017750),
  end = c(86019361)
)
GLI2_enhancers_df <- data.frame(
  chr = c("chr7"),
  start = c(26022718),
  end = c(26024968)
)
LMX1A_enhancers_df <- data.frame(
  chr = c("chr8"),
  start = c(5570901),
  end = c(5572297)
)
LMX1B_enhancers_df <- data.frame(
  chr = c("chr17"),
  start = c(10470326),
  end = c(10472156)
)
MSX1_enhancers_df <- data.frame(
  chr = c("chr4", "chr4", "chr4", "chr4"),
  start = c(78295854, 78663629, 78711612, 79799413),
  end = c(78297772, 78665953, 78713525, 79801348)
)
OTX2_enhancers_df <- data.frame(
  chr = c("chr5"),
  start = c(56120990),
  end = c(56122619)
)
PAX2_enhancers_df <- data.frame(
  chr = c("chr6"),
  start = c(17836767),
  end = c(17837601)
)
PAX7_enhancers_df <- data.frame(
  chr = c("chr21", "chr21", "chr21", "chr21"),
  start = c(3826990, 4238685, 4523315, 4767796),
  end = c(3828981, 4240684, 4525246, 4769743)
)
SETD2_enhancers_df <- data.frame(
  chr = c("chr2"),
  start = c(3958727),
  end = c(3961289)
)
SIX1_enhancers_df <- data.frame(
  chr = c("chr5", "chr5", "chr5", "chr5"),
  start = c(54587474, 54589267, 54590423, 54596822),
  end = c(54587537, 54589478, 54590562, 54597223)
)
SIX3_enhancers_df <- data.frame(
  chr = c("chr3"),
  start = c(25787500),
  end = c(25789339)
)
SOX10_enhancers_df <- data.frame(
  chr = c("chr1"),
  start = c(51065028),
  end = c(51068292)
)
SOX13_enhancers_df <- data.frame(
  chr = c("chr26", "chr26"),
  start = c(1699737, 2397615),
  end = c(1700419, 2399582)
)
SOX2_enhancers_df <- data.frame(
  chr = c("chr9", "chr9", "chr9"),
  start = c(17029120, 17041213, 16883843),
  end = c(17029653, 17041793, 16885306)
)
SOX8_enhancers_df <- data.frame(
  chr = c("chr14"),
  start = c(6235716),
  end = c(6237441)
)
SOX9_enhancers_df <- data.frame(
  chr = c("chr18"),
  start = c(8486733),
  end = c(8488450)
)
SP5_enhancers_df <- data.frame(
  chr = c("chr7"),
  start = c(18346322),
  end = c(18347947)
)
TFAP2A_enhancers_df <- data.frame(
  chr = c("chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2", "chr2"),
  start = c(63311515, 63057851, 63365989, 63389615, 63453876, 63714903, 63772784, 63878567, 64233006, 64421924),
  end = c(63312508, 63059773, 63367907, 63391605, 63455810, 63716810, 63774778, 63880501, 64234932, 64423822)
)
TFAP2B_enhancers_df <- data.frame(
  chr = c("chr3", "chr3", "chr3", "chr3", "chr3"),
  start = c(107735860, 107872557, 107894338, 107927002, 108072186),
  end = c(107737597, 107874552, 107896273, 107928971, 108074077)
)
ZBTB16_enhancers_df <- data.frame(
  chr = c("chr24"),
  start = c(4557929),
  end = c(4558406)
)
ZIC1_enhancers_df <- data.frame(
  chr = c("chr9"),
  start = c(12146874),
  end = c(12148841)
)
ZIC2_enhancers_df <- data.frame(
  chr = c("chr1"),
  start = c(145830051),
  end = c(145831958)
)
ZNF385C_enhancers_df <- data.frame(
  chr = c("chr27"),
  start = c(7610684),
  end = c(7611071)
)

enhancers_df_list <- list(
  AXUD1_enhancers_df, CHAMP1_enhancers_df, DACT2_enhancers_df, DRAXIN_enhancers_df, EPHA4_enhancers_df,
  ETS1_enhancers_df, FOXD3_enhancers_df, FOXI1_enhancers_df, FOXI3_enhancers_df,
  GLI2_enhancers_df, LMX1A_enhancers_df, LMX1B_enhancers_df, MSX1_enhancers_df,
  OTX2_enhancers_df, PAX2_enhancers_df, PAX7_enhancers_df, SETD2_enhancers_df,
  SIX1_enhancers_df, SIX3_enhancers_df, 
  SOX10_enhancers_df, SOX13_enhancers_df, SOX2_enhancers_df, SOX8_enhancers_df, SOX9_enhancers_df,
  SP5_enhancers_df, TFAP2A_enhancers_df, TFAP2B_enhancers_df,
  ZBTB16_enhancers_df, ZIC1_enhancers_df, ZIC2_enhancers_df, ZNF385C_enhancers_df
)
names(enhancers_df_list) <- c("AXUD1", "CHAMP1", "DACT2", "DRAXIN", "EPHA4", 
    "ETS1", "FOXD3", "FOXI1", "FOXI3", 
    "GLI2", "LMX1A", "LMX1B", "MSX1",
    "OTX2", "PAX2", "PAX7", "SETD2",
    "SIX1", "SIX3",
    "SOX10", "SOX13", "SOX2", "SOX8", "SOX9",
    "SP5", "TFAP2A", "TFAP2B", 
    "ZBTB16", "ZIC1", "ZIC2", "ZNF385C"
    )

################################################################################################################
############################## Plot tracks for single cell labels broad ########################################

plot_path = "./output/NF-downstream_analysis/Processing/ss8/Lit_enhancers_genome_browser_plots/plots/single_cell_labels_broad/"
#plot_path = "./plots/single_cell_labels_broad/"
dir.create(plot_path, recursive = T)

pal <- scHelper_cell_type_colours[unique(ArchR$transferred_scHelper_cell_type_broad)]

# loop through genes and then loop through enhancers within each gene
for (gene in names(enhancers_df_list)){
  enhancers_df <- enhancers_df_list[[gene]]
  
  for (row in 1:nrow(enhancers_df)){
    enhancer <- enhancers_df[row,]
    print(enhancer)
    
    p <- ArchR_genome_browser_plot(ArchR, enhancer, group_by = "transferred_scHelper_cell_type_broad", 
                                   pal, extend_by = 5000, show_peaks = FALSE)
    
    # plot
    grid::grid.newpage()
    png(paste0(plot_path, gene, '_enhancer_', row, '_browser_plot.png'), height = 15, width = 18, units = 'cm', res = 400)
    grid::grid.draw(p)
    graphics.off()
  }
  
}

################################################################################################################
############################## Plot tracks for single cell labels ########################################

plot_path = "./output/NF-downstream_analysis/Processing/ss8/Lit_enhancers_genome_browser_plots/plots/single_cell_labels/"
plot_path = "./plots/single_cell_labels/"
dir.create(plot_path, recursive = T)

pal <- scHelper_cell_type_colours[unique(ArchR$transferred_scHelper_cell_type)]

# loop through genes and then loop through enhancers within each gene
for (gene in names(enhancers_df_list)){
  enhancers_df <- enhancers_df_list[[gene]]
  
  for (row in 1:nrow(enhancers_df)){
    enhancer <- enhancers_df[row,]
    print(enhancer)
    
    p <- ArchR_genome_browser_plot(ArchR, enhancer, group_by = "transferred_scHelper_cell_type", 
                                   pal, extend_by = 5000, show_peaks = FALSE)
    
    # plot
    grid::grid.newpage()
    png(paste0(plot_path, gene, '_enhancer_', row, '_browser_plot.png'), height = 15, width = 18, units = 'cm', res = 400)
    grid::grid.draw(p)
    graphics.off()
  }
  
}

################################################################################################################
############################## Plot tracks for metacell cell labels ########################################
plot_path = "./output/NF-downstream_analysis/Processing/ss8/Lit_enhancers_genome_browser_plots/plots/SEACell_scHelper_cell_type/"

plot_path = "./plots/SEACell_scHelper_cell_type/"
dir.create(plot_path, recursive = T)

pal <- scHelper_cell_type_colours[unique(ArchR$SEACell_scHelper_cell_type)]

# loop through genes and then loop through enhancers within each gene
for (gene in names(enhancers_df_list)){
  enhancers_df <- enhancers_df_list[[gene]]
  
  for (row in 1:nrow(enhancers_df)){
    enhancer <- enhancers_df[row,]
    print(enhancer)
    
    p <- ArchR_genome_browser_plot(ArchR, enhancer, group_by = "SEACell_scHelper_cell_type", 
                                   pal, extend_by = 5000, show_peaks = FALSE)
    
    # plot
    grid::grid.newpage()
    png(paste0(plot_path, gene, '_enhancer_', row, '_browser_plot.png'), height = 15, width = 18, units = 'cm', res = 400)
    grid::grid.draw(p)
    graphics.off()
  }
  
}


################################################################################################################
############################## Plot tracks for metacell cell labels broad ########################################

plot_path = "./plots/SEACell_scHelper_cell_type_broad/"
plot_path = "./output/NF-downstream_analysis/Processing/ss8/Lit_enhancers_genome_browser_plots/plots/SEACell_scHelper_cell_type_broad/"

dir.create(plot_path, recursive = T)

pal <- scHelper_cell_type_colours[unique(ArchR$SEACell_scHelper_cell_type_broad)]

# loop through genes and then loop through enhancers within each gene
for (gene in names(enhancers_df_list)){
  enhancers_df <- enhancers_df_list[[gene]]
  
  for (row in 1:nrow(enhancers_df)){
    enhancer <- enhancers_df[row,]
    print(enhancer)
    
    p <- ArchR_genome_browser_plot(ArchR, enhancer, group_by = "SEACell_scHelper_cell_type_broad", 
                                   pal, extend_by = 5000, show_peaks = FALSE)
    
    # plot
    grid::grid.newpage()
    png(paste0(plot_path, gene, '_enhancer_', row, '_browser_plot.png'), height = 15, width = 18, units = 'cm', res = 400)
    grid::grid.draw(p)
    graphics.off()
  }
  
}

####################################################################################
################# New enhancers to test ############################################

PPR_enhancers_df <- data.frame(
  chr = c("chr2", "chr2", "chr9", "chr4", "chr1", "chr2"),
  start = c(10216477, 11936117, 16689582, 69910791, 88390937, 149406064),
  end = c(10216977, 11936617, 16690082, 69911291, 88391437, 149406564)
)

NC_enhancers_df <- data.frame(
  chr = c("chr1", "chr2", "chr8", "chr1", "chr4"),
  start = c(4741943, 42661874, 7216956, 61674271, 39407032),
  end = c(4742443, 42662374, 7217456, 61674771, 39407532)
)


enhancers_df_list <- list(PPR_enhancers_df, NC_enhancers_df)
names(enhancers_df_list) <- c("PPR", "NC")


plot_path = "./output/NF-downstream_analysis/Processing/ss8/New_enhancers/plots/single_cell_labels_broad/"
#plot_path = "./plots/single_cell_labels_broad/"
dir.create(plot_path, recursive = T)

pal <- scHelper_cell_type_colours[unique(ArchR$transferred_scHelper_cell_type_broad)]

# loop through genes and then loop through enhancers within each gene
for (gene in names(enhancers_df_list)){
  enhancers_df <- enhancers_df_list[[gene]]
  
  for (row in 1:nrow(enhancers_df)){
    enhancer <- enhancers_df[row,]
    print(enhancer)
    
    p <- ArchR_genome_browser_plot(ArchR, enhancer, group_by = "transferred_scHelper_cell_type_broad", 
                                   pal, extend_by = 5000, show_peaks = FALSE)
    
    # plot
    grid::grid.newpage()
    png(paste0(plot_path, gene, '_enhancer_', row, '_browser_plot.png'), height = 15, width = 18, units = 'cm', res = 400)
    grid::grid.draw(p)
    graphics.off()
  }
  
}


plot_path = "./output/NF-downstream_analysis/Processing/ss8/New_enhancers/plots/single_cell_labels/"
dir.create(plot_path, recursive = T)

pal <- scHelper_cell_type_colours[unique(ArchR$transferred_scHelper_cell_type)]

# loop through genes and then loop through enhancers within each gene
for (gene in names(enhancers_df_list)){
  enhancers_df <- enhancers_df_list[[gene]]
  
  for (row in 1:nrow(enhancers_df)){
    enhancer <- enhancers_df[row,]
    print(enhancer)
    
    p <- ArchR_genome_browser_plot(ArchR, enhancer, group_by = "transferred_scHelper_cell_type", 
                                   pal, extend_by = 5000, show_peaks = FALSE)
    
    # plot
    grid::grid.newpage()
    png(paste0(plot_path, gene, '_enhancer_', row, '_browser_plot.png'), height = 15, width = 18, units = 'cm', res = 400)
    grid::grid.draw(p)
    graphics.off()
  }
  
}

plot_path = "./output/NF-downstream_analysis/Processing/ss8/New_enhancers/plots/SEACell_scHelper_cell_type/"
dir.create(plot_path, recursive = T)

pal <- scHelper_cell_type_colours[unique(ArchR$SEACell_scHelper_cell_type)]

# loop through genes and then loop through enhancers within each gene
for (gene in names(enhancers_df_list)){
  enhancers_df <- enhancers_df_list[[gene]]
  
  for (row in 1:nrow(enhancers_df)){
    enhancer <- enhancers_df[row,]
    print(enhancer)
    
    p <- ArchR_genome_browser_plot(ArchR, enhancer, group_by = "SEACell_scHelper_cell_type", 
                                   pal, extend_by = 5000, show_peaks = FALSE)
    
    # plot
    grid::grid.newpage()
    png(paste0(plot_path, gene, '_enhancer_', row, '_browser_plot.png'), height = 15, width = 18, units = 'cm', res = 400)
    grid::grid.draw(p)
    graphics.off()
  }
  
}

plot_path = "./output/NF-downstream_analysis/Processing/ss8/New_enhancers/plots/SEACell_scHelper_cell_type_broad/"

dir.create(plot_path, recursive = T)

pal <- scHelper_cell_type_colours[unique(ArchR$SEACell_scHelper_cell_type_broad)]

# loop through genes and then loop through enhancers within each gene
for (gene in names(enhancers_df_list)){
  enhancers_df <- enhancers_df_list[[gene]]
  
  for (row in 1:nrow(enhancers_df)){
    enhancer <- enhancers_df[row,]
    print(enhancer)
    
    p <- ArchR_genome_browser_plot(ArchR, enhancer, group_by = "SEACell_scHelper_cell_type_broad", 
                                   pal, extend_by = 5000, show_peaks = FALSE)
    
    # plot
    grid::grid.newpage()
    png(paste0(plot_path, gene, '_enhancer_', row, '_browser_plot.png'), height = 15, width = 18, units = 'cm', res = 400)
    grid::grid.draw(p)
    graphics.off()
  }
  
}




