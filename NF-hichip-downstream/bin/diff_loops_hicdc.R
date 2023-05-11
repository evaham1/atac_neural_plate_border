#!/usr/bin/env Rscript

### script to use HiCDCPlus to find differential significant loops from HiChip data
print("script to use HiCDCPlus to find differential significant loops from HiChip data")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(parallel)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(HiCDCPlus)
library(DESeq2)
library(VennDiagram)
library(gridExtra)
library(grid)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-m", "--Dmin"), action = "store", type = "integer", help = "Minimum distance (included) to check for significant interactions"),
  make_option(c("-n", "--Dmax"), action = "store", type = "integer", help = "Maximum distance (included) to check for significant interactions"),
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
    
    plot_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_diff_interactions/plots/"
    rds_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_diff_interactions/rds_files/"
    
    ## Created folder with all HiCDC+ outputs to test interactively
    data_path = "./local_test_data/all_HiCDC_outputs/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################################################################################################
################################## Find Differential loops #################################################
############################################################################################################

######################################   Extract file names   ####################################################

print("reading in loops...")

# read in all file paths
print("All samples:")
all_files <- list.files(data_path, full.names = TRUE)
files <-  grep("_HiCDC_output.txt.gz", all_files, value = TRUE)
print(files)

# read in loops for all WE samples
print("WE samples:")
WE_samples <-  grep("WE", files, value = TRUE)
print(WE_samples)

# read in loops for all NF samples
print("NF samples:")
NF_samples <-  grep("NF", files, value = TRUE)
print(NF_samples)


##############################   Read in data and create index file for hicdcdiff  ##############################################
# index file = union of significant interactions, chr, startI, startJ

indexfile <- data.frame()
hidc_outputs <- list()
for (file in files) {
  # extract sample name
  sample_name <- gsub(pattern = "_HiCDC_output_filtered.txt", replacement = "", x = basename(file))
  # read in hicDC+ output
  output <- data.table::fread(file)
  # add unique interactions to indexfile
  indexfile <- unique(rbind(indexfile, output[,c('chrI','startI','startJ')]))
  # add an interaction name to the hicDC+ output
  output <- output %>% mutate(interaction_ID = paste0(chrI, "-", startI, "-", startJ))
  # add the hicDC+ output to the list of outputs
  hidc_outputs[[sample_name]] <- output
}

print(head(indexfile))
dim(indexfile)

#save index file
colnames(indexfile) <- c('chr','startI','startJ')
data.table::fwrite(indexfile,
                   paste0(rds_path,'/Indexfile.txt.gz'),
                   sep='\t',row.names=FALSE,quote=FALSE)

####################################   Run hicdcdiff   ###################################################
# Differential analysis using modified DESeq2 (see ?hicdcdiff)

hicdcdiff(input_paths = list(WE = WE_samples, NF = NF_samples),
          filter_file = paste0(rds_path,'/Indexfile.txt.gz'),
          output_path=paste0(rds_path,'/HicDCDiff_output/'),
          # Bins
          bin_type = 'Bins-uniform',
          binsize = 5000,
          # Interaction sizes
          Dmin = opt$Dmin,
          Dmax = opt$Dmax,
          # fitType in DESeq2::estimateDispersions
          fitType = 'local', # 'parametric' (parametric regression),'local' (local regression), and 'mean' (constant across interaction bins). Default is 'local'.
          # What objects are made
          diagnostics=TRUE, # generates diagnostic plots of normalisation factors and MA plots
          DESeq.save = TRUE, # saves DESEq objects for each chromosome
          chrs = c('chr10')
          )

#### LOOK FOR CONSISTENCY BETWEEN SAMPLES, PCA PLOT??? - venn diagram of all interactions across 6 samples

#### LOOK IF DIFF INTERACTIONS ARE FOUND IN AT LEAST 2 OF 3 REPLICATES

######################   Look at how much overlap there is between samples   ######################################

names(hidc_outputs)

# Table of how many interactions found in each sample

interaction_numbers <- data.frame(
  NF_HiChip_r1 = length(hidc_outputs[[1]]$interaction_ID),
  NF_HiChip_r2 = length(hidc_outputs[[2]]$interaction_ID),
  NF_HiChip_r3 = length(hidc_outputs[[3]]$interaction_ID),
  WE_HiChip_r1 = length(hidc_outputs[[4]]$interaction_ID),
  WE_HiChip_r2 = length(hidc_outputs[[5]]$interaction_ID),
  WE_HiChip_r3 = length(hidc_outputs[[6]]$interaction_ID)
)

png(paste0(plot_path, 'how_many_interaction_in_each_sample_table.png'), height = 10, width = 30, units = 'cm', res = 400)
grid.arrange(tableGrob(interaction_numbers, rows=NULL, theme = ttheme_minimal()))
graphics.off()


# Venn diagram of interactions between WE samples
venn.diagram(
  x = list(hidc_outputs[[4]]$interaction_ID, hidc_outputs[[5]]$interaction_ID, hidc_outputs[[6]]$interaction_ID),
  category.names = c("WE_HiChip_r1", "WE_HiChip_r2", "WE_HiChip_r3"),
  filename = paste0(plot_path, 'Venn_WE_interactions.png'),
  output=TRUE, disable.logging = TRUE,
  # Output features
  imagetype="png",
  height = 600, 
  width = 600, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

# Venn diagram of interactions between NF samples
venn.diagram(
  x = list(hidc_outputs[[1]]$interaction_ID, hidc_outputs[[2]]$interaction_ID, hidc_outputs[[3]]$interaction_ID),
  category.names = c("NF_HiChip_r1", "NF_HiChip_r2", "NF_HiChip_r3"),
  filename = paste0(plot_path, 'Venn_NF_interactions.png'),
  output=TRUE, disable.logging = TRUE,
  # Output features
  imagetype="png",
  height = 600, 
  width = 600, 
  resolution = 300,
  compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)







