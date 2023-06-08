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
    chrs = c("chr21", "chr22")
    
    plot_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_diff_interactions/plots/"
    rds_path = "./output/NF-hichip-downstream/AllSamples/HicDCPlus_diff_interactions/rds_files/"
    
    ## Created folder with all HiCDC+ outputs to test interactively
    data_path = "./local_test_data/all_HiCDC_outputs/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # chrs = NULL
    # chrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr7", "chr9",
    #          "chr11", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
    #          "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28",
    #          "chr31", "chr32", "chr33", "chrZ", "chrW")
    # chr10, chr6, chr8, chr12 all fail with error:
    #  fitType failed to fit for chr8. Overriding with gene estimates
    #   Error in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  : 
    #     newsplit: out of vertex space
    #   Calls: hicdcdiff ... estimateDispersionsFit -> localDispersionFit -> locfit -> lfproc
    #   Execution halted
    # now gives this error:
    # Error in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  : 
    # newsplit: out of vertex space
    # Calls: hicdcdiff ... estimateDispersionsFit -> localDispersionFit -> locfit -> lfproc
    # Execution halted
    chrs = c("chr1", "chr2", "chr3", "chr4", "chr5")
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}


######################################   Extract file names   ####################################################

print("reading in loops...")

# read in all file paths
all_files <- list.files(data_path, full.names = TRUE)

##### read in unfiltered loops #####
files <-  grep("_HiCDC_output.txt.gz", all_files, value = TRUE)
print("Unfiltered outputs:")
print(files)

# read in loops for all WE samples
print("WE samples:")
WE_samples_unfiltered <-  grep("WE", files, value = TRUE)
print(WE_samples_unfiltered)

# read in loops for all NF samples
print("NF samples:")
NF_samples_unfiltered <-  grep("NF", files, value = TRUE)
print(NF_samples_unfiltered)

##### read in filtered loops #####
filtered_files <- grep("_HiCDC_output_filtered.txt", all_files, value = TRUE)
print("Filtered files:")
print(filtered_files)

# read in loops for filtered WE samples
print("WE samples:")
WE_samples <-  grep("WE", filtered_files, value = TRUE)
print(WE_samples)

# read in loops for filtered NF samples
print("NF samples:")
NF_samples <-  grep("NF", filtered_files, value = TRUE)
print(NF_samples)


############################################################################################################
################################## Find consensus loops #################################################
############################################################################################################

print("Finding consensus loops...")

######################   Read in filtered outputs   ######################################

# read in filtered loops
hidc_filtered_outputs <- list()
for (file in filtered_files) {
  # extract sample name
  sample_name <- gsub(pattern = "_HiCDC_output_filtered.txt", replacement = "", x = basename(file))
  # read in hicDC+ output
  output <- data.table::fread(file)
  # add an interaction name to the hicDC+ output
  output <- output %>% mutate(interaction_ID = paste0(chrI, "-", startI, "-", startJ))
  # add the hicDC+ output to the list of outputs
  hidc_filtered_outputs[[sample_name]] <- output
}

print("filtered loops read in!")

######################   How many sig interactions in each sample   ######################################

# Table of how many interactions found in each sample
interaction_numbers <- data.frame(
  NF_HiChip_r1 = length(hidc_filtered_outputs[[1]]$interaction_ID),
  NF_HiChip_r2 = length(hidc_filtered_outputs[[2]]$interaction_ID),
  NF_HiChip_r3 = length(hidc_filtered_outputs[[3]]$interaction_ID),
  WE_HiChip_r1 = length(hidc_filtered_outputs[[4]]$interaction_ID),
  WE_HiChip_r2 = length(hidc_filtered_outputs[[5]]$interaction_ID),
  WE_HiChip_r3 = length(hidc_filtered_outputs[[6]]$interaction_ID)
)

png(paste0(plot_path, 'how_many_interaction_in_each_sample_table.png'), height = 10, width = 30, units = 'cm', res = 400)
grid.arrange(tableGrob(interaction_numbers, rows=NULL, theme = ttheme_minimal()))
graphics.off()

print("Table of how many interactions found in each sample saved!")

######################   Number of overlapping interactions in replicates   ######################################

# Venn diagram of interactions between WE samples
venn.diagram(
  x = list(hidc_filtered_outputs[[4]]$interaction_ID, hidc_filtered_outputs[[5]]$interaction_ID, hidc_filtered_outputs[[6]]$interaction_ID),
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
  x = list(hidc_filtered_outputs[[1]]$interaction_ID, hidc_filtered_outputs[[2]]$interaction_ID, hidc_filtered_outputs[[3]]$interaction_ID),
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

print("Venn diagrams of interactions in WE and NF samples saved!")

######################   Extract consistent interactions from NF samples   ######################################
# consistent interaction = present in at least 2 of the 3 replicates
r1 <- hidc_filtered_outputs[[1]]$interaction_ID
r2 <- hidc_filtered_outputs[[2]]$interaction_ID
r3 <- hidc_filtered_outputs[[3]]$interaction_ID

common_interactions <- unique(c(intersect(r1, r2), intersect(r1, r3), intersect(r2, r3)))
length(common_interactions)

# extract table of these interactions (p vals not important but chrom start end useful)
r1_subset <- hidc_filtered_outputs[[1]] %>% filter(interaction_ID %in% common_interactions)
remaining_interactions <- common_interactions[!common_interactions %in% r1]
length(remaining_interactions) == length(common_interactions) - sum(common_interactions %in% r1)
r2_subset <- hidc_filtered_outputs[[2]] %>% filter(interaction_ID %in% remaining_interactions)
nrow(r2_subset) == length(remaining_interactions)

common_interactions_table <- rbind(r1_subset, r2_subset)
nrow(common_interactions_table) == length(common_interactions)

# save common interactions table
write.csv(common_interactions_table, file = paste0(rds_path, 'Consensus_interactions.txt'))

print("Table of common interactions saved!")

############################################################################################################
################################## Find differential loops #################################################
############################################################################################################

print("Finding differential loops...")

######################   Create index file of interactions to test   ######################################
# index file = union of significant interactions, chr, startI, startJ
# use the filtered file as these interactions have a qval < 0.05

indexfile <- data.frame()
hidc_outputs <- list()

for (file in filtered_files) {
  # extract sample name
  sample_name <- gsub(pattern = "_HiCDC_output_filtered.txt", replacement = "", x = basename(file))
  print(sample_name)
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

## try running with all interactions - error:
# Error in data.table::fread(filepath) : 
#   R character strings are limited to 2^31-1 bytes
# Calls: hicdcdiff ... .readDiffInputFiles -> gi_list_read -> input.file.read -> <Anonymous> 
# hicdcdiff(input_paths = list(WE = WE_samples_unfiltered, NF = NF_samples_unfiltered),
#           filter_file = paste0(rds_path,'/Indexfile.txt.gz'),
#           output_path=paste0(rds_path,'/HicDCDiff_output/'),
#           # Bins
#           bin_type = 'Bins-uniform',
#           binsize = 5000,
#           # Interaction sizes
#           Dmin = opt$Dmin,
#           Dmax = opt$Dmax,
#           # fitType in DESeq2::estimateDispersions
#           fitType = 'local', # 'parametric' (parametric regression),'local' (local regression), and 'mean' (constant across interaction bins). Default is 'local'.
#           # What objects are made
#           diagnostics=TRUE, # generates diagnostic plots of normalisation factors and MA plots
#           DESeq.save = TRUE, # saves DESEq objects for each chromosome
#           chrs = chrs
#           )

# try running with filtered samples
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
          chrs = chrs
          )
