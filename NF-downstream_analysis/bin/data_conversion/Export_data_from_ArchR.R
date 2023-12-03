#!/usr/bin/env Rscript

print("extracting data from ArchR to run SEACells")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(parallel)

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
    addArchRThreads(threads = 1) 
    
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Split_stages/rds_files/" # has all stages
    rds_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_ATAC_WF/0_ATAC_exported_data/exported_ArchR_data/" # for ss8

  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    data_path = "./input/rds_files/"
    rds_path = "./exported_ArchR_data/"
    
    ncores = opt$cores
    #addArchRThreads(threads = ncores)
    addArchRThreads(threads = 1)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(rds_path, recursive = T)
}

set.seed(42)


#################### R script to extract peak counts from ArchR object #########################
# https://github.com/dpeerlab/SEACells/blob/main/notebooks/ArchR/ArchR-preprocessing-nfr-peaks.R

############################## Read in ArchR project #######################################

# If files are not in rds_files subdirectory look in input dir 
label <- sub('_.*', '', list.files(data_path))
print(label) 

# Interactive:
# label <- label[5]

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

############################## Check data #######################################

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getPeakSet(ArchR)
getAvailableMatrices(ArchR)

############################## Export data #######################################

# Dim reduction
write.csv(getReducedDims(ArchR), paste0(rds_path, "svd.csv"), quote=FALSE)
print("Dim reduction exported")

# Cell metadata
write.csv(getCellColData(ArchR), paste0(rds_path, "cell_metadata.csv"), quote=FALSE)
print("Cell metadata exported")

# Gene scores
gene.scores <- getMatrixFromProject(ArchR)
scores <- assays(gene.scores)[['GeneScoreMatrix']]
scores <- as.matrix(scores)
rownames(scores) <- rowData(gene.scores)$name
write.csv(scores, paste0(rds_path, "gene_scores.csv"), quote=FALSE)
print("Gene scores exported")

# Peak counts
dir.create(paste0(rds_path, "peak_counts/"), recursive = T)
peaks <- getPeakSet(ArchR)
peak.counts <- getMatrixFromProject(ArchR, 'PeakMatrix')

#     Reorder peaks by chromosome
chr_order <- sort(seqlevels(peaks))
reordered_features <- list()
for(chr in chr_order) {
  reordered_features[[chr]] = peaks[seqnames(peaks) == chr] }
reordered_features <- Reduce("c", reordered_features)    

#     Export counts
counts <- assays(peak.counts)[['PeakMatrix']]
writeMM(counts, paste0(rds_path, "peak_counts/counts.mtx"))
print("Peak counts matrix exported")
write.csv(colnames(peak.counts), paste0(rds_path, "peak_counts/cells.csv"))
print("Cell names exported")
names(reordered_features) <- paste0("PeakId_", 1:length(reordered_features))
write.csv(as.data.frame(reordered_features), paste0(rds_path, 'peak_counts/peaks.csv'))
print("Peak names exported")