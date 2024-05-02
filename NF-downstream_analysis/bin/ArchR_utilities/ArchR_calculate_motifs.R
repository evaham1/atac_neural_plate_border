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

plot_path = "./plots/ChromVar/"
dir.create(plot_path, recursive = T)

# background peaks used to compute motif deviations
ArchR <- addBgdPeaks(ArchR)

# compute per-cell deviations across all motif annotations
ArchR <- addDeviationsMatrix(ArchR, peakAnnotation = "Motif", force = TRUE)

# save results as df
df <- getVarDeviations(ArchR, name = "MotifMatrix", plot = FALSE)
write.csv(df, file = paste0(csv_path, "Chromvar_df.csv"), row.names = FALSE)

# plot results
plotVarDev <- getVarDeviations(ArchR, name = "MotifMatrix", plot = TRUE)
png(paste0(plot_path, 'ChromVar_results.png'), height = 15, width = 15, units = 'cm', res = 400)
plotVarDev
graphics.off()

print("Chromvar scores calculated!")

############################## Save object #######################################
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)
print("ArchR object saved")

############################## Extract target peaks of key TFs #######################################

# can see where the motif x peak matrix is stored here
anno <- getPeakAnnotation(ArchR, name = "Motif")
anno$Matches

# need to read in this rds file to extract peak annotations with motifs - dont know if this will work in nextflow...
#motifs_peaks <- readRDS("~/output/NF-downstream_analysis/Processing/ss8/ARCHR_INTEGRATING_WF/Single_cell_integration_cluster_identification/rds_files/ss8_Save-ArchR/Annotations/Motif-Matches-In-Peaks.rds")
motifs_peaks <- readRDS(anno$Matches)
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

# now write out the targets into a txt file
for (TF in names(TF_targets)) {
  print(TF)
  write(paste0(TF, "; ", paste0(TF_targets[[TF]], collapse = ", ")), 
          file = paste0(csv_path, "/", "TF_targets.txt"), 
          append = TRUE)
}

print("TF peak targets written to file!")