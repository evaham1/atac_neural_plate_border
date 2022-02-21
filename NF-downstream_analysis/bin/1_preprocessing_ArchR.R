#!/usr/bin/env Rscript

### Script to preprocess in ArchR
print("1_preprocessing_ArchR")

############################## Load libraries #######################################
library(getopt)
library(Signac)
library(Seurat)
library(future)
library(tidyverse)
library(grid)
library(gridExtra)
library(clustree)
library(GenomeInfoDb)
library(ggplot2)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()

#BiocManager::install("GenomicFeatures")
library(GenomicFeatures)

############################## Set up script options #######################################
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    setwd("~/NF-downstream_analysis")
    ncores = 8
    
    plot_path = "../output/NF-downstream_analysis/1_preprocessing_ArchR/plots/"
    rds_path = "../output/NF-downstream_analysis/1_preprocessing_ArchR/rds_files/"
    data_path = "../output/NF-luslab_sc_multiomic/full/cellranger_atac_output/"
    ref_path = "../output/NF-luslab_sc_multiomic/reference/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/cellranger_atac_output/"
    ref_path = "./input/galgal/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 155* 1024^3)
    addArchRThreads(threads = ncores) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Read in data and set up ArchR object #######################################

###   Make gene annotation
# make TxDb file from gtf
# https://seandavi.github.io/ITR/transcriptdb.html
txdb <- makeTxDbFromGFF(paste0(ref_path, "genes.gtf.gz"))
print("txdb made")
genes(txdb)
# download OrgDb for chick from Bioconductor (how do I know this is right?)
if (!requireNamespace("org.Gg.eg.db", quietly = TRUE)){
  BiocManager::install("org.Gg.eg.db")
}
library(org.Gg.eg.db)
# combine TxDB and OrgDB to make gene annotation
geneAnnotation <- createGeneAnnotation(TxDb = txdb, OrgDb = org.Gg.eg.db)
print("gene annotation:")
geneAnnotation

###   Make genome annotation


print("genome annotation:")
genomeAnnotation

###   Read in fragment files to create Arrow files
# Make dataframe with stage and replicate info extracted from path
paths <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)
input <- data.frame(sample = sub('.*/', '', paths), 
                   matrix_path = paste0(paths, "/outs/filtered_peak_bc_matrix.h5"),
                   metadata_path = paste0(paths, "/outs/singlecell.csv"),
                   fragments_path = paste0(paths, "/outs/fragments.tsv.gz"))
fragments_list <- input$fragments_path
names(fragments_list) <- input$sample
print("path df made")

# create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = fragments_list,
  sampleNames = names(fragments_list),
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  filterTSS = 2, #Dont set this too high because you can always increase later
  filterFrags = 500, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
print("Arrow files:")
ArrowFiles


# create ArchR Project
ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ArchR",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
print("ArchR Project:")
ArchR

# save ArchR project
saveArchRProject(ArchRProj = ArchR, outputDirectory = "Save-ArchR", load = FALSE)
