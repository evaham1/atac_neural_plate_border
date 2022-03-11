#!/usr/bin/env Rscript

### Script to preprocess in ArchR
print("1_preprocessing_ArchR")

############################## Load libraries #######################################
library(getopt)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(future)

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
    
    plot_path = "../output/NF-downstream_analysis/1_ArchR_preprocessing/plots/"
    rds_path = "../output/NF-downstream_analysis/1_ArchR_preprocessing/rds_files/"
    data_path = "../output/NF-luslab_sc_multiomic/test/cellranger_atac_output/"
    ref_path = "../output/NF-luslab_sc_multiomic/reference/"

    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/cellranger_atac_output/"
    ref_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    #plan("multicore", workers = ncores)
    #options(future.globals.maxSize = 155* 1024^3)
    addArchRThreads(threads = 1) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############# Get parallel computing to work
# library(datasets)

# # without parallel 
# x <- iris[which(iris[,5] != "setosa"), c(1,5)]
# trials <- 10000
# res <- data.frame()
# system.time({
#   trial <- 1
#   while(trial <= trials) {
#     ind <- sample(100, 100, replace=TRUE)
#     result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
#     r <- coefficients(result1)
#     res <- rbind(res, r)
#     trial <- trial + 1
#   }
# })

# # with parallel 
# x <- iris[which(iris[,5] != "setosa"), c(1,5)]
# trials <- seq(1, 10000)
# boot_fx <- function(trial) {
#   ind <- sample(100, 100, replace=TRUE)
#   result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
#   r <- coefficients(result1)
#   res <- rbind(data.frame(), r)
# }
# system.time({
#   results <- lapply(trials, boot_fx)
# })

############################## Set up Annotation files - will need to revisit #######################################

###   Make gene annotation
# make TxDb file from gtf
# https://seandavi.github.io/ITR/transcriptdb.html

# changed gtf to include 'chr' in each chromosome name
txdb <- makeTxDbFromGFF(paste0(ref_path, "temp.gtf"))
print("txdb made")
genes(txdb)
# download OrgDb for chick from Bioconductor (how do I know this is right?)
# if (!requireNamespace("org.Gg.eg.db", quietly = TRUE)){
#   BiocManager::install("org.Gg.eg.db")
# }
library(org.Gg.eg.db)
# combine TxDB and OrgDB to make gene annotation
geneAnnotation <- createGeneAnnotation(TxDb = txdb, OrgDb = org.Gg.eg.db)
print("gene annotation:")
geneAnnotation

###   Make genome annotation - THIS IS USING UCSC GENOME WILL NEED TO CHANGE
# if (!requireNamespace("BSgenome.Ggallus.UCSC.galGal6", quietly = TRUE)){
#   BiocManager::install("BSgenome.Ggallus.UCSC.galGal6")
# }
library(BSgenome.Ggallus.UCSC.galGal6)
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Ggallus.UCSC.galGal6)

print("genome annotation:")
genomeAnnotation

############################## Fragment files -> Arrow Files #######################################

# Make dataframe with stage and replicate info extracted from path
paths <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)
input <- data.frame(sample = sub('.*/', '', paths), 
                   matrix_path = paste0(paths, "/outs/filtered_peak_bc_matrix.h5"),
                   metadata_path = paste0(paths, "/outs/singlecell.csv"),
                   fragments_path = paste0(paths, "/outs/fragments.tsv.gz"))
fragments_list <- input$fragments_path
names(fragments_list) <- input$sample
print("path df made")

print(fragments_list)

## try creating Arrow files in parallel
CreateOneArrowFile <- function(data){
  createArrowFiles(
    inputFiles = data,
    sampleNames = "ArrowFile",
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    logFile = paste0(plot_path, "createArrows",
    minTSS = 4,
    minFrags = 1000,
    maxFrags = 1e+06))
}

ArrowFiles <- lapply(setNames(fragments_list, fragments_list), CreateOneArrowFile)


# create arrow files - keep thresholds as unrestrictive as possible at this point
ArrowFiles <- createArrowFiles(
  inputFiles = fragments_list,
  sampleNames = names(fragments_list),
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  logFile = paste0(plot_path, "createArrows",
  minTSS = 4,
  minFrags = 1000,
  maxFrags = 1e+06)
)
print("arrow files made")
print("Arrow files:")
ArrowFiles

############################## Create ArchR Project and save #######################################

ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  outputDirectory = paste0(plot_path, "ArchR"),
  copyArrows = FALSE
)
print("ArchR Project:")
ArchR
getAvailableMatrices(ArchR)

# check how much memory
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# add stage to metadata
stage <- substr(ArchR$Sample, 1, 3)
ArchR$stage <- stage

# save ArchR project
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, "Save-ArchR"), load = FALSE)
