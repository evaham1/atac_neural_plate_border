#!/usr/bin/env Rscript

### Script to create seurat object, filter data, remove poor quality clusters and integrate across batches

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

# add this to dockerfile
#devtools::install_github('alexthiery/scHelper@v0.2.4', dependencies = TRUE, force = TRUE)
#library(scHelper)


############################## Set up script options #######################################
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)
test = TRUE

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    setwd("~/NF-downstream_analysis")
    ncores = 8
    ref_path = "../output/NF-luslab_sc_multiomic/reference/"
    
    if(test == TRUE){
      plot_path = "../output/NF-downstream_analysis/2_filtering/TEST/plots/"
      rds_path = "../output/NF-downstream_analysis/2_filtering/TEST/rds_files/"
      data_path = "../output/NF-downstream_analysis/1_preprocessing/TEST/rds_files/"
      }else{
      plot_path = "../output/NF-downstream_analysis/2_filtering/plots/"
      rds_path = "../output/NF-downstream_analysis/2_filtering/rds_files/"
      data_path = "../output/NF-downstream_analysis/1_preprocessing/rds_files/"}
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 75* 1024^3)
    plan()
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Read merged Seurat RDS object and make plots #######################################

seurat_all <- readRDS(paste0(data_path, "seurat_all.RDS"))

seurat_all

seurat_all@meta.data

#######   TEST CODE FROM PREVIOUS SCRIPT TO CHECK IT WORKS AND MAKES PLOTS
filter_thresholds <- data.frame(pct_reads_in_peaks = c(0, 40, 50, 60), TSS.enrichment = c(0, 0.3, 2.5, 3), nucleosome_signal = c(Inf, 1.5, 0.8, 0.6), row.names = c("unfilt", "low", "med", "high"))

# Calculate remaining cells following different filter thresholds
filter_qc <- lapply(rownames(filter_thresholds), function(condition){
  seurat_all@meta.data %>%
    filter(pct_reads_in_peaks > filter_thresholds[condition,'pct_reads_in_peaks']) %>%
    filter(TSS.enrichment > filter_thresholds[condition,'TSS.enrichment']) %>%
    filter(nucleosome_signal < filter_thresholds[condition,'nucleosome_signal']) %>%
    group_by(orig.ident) %>%
    tally() %>%
    dplyr::rename(!!condition := n)
})

filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc)

# Plot remaining cell counts
filter_qc <-  filter_qc %>% column_to_rownames('orig.ident')
filter_qc <- rbind(filter_qc, Total = colSums(filter_qc)) %>% rownames_to_column("orig.ident")

png(paste0(plot_path, 'remaining_cell_table.png'), height = 10, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(filter_qc, rows=NULL, theme = ttheme_minimal()))
graphics.off()

png(paste0(plot_path, 'remaining_cell_bar.png'), height = 15, width = 21, units = 'cm', res = 400)
ggplot(filter_qc[filter_qc$orig.ident != "Total",] %>% reshape2::melt(), aes(x=variable, y=value, fill=orig.ident)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  xlab("Filter Condition") +
  ylab("Cell Count") +
  ggtitle("Cell count after filtering") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()
