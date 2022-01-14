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

# add this to dockerfile
#devtools::install_github('alexthiery/scHelper@v0.2.4', dependencies = TRUE, force = TRUE)
#library(scHelper)
install.packages("rtracklayer")
library(rtracklayer)

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
      plot_path = "../output/NF-downstream_analysis/TEST/plots/"
      rds_path = "../output/NF-downstream_analysis/TEST/rds_files/"
      data_path = "../output/NF-luslab_sc_multiomic/cellranger_atac_output_test/"
      }else{
      plot_path = "../output/NF-downstream_analysis/1_preprocessing/plots/"
      rds_path = "../output/NF-downstream_analysis/1_preprocessing/rds_files/"
      data_path = "../output/NF-luslab_sc_multiomic/cellranger_atac_output/"}
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/cellranger_atac_output/"
    ref_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 16gb
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Read in data and set up Signac object #######################################

# Make dataframe with stage and replicate info extracted from path
paths <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)

input <- data.frame(sample = sub('.*/', '', paths), 
                   matrix_path = paste0(paths, "/outs/filtered_peak_bc_matrix.h5"),
                   metadata_path = paste0(paths, "/outs/singlecell.csv"),
                   fragments_path = paste0(paths, "/outs/fragments.tsv.gz"))

# Read in the 3 files needed in list format
counts_list <- apply(input, 1, function(x) Read10X_h5(filename = x[["matrix_path"]]))
metadata_list <- apply(input, 1, function(x) read.csv(file = x[["metadata_path"]], header = TRUE, row.names = 1))
fragments_list <- as.list(input$fragments_path)

# Build list of assays using these files
chrom_assays <- lapply(1:nrow(input), function(x) CreateChromatinAssay(
                                                    counts = counts_list[[x]],
                                                    sep = c(":", "-"),
                                                    fragments = fragments_list[[x]],
                                                    min.cells = 10,
                                                    min.features = 200))
signac_datas <- lapply(1:nrow(input), function(x) CreateSeuratObject(
  counts = chrom_assays[[x]],
  project = input$sample[x],
  assay = "peaks",
  meta.data = metadata_list[[x]]))

# add annotations using chick gtf
gtf <- rtracklayer::import(paste0(ref_path, "genes.gtf.gz"))
gene.coords <- gtf[gtf$type == 'gene']

signac_list <- lapply(signac_datas, function(x) SetAssayData(x, slot = "annotation", new.data = gene.coords))

# Init list of signac objects then merge
names(signac_list) <- input$sample
seurat_all <- merge(x = signac_list[[1]], y=signac_list[-1], add.cell.ids = names(signac_list), project = "chick.10x.atac")

# Add metadata col for stage and flow cell
seurat_all@meta.data[["stage"]] <- substr(seurat_all@meta.data$orig.ident, 1, 3)
seurat_all@meta.data[["flow_cell"]] <- substr(seurat_all@meta.data$orig.ident, 5, 5)

# Convert metadata character cols to factors
seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)] <- lapply(seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)], as.factor)

# to test: save RDS
saveRDS(seurat_all, paste0(rds_path, "seurat_all.RDS"), compress = FALSE)

############################## Try low/med/high filtering thresholds to QC #######################################
seurat_all <- NucleosomeSignal(object = seurat_all)
seurat_all <- TSSEnrichment(object = seurat_all, fast = FALSE)
seurat_all$pct_reads_in_peaks <- seurat_all$peak_region_fragments / seurat_all$passed_filters * 100

# make dataframe with different filtering parameters which can be put into a loop for carrying out downstream analysis
# this doesnt work if strictest threshold has no cells in it? how do I choose these numbers?
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
