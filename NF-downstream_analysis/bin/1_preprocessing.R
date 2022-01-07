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


############################## Set up script options #######################################
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data - NEED TO CHANGE TO WORK DIRS FOR INPUT PATHS
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')

    plot_path = "../output/NF-downstream_analysis/1_preprocessing/plots/"
    rds_path = "../output/NF-downstream_analysis/1_preprocessing/rds_files/"
    data_path = "./input/cellranger_atac_output/"
    ncores = 8

  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')

    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/cellranger_atac_output/"
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
                   matrix_path = paste0(paths, "/filtered_peak_bc_matrix.h5"),
                   metadata_path = paste0(paths, "/singlecell.csv"),
                   fragments_path = paste0(paths, "/fragments.tsv.gz"))

# # Init list of seurat objects then merge
# seurat_list <- apply(input, 1, function(x) CreateSeuratObject(counts= Read10X(data.dir = x[["path"]]), project = x[["sample"]]))
# names(seurat_list) <- input$sample
# seurat_all <- merge(x = seurat_list[[1]], y=seurat_list[-1], add.cell.ids = names(seurat_list), project = "chick.10x")

# # Add metadata col for seq run
# seurat_all@meta.data[["run"]] <- gsub(".*-", "", as.character(seurat_all@meta.data$orig.ident))
# seurat_all@meta.data[["stage"]] <- gsub("-.*", "", as.character(seurat_all@meta.data$orig.ident))

# # Convert metadata character cols to factors
# seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)] <- lapply(seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)], as.factor)

# # Make seurat gene annotation dataframe and save
# annotations <- read.table(paste0(input[1,'path'], '/features.tsv.gz'), col.names = c('Accession', 'Gene', 'V3', 'V4'))[,1:2]
# # make gene names unique in annotations dataframe in order to match seurat annotations
# annotations$Gene <- make.unique(annotations$Gene)
# # Save annnotation dataframe
# write.table(annotations, 'seurat_annotations.csv', row.names=FALSE, quote=FALSE, sep=',')

## hh7_1 test data
#filtered_peak_bc_matrix = "../work/98/4887a345fb1f0935b6e5a404f3683c/hh7_1_cellranger_atac/outs/filtered_peak_bc_matrix.h5"
#singlecell = "../work/98/4887a345fb1f0935b6e5a404f3683c/hh7_1_cellranger_atac/outs/singlecell.csv"
#fragments = "../work/98/4887a345fb1f0935b6e5a404f3683c/hh7_1_cellranger_atac/outs/fragments.tsv.gz"

#counts <- Read10X_h5(filename = filtered_peak_bc_matrix)
# metadata <- read.csv(
#   file = singlecell,
#   header = TRUE,
#   row.names = 1
# )

# chrom_assay <- CreateChromatinAssay(
#   counts = counts,
#   sep = c(":", "-"),
#   fragments = fragments,
#   min.cells = 10,
#   min.features = 200
# )

# signac_data <- CreateSeuratObject(
#   counts = chrom_assay,
#   assay = "peaks",
#   meta.data = metadata
# )

# add annotations
# gtf <- rtracklayer::import('../work/1b/b33a424af93ffb7603d5a968171b87/REFERENCE/genes/genes.gtf.gz')
# gene.coords <- gtf[gtf$type == 'gene']
# Annotation(signac_data) <- gene.coords

# to test: save RDS
saveRDS(metadata, paste0(rds_path, "TEST.RDS"), compress = FALSE)

# ############################## Try low/med/high filtering thresholds to QC #######################################
# signac_data <- NucleosomeSignal(object = signac_data)
# signac_data <- TSSEnrichment(object = signac_data, fast = FALSE)
# signac_data$pct_reads_in_peaks <- signac_data$peak_region_fragments / signac_data$passed_filters * 100

# # make dataframe with different filtering parameters which can be put into a loop for carrying out downstream analysis
# # this doesnt work if strictest threshold has no cells in it? how do I choose these numbers?
# filter_thresholds <- data.frame(pct_reads_in_peaks = c(0, 40, 50, 60), TSS.enrichment = c(0, 0.3, 2.5, 3), nucleosome_signal = c(Inf, 1.5, 0.8, 0.6), row.names = c("unfilt", "low", "med", "high"))

# # Calculate remaining cells following different filter thresholds
# filter_qc <- lapply(rownames(filter_thresholds), function(condition){
#   signac_data@meta.data %>%
#     filter(pct_reads_in_peaks > filter_thresholds[condition,'pct_reads_in_peaks']) %>%
#     filter(TSS.enrichment > filter_thresholds[condition,'TSS.enrichment']) %>%
#     filter(nucleosome_signal < filter_thresholds[condition,'nucleosome_signal']) %>%
#     group_by(orig.ident) %>%
#     tally() %>%
#     dplyr::rename(!!condition := n)
# })

# filter_qc <- Reduce(function(x, y) merge(x, y), filter_qc)

# # Plot remaining cell counts
# filter_qc <-  filter_qc %>% column_to_rownames('orig.ident')
# filter_qc <- rbind(filter_qc, Total = colSums(filter_qc)) %>% rownames_to_column("orig.ident")

# png(paste0(plot_path, 'remaining_cell_table.png'), height = 10, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(filter_qc, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# png(paste0(plot_path, 'remaining_cell_bar.png'), height = 15, width = 21, units = 'cm', res = 400)
# ggplot(filter_qc[filter_qc$orig.ident != "Total",] %>% reshape2::melt(), aes(x=variable, y=value, fill=orig.ident)) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   xlab("Filter Condition") +
#   ylab("Cell Count") +
#   ggtitle("Cell count after filtering") +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5))
# graphics.off()
