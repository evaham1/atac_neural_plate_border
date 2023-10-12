#!/usr/bin/env Rscript

print("Pair RNA and ATAC cells to make multiome data, then run ChromVar")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(Seurat)
library(Signac)
library(scMEGA)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Ggallus.UCSC.galGal6)
library(SummarizedExperiment)
library(igraph)
library(ggraph)
library(BiocParallel)

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
    
    # ss8 for faster testing
    rds_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA_GRNi/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA_GRNi/plots/"
    data_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA_integrated/rds_files/"
    
    # full data
    data_path = "./output/NF-downstream_analysis/Processing/FullData/scMEGA/ArchR-to_seurat/" # ATAC seurat object + cell pairings
    data_path = "./output/NF-downstream_analysis/rna_objects/" # RNA seurat objec
    
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/"
    rds_path = "./rds_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

set.seed(42)

############################################################################################
#########################             FUNCTIONS             ################################

AddMotifs_to_Seurat <- function(object, assay, pfm, genome, cutOff = 5e-05){
  # added this myself
  DefaultAssay(object) <- assay
  peak_ranges <- granges(object)
  # this is from github
  motif.matrix <- CreateMotifMatrix( 
    features = peak_ranges, 
    pwm = pfm, 
    genome = genome, 
    use.counts = FALSE,
    p.cutoff = cutOff
  ) 
  message("Finding motif positions") 
  
  # for positions, a list of granges is returned 
  # each element of list is a PFM name 
  # each entry in granges is the position within a feature that matches motif 
  peak_ranges_keep <- as.character(seqnames(x = peak_ranges)) %in% seqlevels(x = genome) 
  motif.positions <- motifmatchr::matchMotifs( 
    pwms = pfm, 
    subject = peak_ranges[peak_ranges_keep], 
    out = 'positions', 
    genome = genome,
    p.cutoff = cutOff
  ) 
  message("Creating Motif object") 
  
  motif <- CreateMotifObject( 
    data = motif.matrix, 
    positions = motif.positions, 
    pwm = pfm 
  ) 
  object <- SetAssayData(object = object, slot = "motifs", 
                         new.data = motif)
  return(object)
}

############################## Read in seurat objects #######################################

print("reading in data...")

# if reading in transfer labels object
obj.rna <- readRDS(paste0(data_path, "seurat_label_transfer_minus_HH4.RDS"))

# read in atac data - in input/rds_files folder
obj.atac <- readRDS(paste0(data_path, "rds_files/ATAC_seurat.RDS"))

# read in cell pairings from archr integration - in input/rds_files folder
df.pair <- read.csv(paste0(data_path, "rds_files/archr_cell_pairings.csv"))
head(df.pair)

# read in gene activity matrix - in input/rds_files folder
gene.activity <- readRDS(paste0(data_path, "rds_files/gene_score_matrix.RDS"))

print("data read in!")

############################## Plot UMAPs of RNA and ATAC data #######################################

print("plotting UMAPs..")

## set cols
# schelper cell type colours
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                "#A5718D", "#3F918C", "#ed5e5f", "9792A3")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                       'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo',
                                       'Neural', 'Placodal', 'Non-neural', 'Contam')

# stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_cols = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_cols) <- stage_order

# set colour palettes for UMAPs
rna_cols <- scHelper_cell_type_colours[as.character(unique(obj.rna$scHelper_cell_type))]
atac_cols <- scHelper_cell_type_colours[as.character(unique(obj.atac$scHelper_cell_type))]

## UMAPs
p1 <- DimPlot(obj.rna, pt.size = 1, reduction = "umap", group.by = "scHelper_cell_type", cols = rna_cols, shuffle = TRUE) +
  ggtitle("scRNA-seq")
p2 <- DimPlot(obj.atac, pt.size = 1, reduction = "umap_iLSI", group.by = "scHelper_cell_type", cols = atac_cols, shuffle = TRUE) +
  ggtitle("snATAC-seq")
png(paste0(plot_path, '0_UMAPs_scHelper_cell_type.png'), height = 10, width = 30, units = 'cm', res = 400)
p1 + p2
graphics.off()

p1 <- DimPlot(obj.rna, pt.size = 1, reduction = "umap", group.by = "stage", cols = stage_cols, shuffle = TRUE) +
  ggtitle("scRNA-seq")
p2 <- DimPlot(obj.atac, pt.size = 1, reduction = "umap_iLSI", group.by = "stage", cols = stage_cols, shuffle = TRUE) +
  ggtitle("snATAC-seq")
png(paste0(plot_path, '0_UMAPs_stage.png'), height = 10, width = 24, units = 'cm', res = 400)
p1 + p2
graphics.off()

print("UMAPs plotted!")

############################## Coembed RNA and ATAC #######################################

print("Coembedding modalities...")

# set pca dim reduction of rna seurat as 'dr' slot so slots are named the same for atac and rna
obj.rna[["dr"]] <- obj.rna[["pca"]]

## coembedding
obj.coembed <- CoembedData(
  obj.rna,
  obj.atac, 
  gene.activity, 
  weight.reduction = "dr", # 'dr' slot holds the iterative LSI for atac and pca for rna
  verbose = FALSE
)

print("Coembedding complete!")
print(obj.coembed)

# don't bother re-run dimensionality reduction or plot UMAPs as they will be crap

############################## Clean up cell pairings #######################################
# use previously identified (by ArchR) RNA-ATAC single cell pairings

print("pairing cells...")

# how many unique ATAC and RNA cells left in paired object
print("cell numbers:")
length(unique(df.pair$ATAC)) # 86,217
length(unique(df.pair$RNA)) # 1,895
dim(df.pair) # 86217     3

# the RNA transfer labels object doesn't include contamination, so 191 RNA cells which are in the df.pair are missing
# so first lets remove these cells from the df.pair object, so there should be no contam cells in the final paired seurat object
RNA_cells_to_remove <- setdiff(unique(df.pair$RNA), Cells(obj.coembed))
length(RNA_cells_to_remove) # 191
filtered_df_pair <- df.pair %>% filter(!RNA %in% RNA_cells_to_remove)
head(filtered_df_pair)
dim(filtered_df_pair) # 81,390   3

############################## Pair object #######################################

# how many unique ATAC and RNA cells left in paired object
print("RNA data")
length(unique(filtered_df_pair$RNA)) # 1,704
sum(unique(filtered_df_pair$RNA) %in% Cells(obj.coembed)) # 1,704
print("ATAC data")
length(unique(filtered_df_pair$ATAC)) # 81,390
sum(unique(filtered_df_pair$ATAC) %in% Cells(obj.coembed)) # 81,390

# only keep paired cells in the seurat object
sel_cells <- c(filtered_df_pair$ATAC, filtered_df_pair$RNA)
coembed.sub2 <- obj.coembed[, sel_cells]
print(coembed.sub2) # 83,094

## create paired object
obj.pair <- CreatePairedObject(df.pair = filtered_df_pair,
                               object = coembed.sub2,
                               use.assay1 = "RNA",
                               use.assay2 = "ATAC")

# see how many cells are left in the paired object
cell_counts <- data.frame(dim(obj.coembed)[2], dim(coembed.sub2)[2], dim(obj.pair)[2])
colnames(cell_counts) <- c("Before removed contam", "Before paired obj", "After paired obj")

png(paste0(plot_path, 'Cell_counts_after_creating_paired_object.png'), height = 10, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

print("cells paired!")

############################## Check object #######################################

print("Checking metadata...")

# check final object has all the correct metadata (including the latent time values)
print(obj.pair)
print(head(obj.pair@meta.data))

############################## Add motif information #######################################

print("Adding motif information...")

# download motif database
motifList <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))

# add motif information to ATAC data

# in-built signac function where I can't change the p val:
# obj.pair <- AddMotifs(
#   object = obj.pair,
#   genome = BSgenome.Ggallus.UCSC.galGal6,
#   pfm = motifList,
#   assay = "ATAC"
# )

# edited function where I can edit p-val cutoff:
obj.pair <- AddMotifs_to_Seurat(obj.pair, assay = "ATAC", pfm = motifList, genome = BSgenome.Ggallus.UCSC.galGal6, cutOff = 1e-05)

############################## Check motif-peak annotations #######################################

print("Checking motif information...")

motif.matching <- obj.pair@assays$ATAC@motifs@data
colnames(motif.matching) <- obj.pair@assays$ATAC@motifs@motif.names
dim(motif.matching) # 271391    746
print(paste0("Total motif-peak hits: ", sum(motif.matching)))

# distribution of motifs by TF
n_hits_per_TF <- colSums(motif.matching)
png(paste0(plot_path, 'Motif_hits_per_TF.png'), height = 8, width = 10, units = 'cm', res = 400)
hist(n_hits_per_TF, breaks = 100)
graphics.off()
print("Number of hits per TF summary stats:")
summary(n_hits_per_TF)

# distribution of peaks by TF
n_hits_per_peak <- rowSums(motif.matching)
png(paste0(plot_path, 'Motif_hits_per_peak.png'), height = 8, width = 10, units = 'cm', res = 400)
hist(n_hits_per_peak, breaks = 100)
graphics.off()
print("Number of hits per peak summary stats:")
summary(n_hits_per_peak)

############################## Run chromvar #######################################

## NB if I try to rename motifList by TF names like I do when running chromvar in ArchR this fails

print("Running chromVar...")

BiocParallel::register(SerialParam())

# run chromvar
obj.pair <- RunChromVAR(
  object = obj.pair,
  genome = BSgenome.Ggallus.UCSC.galGal6,
  assay = 'ATAC'
)

############################## Save #######################################

saveRDS(obj.pair, paste0(rds_path, "paired_object_chromvar.RDS"), compress = FALSE)