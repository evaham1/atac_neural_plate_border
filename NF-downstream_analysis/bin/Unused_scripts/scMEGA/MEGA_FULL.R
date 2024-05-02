#!/usr/bin/env Rscript

print("Run scMEGA on integrated RNA/ATAC object to infer GRN")

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
    data_path = "./input/rds_files/"
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

############################## Read in seurat objects #######################################

print("reading in data...")

# if reading in transfer labels object
obj.rna <- readRDS(paste0(data_path, "seurat_label_transfer_minus_HH4.RDS"))

# read in atac data - in input/rds_files folder
obj.atac <- readRDS(paste0(data_path, "rds_files/ATAC_seurat.RDS"))

# read in cell pairings from archr integration - in input/rds_files folder
df.pair <- read.csv(paste0(data_path, "./rds_files/archr_cell_pairings.csv"))
head(df.pair)

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

############################## Transfer over RNA trajectories #######################################
# when tried recalculating trajectories was crap so instead use RNA velocity latent time and
# lineage probabilities which were calculated on RNA data

## NEED TO SEE WHERE THIS TRAJECTORY COMES FROM IN SEURAT OBJECT

############################## Run chromvar #######################################

## NB if I try to rename motifList by TF names like I do when running chromvar in ArchR this fails

BiocParallel::register(SerialParam())

# run chromvar
obj_chromvar <- RunChromVAR(
  object = obj.motifs,
  genome = BSgenome.Ggallus.UCSC.galGal6,
  assay = 'ATAC'
)

## NEED TO SEE WHERE THIS IS SAVED TO SEE IF I CAN JUST USE CHROMVAR FROM ARCHR

############################## Filter genes and TFs #######################################

# select the TFs that correlate in 'activity' and expression
res <- SelectTFs(object = obj.trajectory, trajectory.name = "Trajectory", return.heatmap = TRUE,
                 cor.cutoff = 0.1)

# save the (top 100?) TFs that correlate
df.cor <- res$tfs
write.csv(df.cor, file = paste0(temp_rds_path, "TF_correlations.csv"), row.names = FALSE)

# plot TF activity dynamics across trajectory
ht <- res$heatmap
png(paste0(temp_plot_path, 'TFs_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
draw(ht)
graphics.off()

# select top 10% genes that vary with the trajectory and correlate with peaks
res <- SelectGenes(object = obj.trajectory,
                   labelTop1 = 0,
                   labelTop2 = 0)

# save the most variable genes
df.p2g <- res$p2g
write.csv(df.p2g, file = paste0(temp_rds_path, "Variable_genes_and_matched_enhancers.csv"), row.names = FALSE)

# plot the dynamics of these genes across trajectory
ht <- res$heatmap
png(paste0(temp_plot_path, 'Genes_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
draw(ht)
graphics.off()

############################## GRN inference #######################################

# GRN inference
tf.gene.cor <- GetTFGeneCorrelation(object = obj.trajectory, 
                                    tf.use = df.cor$tfs, 
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar", 
                                    gene.assay = "RNA",
                                    trajectory.name = "Trajectory")

## plot TF-gene correlation heatmap
ht <- GRNHeatmap(tf.gene.cor, tf.timepoint = df.cor$time_point)
png(paste0(temp_plot_path, 'TF_gene_corr_heatmap.png'), height = 30, width = 45, units = 'cm', res = 400)
ht
graphics.off()

# associate genes with TFs by looking for TF binding sites in their respective linked enhancers
motif.matching <- obj.trajectory@assays$ATAC@motifs@data
colnames(motif.matching) <- obj.trajectory@assays$ATAC@motifs@motif.names
motif.matching <- motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]

# take the TF-gene correlation, peak-TF binding prediction and peaks-to-genes linkage to build network
df.grn <- GetGRN(motif.matching = motif.matching, 
                 df.cor = tf.gene.cor, 
                 df.p2g = df.p2g)

# save network
write.csv(df.grn, file = paste0(temp_rds_path, "GRN_data.csv"), row.names = FALSE)

############################## GRN visualisation #######################################

# define colors for nodes representing TFs (i.e., regulators)
df.cor <- df.cor[order(df.cor$time_point), ]
tfs.timepoint <- df.cor$time_point
names(tfs.timepoint) <- df.cor$tfs

# filter the grn before plotting
df.grn2 <- df.grn %>%
  subset(correlation > 0.8) %>%
  dplyr::select(c(tf, gene, correlation)) %>%
  rename(weights = correlation)

# which TFs are left in the GRN
unique(df.grn2$tf)

# plot GRN, specifying set of TFS
p <- GRNPlot(df.grn2, 
             tfs.use = c("SIX1", "DLX5", "DLX6", "GATA2", "GATA3", "SNAI2", "SOX13", "SOX2", "SOX21",
                         "TFAP2A", "TFAP2B", "TFAP2C", "ZIC1", "ZIC3"),
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 42, 
             plot.importance = TRUE,
             min.importance = 2,
             remove.isolated = FALSE)

options(repr.plot.height = 20, repr.plot.width = 20)

png(paste0(temp_plot_path, 'Network_filtered_TFs.png'), height = 30, width = 45, units = 'cm', res = 400)
print(p)
graphics.off()

# plot all TFs 
p <- GRNPlot(df.grn2,
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 42, 
             plot.importance = TRUE,
             min.importance = 2,
             remove.isolated = FALSE)

png(paste0(temp_plot_path, 'Network_all.png'), height = 30, width = 45, units = 'cm', res = 400)
print(p)
graphics.off()