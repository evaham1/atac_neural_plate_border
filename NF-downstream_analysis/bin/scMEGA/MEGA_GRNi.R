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
    #data_path = "./output/NF-downstream_analysis/Processing/FullData/TransferLabels/scMEGA_integrated/rds_files/"
    plot_path = "./output/NF-downstream_analysis/Processing/FullData/scMEGA/MEGA_GRNi/plots/"
    
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

## read in paired seurat object
obj.pair <- readRDS(paste0(data_path, "paired_object_chromvar.RDS"))
obj.pair
getAvailableMatrices(ArchRProj = obj.pair)

############################## Plotting colours #######################################

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
cols <- scHelper_cell_type_colours[as.character(unique(obj.pair$scHelper_cell_type))]

# stage cols
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_cols = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_cols) <- stage_order

###########################################################################################################
############################## Create cell groupings for trajectory #######################################
###########################################################################################################

###### add trajectory groupings
## HH5 and HH6, then split HH7, ss4 and ss8 by neural vs non-neural
stages <- as.data.frame(obj.pair$stage)
cell_types <- as.data.frame(obj.pair$scHelper_cell_type_broad)
metadata <- stages %>% mutate(cell_types = cell_types$`obj.pair$scHelper_cell_type_broad`)
colnames(metadata) <- c("stage", "cell_type")
head(metadata)

## first merge NC and neural // placodal and non-neural
unique(metadata$cell_type)
metadata <- metadata %>% mutate(broad = cell_type)
metadata <- metadata %>% mutate(
  broad = case_when(
    broad == "Neural" ~ cell_type,
    broad == "Non-neural" ~ cell_type,
    broad == "NC" ~ "Neural",
    broad == "Placodal" ~ "Non-neural",
  )
)
head(metadata)
unique(metadata$broad)

## add new metadata called 'Order'
# all HH5 -> HH5
# HH6 split: non-neural and placodal cells = HH6_NN, neural and NC cells = HH6_Neural
# HH7 split: non-neural and placodal cells = HH7_NN, neural and NC cells = HH7_Neural
# ss4 split: non-neural and placodal cells = ss4_NN, neural and NC cells = ss4_Neural
# ss8 split: non-neural and placodal cells = ss8_NN, neural and NC cells = ss8_Neural
metadata <- metadata %>% mutate(
  order = case_when(
    stage == "HH5" ~ stage,
    stage == "HH6" ~ paste0(stage, "_", broad),
    stage == "HH7" ~ paste0(stage, "_", broad),
    stage == "ss4" ~ paste0(stage, "_", broad),
    stage == "ss8" ~ paste0(stage, "_", broad),
  )
)

## add new cell groupings to seurat object
obj.pair$broad <- metadata$broad
unique(obj.pair$broad)

obj.pair$order <- metadata$order
unique(obj.pair$order)

############################## Plot cell groupings for trajectory #######################################

## colours for broad neural/non-neural division
broad <- c("Neural", "Non-neural")
broad_cols = c("#4D4799", "#EE5E5F")
names(broad_cols) <- broad

## colours for specific cell ordering
noneural_order <- c("HH5", "HH6_Non-neural", "HH7_Non-neural", "ss4_Non-neural", "ss8_Non-neural")
noneural_order_cols = c("#8DA0CB", "#D0AAAA", "#DD8889", "#E77071", "#EE5E5F")
names(noneural_order_cols) <- noneural_order

neural_order <- c("HH5", "HH6_Neural", "HH7_Neural", "ss4_Neural", "ss8_Neural")
neural_order_cols = c("#8DA0CB", "#A6A4BA", "#8885AF", "#6964A3", "#4D4799")
names(neural_order_cols) <- neural_order

order <- c("HH5",
           "HH6_Non-neural", "HH7_Non-neural", "ss4_Non-neural", "ss8_Non-neural",
           "HH6_Neural", "HH7_Neural", "ss4_Neural", "ss8_Neural")
order_cols = c("#8DA0CB", 
               "#D0AAAA", "#D3A3A3", "#E37B7B", "#EE5E5F",
               "#A6A4BA", "#9D9BB6", "#6E6AA5", "#4D4799")
names(order_cols) <- order

# plot broad neural/non-neural division
p2 <- DimPlot(obj.pair, pt.size = 2, group.by = "stage", shuffle = TRUE, label = FALSE, reduction = "umap_iLSI", cols = stage_cols)

p1 <- DimPlot(obj.pair, group.by = "broad", shuffle = TRUE, label = FALSE, reduction = "umap_iLSI", pt.size = 2, cols = broad_cols)
png(paste0(plot_path, 'Broad_grouping.png'), height = 10, width = 24, units = 'cm', res = 400)
p1 + p2
graphics.off()

# plot all ordering
p1 <- DimPlot(obj.pair, group.by = "order", shuffle = TRUE, label = FALSE, reduction = "umap_iLSI", pt.size = 2, cols = order_cols)
png(paste0(plot_path, 'Order_grouping.png'), height = 10, width = 14, units = 'cm', res = 400)
p1
graphics.off()

####################################################################################################
############################## NON-NEURAL: Create trajectory #######################################
####################################################################################################

# set up params to swap between non-neural and neural trajectories
temp_plot_path = paste0(plot_path, "non-neural-trajectory/")
dir.create(temp_plot_path, recursive = T)
temp_rds_path = paste0(rds_path, "non-neural-trajectory/")
dir.create(temp_rds_path, recursive = T)

# set ordering and cols
trajectory = c("HH5", "HH6_Non-neural", "HH7_Non-neural", "ss4_Non-neural", "ss8_Non-neural")
specific_cols = noneural_order_cols

# plot the non-neural ordering
p1 <- DimPlot(obj.pair, group.by = "order", shuffle = TRUE, label = FALSE, reduction = "umap_iLSI", pt.size = 2, cols = specific_cols)
png(paste0(temp_plot_path, 'Order_grouping.png'), height = 10, width = 14, units = 'cm', res = 400)
p1
graphics.off()

## create trajectory
obj.pair <- AddTrajectory(object = obj.pair, 
                          trajectory = trajectory,
                          group.by = "order", 
                          reduction = "umap_iLSI",
                          dims = 1:2, 
                          use.all = TRUE)

# visualise trajectory on full data
p2 <- TrajectoryPlot(object = obj.pair, 
                     reduction = "umap_iLSI",
                     continuousSet = "blueYellow",
                     size = 2,
                     addArrow = FALSE, randomize = T)
png(paste0(temp_plot_path, 'Trajectory.png'), height = 10, width = 12, units = 'cm', res = 400)
p2
graphics.off()

# visualise ordering and trajectory together
png(paste0(temp_plot_path, 'Trajectory_UMAPs.png'), height = 10, width = 22, units = 'cm', res = 400)
p1 + p2
graphics.off()

# subset to remove cells that are not part of trajectory
obj.trajectory <- obj.pair[, !is.na(obj.pair$Trajectory)]

# see how many cells left after filtering
cell_counts <- data.frame(dim(obj.pair)[2], dim(obj.trajectory)[2])
colnames(cell_counts) <- c("Before trajectory", "After trajectory total")

png(paste0(temp_plot_path, 'cell_counts_after_trajectory.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# visualise trajectory on subsetted data
cols <- scHelper_cell_type_colours[as.character(unique(obj.trajectory$scHelper_cell_type_broad))]
p1 <- DimPlot(obj.trajectory, reduction = "umap_iLSI", 
              group.by = "scHelper_cell_type_broad", pt.size = 2.5,
              cols = cols, shuffle = T)
p2 <- DimPlot(obj.trajectory, reduction = "umap_iLSI", 
              group.by = "stage", pt.size = 2.5,
              cols = stage_cols, shuffle = T)
p3 <- TrajectoryPlot(object = obj.trajectory, 
                     reduction = "umap_iLSI",
                     continuousSet = "blueYellow",
                     size = 2,
                     addArrow = FALSE, randomize = T)
png(paste0(temp_plot_path, 'Subsetted_trajectory_UMAPs.png'), height = 10, width = 32, units = 'cm', res = 400)
p1 + p2 + p3
graphics.off()

############################## Subset TFs and genes #######################################

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

# select 1op 10% genes that vary with the trajectory and correlate with peaks
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

####################################################################################################
############################## NEURAL: Create trajectory #######################################
####################################################################################################

# set up params to swap between non-neural and neural trajectories
temp_plot_path = paste0(plot_path, "neural-trajectory/")
dir.create(temp_plot_path, recursive = T)
temp_rds_path = paste0(rds_path, "neural-trajectory/")
dir.create(temp_rds_path, recursive = T)

# set ordering and cols
trajectory = c("HH5", "HH6_Neural", "HH7_Neural", "ss4_Neural", "ss8_Neural")
specific_cols = neural_order_cols

# plot the ordering
p1 <- DimPlot(obj.pair, group.by = "order", shuffle = TRUE, label = FALSE, reduction = "umap_iLSI", pt.size = 2, cols = specific_cols)
png(paste0(temp_plot_path, 'Order_grouping.png'), height = 10, width = 14, units = 'cm', res = 400)
p1
graphics.off()

## create trajectory
obj.pair <- AddTrajectory(object = obj.pair, 
                          trajectory = trajectory,
                          group.by = "order", 
                          reduction = "umap_iLSI",
                          dims = 1:2, 
                          use.all = TRUE)

# visualise trajectory on full data
p2 <- TrajectoryPlot(object = obj.pair, 
                     reduction = "umap_iLSI",
                     continuousSet = "blueYellow",
                     size = 2,
                     addArrow = FALSE, randomize = T)
png(paste0(temp_plot_path, 'Trajectory.png'), height = 10, width = 12, units = 'cm', res = 400)
p2
graphics.off()

# visualise ordering and trajectory together
png(paste0(temp_plot_path, 'Trajectory_UMAPs.png'), height = 10, width = 22, units = 'cm', res = 400)
p1 + p2
graphics.off()

# subset to remove cells that are not part of trajectory
obj.trajectory <- obj.pair[, !is.na(obj.pair$Trajectory)]

# see how many cells left after filtering
cell_counts <- data.frame(dim(obj.pair)[2], dim(obj.trajectory)[2])
colnames(cell_counts) <- c("Before trajectory", "After trajectory total")

png(paste0(temp_plot_path, 'cell_counts_after_trajectory.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("Remaining Cell Count", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# visualise trajectory on subsetted data
cols <- scHelper_cell_type_colours[as.character(unique(obj.trajectory$scHelper_cell_type_broad))]
p1 <- DimPlot(obj.trajectory, reduction = "umap_iLSI", 
              group.by = "scHelper_cell_type_broad", pt.size = 2.5,
              cols = cols, shuffle = T)
p2 <- DimPlot(obj.trajectory, reduction = "umap_iLSI", 
              group.by = "stage", pt.size = 2.5,
              cols = stage_cols, shuffle = T)
p3 <- TrajectoryPlot(object = obj.trajectory, 
                     reduction = "umap_iLSI",
                     continuousSet = "blueYellow",
                     size = 2,
                     addArrow = FALSE, randomize = T)
png(paste0(temp_plot_path, 'Subsetted_trajectory_UMAPs.png'), height = 10, width = 32, units = 'cm', res = 400)
p1 + p2 + p3
graphics.off()

############################## Subset TFs and genes #######################################

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

# select 1op 10% genes that vary with the trajectory and correlate with peaks
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