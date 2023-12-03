#!/usr/bin/env Rscript

print("ArchR constrained integration, assigning cluster labels and peak to gene linkage")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)
library(presto)
library(Seurat)
library(plyr)
library(gtools)
library(scHelper)
library(ggrepel)

############################## Set up script options #######################################

option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-m", "--min_threshold"), action = "store", type = "integer", help = "minimum percentage of cells a label must contain", default = 40),
  make_option(c("-l", "--max_label"), action = "store", type = "integer", help = "maximum number of labels a cluster can have", default = 3),
  make_option(c("", "--verbose"), action = "store", type = "logical", help = "Verbose", default = TRUE)
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
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    csv_path = "./csv_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  dir.create(csv_path, recursive = T)
}

set.seed(42)

###### FUNCTIONS TO UPDATE IN SCHELPER:
ArchRCellCountsHeatmap <- function (ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters") 
{
  group1_data <- getCellColData(ArchR, select = group1)[, 1]
  group2_data <- getCellColData(ArchR, select = group2)[, 1]
  cM <- confusionMatrix(paste0(group2_data), paste0(group1_data))
  cM <- cM/Matrix::rowSums(cM)
  p <- pheatmap::pheatmap(mat = cM, color = paletteContinuous("whiteBlue"), 
                          border_color = "black", fontsize = 25)
}

############################## Read in ArchR project and seurat object #######################################

# Load atac data in rds_files
ArchR <- loadArchRProject(path = paste0(data_path, "rds_files/", "FullData_Save-ArchR"), force = FALSE, showLogo = TRUE)

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getPeakSet(ArchR)
getAvailableMatrices(ArchR)

# set output directory to make sure new files are saved in the right place
output_directory <- paste0(rds_path, "FullData_Save-ArchR")
getOutputDirectory(ArchR)

# load seurat object
seurat_data <- readRDS(paste0(data_path, "seurat_label_transfer.RDS"))

############################################################################################
############################## Variance plots #######################################

print("Making plots to explore stage variance...")

stage_cols = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
temp_plot_path = paste0(plot_path, "variance_plots/")
dir.create(temp_plot_path, recursive = T)

# Plot PCs for RNA data
for (i in 1:50){
  pc_oi <- paste0('PC_', i)
  print(pc_oi)
  
  # Subset embeddings for PC of interest
  pc_embeddings <- seurat_data@reductions$pca@cell.embeddings[,pc_oi]
  
  # Prepare plot data
  plot_data <- data.frame(cell_type = as.character(seurat_data@meta.data[names(pc_embeddings), 'stage']), pc_embeddings = pc_embeddings)
  
  # Plot histogram of cell state distributions along PC of interest
  png(paste0(temp_plot_path, "RNA_", pc_oi, "_stage.png"), width = 11, height = 6, units = 'cm', res = 720)
  print(
    ggplot(plot_data, aes(x = pc_embeddings, colour = cell_type, fill = cell_type)) +
      geom_density(alpha = 0.2) +
      scale_fill_manual(values = stage_cols) +
      scale_colour_manual(values = stage_cols) +
      xlab(sub('_', ' ', pc_oi)) +
      ylab('Density') +
      theme_classic() +
      theme(legend.title=element_blank()) +
      scale_x_reverse()
  )
  graphics.off()
}

# Extract total variance explained for each PC and plot
mat <- Seurat::GetAssayData(seurat_data, assay = "RNA", slot = "scale.data")
pca <- seurat_data[["pca"]]
total_variance = seurat_data@reductions$pca@misc$total.variance
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = (eigValues / total_variance)*100
totalvarExplained = sum(varExplained)
data <- data.frame(PC = c(colnames(seurat_data@reductions$pca@cell.embeddings), "Total"), 
                   VarianceExplained = c(varExplained, totalvarExplained))
data <- data %>% mutate(RelativeVarianceExplained = (VarianceExplained / totalvarExplained)*100)

png(paste0(temp_plot_path, 'RNA_PCs_variance_explained.png'), height = 45, width = 15, units = 'cm', res = 400)
grid.arrange(top=textGrob("PCs", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(data, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# Plot SVDs for ATAC data
for (i in 1:30){
  svd_oi <- paste0('LSI', i)
  print(svd_oi)
  
  # Subset embeddings for PC of interest
  embeddings <- ArchR@reducedDims$IterativeLSI@listData$matSVD[,svd_oi]
  
  # Prepare plot data
  plot_data <- as.data.frame(merge(getCellColData(ArchR, select = "stage"), as.data.frame(embeddings), by = "row.names"))
  
  # Plot histogram of cell state distributions along PC of interest
  png(paste0(temp_plot_path, "ATAC_", svd_oi, "_stage.png"), width = 11, height = 6, units = 'cm', res = 720)
  print(
    ggplot(plot_data, aes(x = embeddings, colour = stage, fill = stage)) +
      geom_density(alpha = 0.2) +
      scale_fill_manual(values = stage_cols) +
      scale_colour_manual(values = stage_cols) +
      xlab(sub('_', ' ', svd_oi)) +
      ylab('Density') +
      theme_classic() +
      theme(legend.title=element_blank()) +
      scale_x_reverse()
  )
  graphics.off()
  
}

# Extract the d for each SVD component, but maybe shouldnt actually use this
data <- data.frame(SVD = colnames(ArchR@reducedDims$IterativeLSI@listData$matSVD),
                   d = ArchR@reducedDims$IterativeLSI@listData$svd[["d"]])
png(paste0(temp_plot_path, 'ATAC_SVDss_d_value.png'), height = 45, width = 15, units = 'cm', res = 400)
grid.arrange(top=textGrob("SVDss", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(data, rows=NULL, theme = ttheme_minimal()))
graphics.off()

############################################################################################
############################## Pre-Integration Plots #######################################

print("Making pre-integration plots...")

#### Prepare RNA labels and colours #######

#### these labels were recalculated on rna data with contam, save them but wont use them
seurat_data@meta.data$scHelper_cell_type_new <- seurat_data@meta.data$scHelper_cell_type

#### combine contamination and old scHelper_cell_state labels --> call this scHelper_cell_type
contam <- seurat_data@meta.data$contamination
original <- seurat_data@meta.data$scHelper_cell_type_original
old_labels <- contam
old_labels[is.na(contam)] <- as.character(original[is.na(contam)])
seurat_data@meta.data$scHelper_cell_type <- old_labels

#### add broad labels to RNA data
scHelper_cell_types <- data.frame(seurat_data$scHelper_cell_type)
broad <- scHelper_cell_types %>% mutate(broad = mapvalues(seurat_data.scHelper_cell_type, 
                                                          from=c("NP", "aNP", "iNP", "pNP", "eN", "vFB", "FB", "MB", "HB", "eCN", "eN",
                                                                 'PPR', 'aPPR', 'pPPR',
                                                                 'eNPB', 'NPB', 'aNPB', 'pNPB',
                                                                 'NC', 'dNC',
                                                                 'NNE', 'pEpi',
                                                                 'EE', 'meso', 'endo', 'BI', 'PGC'),
                                                          to=c(
                                                            rep("Neural", 11),
                                                            rep("Placodal", 3),
                                                            rep("NPB", 4),
                                                            rep("NC", 2),
                                                            rep("Non-neural", 2),
                                                            rep("Contam", 5)
                                                          )
))
seurat_data$scHelper_cell_type_broad <- broad$broad

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order
stage_cols <- stage_colours[levels(droplevels(seurat_data@meta.data$stage))]

###### schelper cell type colours
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak', 
                              'PGC', 'BI', 'meso', 'endo',
                              'Neural', 'Placodal', 'Non-neural', 'Contam')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                "#A5718D", "#3F918C", "#ed5e5f", "#9792A3")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                       'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo',
                                       'Neural', 'Placodal', 'Non-neural', 'Contam')
seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = scHelper_cell_type_order)
scHelper_new_cols <- scHelper_cell_type_colours[levels(droplevels(seurat_data@meta.data$scHelper_cell_type_new))]
scHelper_old_cols <- scHelper_cell_type_colours[levels(droplevels(seurat_data@meta.data$scHelper_cell_type))]

###### UMAPs
plot_path = "./plots/before_integration/"
dir.create(plot_path, recursive = T)

umap_rna_new <- DimPlot(seurat_data, group.by = 'scHelper_cell_type_new', label = TRUE, 
                    label.size = 12,
                    label.box = TRUE, repel = TRUE,
                    pt.size = 1.2, 
                    cols = scHelper_new_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
umap_rna_old <- DimPlot(seurat_data, group.by = 'scHelper_cell_type', label = TRUE, 
                    label.size = 12,
                    label.box = TRUE, repel = TRUE,
                    pt.size = 1.2, 
                    cols = scHelper_old_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())

png(paste0(plot_path, 'RNA_UMAPs_old_vs_new.png'), height = 20, width = 40, units = 'cm', res = 400)
print(umap_rna_old + umap_rna_new)
graphics.off()

scHelper_broad_cols <- scHelper_cell_type_colours[seurat_data@meta.data$scHelper_cell_type_broad]
umap_broad <- DimPlot(seurat_data, group.by = 'scHelper_cell_type_broad', label = TRUE, 
                        label.size = 12,
                        label.box = TRUE, repel = TRUE,
                        pt.size = 1.2, 
                        cols = scHelper_broad_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())

png(paste0(plot_path, 'RNA_UMAP_broad.png'), height = 20, width = 20, units = 'cm', res = 400)
print(umap_broad)
graphics.off()

############################## UMAPs before integration #######################################
# UMAPs of RNA and ATAC data, with RNA coloured by cell state and ATAC by clusters
umap_atac <- plotEmbedding(ArchR, name = "clusters", plotAs = "points", size = 1.8, baseSize = 0, labelSize = 8, legendSize = 0)

png(paste0(plot_path, 'UMAPs_before_integration_scHelper_cell_states.png'), height = 20, width = 40, units = 'cm', res = 400)
print(umap_rna_old + umap_atac)
graphics.off()

umap_atac <- plotEmbedding(ArchR, name = "stage", plotAs = "points", size = 1.8, baseSize = 0, 
labelSize = 8, legendSize = 0, cols = stage_colours)
umap_rna <- DimPlot(seurat_data, group.by = 'stage', label = FALSE, 
                    pt.size = 1.8, 
                    cols = stage_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())

png(paste0(plot_path, 'UMAPs_before_integration_stage.png'), height = 20, width = 40, units = 'cm', res = 400)
print(umap_rna + umap_atac)
graphics.off()

print("Pre-integration plots made.")

################################################################################################
############################## Constrained integration #######################################

Idents(object = seurat_data) <- "stage"

groupList <- SimpleList(
  GroupHH5 = SimpleList(
    ATAC = ArchR$cellNames[ArchR$stage %in% "HH5"],
    RNA = SeuratObject::WhichCells(seurat_data, idents = "HH5")
  ),
  GroupHH6 = SimpleList(
    ATAC = ArchR$cellNames[ArchR$stage %in% "HH6"],
    RNA = SeuratObject::WhichCells(seurat_data, idents = "HH6")
  ),
  GroupHH7 = SimpleList(
    ATAC = ArchR$cellNames[ArchR$stage %in% "HH7"],
    RNA = SeuratObject::WhichCells(seurat_data, idents = "HH7")
  ),
  Groupss4 = SimpleList(
    ATAC = ArchR$cellNames[ArchR$stage %in% "ss4"],
    RNA = SeuratObject::WhichCells(seurat_data, idents = "ss4")
  ),
  Groupss8 = SimpleList(
    ATAC = ArchR$cellNames[ArchR$stage %in% "ss8"],
    RNA = SeuratObject::WhichCells(seurat_data, idents = "ss8")
  )
)
print(groupList)

print("Starting constrained integration...")

ArchR <- addGeneIntegrationMatrix(
  ArchRProj = ArchR, seRNA = seurat_data,
  useMatrix = "GeneScoreMatrix", matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  groupList = groupList,
  groupRNA = "scHelper_cell_type", nameCell = "predictedCell",
  nameGroup = "predictedGroup", nameScore = "predictedScore",
  force = TRUE, addToArrow = TRUE, threads = 1,
  sampleCellsATAC = 85000, sampleCellsRNA = 15000
)
print("integration completed")

# use matched RNA cells to add new and old labels to ATAC cells
extracted_rna_labels <- seurat_data@meta.data[ArchR$predictedCell, c("scHelper_cell_type_new", "scHelper_cell_type", "scHelper_cell_type_broad")]
ArchR$scHelper_cell_type_new <- extracted_rna_labels[, "scHelper_cell_type_new"]
ArchR$scHelper_cell_type <- extracted_rna_labels[, "scHelper_cell_type"]
ArchR$scHelper_cell_type_broad <- extracted_rna_labels[, "scHelper_cell_type_broad"]
print("scHelper cell type labels added")

# use matched RNA cells to add rna metadata to ATAC cells
extracted_rna_metadata <- seurat_data@meta.data[ArchR$predictedCell, c("run", "stage", "seurat_clusters")]
ArchR$rna_stage <- extracted_rna_metadata[, "stage"]
ArchR$rna_run <- extracted_rna_metadata[, "run"]
ArchR$rna_clusters <- extracted_rna_metadata[, "seurat_clusters"]
print("RNA metadata added")

print(head(getCellColData(ArchR)))

#############################################################################################
################################ QC of Integration ##########################################

print("Making integration QC plots...")

plot_path = "./plots/after_integration/integration_QC/"
dir.create(plot_path, recursive = T)

# highlight which RNA cells were used for integration
umap_rna_used <- DimPlot(seurat_data,
                    pt.size = 1.2, 
                    shuffle = TRUE,
                    cells.highlight = unique(ArchR$predictedCell), cols.highlight = "red", cols = "gray") +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())

png(paste0(plot_path, 'RNA_cells_used_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
print(umap_rna_used)
graphics.off()

# count how many RNA cells used for integration
df <- data.frame(
  RNA_counts = length(Cells(seurat_data)),
  ATAC_counts = length(ArchR$cellNames),
  Transferred_RNA_counts = length(unique(ArchR$predictedCell))
)

png(paste0(plot_path, 'RNA_cell_counts_used.png'), height = 10, width = 20, units = 'cm', res = 400)
grid.arrange(top=textGrob("Cell counts", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

#### distribution of broad labels in RNA and ATAC data

# RNA: extract cell state proportions and order and colour them -> plot pie chart
counts <- as.data.frame(table(seurat_data$scHelper_cell_type))
counts$Freq <- as.numeric(counts$Freq)
colnames(counts) <- c("Cell state", "Frequency")
counts <- counts %>% mutate('Percentage' = round(100*(Frequency/sum(Frequency))))

counts2 <- counts %>% 
  mutate(csum = rev(cumsum(rev(Frequency))),
         pos = Frequency/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Frequency/2, pos))

png(paste0(plot_path, 'RNA_cell_state_distribution_piechart.png'), height = 20, width = 30, units = 'cm', res = 400)
ggplot(counts2, aes(x = "" , y = Frequency, fill = fct_inorder(`Cell state`))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y", direction=-1) +
  scale_fill_manual(values = scHelper_cell_type_colours) +
  ggrepel::geom_label_repel(data = counts2,
                   aes(y = pos, label = `Cell state`),
                   size = 6, nudge_x = 1, show.legend = FALSE) +
  theme_void() +
  theme(text = element_text(size = 30), legend.position = "none")
graphics.off()

png(paste0(plot_path, 'RNA_cell_state_distribution_table.png'), height = 50, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("RNA data", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# RNA: extract cell state proportions and order and colour them -> plot pie chart
counts <- as.data.frame(table(seurat_data$scHelper_cell_type_broad))
counts$Freq <- as.numeric(counts$Freq)
colnames(counts) <- c("Cell state", "Frequency")
counts <- counts %>% mutate('Percentage' = round(100*(Frequency/sum(Frequency))))

counts2 <- counts %>% 
  mutate(csum = rev(cumsum(rev(Frequency))),
         pos = Frequency/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Frequency/2, pos))

png(paste0(plot_path, 'RNA_cell_state_distribution_broad_piechart.png'), height = 20, width = 30, units = 'cm', res = 400)
ggplot(counts2, aes(x = "" , y = Frequency, fill = fct_inorder(`Cell state`))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y", direction=-1) +
  scale_fill_manual(values = scHelper_cell_type_colours) +
  ggrepel::geom_label_repel(data = counts2,
                   aes(y = pos, label = `Cell state`),
                   size = 6, nudge_x = 1, show.legend = FALSE) +
  theme_void() +
  theme(text = element_text(size = 30), legend.position = "none")
graphics.off()

png(paste0(plot_path, 'RNA_cell_state_distribution_broad_table.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("RNA data", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()


# ATAC: extract cell state proportions and order and colour them -> plot pie chart
counts <- as.data.frame(table(ArchR$scHelper_cell_type))
counts$Freq <- as.numeric(counts$Freq)
colnames(counts) <- c("Cell state", "Frequency")
counts <- counts %>% mutate('Percentage' = round(100*(Frequency/sum(Frequency))))

counts2 <- counts %>% 
  mutate(csum = rev(cumsum(rev(Frequency))),
         pos = Frequency/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Frequency/2, pos))

png(paste0(plot_path, 'ATAC_cell_state_distribution_piechart.png'), height = 20, width = 30, units = 'cm', res = 400)
ggplot(counts2, aes(x = "" , y = Frequency, fill = fct_inorder(`Cell state`))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y", direction=-1) +
  scale_fill_manual(values = scHelper_cell_type_colours) +
  ggrepel::geom_label_repel(data = counts2,
                   aes(y = pos, label = `Cell state`),
                   size = 6, nudge_x = 1, show.legend = FALSE) +
  theme_void() +
  theme(text = element_text(size = 30), legend.position = "none")
graphics.off()

png(paste0(plot_path, 'ATAC_cell_state_distribution_table.png'), height = 20, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("ATAC data", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# ATAC: extract cell state proportions and order and colour them -> plot pie chart
counts <- as.data.frame(table(ArchR$scHelper_cell_type_broad))
counts$Freq <- as.numeric(counts$Freq)
colnames(counts) <- c("Cell state", "Frequency")
counts <- counts %>% mutate('Percentage' = round(100*(Frequency/sum(Frequency))))

counts2 <- counts %>% 
  mutate(csum = rev(cumsum(rev(Frequency))),
         pos = Frequency/2 + lead(csum, 1),
         pos = if_else(is.na(pos), Frequency/2, pos))

png(paste0(plot_path, 'ATAC_cell_state_distribution_broad_piechart.png'), height = 20, width = 20, units = 'cm', res = 400)
ggplot(counts2, aes(x = "" , y = Frequency, fill = fct_inorder(`Cell state`))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y", direction=-1) +
  scale_fill_manual(values = scHelper_cell_type_colours) +
  ggrepel::geom_label_repel(data = counts2,
                   aes(y = pos, label = `Cell state`),
                   size = 6, nudge_x = 1, show.legend = FALSE) +
  theme_void() +
  theme(text = element_text(size = 30), legend.position = "none")
graphics.off()

png(paste0(plot_path, 'ATAC_cell_state_distribution_broad_table.png'), height = 10, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("ATAC data", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

#############################################################################################
############################## Post-Integration Plots #######################################

print("Making post-integration plots...")

# set colour palettes for UMAPs
atac_scHelper_new_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_new)]
atac_scHelper_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type)]
atac_scHelper_broad <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_broad)]

############################## RNA cell labels on ATAC data #######################################

### New labels
plot_path = "./plots/after_integration/new_labels/"
dir.create(plot_path, recursive = T)

png(paste0(plot_path, 'UMAP_integrated.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type_new", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_new_cols, labelAsFactors = FALSE)
graphics.off()

png(paste0(plot_path, 'UMAP_integrated_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type_new", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_new_cols)
graphics.off()

### Old labels / ones I will use
plot_path = "./plots/after_integration/old_labels_to_use/"
dir.create(plot_path, recursive = T)

png(paste0(plot_path, 'UMAP_integrated.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_cols, labelAsFactors = FALSE)
graphics.off()

png(paste0(plot_path, 'UMAP_integrated_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_cols)
graphics.off()

### Broad labels
plot_path = "./plots/after_integration/broad_labels/"
dir.create(plot_path, recursive = T)

png(paste0(plot_path, 'UMAP_integrated.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type_broad", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_broad, labelAsFactors = FALSE)
graphics.off()

png(paste0(plot_path, 'UMAP_integrated_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type_broad", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_broad)
graphics.off()

############################## Integration scores plots #######################################

plot_path = "./plots/after_integration/integration_scores/"
dir.create(plot_path, recursive = T)

png(paste0(plot_path, 'Integration_Scores_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "predictedScore", plotAs = "points", size = 1.8, baseSize = 0, 
              legendSize = 10)
graphics.off()

png(paste0(plot_path, "Integration_Scores_Vln.png"), width=50, height=20, units = 'cm', res = 200)
plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", 
  name = "predictedScore", plotAs = "Violin", alpha = 0.4)
graphics.off()

print("Post-integration plots made.")

############################## Gene Integration and Score Plots for key TFs #######################################

plot_path = "./plots/plots_of_key_TFs/"
dir.create(plot_path, recursive = T)

ArchR <- addImputeWeights(ArchR, seed = 1)

# set genes of interest
TFs <- c("SIX1", "EYA2", "IRF6", "DLX5", "DLX6", "GATA2", "GATA3", 
         "TFAP2A", "TFAP2B", "TFAP2C", "TFAP2E",
         "PITX1", "PITX2",
         "PAX7", "MSX1", "ETS1", 
         "SOX2", "SOX9", "SOX8", "SOX10", "SOX5", "SOX21", "SOX3",
         "NKX6-2", "CSRNP1", "SNAI2", "LMX1B", "ZEB1", "ZEB2",
         "EPAS1", "BMP4", "YEATS4", "HOXB1", "EOMES", "ADMP",
         "TCF3", "LEF1", "TEAD3", "TEAD4", "FOXK1", "FOXK2", # no TCF4
         "NKX2-3", "NKX2-5", "NKX6-2", "RFX2", "ZIC1", "ZIC3", "TBXT",
         "NR6A1", "HIC2", "HES1", "HES5", "FOS", "JUND")

# ArchR <- addImputeWeights(ArchR) # shouldnt really do this for gene counts

# Plot UMAP for expression of each TF
for (TF in TFs){
  print(TF)
  
  # Plot GeneIntegration values on UMAP
  png(paste0(plot_path, TF, '_gene_integration_UMAP.png'), height = 12, width = 10, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = TF,
                plotAs = "points", size = 1.8,
                colorBy = "GeneIntegrationMatrix", continuousSet = "blueYellow",
                imputeWeights = NULL) + 
    theme_ArchR(legendTextSize = 17, baseSize = 17, plotMarginCm = 0.5))
  graphics.off()
  
  # Plot GeneScore values on UMAP
  png(paste0(plot_path, TF, '_gene_score_UMAP.png'), height = 12, width = 14, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = TF,
                plotAs = "points", size = 1.8,
                colorBy = "GeneScoreMatrix", continuousSet = "horizon",
                imputeWeights = getImputeWeights(ArchR)) + 
    theme_ArchR(legendTextSize = 12, baseSize = 16, plotMarginCm = 0.5))
  graphics.off()
  
}

##################### Distribution of labels across clusters ##################################

plot_path = "./plots/label_by_cluster_distribution/"
dir.create(plot_path, recursive = T)

# visualise distribution across clusters: table of cell counts
png(paste0(plot_path, 'label_by_cluster_cell_number_table.png'), height = 25, width = 40, units = 'cm', res = 400)
scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters", print_table = TRUE, scHelper_cell_type_order = scHelper_cell_type_order)
graphics.off()

# visualise distribution across clusters: confusion matrix
png(paste0(plot_path, "label_by_cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
ArchRCellCountsHeatmap(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters")
graphics.off()

# visualise distribution across clusters: table of cell percentages
cell_counts <- scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
percentage_counts <- as.data.frame(lapply(cell_counts, function(x) (x / sum(x))*100))
rownames(percentage_counts) <- rownames(cell_counts)

png(paste0(plot_path, 'label_by_cluster_cell_percentage_table.png'), height = 25, width = 40, units = 'cm', res = 400)
grid.arrange(tableGrob(round(percentage_counts, 2), theme = ttheme_minimal()))
graphics.off()

# visualise distribution across clusters: piecharts
counts <- scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
png(paste0(plot_path, "label_by_cluster_piecharts.png"), width=50, height=40, units = 'cm', res = 200)
scHelper::CellLabelPieCharts(counts, col = scHelper_cell_type_colours)
graphics.off()


###############################################################################################
############################## CO-ACCESSIBILITY BETWEEN PEAKS #################################

print("Calculating coaccessibility...")

# calculate co-accessibility between all peaks
ArchR <- addCoAccessibility(ArchR)

# extract interactions - returns indexes of queryHits and subjectHits
cA <- getCoAccessibility(ArchR, corCutOff = 0, returnLoops = FALSE)
cA
  # DataFrame with 120270 rows and 11 columns
  # queryHits subjectHits seqnames correlation Variability1 Variability2     TStat        Pval         FDR VarQuantile1 VarQuantile2
  # <integer>   <integer>    <Rle>   <numeric>    <numeric>    <numeric> <numeric>   <numeric>   <numeric>    <numeric>    <numeric>
  #   1              3           4     chr1    0.548725   0.00437754   0.00683964   14.5441 4.15759e-40 4.52151e-38     0.911185     0.965430
  # 2              4           3     chr1    0.548725   0.00683964   0.00437754   14.5441 4.15759e-40 4.52151e-38     0.965430     0.911185
  # 3              4           5     chr1    0.517190   0.00683964   0.00356568   13.3901 4.49249e-35 3.64027e-33     0.965430     0.870967
  # 4              5           4     chr1    0.517190   0.00356568   0.00683964   13.3901 4.49249e-35 3.64027e-33     0.870967     0.965430
  # 5             27          40     chr1    0.761607   0.01690577   0.00855042   26.0418 1.47916e-94 2.12498e-91     0.995825     0.978303
coacessibility_df <- as.data.frame(cA)

# Need to use indices from df to extract granges and therefore informative peak IDs
coacessibility_df <- coacessibility_df %>% 
  mutate(query_PeakID = paste0(seqnames(metadata(cA)[[1]][queryHits]), "-", start(metadata(cA)[[1]][queryHits]), "-", end(metadata(cA)[[1]][queryHits]))) %>%
  mutate(subject_PeakID = paste0(seqnames(metadata(cA)[[1]][subjectHits]), "-", start(metadata(cA)[[1]][subjectHits]), "-", end(metadata(cA)[[1]][subjectHits])))

head(coacessibility_df)

# sanity check that all interaction Peak IDs are in the ArchR peakset
table(coacessibility_df$subject_PeakID %in% getPeakSet(ArchR)$name)

# save df
write.csv(coacessibility_df, file = paste0(csv_path, "Peak_coaccessibility_df.csv"), row.names = FALSE)

# extract resulting interactions as granges object
granges <- getCoAccessibility(ArchR, corCutOff = 0, returnLoops = TRUE)[[1]]
saveRDS(granges, paste0(csv_path, "Peak_coaccessibility.RDS"))

print("Coaccessibility calculated and saved.")

####################################################################################################################
############################## CO-ACCESSIBILITY BETWEEN PEAKS AND GENES - FIXED DIST #################################

print("Calculating coaccessibility between peaks and genes within a 250000 distance...")

# calculate gene-to-peak co-accessibility using GeneIntegrationMatrix
ArchR <- addPeak2GeneLinks(ArchR, maxDist = 250000)
# biggest chrom chrom 1 size: 197608386 (200000000), default is 250000
# tried running at max distance but then found very few interactions very far away...

# extract all positively correlated interactions as df
p2g <- getPeak2GeneLinks(ArchR, corCutOff = 0, FDRCutOff = 1, returnLoops = FALSE)
p2g_df <- as.data.frame(p2g)

# add correct Peak IDs and gene names to df
p2g_df <- p2g_df %>% 
  mutate(PeakID = paste0(seqnames(metadata(p2g)$peakSet[idxATAC]), "-", start(metadata(p2g)$peakSet[idxATAC]), "-", end(metadata(p2g)$peakSet[idxATAC]))) %>%
  mutate(gene_name = metadata(p2g)$geneSet[idxRNA]$name)
head(p2g_df)
print(paste0(nrow(p2g_df), " interactions identified by coaccessibility!"))

# sanity check that all interaction Peak IDs are in the ArchR peakset
if(sum(p2g_df$PeakID %in% getPeakSet(ArchR)$name) != nrow(p2g_df)){stop("Issue with peak IDs in interactions!")}

# check size of df
print("Extracting P2G df:")
print(dim(p2g_df))
head(p2g_df)

# save df
write.csv(p2g_df, file = paste0(csv_path, "Peak_to_gene_linkage_df_250000_distance.csv"), row.names = FALSE)

# extract all positively correlated interactions as granges object
print("Extracting P2G granges:")
granges <- getPeak2GeneLinks(ArchR, returnLoops = TRUE,
                             corCutOff = 0, FDRCutOff = 1)[[1]]
print(granges)

# check that interactions df and interactions granges are the same size
if(nrow(p2g_df) != length(granges)){stop("Df and granges extraction resulted in different interaction numebers!")}

# save the granges
saveRDS(granges, paste0(csv_path, "Peak_to_gene_linkage_250000_distance.RDS"))

print("Coaccessibility between peaks and genes (250000 dist) calculated and saved.")

####################################################################################################################
############################## EXTRACT GENE LOCATIONS #################################

# extract gene locations
gene_metadata <- metadata(getPeak2GeneLinks(ArchR, corCutOff = 0, FDRCutOff = 1, returnLoops = FALSE))$geneSet
print("Extracting gene metadata:")
print(gene_metadata)

# check that every gene from interaction has a corresponding location
if(sum(unique(p2g_df$gene_name) %in% gene_metadata$name) != length(unique(p2g_df$gene_name))){stop("Some genes from interactions don't have a genomic location!")}

# save the gene locations
saveRDS(gene_metadata, paste0(csv_path, "Gene_locations.RDS"))

##################################################################################
############################## SAVE ARCHR OBJECT #################################

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getAvailableMatrices(ArchR)

# save integrated ArchR project
print("Saving integrated ArchR project...")
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
print(paste0("Output filename = ", rds_path, "FullData_Save-ArchR"))
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, "FullData_Save-ArchR"), load = FALSE)
print("Integrated ArchR project saved.")

################################################################################
############################## HEATMAPS OF P2L #################################

plot_path = "./plots/peak2gene_250000_dist/"
dir.create(plot_path, recursive = T)

# This sporadically fails so comment out for now, not v useful plot anyway
## Heatmap of linkage across clusters
p <- plotPeak2GeneHeatmap(ArchRProj = ArchR, groupBy = "clusters")
png(paste0(plot_path, 'Peak_to_gene_linkage_clusters_heatmap.png'), height = 80, width = 60, units = 'cm', res = 400)
print(p)
graphics.off()

p <- plotPeak2GeneHeatmap(ArchRProj = ArchR, groupBy = "stage")
png(paste0(plot_path, 'Peak_to_gene_linkage_stage_heatmap.png'), height = 80, width = 60, units = 'cm', res = 400)
print(p)
graphics.off()