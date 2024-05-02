#!/usr/bin/env Rscript

print("ArchR integration ")

############################## Load libraries #######################################
library(getopt)
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

# check that ArchR version is 1.0.3
sessionInfo()

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
    
    ncores = 8
    
    plot_path = "./output/NF-downstream_analysis/ArchR_integration/plots/"
    rds_path = "./output/NF-downstream_analysis/ArchR_integration/rds_files/"
    
    # temp folder with ss8 ATAC and ss8 RNA input
    data_path = "./output/NF-downstream_analysis/temp_archr_integration_inputs/"

    # already integrated
    # data_path = "./output/NF-downstream_analysis/ArchR_integration//ss8/1_unconstrained_integration/rds_files/"
    # data_path = "./output/NF-downstream_analysis/ArchR_integration/HH5/1_unconstrained_integration/rds_files/"
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}



############################## Read in ArchR project and seurat object #######################################

# Pull out label, input folder should look like: rds_files/HH5-SaveArchR, HH5_clustered_data.RDS

files <- list.files(data_path)
print("Files: ")
print(files)

files <- files[!files %in% "rds_files"]
label <- unique(sub('_.*', '', files))
print("label: ")
print(label)

# Load atac data in rds_files
ArchR <- loadArchRProject(path = paste0(data_path, "rds_files/", label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)

# set output directory to make sure new files are saved in the right place
output_directory <- paste0(rds_path, label[1], "_Save-ArchR")
getOutputDirectory(ArchR)

# load seurat object by reading in any rds object
rna_path <- list.files(path = data_path, pattern = "*.RDS", full.names = TRUE)
seurat_data <- readRDS(rna_path)


############################################################################################
############################## Pre-Integration Plots #######################################

print("Making pre-integration plots...")

#### Prepare RNA labels and colours #######

#### combine contamination and old scHelper_cell_state labels
contam <- seurat_data@meta.data$contamination
original <- seurat_data@meta.data$scHelper_cell_type_original
old_labels <- contam
old_labels[is.na(contam)] <- as.character(original[is.na(contam)])
seurat_data@meta.data$old_labels <- old_labels

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
                              'PGC', 'BI', 'meso', 'endo')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                                "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                                "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                       'vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo')
seurat_data@meta.data$old_labels <- factor(seurat_data@meta.data$old_labels, levels = scHelper_cell_type_order)
scHelper_new_cols <- scHelper_cell_type_colours[levels(droplevels(seurat_data@meta.data$scHelper_cell_type))]
scHelper_old_cols <- scHelper_cell_type_colours[levels(droplevels(seurat_data@meta.data$old_labels))]

###### UMAPs
plot_path = "./plots/before_integration/"
dir.create(plot_path, recursive = T)

umap_rna_new <- DimPlot(seurat_data, group.by = 'scHelper_cell_type', label = TRUE, 
                    label.size = ifelse(length(unique(seurat_data$stage)) == 1, 9, 3),
                    label.box = TRUE, repel = TRUE,
                    pt.size = ifelse(length(unique(seurat_data$stage)) == 1, 1.2, 1), 
                    cols = scHelper_new_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
umap_rna_old <- DimPlot(seurat_data, group.by = 'old_labels', label = TRUE, 
                    label.size = ifelse(length(unique(seurat_data$stage)) == 1, 9, 3),
                    label.box = TRUE, repel = TRUE,
                    pt.size = ifelse(length(unique(seurat_data$stage)) == 1, 1.2, 1), 
                    cols = scHelper_old_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())

png(paste0(plot_path, 'RNA_UMAPs_old_vs_new.png'), height = 20, width = 40, units = 'cm', res = 400)
print(umap_rna_new + umap_rna_old)
graphics.off()

############################## UMAPs before integration #######################################
# UMAPs of RNA and ATAC data, with RNA coloured by cell state and ATAC by clusters
umap_atac <- plotEmbedding(ArchR, name = "clusters", plotAs = "points", size = 1.8, baseSize = 0, labelSize = 8, legendSize = 0)

png(paste0(plot_path, 'UMAPs_before_integration_new_scHelper_cell_states.png'), height = 20, width = 40, units = 'cm', res = 400)
print(umap_rna_new + umap_atac)
graphics.off()

png(paste0(plot_path, 'UMAPs_before_integration_old_scHelper_cell_states.png'), height = 20, width = 40, units = 'cm', res = 400)
print(umap_rna_old + umap_atac)
graphics.off()

print("Pre-integration plots made.")

################################################################################################
############################## Unconstrained integration #######################################

print("Starting unconstrained integration...")

ArchR <- addGeneIntegrationMatrix(
  ArchRProj = ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seurat_data,
  addToArrow = TRUE,
  groupRNA = "scHelper_cell_type",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  force = TRUE,
  threads = 1
)
print("integration completed")

# use matched RNA cells to add new and old labels to ATAC cells
extracted_rna_labels <- seurat_data@meta.data[ArchR$predictedCell_Un, c("scHelper_cell_type", "old_labels")]
ArchR$scHelper_cell_type_new <- extracted_rna_labels[, "scHelper_cell_type"]
ArchR$scHelper_cell_type <- extracted_rna_labels[, "old_labels"]
print("scHelper cell type labels added")

print(head(extracted_rna_labels[, "old_labels"]))

# use matched RNA cells to add rna metadata to ATAC cells
extracted_rna_metadata <- seurat_data@meta.data[ArchR$predictedCell_Un, c("run", "stage", "seurat_clusters")]
ArchR$rna_stage <- extracted_rna_metadata[, "stage"]
ArchR$rna_run <- extracted_rna_metadata[, "run"]
ArchR$rna_clusters <- extracted_rna_metadata[, "seurat_clusters"]
print("RNA metadata added")

# add more general scHelpercelltype labels (ie group together all PPRs, neurals, etc)
scHelper_cell_types <- data.frame(getCellColData(ArchR, select = "scHelper_cell_type"))
broad <- scHelper_cell_types %>% mutate(broad = mapvalues(scHelper_cell_type, 
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
ArchR$scHelper_cell_type_broad <- broad$broad
print("Broad scHelper cell type labels added")

print(head(getCellColData(ArchR)))

# # save integrated ArchR project
# print("Saving integrated ArchR project...")
# paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
# print(paste0("Output filename = ", rds_path, label[1], "_Save-ArchR"))
# saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label[1], "_Save-ArchR"), load = FALSE)
# print("Integrated ArchR project saved.")

#############################################################################################
############################## Post-Integration Plots #######################################

print("Making post-integration plots...")

# set colour palettes for UMAPs
atac_scHelper_new_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_new)]
atac_scHelper_old_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_old)]

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
plot_path = "./plots/after_integration/old_labels/"
dir.create(plot_path, recursive = T)

png(paste0(plot_path, 'UMAP_integrated.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE)
graphics.off()

png(paste0(plot_path, 'UMAP_integrated_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_old_cols)
graphics.off()

### Broad labels
plot_path = "./plots/after_integration/broad_labels/"
dir.create(plot_path, recursive = T)

png(paste0(plot_path, 'UMAP_integrated.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type_broad", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE)
graphics.off()

png(paste0(plot_path, 'UMAP_integrated_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "scHelper_cell_type_broad", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_old_cols)
graphics.off()

############################## Integration scores plots #######################################

plot_path = "./plots/after_integration/"

png(paste0(plot_path, 'Integration_Scores_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "predictedScore_Un", plotAs = "points", size = 1.8, baseSize = 0, 
              legendSize = 10)
graphics.off()

png(paste0(plot_path, "Integration_Scores_Vln.png"), width=50, height=20, units = 'cm', res = 200)
plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", 
  name = "predictedScore_Un", plotAs = "Violin", alpha = 0.4)
graphics.off()

print("Post-integration plots made.")

#############################################################################################
#############################   IDENTIFY  CLUSTERS    #######################################
#############################################################################################

print("ArchR cell state plots and assign identities based on label proportions")


############################## Gene scores plots #######################################
#### compare gene scores with integrated gene exp values

plot_path = "./plots/gene_scores_vs_integrated_gex/"
dir.create(plot_path, recursive = T)

ArchR <- addImputeWeights(ArchR)
print("Impute weights added")

# look for late marker genes
late_markers <- c(
  "GATA3", "DLX5", "SIX1", "EYA2", #PPR
  "MSX1", "TFAP2A", "TFAP2B", #mix
  "PAX7", "CSRNP1", "SNAI2", "SOX10", #NC
  "SOX2", "SOX21" # neural
)

png(paste0(plot_path, 'late_markers_GeneScoreMatrix.png'), height = 40, width = 25, units = 'cm', res = 400)
scHelper::ArchR_FeaturePlotGrid(ArchR, matrix = "GeneScoreMatrix", late_markers)
graphics.off()

png(paste0(plot_path, 'late_markers_GeneIntegrationMatrix.png'), height = 40, width = 25, units = 'cm', res = 400)
scHelper::ArchR_FeaturePlotGrid(ArchR, matrix = "GeneIntegrationMatrix", late_markers)
graphics.off()

##################### Distribution of labels across clusters ##################################

plot_path = "./plots/label_by_cluster_distribution/"
dir.create(plot_path, recursive = T)

# visualise distribution across clusters: table of cell counts
png(paste0(plot_path, 'label_by_cluster_cell_number_table.png'), height = 25, width = 40, units = 'cm', res = 400)
scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters", print_table = TRUE, scHelper_cell_type_order = scHelper_cell_type_order)
graphics.off()

# visualise distribution across clusters: confusion matrix
png(paste0(plot_path, "label_by_cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
scHelper::ArchRCellCountsHeatmap(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters")
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

##################### Label clusters based on thresholds ##################################

plot_path = "./plots/label_by_cluster_distribution/assigned_cluster_labels/"
dir.create(plot_path, recursive = T)

identities <- c()

for(i in 1:ncol(percentage_counts)) {       # for-loop over columns (ie clusters)
  column <- percentage_counts[ , i]
  names(column) <- rownames(percentage_counts)
  
  identity <- ifelse(
    # do any labels pass the minimum threshold? 
    # if NO labels label more than the minimum threshold, OR MORE THAN max_label labels pass the minimum threshold -> mixed
    any(sum(column > opt$min_threshold) == 0 | sum(column > opt$min_threshold) > opt$max_label), "MIXED", # condition - mixed
                     ifelse(
                       sum(column > opt$min_threshold) == 1, names(column[column > opt$min_threshold]), # condition 2 - monolabel
                            paste(names( sort(column[column > opt$min_threshold], decreasing = TRUE) ), collapse='/') # condition 3 - multilabel
                     ) 
  )
  identities <- c(identities, identity)
  
}

# assign cluster labels
cluster_idents <- cbind(colnames(percentage_counts), identities)

png(paste0(plot_path, 'assigned_cluster_idents_table.png'), height = 20, width = 10, units = 'cm', res = 400)
grid.arrange(tableGrob(cluster_idents, rows=NULL, theme = ttheme_minimal()))
graphics.off()

new_labels <- cluster_idents[,2]
names(new_labels) <- cluster_idents[,1]
ArchR$cluster_labels <- mapLabels(ArchR$clusters, newLabels = new_labels)

# # save ArchR object with cluster labels
# paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
# saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label[1], "_Save-ArchR"), load = FALSE)

# plot cluster labels on UMAPs
p1 <- plotEmbedding(ArchR, name = "cluster_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, labelAsFactors = FALSE)
p2 <- plotEmbedding(ArchR, name = "scHelper_cell_type", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE)

png(paste0(plot_path, 'assigned_cluster_idents_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
print(p1)
graphics.off()

png(paste0(plot_path, 'assigned_cluster_idents_UMAP_comparison.png'), height = 20, width = 40, units = 'cm', res = 400)
print(p1 + p2)
graphics.off()