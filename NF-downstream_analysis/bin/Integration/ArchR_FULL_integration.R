#!/usr/bin/env Rscript

print("ArchR integration, assigning cluster labels and peak to gene linkage")

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

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getPeakSet(ArchR)
getAvailableMatrices(ArchR)

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

# ############################## Integration scores plots #######################################

# plot_path = "./plots/after_integration/integration_scores/"
# dir.create(plot_path, recursive = T)

# png(paste0(plot_path, 'Integration_Scores_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
# plotEmbedding(ArchR, name = "predictedScore_Un", plotAs = "points", size = 1.8, baseSize = 0, 
#               legendSize = 10)
# graphics.off()

# png(paste0(plot_path, "Integration_Scores_Vln.png"), width=50, height=20, units = 'cm', res = 200)
# plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", 
#   name = "predictedScore_Un", plotAs = "Violin", alpha = 0.4)
# graphics.off()

# print("Post-integration plots made.")

# ############################## Gene Integration Count Plots for key TFs #######################################

# plot_path = "./plots/gene_integration_plots_of_key_TFs/"
# dir.create(plot_path, recursive = T)

# # set genes of interest
# TFs <- c("SIX1", "IRF6", "DLX5", "DLX6", "GATA2", "GATA3", "TFAP2A", "TFAP2B", "TFAP2C", "PITX1", "PITX2",
#          "PAX7", "MSX1", "ETS1", "SOX9", "SOX8", "SOX10", "SOX5", "SOX21", "NKX6-2")
# # CTNRP and LMX1B and ZEB2 not found

# ArchR <- addImputeWeights(ArchR)

# # Plot ridge plot of each TF deviation
# for (TF in TFs){
#   print(TF)
  
#   p <- plotGroups(ArchR, groupBy = "clusters", 
#                   colorBy = "GeneIntegrationMatrix", 
#                   name = TF,
#                   imputeWeights = getImputeWeights(ArchR))
  
#   # Plot distribution of GeneIntegration values for each cluster
#   png(paste0(plot_path, TF, '_gene_integration_ridge_plot.png'), height = 12, width = 10, units = 'cm', res = 400)
#   print(p)
#   graphics.off()
  
#   # Plot GeneIntegration values on UMAP
#   p <- plotEmbedding(ArchR, colorBy = "GeneIntegrationMatrix", name = TF, embedding = "UMAP", continuousSet = "blueYellow", 
#                      imputeWeights = getImputeWeights(ArchR), plotAs = "points", size = 1.8,)
#   png(paste0(plot_path, TF, '_gene_integration_UMAP.png'), height = 12, width = 10, units = 'cm', res = 400)
#   print(p)
#   graphics.off()
  
# }

# ############################## Gene scores plots #######################################
# #### compare gene scores with integrated gene exp values

# plot_path = "./plots/gene_scores_vs_integrated_gex_marker_genes/"
# dir.create(plot_path, recursive = T)

# ArchR <- addImputeWeights(ArchR)
# print("Impute weights added")

# # look for late marker genes
# late_markers <- c(
#   "GATA3", "DLX5", "SIX1", "EYA2", #PPR
#   "MSX1", "TFAP2A", "TFAP2B", #mix
#   "PAX7", "CSRNP1", "SNAI2", "SOX10", #NC
#   "SOX2", "SOX21" # neural
# )

# png(paste0(plot_path, 'late_markers_GeneScoreMatrix.png'), height = 40, width = 25, units = 'cm', res = 400)
# scHelper::ArchR_FeaturePlotGrid(ArchR, matrix = "GeneScoreMatrix", late_markers)
# graphics.off()

# png(paste0(plot_path, 'late_markers_GeneIntegrationMatrix.png'), height = 40, width = 25, units = 'cm', res = 400)
# scHelper::ArchR_FeaturePlotGrid(ArchR, matrix = "GeneIntegrationMatrix", late_markers)
# graphics.off()

# # look for early marker genes
# early_markers <- c(
#   "EPAS1", "BMP4", "YEATS4", "SOX3",
#   "HOXB1", "EOMES", "ADMP"
# )

# png(paste0(plot_path, 'early_markers_GeneScoreMatrix.png'), height = 20, width = 25, units = 'cm', res = 400)
# scHelper::ArchR_FeaturePlotGrid(ArchR, matrix = "GeneScoreMatrix", early_markers)
# graphics.off()

# png(paste0(plot_path, 'early_markers_GeneIntegrationMatrix.png'), height = 20, width = 25, units = 'cm', res = 400)
# scHelper::ArchR_FeaturePlotGrid(ArchR, matrix = "GeneIntegrationMatrix", early_markers)
# graphics.off()

# ##################### Distribution of labels across clusters ##################################

# plot_path = "./plots/label_by_cluster_distribution/"
# dir.create(plot_path, recursive = T)

# # visualise distribution across clusters: table of cell counts
# png(paste0(plot_path, 'label_by_cluster_cell_number_table.png'), height = 25, width = 40, units = 'cm', res = 400)
# scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters", print_table = TRUE, scHelper_cell_type_order = scHelper_cell_type_order)
# graphics.off()

# # visualise distribution across clusters: confusion matrix
# png(paste0(plot_path, "label_by_cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
# scHelper::ArchRCellCountsHeatmap(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters")
# graphics.off()

# # visualise distribution across clusters: table of cell percentages
# cell_counts <- scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
# percentage_counts <- as.data.frame(lapply(cell_counts, function(x) (x / sum(x))*100))
# rownames(percentage_counts) <- rownames(cell_counts)

# png(paste0(plot_path, 'label_by_cluster_cell_percentage_table.png'), height = 25, width = 40, units = 'cm', res = 400)
# grid.arrange(tableGrob(round(percentage_counts, 2), theme = ttheme_minimal()))
# graphics.off()

# # visualise distribution across clusters: piecharts
# counts <- scHelper::ArchRCellCounting(ArchR = ArchR, group1 = "scHelper_cell_type", group2 = "clusters", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
# png(paste0(plot_path, "label_by_cluster_piecharts.png"), width=50, height=40, units = 'cm', res = 200)
# scHelper::CellLabelPieCharts(counts, col = scHelper_cell_type_colours)
# graphics.off()


# ###############################################################################################
# ############################## CO-ACCESSIBILITY BETWEEN PEAKS #################################

# print("Calculating coaccessibility...")

# # calculate co-accessibility between all peaks
# ArchR <- addCoAccessibility(ArchR)

# # extract interactions - returns indexes of queryHits and subjectHits
# cA <- getCoAccessibility(ArchR, corCutOff = 0.5, returnLoops = FALSE)
# cA
#   # DataFrame with 120270 rows and 11 columns
#   # queryHits subjectHits seqnames correlation Variability1 Variability2     TStat        Pval         FDR VarQuantile1 VarQuantile2
#   # <integer>   <integer>    <Rle>   <numeric>    <numeric>    <numeric> <numeric>   <numeric>   <numeric>    <numeric>    <numeric>
#   #   1              3           4     chr1    0.548725   0.00437754   0.00683964   14.5441 4.15759e-40 4.52151e-38     0.911185     0.965430
#   # 2              4           3     chr1    0.548725   0.00683964   0.00437754   14.5441 4.15759e-40 4.52151e-38     0.965430     0.911185
#   # 3              4           5     chr1    0.517190   0.00683964   0.00356568   13.3901 4.49249e-35 3.64027e-33     0.965430     0.870967
#   # 4              5           4     chr1    0.517190   0.00356568   0.00683964   13.3901 4.49249e-35 3.64027e-33     0.870967     0.965430
#   # 5             27          40     chr1    0.761607   0.01690577   0.00855042   26.0418 1.47916e-94 2.12498e-91     0.995825     0.978303
# coacessibility_df <- as.data.frame(cA)

# # Need to use indices from df to extract granges and therefore informative peak IDs
# coacessibility_df <- coacessibility_df %>% 
#   mutate(query_PeakID = paste0(seqnames(metadata(cA)[[1]][queryHits]), "-", start(metadata(cA)[[1]][queryHits]), "-", end(metadata(cA)[[1]][queryHits]))) %>%
#   mutate(subject_PeakID = paste0(seqnames(metadata(cA)[[1]][subjectHits]), "-", start(metadata(cA)[[1]][subjectHits]), "-", end(metadata(cA)[[1]][subjectHits])))

# head(coacessibility_df)

# # sanity check that all interaction Peak IDs are in the ArchR peakset
# table(coacessibility_df$subject_PeakID %in% getPeakSet(ArchR)$name)

# # save df
# write.csv(coacessibility_df, file = paste0(csv_path, label[1], "_peak_coaccessibility_df.csv"), row.names = FALSE)

# print("Coaccessibility calculated and saved.")

# ####################################################################################################################
# ############################## CO-ACCESSIBILITY BETWEEN PEAKS AND GENES - FIXED DIST #################################

# print("Calculating coaccessibility between peaks and genes within a 250000 distance...")

# # calculate gene-to-peak co-accessibility using GeneIntegrationMatrix
# ArchR <- addPeak2GeneLinks(ArchR, maxDist = 250000)
# # biggest chrom chrom 1 size: 197608386 (200000000), default is 250000
# # tried running at max distance but then found very few interactions very far away...

# # extract resulting interactions
# p2g <- getPeak2GeneLinks(ArchR, corCutOff = 0.5, returnLoops = FALSE)
# p2g_df <- as.data.frame(p2g)

# # add correct Peak IDs and gene names to df
# p2g_df <- p2g_df %>% 
#   mutate(PeakID = paste0(seqnames(metadata(p2g)$peakSet[idxATAC]), "-", start(metadata(p2g)$peakSet[idxATAC]), "-", end(metadata(p2g)$peakSet[idxATAC]))) %>%
#   mutate(gene_name = metadata(p2g)$geneSet[idxRNA]$name)
# head(p2g_df)
# print(paste0(nrow(p2g_df), " interactions identified by coaccessibility!"))

# # sanity check that all interaction Peak IDs are in the ArchR peakset
# if(sum(p2g_df$PeakID %in% getPeakSet(ArchR)$name) != nrow(p2g_df)){stop("Issue with peak IDs in interactions!")}

# # save df
# write.csv(p2g_df, file = paste0(csv_path, label[1], "_peak_to_gene_linkage_df_250000_distance.csv"), row.names = FALSE)

# print("Coaccessibility between peaks and genes (250000 dist) calculated and saved.")

##################################################################################
############################## SAVE ARCHR OBJECT #################################

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getAvailableMatrices(ArchR)

# save integrated ArchR project
print("Saving integrated ArchR project...")
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
print(paste0("Output filename = ", rds_path, label[1], "_Save-ArchR"))
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label[1], "_Save-ArchR"), load = FALSE)
print("Integrated ArchR project saved.")

# ################################################################################
# ############################## HEATMAPS OF P2L #################################

# plot_path = "./plots/peak2gene_250000_dist/"
# dir.create(plot_path, recursive = T)

# # This sporadically fails so comment out for now, not v useful plot anyway
# ## Heatmap of linkage across clusters
# p <- plotPeak2GeneHeatmap(ArchRProj = ArchR, groupBy = "clusters")
# png(paste0(plot_path, 'Peak_to_gene_linkage_clusters_heatmap.png'), height = 80, width = 60, units = 'cm', res = 400)
# print(p)
# graphics.off()

# p <- plotPeak2GeneHeatmap(ArchRProj = ArchR, groupBy = "stage")
# png(paste0(plot_path, 'Peak_to_gene_linkage_stage_heatmap.png'), height = 80, width = 60, units = 'cm', res = 400)
# print(p)
# graphics.off()