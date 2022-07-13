#!/usr/bin/env Rscript

print("clustering ArchR")

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
library(clustree)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
    make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
    make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
    make_option(c("", "--stage_clust_res"), action = "store", type = "double", help = "clustering resolution for stage data", default = 1),
    make_option(c("", "--full_clust_res"), action = "store", type = "double", help = "clustering resolution for full data", default = 2),
    make_option(c("", "--clustree"), action = "store", type = "logical", help = "whether to run clustree plot", default = TRUE),
    make_option(c("", "--stage_clustree_by"), action = "store", type = "double", help = "clustering res intervals for clustree for stages", default = 0.1),
    make_option(c("", "--full_clustree_by"), action = "store", type = "double", help = "clustering res intervals for clustree for full data", default = 0.2),
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
    
    plot_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ss8_Save-ArchR/ArchR_clustering/plots/"
    rds_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ss8_Save-ArchR/ArchR_clustering/rds_files/"
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ArchR_split/rds_files/"

    # already clustered
    data_path = "./output/NF-downstream_analysis/ArchR_preprocessing/ss8/ArchR_clustering/rds_files/"

    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################### FUNCTIONS - adapted from scHelper #################################################
ArchR_IdentifyOutliers <- function(ArchR, group_by = 'clusters', metrics, intersect_metrics = TRUE, quantiles){
  outlier <- list()
  if(!length(quantiles) == 2){
    stop('quantiles must be an array of length == 2')
  }
  for(metric in metrics){
    min = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1])
    max = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2])
    
    outlier[[metric]] <- as.tibble(getCellColData(ArchR)) %>%
      group_by((!!as.symbol(group_by))) %>%
      summarise(median = median((!!as.symbol(metric)))) %>%
      filter(median > max | median < min) %>%
      pull(!!as.symbol(group_by))
  }
  
  if(intersect_metrics){
    if(length(Reduce(intersect, outlier)) == 0){
      cat('No outliers detected!')
    } else {
      return(Reduce(intersect, outlier))
    }
  } else{
    if(length(as.character(unique(unlist(outlier)))) == 0){
      cat('No outliers detected!')
    } else {
      return(as.character(unique(unlist(outlier))))
    }
  }
}

ArchR_ClustRes <- function(ArchR, by = 0.1, starting_res = 0){
  plots <- list()
  resolutions <- c(seq(starting_res, starting_res+9*by, by=by))
  cluster_df <- data.frame(ArchR$cellNames)
  
  if(length(ArchR@reducedDims) == 0){stop("Carry out dimensionality reduction before clustering")}
  
  for(res in resolutions[2:length(resolutions)]){
    ArchR_clustered <- addClusters(input = ArchR, name = "clusters", force = TRUE, resolution = res)
    plots[[paste(res)]] <- plotEmbedding(ArchR_clustered, name = "clusters") +
      ggtitle(paste("resolution = ", res))
    title <- paste0("clustering_res_", res)
    cluster_df <- cluster_df %>% mutate(!!title := ArchR_clustered@cellColData$clusters)
  }
  
  plots[["clustree"]] <- clustree(cluster_df, prefix = "clustering_res_")
  lay <- rbind(c(1,1,1,2,3,4),
             c(1,1,1,5,6,7),
             c(1,1,1,8,9,10))
  lay <- rbind(c(10,10,10,1,2,3),
               c(10,10,10,4,5,6),
               c(10,10,10,7,8,9))
  plots2 <- gridExtra::arrangeGrob(grobs = plots, layout_matrix = lay)
  return(gridExtra::grid.arrange(plots2))
}

# function to make heatmap showing contribution of cell groups to other cell groups
cell_counts_heatmap <- function(ArchR = ArchR, group1 = "scHelper_cell_type_new", group2 = "clusters") {
  group1_data <- getCellColData(ArchR, select = group1)[,1]
  group2_data <- getCellColData(ArchR, select = group2)[,1]
  cM <- confusionMatrix(paste0(group2_data), paste0(group1_data))
  cM <- cM / Matrix::rowSums(cM)
  
  p <- pheatmap::pheatmap(
    mat = cM,
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
  )
}

############################### FUNCTIONS - custom made #################################################

# function to count how many cells from each cluster/sample are assigned the same label/cluster
cell_counting <- function(ArchR = ArchR, group1 = "clusters", group2 = "stage", print_table = TRUE, scHelper_cell_type_order = scHelper_cell_type_order) {
  group1_data <- getCellColData(ArchR, select = group1)[,1]
  group1_cell_counts <- as.data.frame(table(group1_data))
  colnames(group1_cell_counts) <- c("ID", "Total_count")
  
  group2_cell_counts <- data.frame()
  group2_data <- getCellColData(ArchR, select = group2)[,1]
  data_group1 <- getCellColData(ArchR, select = group1)[,1]
  for (i in unique(group1_data)) {
    cells <- ArchR$cellNames[BiocGenerics::which(data_group1 == i)]
    if (length(cells) > 1){
      ArchR_subset <- ArchR[cells, ]
      data_group2 <- getCellColData(ArchR_subset, select = group2)[,1]
      group2_cell_counts_i <- as.data.frame(table(data_group2)) %>%
        pivot_wider(names_from = data_group2, values_from = Freq) %>% 
        add_column(ID = !!i)
      group2_cell_counts <- rbind.fill(group2_cell_counts, group2_cell_counts_i)
    }
  }
  
  cell_counts <- merge(group1_cell_counts, group2_cell_counts) %>%
    column_to_rownames(., var = "ID")
  totals <- cell_counts$Total_count
  cell_counts <- cell_counts[, -1]
  cell_counts[is.na(cell_counts)] <- 0
  
  # Order rows and columns - group 1 = order rows
  if (group1 %in% c("clusters", "stage")) {
    cell_counts <- cell_counts[ mixedsort(rownames(cell_counts)) , ] }
  if (group1 == "scHelper_cell_type_old") {
    if (is.null(scHelper_cell_type_order)){
      print("scHelper_cell_type_order not specified!")
    } else {
      order <- scHelper_cell_type_order[scHelper_cell_type_order %in% rownames(cell_counts)]
      cell_counts <- cell_counts[ order , ] }
  }
  
  # Order rows and columns - group 2 = order columns
  if (group2 %in% c("clusters", "stage")){
    cell_counts <- cell_counts[ , mixedsort(colnames(cell_counts)) ]
  }
  if (group2 == "scHelper_cell_type_old") {
    if (is.null(scHelper_cell_type_order)){
      print("scHelper_cell_type_order not specified!")
    } else {
      order <- scHelper_cell_type_order[scHelper_cell_type_order %in% colnames(cell_counts)]
      cell_counts <- cell_counts[ , order ] }
  }
  
  # either return table or print it
  if (print_table == FALSE){
    return(cell_counts)
  } else {
    cell_counts <- cell_counts %>% mutate(., Total = totals)
    grid.arrange(tableGrob(cell_counts, theme = ttheme_minimal()))
  }
}

# function to make heatmap showing contribution of cell groups to other cell groups
cell_counts_heatmap <- function(ArchR = ArchR, group1 = "scHelper_cell_type_new", group2 = "clusters") {
  group1_data <- getCellColData(ArchR, select = group1)[,1]
  group2_data <- getCellColData(ArchR, select = group2)[,1]
  cM <- confusionMatrix(paste0(group2_data), paste0(group1_data))
  cM <- cM / Matrix::rowSums(cM)
  
  p <- pheatmap::pheatmap(
    mat = cM,
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
  )
}

# function to make piecharts/barcharts showing contribution of cell groups to other cell groups
cell_counts_piecharts <- function(counts, cols, scale = FALSE) {
  
  # remove any columns that have no cells in
  if (0 %in% colSums(counts)){
    counts <- counts[,-(which(colSums(counts)==0))] }
  
  # scale by number of cells in each row
  if (scale == TRUE){
    count_data <- t(apply(counts, 1, function(x) x/sum(x)))
  } else { 
    count_data <- counts }
  
  # make piecharts
  plots <- list()
  for (i in colnames(counts)){
    print(i)
    # calculate totals to add to plot titles
    raw_data <- data.frame(group = rownames(counts), value = (counts[,i]))
    total <- sum(raw_data$value)
    # extract either scaled or raw data for plotting
    data <- data.frame(group = rownames(count_data), value = (count_data[,i]))
    plot <- ggplot(data, aes(x="", y=value, fill=group)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() + scale_fill_manual(values= cols) +
      ggtitle(paste0(i, " (cells: ", total, ")")) +
      theme(legend.position="none", plot.title = element_text(size = 20, hjust = 0.5, vjust = 0))
    plots[[i]] <- plot
  }
  do.call(grid.arrange,plots)
}

############################## Read in ArchR project #######################################

# If files are not in rds_files subdirectory look in input dir 
label <- sub('_.*', '', list.files(data_path))
print(label) 

if (length(label) == 0){
  data_path = "./input/"
  label <- sub('_.*', '', list.files(data_path))
  print(label)
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
} else {
  ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
  paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
}

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

#################################################################################
############################## PROCESSING #######################################
# Dimensionality reduction
ArchR <- addIterativeLSI(ArchR, force = TRUE)
print("iterative LSI ran")

# Run UMAP
ArchR <- addUMAP(ArchR, force = TRUE)
print("UMAP added")

# Cluster
if (length(unique(ArchR$stage)) == 1){
  ArchR <- addClusters(ArchR, name = "clusters", resolution = opt$stage_clust_res, force = TRUE)
} else {
  ArchR <- addClusters(ArchR, name = "clusters", resolution = opt$full_clust_res, force = TRUE)
}
print("clustering ran")

####### OLD PROCESSING ########
# ArchR <- addIterativeLSI(
#   ArchRProj = ArchR,
#   useMatrix = "TileMatrix", 
#   name = "IterativeLSI", 
#   iterations = 2, 
#   clusterParams = list( #See Seurat::FindClusters
#     resolution = c(0.2), 
#     sampleCells = 10000, 
#     n.start = 10
#   ), 
#   varFeatures = 25000, 
#   dimsToUse = 1:30,
#   seed = 1,
#   force = TRUE
# )

# ArchR <- addUMAP(
#   ArchRProj = ArchR, 
#   reducedDims = "IterativeLSI", 
#   name = "UMAP", 
#   nNeighbors = 30, 
#   minDist = 0.5, 
#   metric = "cosine",
#   force = TRUE
# )
# # Set clustering for stage and full data
# if (length(unique(ArchR$stage)) == 1){
#   ArchR <- addClusters(
#   input = ArchR,
#   reducedDims = "IterativeLSI",
#   method = "Seurat",
#   name = "clusters",
#   resolution = opt$stage_clust_res,
#   force = TRUE)
# } else {
#   ArchR <- addClusters(
#   input = ArchR,
#   reducedDims = "IterativeLSI",
#   method = "Seurat",
#   name = "clusters",
#   resolution = opt$full_clust_res,
#   force = TRUE)
# }

# Plot Clustree
if (isTRUE(opt$clustree)) {
  print("running clustree plot...")

  if (length(unique(ArchR$stage)) == 1){
    png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
    print(ArchR_ClustRes(ArchR, by = opt$stage_clustree_by, starting_res = -opt$stage_clustree_by))
    graphics.off()
    print("clustree plot ran")
} else {
    png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
    print(ArchR_ClustRes(ArchR, by = opt$full_clustree_by, starting_res = -opt$full_clustree_by))
    graphics.off()
    print("clustree plot ran")
}}

################## Save clustered ArchR project #################################
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)

#######################################################################################
############################ CELL COUNTS PER CLUSTER ##################################

# Plot number of cells in each cluster
cluster_cell_counts <- as.data.frame(table(substr(ArchR$clusters, 2, nchar(ArchR$clusters))))
cluster_cell_counts <- cluster_cell_counts %>% 
  rename(Cell_count = Freq, Cluster_number = Var1) %>%
  mutate(Cluster_number = as.numeric(as.character(Cluster_number))) %>%
  arrange(Cluster_number)

png(paste0(plot_path, 'cluster_cell_counts_table.png'), height = 25, width = 10, units = 'cm', res = 400)
grid.arrange(tableGrob(cluster_cell_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

p<-ggplot(data=cluster_cell_counts, aes(x=`Cluster_number`, y=`Cell_count`)) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks = round(seq(min(cluster_cell_counts$Cluster_number), max(cluster_cell_counts$Cluster_number), by = 1),1))

png(paste0(plot_path, 'cluster_cell_counts_barchart.png'), height = 10, width = 20, units = 'cm', res = 400)
print(p)
graphics.off()

# Plot contribution of each stage to each cluster
if (length(unique(ArchR$stage)) > 1){
  png(paste0(plot_path, "cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
  cell_counts_heatmap(ArchR = ArchR, group1 = "clusters", group2 = "stage")
  graphics.off()
}

#######################################################################################
################################### UMAPS #############################################

p1 <- plotEmbedding(ArchR, 
                    name = "stage",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0, 
                    pal = stage_colours, randomize = TRUE)
p2 <- plotEmbedding(ArchR, 
                    name = "clusters",
                    plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                    baseSize = 0, labelSize = 0, legendSize = 0,
                    randomize = TRUE)

png(paste0(plot_path, "UMAPs.png"), width=60, height=40, units = 'cm', res = 200)
ggAlignPlots(p1, p2, type = "h")
graphics.off()

png(paste0(plot_path, 'UMAP_clusters.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "clusters", plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1), baseSize = 0, 
              labelSize = 10, legendSize = 0, randomize = TRUE, labelAsFactors = FALSE)
graphics.off()

#################################################################################
############################### QC PLOTS ########################################

plot_path <- "./plots/QC_plots/"
dir.create(plot_path, recursive = T)

######################## QC Vioin Plots #######################################

quantiles = c(0.2, 0.8)

##### nFrags
p <- plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", 
  name = "nFrags", plotAs = "Violin", baseSize = 12)

metrics = "nFrags"
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1]), linetype = "dashed", 
                   color = "red")
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2]), linetype = "dashed", 
                   color = "red")
png(paste0(plot_path, "VlnPlot_thresholds_nFrags.png"), width=50, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

# automatically identify outlier clusters using adapted scHelper function
outliers <- ArchR_IdentifyOutliers(ArchR, group_by = 'clusters', metrics = metrics, intersect_metrics = FALSE, quantiles = quantiles)

# highlight outlier clusters on UMAP
if (is.null(outliers) == FALSE){
  idxSample <- BiocGenerics::which(ArchR$clusters %in% outliers)
  cellsSample <- ArchR$cellNames[idxSample]
  png(paste0(plot_path, "UMAP_nFrags_outliers.png"), width=20, height=20, units = 'cm', res = 200)
  print(plotEmbedding(ArchR, name = "clusters", highlightCells = cellsSample,
                      plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
                      baseSize = 20, labelSize = 14, legendSize = 0, randomize = TRUE, labelAsFactors = FALSE))
  graphics.off()
}

#### TSS Enrichment
p <- plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", 
  name = "TSSEnrichment",plotAs = "violin",
  alpha = 0.4, addBoxPlot = TRUE, baseSize = 12)

metrics = "TSSEnrichment"
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1]), linetype = "dashed", 
                   color = "red")
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2]), linetype = "dashed", 
                   color = "red")
png(paste0(plot_path, "VlnPlot_thresholds_TSSEnrichment.png"), width=50, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

# automatically identify outlier clusters using adapted scHelper function
outliers <- ArchR_IdentifyOutliers(ArchR, group_by = 'clusters', metrics = metrics, intersect_metrics = FALSE, quantiles = quantiles)

# highlight outlier clusters on UMAP
if (is.null(outliers) == FALSE){
  idxSample <- BiocGenerics::which(ArchR$clusters %in% outliers)
  cellsSample <- ArchR$cellNames[idxSample]
  png(paste0(plot_path, "UMAP_TSSEnrichment_outliers.png"), width=20, height=20, units = 'cm', res = 200)
  print(plotEmbedding(ArchR, name = "clusters", highlightCells = cellsSample,
      plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
      baseSize = 20, labelSize = 14, legendSize = 0, randomize = TRUE, labelAsFactors = FALSE))
  graphics.off()
}

#### Nucleosome signal
p <- plotGroups(ArchR, 
  groupBy = "clusters", colorBy = "cellColData", 
  name = "NucleosomeRatio", plotAs = "violin",
  alpha = 0.4, addBoxPlot = TRUE, baseSize = 12
)

metrics = "NucleosomeRatio"
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1]), linetype = "dashed", 
                   color = "red")
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2]), linetype = "dashed", 
                   color = "red")
png(paste0(plot_path, "VlnPlot_thresholds_NucleosomeRatio.png"), width=50, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

# automatically identify outlier clusters using adapted scHelper function
outliers <- ArchR_IdentifyOutliers(ArchR, group_by = 'clusters', metrics = metrics, intersect_metrics = FALSE, quantiles = quantiles)

# highlight outlier clusters on UMAP
if (is.null(outliers) == FALSE){
  idxSample <- BiocGenerics::which(ArchR$clusters %in% outliers)
  cellsSample <- ArchR$cellNames[idxSample]
  png(paste0(plot_path, "UMAP_NucleosomeRatio_outliers.png"), width=20, height=20, units = 'cm', res = 200)
  print(plotEmbedding(ArchR, name = "clusters", highlightCells = cellsSample,
      plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
      baseSize = 20, labelSize = 14, legendSize = 0, randomize = TRUE, labelAsFactors = FALSE))
  graphics.off()
}

#################################################################################
############################ CELL LABEL PLOTS ###################################

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
atac_scHelper_old_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_old)]

if ( !(is.null(ArchR$scHelper_cell_type_old)) ) {

  plot_path <- "./plots/RNA_label_plots/"
  dir.create(plot_path, recursive = T)
  
  png(paste0(plot_path, 'UMAP_integrated.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE))
  graphics.off()

  png(paste0(plot_path, 'UMAP_integrated_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_old_cols))
  graphics.off()

  ############################## Integration scores plots #######################################

  plot_path = paste0(plot_path, "integration_scores/")
  dir.create(plot_path, recursive = T)

  png(paste0(plot_path, 'Integration_Scores_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
  plotEmbedding(ArchR, name = "predictedScore_Un", plotAs = "points", size = 1.8, baseSize = 0, 
                legendSize = 10)
  graphics.off()

  png(paste0(plot_path, "Integration_Scores_Vln.png"), width=50, height=20, units = 'cm', res = 200)
  plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", 
            name = "predictedScore_Un", plotAs = "Violin")
  graphics.off()

  ##################### Distribution of labels across clusters ##################################

  plot_path = paste0(plot_path, "label_by_cluster_distribution/")
  dir.create(plot_path, recursive = T)

  # visualise distribution across clusters: table
  png(paste0(plot_path, 'label_by_cluster_table.png'), height = 25, width = 40, units = 'cm', res = 400)
  cell_counting(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "clusters", print_table = TRUE, scHelper_cell_type_order = scHelper_cell_type_order)
  graphics.off()

  # visualise distribution across clusters: confusion matrix
  png(paste0(plot_path, "label_by_cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
  cell_counts_heatmap(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "clusters")
  graphics.off()

  # visualise distribution across clusters: piecharts
  counts <- cell_counting(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "clusters", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
  png(paste0(plot_path, "label_by_cluster_piecharts.png"), width=50, height=40, units = 'cm', res = 200)
  cell_counts_piecharts(counts, col = scHelper_cell_type_colours)
  graphics.off()

  # assign cluster labels
  cM <- as.matrix(confusionMatrix(ArchR$clusters, ArchR$scHelper_cell_type_old))
  scHelper_cell_types <- colnames(cM)[apply(cM, 1 , which.max)]
  cluster_idents <- cbind(scHelper_cell_types, rownames(cM))

  png(paste0(plot_path, 'assigned_cluster_idents_table.png'), height = 20, width = 10, units = 'cm', res = 400)
  grid.arrange(tableGrob(cluster_idents, rows=NULL, theme = ttheme_minimal()))
  graphics.off()

  new_labels <- cluster_idents[,1]
  names(new_labels) <- cluster_idents[,2]
  ArchR$cluster_old_labels <- mapLabels(ArchR$clusters, newLabels = new_labels)

  p1 <- plotEmbedding(ArchR, name = "cluster_old_labels", plotAs = "points", size = 1.8, baseSize = 0, 
                labelSize = 8, legendSize = 0, pal = atac_scHelper_new_cols, labelAsFactors = FALSE)
  p2 <- plotEmbedding(ArchR, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
                labelSize = 8, legendSize = 0, pal = atac_scHelper_new_cols, labelAsFactors = FALSE)

  png(paste0(plot_path, 'assigned_cluster_idents_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(p1)
  graphics.off()

  png(paste0(plot_path, 'assigned_cluster_idents_UMAP_comparison.png'), height = 20, width = 40, units = 'cm', res = 400)
  print(p1 + p2)
  graphics.off()

  ##################### Distribution of labels across stages ##################################

  if (length(unique(ArchR$stage)) > 1){

    plot_path = paste0(plot_path, "labels_by_stage_distribution/")
    dir.create(plot_path, recursive = T)
  
    png(paste0(plot_path, 'counts_by_stage_table.png'), height = 25, width = 40, units = 'cm', res = 400)
    cell_counting(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "stage", scHelper_cell_type_order = scHelper_cell_type_order)
    graphics.off()
  
    # visualise distribution across stages: confusion matrix
    png(paste0(plot_path, "stage_distribution.png"), width=25, height=20, units = 'cm', res = 200)
    cell_counts_heatmap(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "stage")
    graphics.off()
  
    # visualise distribution across stages: piecharts
    counts <- cell_counting(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "stage", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
    png(paste0(plot_path, "label_by_stage_piecharts_unscaled.png"), width=50, height=40, units = 'cm', res = 200)
    cell_counts_piecharts(counts, col = scHelper_cell_type_colours)
    graphics.off()
  
    png(paste0(plot_path, "label_by_stage_piecharts_scaled.png"), width=50, height=40, units = 'cm', res = 200)
    cell_counts_piecharts(counts, col = scHelper_cell_type_colours, scale = TRUE)
    graphics.off()
  
    ##################### Distribution of stages across labels ##################################
  
    # visualise distribution across stages: piecharts
    counts <- cell_counting(ArchR = ArchR, group1 = "stage", group2 = "scHelper_cell_type_old", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
    png(paste0(plot_path, "stage_by_label_piecharts_unscaled.png"), width=50, height=40, units = 'cm', res = 200)
    cell_counts_piecharts(counts, col = stage_colours)
    graphics.off()
  
    png(paste0(plot_path, "stage_by_label_piecharts_scaled.png"), width=50, height=40, units = 'cm', res = 200)
    cell_counts_piecharts(counts, col = stage_colours, scale = TRUE)
    graphics.off()
  
    ##################### Distribution of rna stages across atac stages ##################################

    rna_stages <- plotEmbedding(ArchR, name = "rna_stage", plotAs = "points", size = 1.8, baseSize = 0, 
                              labelSize = 8, legendSize = 0, labelAsFactors = FALSE, pal = stage_colours)
    atac_stages <- plotEmbedding(ArchR, name = "stage", plotAs = "points", size = 1.8, baseSize = 0, 
                               labelSize = 8, legendSize = 0, labelAsFactors = FALSE, pal = stage_colours)
    png(paste0(plot_path, 'UMAPs_rna_stages_VS_atac_stages.png'), height = 20, width = 40, units = 'cm', res = 400)
    print(rna_stages + atac_stages)
    graphics.off()

    png(paste0(plot_path, "rna_atac_stage_distribution.png"), width=25, height=20, units = 'cm', res = 200)
    cell_counts_heatmap(ArchR = ArchR, group1 = "rna_stage", group2 = "stage")
    graphics.off()

  }


}