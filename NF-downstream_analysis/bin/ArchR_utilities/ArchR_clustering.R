#!/usr/bin/env Rscript

print("clustering ArchR")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(hexbin)
library(pheatmap)
library(gridExtra)
library(grid)
library(parallel)
library(clustree)
library(presto)
library(Seurat)
library(gtools)
library(ComplexHeatmap)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
    make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
    make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
    make_option(c("", "--stage_clust_res"), action = "store", type = "double", help = "clustering resolution for stage data", default = 1),
    make_option(c("", "--full_clust_res"), action = "store", type = "double", help = "clustering resolution for full data", default = 2),
    make_option(c("", "--clustree_stage"), action = "store", type = "logical", help = "whether to run clustree plot on stage data", default = FALSE),
    make_option(c("", "--clustree_full"), action = "store", type = "logical", help = "whether to run clustree plot on full data", default = FALSE),
    make_option(c("", "--stage_clustree_by"), action = "store", type = "double", help = "clustering res intervals for clustree for stages", default = 0.1),
    make_option(c("", "--full_clustree_by"), action = "store", type = "double", help = "clustering res intervals for clustree for full data", default = 0.2),
    make_option(c("", "--GeneScore_heatmaps_stage"), action = "store", type = "logical", help = "whether to run gene score heatmaps on stage data", default = FALSE),
    make_option(c("", "--GeneScore_heatmaps_full"), action = "store", type = "logical", help = "whether to run gene score heatmaps on full data", default = FALSE),
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
    
    # already clustered
    data_path = "./output/NF-downstream_analysis/Processing/ss8/1_clustered/rds_files/"
    
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

ArchR_ClustRes <- function(ArchR, starting_res = 0, by = 0.1){
  plots <- list()
  resolutions <- c(seq(starting_res, 8*by + starting_res, by=by))
  print(paste0("resolutions: ", resolutions))
  cluster_df <- data.frame(ArchR$cellNames)
  
  if(length(ArchR@reducedDims) == 0){stop("Carry out dimensionality reduction before clustering")}
  
  for(res in resolutions[1:length(resolutions)]){
    print(paste0("resolution: ", res))
    ArchR_clustered <- addClusters(input = ArchR, name = "clusters", force = TRUE, resolution = res)
    plots[[paste(res)]] <- plotEmbedding(ArchR_clustered, name = "clusters") +
      ggtitle(paste("resolution = ", res))
    title <- paste0("clustering_res_", res)
    cluster_df[, title] <- ArchR_clustered@cellColData$clusters
  }
  plots[["clustree"]] <- clustree(cluster_df, prefix = "clustering_res_")
  
  lay <- rbind(c(10,10,10,1,2,3),
               c(10,10,10,4,5,6),
               c(10,10,10,7,8,9))
  plots2 <- gridExtra::arrangeGrob(grobs = plots, layout_matrix = lay)
  
  return(plots2)
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
        tidyr::pivot_wider(names_from = data_group2, values_from = Freq) %>% 
        tibble::add_column(ID = !!i)
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

# Feature plot function to create grid of feature plots
feature_plot_grid <- function(ArchRProj = ArchR, matrix = "GeneScoreMatrix", gene_list) {
  p <- plotEmbedding(ArchRProj, colorBy = matrix, name = gene_list, 
                     plotAs = "points", size = 1.8, baseSize = 0, labelSize = 8, legendSize = 10)
  p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
      )
  })
  do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
}

### Function to add unique ids to se peak object so can subset properly
add_unique_ids_to_se <- function(seMarker, ArchR, matrix_type) {
  
  if (matrix_type == "PeakMatrix") {
    tmp_peaks = data.frame(ArchR@peakSet)
    tmp_diff_peaks = data.frame(rowData(seMarker))
    diff_peaks_join_peakset = left_join(tmp_diff_peaks, tmp_peaks, 
                                        by = c("seqnames" = "seqnames", "start" = "start", "end" = "end"))
    diff_peaks_join_peakset$gene_name = paste(diff_peaks_join_peakset$nearestGene, diff_peaks_join_peakset$distToTSS,sep="_")
    diff_peaks_join_peakset$unique_id = paste0(diff_peaks_join_peakset$seqnames, ":", diff_peaks_join_peakset$start, "-", diff_peaks_join_peakset$end)
    rowData(seMarker) = diff_peaks_join_peakset
  }
  
  if (matrix_type == "GeneScoreMatrix") {
    rowData <- as.data.frame(rowData(seMarker))
    
    duplicated_gene_names <- rowData$name[duplicated(rowData$name)]
    duplicated_genes <- rowData[which(rowData$name %in% duplicated_gene_names), ]
    duplicated_genes <- duplicated_genes %>% group_by(name) %>% 
      arrange(name) %>% mutate(unique_id = paste0(name, "-", rowid(name)))
    join_df = left_join(rowData, duplicated_genes,
                        by = c("seqnames" = "seqnames", "start" = "start", "end" = "end", "name" = "name", "idx" = "idx", "strand" = "strand"))
    
    join_df <- join_df %>% mutate(unique_id = coalesce(unique_id, name))
    
    rowData(seMarker) = join_df
  }
  
  return(seMarker)
}

### Function to extract means from se object into matrix for plotting
extract_means_from_se <- function(seMarker) {
  mat <- as.data.frame(SummarizedExperiment::assays(seMarker)[["Mean"]])
  rownames(mat) <- rowData(seMarker)$unique_id
  
  return(mat)
}

### Function to Log2normalise matrix (must be run before subsetting!)
Log2norm <- function(mat, scaleTo = 10^4) {
  mat <- log2(t(t(mat)/colSums(mat)) * scaleTo + 1) # normalising means for depth of cluster
  return(mat)
}

### Function to extract IDs to be plotted, either by cut off or cut off + top n features
extract_ids <- function(seMarker, cutOff = "FDR <= 1 & Log2FC >= 0", top_n = TRUE, n = 10, group_name = "clusters") {
  
  markerList <- getMarkers(seMarker, cutOff = cutOff) # extract features that pass threshold
  
  df <- data.frame() # merged all features into a df
  for (i in 1:length(names(markerList))) {
    print(i)
    df_i <- as.data.frame(markerList[i])
    df <- rbind(df, df_i)
  }
  
  if (top_n == FALSE){
    ids <- df$unique_id
  } else {
    df <- df %>%
      group_by(group_name) %>%
      top_n(n, Log2FC) %>%
      dplyr::arrange(Log2FC, .by_group = TRUE)
    ids <- unique(df$unique_id)
  }
  
  return(ids)
}

### Function to subset normalised matrix using IDs
subset_matrix <- function(mat, ids) {
  subsetted_matrix <- mat[ids, ]
  return(subsetted_matrix)
}

### Function to plot marker heatmap
marker_heatmap <- function(mat, pal = NULL, 
                           labelRows = FALSE, fontSizeRows = 12,
                           labelCols = TRUE, fontSizeCols = 12,
                           cluster_columns = TRUE, cluster_rows = TRUE) {
  
  # scale each feature independently and add min/max limits
  limits <- c(-2, 2) # could make this user-defined
  mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), 
               `/`)
  mat[mat > max(limits)] <- max(limits)
  mat[mat < min(limits)] <- min(limits)
  
  # colours - set default if NULL
  if (is.null(pal) == TRUE) {
    pal <- paletteContinuous(set = "solarExtra", n = 100)
  }
  
  # # order rows by eucladian distance
  # dist_mat <- dist(mat, method = 'euclidean')
  # hclust_avg <- hclust(dist_mat, method = 'average')
  # ordered_features <- hclust_avg$labels[c(hclust_avg$order)]
  # mat <- mat[match(ordered_features, rownames(mat)), ]
  # 
  # # order columns by eucladian distance
  # dist_mat <- dist(t(mat), method = 'euclidean')
  # hclust_avg <- hclust(dist_mat, method = 'average')
  # ordered_cell_groups <- hclust_avg$labels[c(hclust_avg$order)]
  # mat <- mat[ , match(ordered_cell_groups, colnames(mat))]
  
  Heatmap(
    matrix = mat,
    col = pal,
    heatmap_legend_param = list(title = "z-scores"),
    #top_annotation = topAnno, 
    # add raster stuff?
    
    #Column Options
    cluster_columns = cluster_columns,
    show_column_names = labelCols,
    column_names_gp = gpar(fontsize = fontSizeCols),
    column_names_max_height = unit(100, "mm"),
    #column_split = colData$stage,
    
    #Row Options
    cluster_rows = cluster_rows,
    show_row_names = labelRows,
    row_names_gp = gpar(fontsize = fontSizeRows)
    #row_split = row_split_params
  )
  
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

# Plot Clustree (optional)
if (length(unique(ArchR$stage)) == 1){
  if (isTRUE(opt$clustree_stage)) {
    print("running clustree plot...")
    png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
      plot(ArchR_ClustRes(ArchR, by = opt$stage_clustree_by))
      graphics.off()
    print("clustree plot ran") }

} else {
  if (isTRUE(opt$clustree_full)) {
    print("running clustree plot...")
    png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
      plot(ArchR_ClustRes(ArchR, by = opt$full_clustree_by))
      graphics.off()
    print("clustree plot ran") }
}

################## Save clustered ArchR project #################################
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label, "_Save-ArchR"), load = FALSE)
print("ArchR object saved")

#######################################################################################
############################ CELL COUNTS PER CLUSTER ##################################
plot_path_temp = "./plots/cell_counts/"
dir.create(plot_path_temp, recursive = T)

# Plot number of cells in each cluster
cluster_cell_counts <- as.data.frame(table(substr(ArchR$clusters, 2, nchar(ArchR$clusters))))
cluster_cell_counts <- cluster_cell_counts %>% 
  dplyr::rename(Cell_count = Freq, Cluster_number = Var1) %>%
  dplyr::mutate(Cluster_number = as.numeric(as.character(Cluster_number))) %>%
  dplyr::arrange(Cluster_number)
print("Cluster cell counts: ")
print(cluster_cell_counts)

cluster_cell_counts_totals <- cluster_cell_counts %>%
  rbind(c("Total", sum(cluster_cell_counts$Cell_count)))
print("Cluster cell counts with totals: ")
print(cluster_cell_counts_totals)
png(paste0(plot_path_temp, 'cluster_cell_counts_table.png'), height = 25, width = 10, units = 'cm', res = 400)
grid.arrange(tableGrob(cluster_cell_counts_totals, rows=NULL, theme = ttheme_minimal()))
graphics.off()

p<-ggplot(data=cluster_cell_counts, aes(x=`Cluster_number`, y=`Cell_count`)) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks = round(seq(min(cluster_cell_counts$Cluster_number), max(cluster_cell_counts$Cluster_number), by = 1),1))

png(paste0(plot_path_temp, 'cluster_cell_counts_barchart.png'), height = 10, width = 20, units = 'cm', res = 400)
print(p)
graphics.off()

# Plot contribution of each stage to each cluster
if (length(unique(ArchR$stage)) > 1){
  png(paste0(plot_path_temp, "cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
  cell_counts_heatmap(ArchR = ArchR, group1 = "clusters", group2 = "stage")
  graphics.off()
}

print("cell counts calculated")

###############################################################################################
################################### cluster UMAPS #############################################
plot_path_temp = "./plots/UMAPs/"
dir.create(plot_path_temp, recursive = T)

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

png(paste0(plot_path_temp, "UMAPs.png"), width=60, height=40, units = 'cm', res = 200)
ggAlignPlots(p1, p2, type = "h")
graphics.off()

png(paste0(plot_path_temp, 'UMAP_clusters.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "clusters", plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1), baseSize = 0, 
              labelSize = 10, legendSize = 0, randomize = TRUE, labelAsFactors = FALSE)
graphics.off()

if ( !(is.null(ArchR$stage_clusters)) ) {
  
  png(paste0(plot_path_temp, 'UMAP_stage_clusters.png'), height = 20, width = 20, units = 'cm', res = 400)
  plotEmbedding(ArchR, name = "stage_clusters", plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1), baseSize = 0, 
              labelSize = 10, legendSize = 0, randomize = TRUE, labelAsFactors = FALSE)
  graphics.off()

}

#################################################################################
############################### QC PLOTS ########################################
plot_path_temp = "./plots/QC_plots/"
dir.create(plot_path_temp, recursive = T)

quantiles = c(0.2, 0.8)

##### nFrags
p <- plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", alpha = 0.4,
  name = "nFrags", plotAs = "Violin", baseSize = 12)
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = "nFrags")[,1], probs = quantiles[1]), linetype = "dashed",
                   color = "red")
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = "nFrags")[,1], probs = quantiles[2]), linetype = "dashed",
                   color = "red")
png(paste0(plot_path_temp, "VlnPlot_thresholds_nFrags.png"), width=50, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

#### TSS Enrichment
p <- plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", alpha = 0.4,
                name = "TSSEnrichment", plotAs = "Violin", baseSize = 12)
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = "TSSEnrichment")[,1], probs = quantiles[1]), linetype = "dashed",
                   color = "red")
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = "TSSEnrichment")[,1], probs = quantiles[2]), linetype = "dashed",
                   color = "red")
png(paste0(plot_path_temp, "VlnPlot_thresholds_TSSEnrichment.png"), width=50, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

#### Nucleosome signal
p <- plotGroups(ArchR, groupBy = "clusters", colorBy = "cellColData", alpha = 0.4,
                name = "NucleosomeRatio", plotAs = "Violin", baseSize = 12)
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = "NucleosomeRatio")[,1], probs = quantiles[1]), linetype = "dashed",
                   color = "red")
p = p + geom_hline(yintercept = quantile(getCellColData(ArchR, select = "NucleosomeRatio")[,1], probs = quantiles[2]), linetype = "dashed",
                   color = "red")
png(paste0(plot_path_temp, "VlnPlot_thresholds_NucleosomeRatio.png"), width=50, height=20, units = 'cm', res = 200)
print(p)
graphics.off()

print("QC plots done")

#################################################################################
############################ GENE SCORE PLOTS ###################################

plot_path_temp = "./plots/Gene_score_plots/"
dir.create(plot_path_temp, recursive = T)

##########    Feature plots

ArchR <- addImputeWeights(ArchR)

# Contaminating markers
contaminating_markers <- c(
  'DAZL', #PGC
  'CDH5', 'TAL1', 'HBZ', # Blood island
  'CDX2', 'GATA6', 'ALX1', 'PITX2', 'TWIST1', 'TBXT', 'MESP1', #mesoderm
  'SOX17', 'CXCR4', 'FOXA2', 'NKX2-2', 'GATA6' #endoderm
)
# Late marker genes
late_markers <- c(
  "GATA3", "DLX5", "SIX1", "EYA2", #PPR
  "MSX1", "TFAP2A", "TFAP2B", #mix
  "PAX7", "CSRNP1", "SNAI2", "SOX10", #NC
  "SOX2", "SOX21" # neural
)
# look for ap marker genes
ap_markers <- c(
  "PAX2", "WNT4", "SIX3", "SHH" # no GBX2 in matrix
)
# look for early markers
early_markers <- c(
  "EPAS1", "BMP4", "YEATS4", "SOX3", "HOXB1", "ADMP", "EOMES"
)
feature_plot_genes <- c("SIX1", "PAX7", "DLX5", "CSRNP1", "SOX10",
                        "SOX21", "SOX2", "BMP4", "HOXB1")

png(paste0(plot_path_temp, 'Contaminating_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = contaminating_markers)
graphics.off()

png(paste0(plot_path_temp, 'Late_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = late_markers)
graphics.off()

png(paste0(plot_path_temp, 'AP_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = ap_markers)
graphics.off()

png(paste0(plot_path_temp, 'Early_markers_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = early_markers)
graphics.off()

png(paste0(plot_path_temp, 'Useful_FeaturePlots.png'), height = 25, width = 25, units = 'cm', res = 400)
feature_plot_grid(ArchR, gene_list = feature_plot_genes)
graphics.off()

print("Feature plots done")

##########    Heatmaps (optional)

run_heatmaps <- ifelse(length(unique(ArchR$stage)) == 1 & isTRUE(opt$GeneScore_heatmaps_stage) | length(unique(ArchR$stage)) > 1 & isTRUE(opt$GeneScore_heatmaps_full),
       TRUE, FALSE)

if (isTRUE(run_heatmaps)) {

  seMarker <- getMarkerFeatures(
    ArchRProj = ArchR, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "clusters")
  seMarker <- add_unique_ids_to_se(seMarker, ArchR, matrix_type = "GeneScoreMatrix")

  # prepare for plotting
  matrix <- extract_means_from_se(seMarker) # extract means df from se object
  normalised_matrix <- Log2norm(matrix) # log2norm across all features in each cell group
  
  pal = viridis::magma(100)

  # Heatmap of positive markers which pass cutoff thresholds
  ids <- extract_ids(seMarker, cutOff = "FDR <= 0.01 & Log2FC >= 1", top_n = FALSE) # extract ids
  if (length(ids) < 2){
    print(paste0(length(ids), " features passed cutoff - not enough to make heatmap"))
  } else {
    print(paste0(length(ids), " features passed cutoff - now plotting heatmap"))
    subsetted_matrix <- subset_matrix(normalised_matrix, ids) # subset matrix to only include features of interest
    
    png(paste0(plot_path_temp, 'diff_cutoff_heatmap.png'), height = 40, width = 20, units = 'cm', res = 400)
    print(marker_heatmap(subsetted_matrix, pal = pal))
    graphics.off()
  }

  # Heatmap of positive markers top 10 per cell group
  ids <- extract_ids(seMarker, cutOff = "FDR <= 0.05 & Log2FC >= 0", top_n = TRUE, n = 10) # extract ids
  subsetted_matrix <- subset_matrix(normalised_matrix, ids) # subset matrix to only include features of interest
  
  png(paste0(plot_path_temp, 'diff_top10_heatmap.png'), height = 40, width = 20, units = 'cm', res = 400)
  marker_heatmap(subsetted_matrix, labelRows = TRUE, pal = pal, cluster_columns = FALSE, cluster_rows = FALSE)
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

if ( !(is.null(ArchR$scHelper_cell_type_old)) ) {

  atac_scHelper_old_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_old)]

  plot_path_temp <- "./plots/RNA_label_plots/"
  dir.create(plot_path_temp, recursive = T)
  
  png(paste0(plot_path_temp, 'UMAP_scHelper_cell_type.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE))
  graphics.off()

  png(paste0(plot_path_temp, 'UMAP_scHelper_cell_type_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_old_cols))
  graphics.off()

}

if ( !(is.null(ArchR$stage_scHelper_cell_type_old)) ) {

  atac_scHelper_old_cols <- scHelper_cell_type_colours[unique(ArchR$stage_scHelper_cell_type_old)]

  plot_path_temp <- "./plots/RNA_label_plots/"
  dir.create(plot_path_temp, recursive = T)
  
  png(paste0(plot_path_temp, 'UMAP_stage_scHelper_cell_type.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "stage_scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE))
  graphics.off()

  png(paste0(plot_path_temp, 'UMAP_stage_scHelper_cell_type_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "stage_scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_old_cols))
  graphics.off()

}

if ( !(is.null(ArchR$cluster_labels)) ) {

  plot_path_temp <- "./plots/RNA_label_plots/"
  dir.create(plot_path_temp, recursive = T)
  
  png(paste0(plot_path_temp, 'UMAP_cluster_labels.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "cluster_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE))
  graphics.off()

  png(paste0(plot_path_temp, 'UMAP_cluster_labels_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "cluster_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_old_cols))
  graphics.off()

}

if ( !(is.null(ArchR$stage_cluster_labels)) ) {

  plot_path_temp <- "./plots/RNA_label_plots/"
  dir.create(plot_path_temp, recursive = T)
  
  png(paste0(plot_path_temp, 'UMAP_stage_cluster_labels.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "stage_cluster_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE))
  graphics.off()

  png(paste0(plot_path_temp, 'UMAP_stage_cluster_labels_nolabel.png'), height = 20, width = 20, units = 'cm', res = 400)
  print(plotEmbedding(ArchR, name = "stage_cluster_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = atac_scHelper_old_cols))
  graphics.off()

}