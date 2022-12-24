#!/usr/bin/env Rscript

print("ArchR cell state plots and assign identities based on label proportions")

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
library(gtools)

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

    data_path = "./output/NF-downstream_analysis/Processing/ss8/2_unconstrained_integration/rds_files/"
    
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


############################### FUNCTIONS ####################################

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
  

############################## Read in ArchR project  #######################################

# Retrieve object label
label <- sub('_.*', '', list.files(data_path))
print(label)

# load ArchR object using its retrieved name
ArchR <- loadArchRProject(path = paste0(data_path, label, "_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")


############################################################################################
############################## COLOURS #######################################

###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

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
# set colour palettes for UMAPs
atac_scHelper_old_cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_old)]


############################## Gene scores plots #######################################
#### compare gene scores with integrated gene exp values

ArchR <- addImputeWeights(ArchR)

# look for late marker genes
late_markers <- c(
  "GATA3", "DLX5", "SIX1", "EYA2", #PPR
  "MSX1", "TFAP2A", "TFAP2B", #mix
  "PAX7", "CSRNP1", "SNAI2", "SOX10", #NC
  "SOX2", "SOX21" # neural
)
p1 <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = late_markers, 
  embedding = "UMAP"
)
p2 <- plotEmbedding(
  ArchRProj = ArchR, 
  colorBy = "GeneIntegrationMatrix", 
  name = late_markers, 
  embedding = "UMAP"
)

p1c <- lapply(p1, function(x){
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

p2c <- lapply(p2, function(x){
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

png(paste0(plot_path, 'late_markers_GeneScoreMatrix.png'), height = 40, width = 25, units = 'cm', res = 400)
do.call(cowplot::plot_grid, c(list(ncol = 3), p1c))
graphics.off()

png(paste0(plot_path, 'late_markers_GeneIntegrationMatrix.png'), height = 40, width = 25, units = 'cm', res = 400)
do.call(cowplot::plot_grid, c(list(ncol = 3), p2c))
graphics.off()


##################### Distribution of labels across clusters ##################################

plot_path = "./plots/label_by_cluster_distribution/"
dir.create(plot_path, recursive = T)

# visualise distribution across clusters: table of cell counts
png(paste0(plot_path, 'label_by_cluster_cell_number_table.png'), height = 25, width = 40, units = 'cm', res = 400)
cell_counting(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "clusters", print_table = TRUE, scHelper_cell_type_order = scHelper_cell_type_order)
graphics.off()

# visualise distribution across clusters: confusion matrix
png(paste0(plot_path, "label_by_cluster_distribution.png"), width=25, height=20, units = 'cm', res = 200)
cell_counts_heatmap(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "clusters")
graphics.off()

# visualise distribution across clusters: table of cell percentages
cell_counts <- cell_counting(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "clusters", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
percentage_counts <- as.data.frame(lapply(cell_counts, function(x) x / sum(x)))
rownames(percentage_counts) <- rownames(cell_counts)

png(paste0(plot_path, 'label_by_cluster_cell_percentage_table.png'), height = 25, width = 40, units = 'cm', res = 400)
grid.arrange(tableGrob(round(percentage_counts, 2), theme = ttheme_minimal()))
graphics.off()

# visualise distribution across clusters: piecharts
counts <- cell_counting(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "clusters", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
png(paste0(plot_path, "label_by_cluster_piecharts.png"), width=50, height=40, units = 'cm', res = 200)
cell_counts_piecharts(counts, col = scHelper_cell_type_colours)
graphics.off()

##################### Label clusters based on thresholds ##################################

plot_path = "./plots/label_by_cluster_distribution/assigned_cluster_labels/"
dir.create(plot_path, recursive = T)

min_threshold = 0.15
max_label = 3
identities <- c()

for(i in 1:ncol(percentage_counts)) {       # for-loop over columns
  column <- percentage_counts[ , i]
  names(column) <- rownames(percentage_counts)
  
  identity <- ifelse(
    any(sum(column > min_threshold) == 0 | sum(column > min_threshold) > max_label), "MIXED", # condition - mixed
                     ifelse(
                       sum(column > min_threshold) == 1, names(column[column > min_threshold]), # condition 2 - monolabel
                            paste(names( sort(column[column > min_threshold], decreasing = TRUE) ), collapse='/') # condition 3 - multilabel
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

# save ArchR object with cluster labels
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, label[1], "_Save-ArchR"), load = FALSE)

# plot cluster labels on UMAPs
p1 <- plotEmbedding(ArchR, name = "cluster_labels", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, labelAsFactors = FALSE)
p2 <- plotEmbedding(ArchR, name = "scHelper_cell_type_old", plotAs = "points", size = 1.8, baseSize = 0, 
              labelSize = 8, legendSize = 0, pal = atac_scHelper_old_cols, labelAsFactors = FALSE)

png(paste0(plot_path, 'assigned_cluster_idents_UMAP.png'), height = 20, width = 20, units = 'cm', res = 400)
print(p1)
graphics.off()

png(paste0(plot_path, 'assigned_cluster_idents_UMAP_comparison.png'), height = 20, width = 40, units = 'cm', res = 400)
print(p1 + p2)
graphics.off()