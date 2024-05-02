#!/usr/bin/env Rscript

# Script to use BNM to classify clusters of cells/metacells based on median gene scores
# Have set the clustering resolution manually to be stable but separate metacells as much as possible using clustree
# Have also simplified the options of cell types because there are less metacells than single cells

# Load packages
library(getopt)
library(optparse)
library(parallel)
library(Seurat)
library(dplyr)
library(tibble)
library(scHelper)
library(ggplot2)
library(future)
library(cowplot)
library(clustree)
library(gridExtra)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-i", "--input"), action = "store", type = "character", help = "Name of seurat input file to process", default = "seacells_seurat_processed.RDS"),
  make_option(c("", "--verbose"), action = "store", type = "logical", help = "Verbose", default = TRUE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

########################       CELL STATE COLOURS    ########################################
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

stage_order <- c("HH4", "HH5", "HH6", "HH7", "ss4", "ss8")
############################################################################################

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    # paths for RNA metacells testing locally
    #plot_path = "./local_test_data/state_classification/plots/"
    #rds_path = "./local_test_data/state_classification/rds_files/"
    #data_path = "./local_test_data/test_inputs/test_input_seacells_classify/"
    
    # paths for ATAC metacells testing locally
    plot_path = "./local_test_data/state_classification_ATAC/plots/"
    rds_path = "./local_test_data/state_classification_ATAC/rds_files/"
    data_path = "./local_test_data/test_inputs/test_input_seacells_classify_ATAC/"
    
    # paths for testing on nemo
    #plot_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_RNA_WF/5_Classify_metacells/plots/"
    #rds_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_RNA_WF/5_Classify_metacells/rds_files/"
    #data_path = "./output/NF-downstream_analysis/Processing/ss8/SEACELLS_RNA_WF/4_Process_metacells/"
    #data_path = "./output/NF-downstream_analysis/Processing/ss4/SEACELLS_RNA_WF/4_Process_metacells/"
    #data_path = "./output/NF-downstream_analysis/Processing/HH7/SEACELLS_RNA_WF/4_Process_metacells/"
    #data_path = "./output/NF-downstream_analysis/Processing/HH6/SEACELLS_RNA_WF/4_Process_metacells/"
    #data_path = "./output/NF-downstream_analysis/Processing/HH5/SEACELLS_RNA_WF/4_Process_metacells/"
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 16* 1024^3) # 32gb
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}



########################################################################################################
#                                      Cell state classification                                    #
########################################################################################################

### For classification at a single cell level the following number of clusters and cell states were found:
# HH5 - 10 clusters, 8 cell states
# HH6 - 9 clusters, 8 cell states
# HH7 - 11 clusters, 10 cell states
# ss4 - 12 clusters, 9 cell states
# ss8 - 10 clusters, 9 cell states

# Read in seurat object using --input arg as file name in `.input/rds_files/`
seurat_data <- readRDS(paste0(data_path, "rds_files/", opt$input))
seurat_data

# Convert knowledge matrix to gene list
BNM <- list.files(path = data_path, pattern = "*.csv", full.names = TRUE)
# interactively:
#BNM = "./NF-downstream_analysis/binary_knowledge_matrix_contam.csv"
print(BNM)

cell_state_markers <- read.csv(BNM, row.names = 1) %>% select(!c(evidence))
cell_state_markers <- apply(cell_state_markers, 2, function(x) rownames(cell_state_markers)[x > 0])

cell_states = list(
  HH4 = c('NNE', 'node', 'streak', 'EE', 'eNPB', 'eN', 'eCN', 
          'PGC', 'BI', 'meso', 'endo'),

  HH5 = c('NNE', 'node', 'streak', 'EE', 'eNPB', 'eN', 'eCN',
          'NPB', 'aNPB', 'pNPB', 'NP', 'pNP', 'iNP', 'aNP', 'PPR', 'aPPR', 'pPPR',
          'PGC', 'BI', 'meso', 'endo'),
  
  HH6 = c('NNE', 'node', 'streak', 'eN', 'eCN',
          'NPB', 'aNPB', 'pNPB', 'NP', 'pNP', 'iNP', 'aNP', 'PPR', 'aPPR', 'pPPR',
          'PGC', 'BI', 'meso', 'endo'),

  HH7 = c('pEpi', 'NPB', 'aNPB', 'pNPB', 'NC', 'dNC', 'NP', 'pNP', 'iNP',
          'aNP', 'HB', 'MB', 'FB', 'vFB', 'PPR', 'aPPR', 'pPPR',
          'PGC', 'BI', 'meso', 'endo'),

  ss4 = c('pEpi', 'NPB', 'aNPB', 'pNPB', 'NC', 'dNC', 'NP', 'pNP', 'iNP',
          'aNP', 'HB', 'MB', 'FB', 'vFB', 'PPR', 'aPPR', 'pPPR',
          'PGC', 'BI', 'meso', 'endo'),

  ss8 = c('pEpi', 'NPB', 'aNPB', 'pNPB', 'NC', 'dNC', 'NP', 'pNP', 'iNP',
          'aNP', 'HB', 'MB', 'FB', 'vFB', 'PPR', 'aPPR', 'pPPR',
          'PGC', 'BI', 'meso', 'endo')
)

# simplified cell_states for SEACells:
# cell_states = list(

#   HH5 = c('node', 'streak', 
#           'NNE', 'EE', 
#           'eNPB', 'NPB',
#           'eN', 'eCN', 'NP', 'pNP', 'iNP', 'aNP', 
#           'PPR',
#           'PGC', 'BI', 'meso', 'endo'),
  
#   HH6 = c('node', 'streak', 
#           'NNE',
#           'NPB',
#           'eN', 'eCN','NP', 'pNP', 'iNP', 'aNP', 
#           'PPR',
#           'PGC', 'BI', 'meso', 'endo'),

#   HH7 = c('pEpi', 
#           'NPB', 
#           'NC', 'dNC', 
#           'NP', 'pNP', 'iNP','aNP', 
#           'HB', 'MB', 'FB', 
#           'PPR', 'aPPR', 'pPPR',
#           'PGC', 'BI', 'meso', 'endo'),

#   ss4 = c('pEpi', 
#           'NPB', 
#           'NC', 'dNC',
#           'HB', 'MB', 'FB', 'vFB', 
#           'PPR', 'aPPR', 'pPPR',
#           'PGC', 'BI', 'meso', 'endo'),

#   ss8 = c('pEpi', 
#           'NPB', 
#           'NC', 'dNC', 
#           'NP', 'pNP', 'iNP','aNP', 
#           'HB', 'MB', 'FB', 'vFB', 
#           'PPR',
#           'PGC', 'BI', 'meso', 'endo')
# )

cell_state_markers <- lapply(cell_states, function(x) cell_state_markers[names(cell_state_markers) %in% x])
print(cell_state_markers)

############################################################################
############################    Cluster data   #############################

# Find optimal cluster resolution - want about 10 clusters
png(paste0(plot_path, "clustree.png"), width=70, height=35, units = 'cm', res = 200)
ClustRes(seurat_object = seurat_data, by = 0.4, prefix = "RNA_snn_res.")
graphics.off()

# Run classification using different resolutions for different stages (from looking at clustree)
stage = unique(seurat_data@meta.data$stage)
print(paste0("Stage: ", stage))

if(length(stage) == 1){
  cell_state_markers = cell_state_markers[[stage]]
  cluster_res = list(HH5 = 5, HH6 = 5, HH7 = 5, ss4 = 5, ss8 = 5)[[stage]]
  # cluster res set so same number cell states found in sc: 11 at HH5 D, 10 at HH6 D, 14 at HH7, 13 at 4ss D and 13 at 8ss D
  metadata = c('scHelper_cell_type')
} else {
  cell_state_markers = flatten(cell_state_markers)
  cell_state_markers = cell_state_markers[!duplicated(cell_state_markers)]
  metadata = c('scHelper_cell_type', 'stage')
}

print(paste0("Cluster_res: ", cluster_res))
DefaultAssay(seurat_data) <- "RNA"
seurat_data <- FindClusters(seurat_data, resolution = cluster_res)

png(paste0(plot_path, "clusters_UMAP.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = 'seurat_clusters', label = TRUE, 
        label.size = ifelse(length(unique(seurat_data$stage)) == 1, 9, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(seurat_data$stage)) == 1, 6, 6), 
        shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

df <- as.data.frame(table(seurat_data@meta.data$seurat_clusters))
colnames(df) <- c("Cluster", "nCells")
png(paste0(plot_path, 'cluster_counts.png'), height = 40, width = 10, units = 'cm', res = 400)
grid.arrange(top=textGrob("", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

####################################################################################
############################    Plot gene expression   #############################

# Set RNA to default assay for plotting expression data
DefaultAssay(seurat_data) <- "RNA"

# Print how many cell state markers are not in the data
print(paste0("Number of cell state markers: ", length(unlist(cell_state_markers))))

# Filter cell state markers to remove those that aren't in seurat object
cell_state_markers <- lapply(cell_state_markers, function(x) x[x %in% rownames(seurat_data)])
print(paste0("Number of cell state markers in seurat object: ", length(unlist(cell_state_markers))))
print(cell_state_markers)

cell_type_df <- lapply(cell_state_markers, function(x) t(GetAssayData(object = seurat_data, assay = 'RNA', slot = 'scale.data'))[,x] %>% as.data.frame(.) %>% rowSums(.)) %>%
  do.call('cbind', .) %>%
  merge(., seurat_data@meta.data[,'seurat_clusters', drop=FALSE], by=0, all=TRUE)

cell_type_df <- cell_type_df %>% column_to_rownames('Row.names') %>%
  pivot_longer(cols = !seurat_clusters) %>%
  group_by(seurat_clusters, name)

# Plot cell classification per cluster
dir.create(paste0(plot_path, "scHelper_log/"))
png(paste0(plot_path, "scHelper_log/classification_boxplots.png"), width = 30, height = length(unique(cell_type_df$seurat_clusters)) * 3, units = 'cm', res = 200)
ggplot(cell_type_df, aes(x = name, y = value, fill = name)) +
  geom_boxplot() +
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(cell_type_df$name)))) +
  facet_wrap(~seurat_clusters, ncol = 2) +
  xlab(NULL) + 
  ylab('Average scaled expression') +
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = 10), 
        axis.text = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 10)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
graphics.off()

# Select top cell type per cluster
cell_type_df <- cell_type_df %>%
  summarise(value = median(value)) %>%
  group_by(seurat_clusters) %>%
  filter(value == max(value))

# add cell_type to seurat metadata
seurat_data@meta.data[['scHelper_cell_type']] <- unlist(apply(seurat_data@meta.data, 
                                                              1, function(x) ifelse(x[['seurat_clusters']] %in% cell_type_df[['seurat_clusters']], 
                                                                                    cell_type_df[cell_type_df[['seurat_clusters']] == x[['seurat_clusters']], "name"], NA)))

#####################   Set levels
seurat_data@meta.data$scHelper_cell_type <- factor(seurat_data@meta.data$scHelper_cell_type, levels = scHelper_cell_type_order)
seurat_data@meta.data$stage <- factor(seurat_data@meta.data$stage, levels = stage_order)

#####################   Set colours
scHelper_cols <- scHelper_cell_type_colours[levels(droplevels(seurat_data@meta.data$scHelper_cell_type))]

# UMAP for cell state
png(paste0(plot_path, "scHelper_celltype_umap.png"), width=12, height=12, units = 'cm', res = 200)
DimPlot(seurat_data, group.by = 'scHelper_cell_type', label = TRUE, 
        label.size = ifelse(length(unique(seurat_data$stage)) == 1, 9, 3),
        label.box = TRUE, repel = TRUE,
        pt.size = ifelse(length(unique(seurat_data$stage)) == 1, 6, 6), 
        cols = scHelper_cols, shuffle = TRUE) +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "none", 
                 plot.title = element_blank())
graphics.off()

saveRDS(seurat_data, paste0(rds_path, "Classified_metacells.RDS"), compress = FALSE)


# Plot stacked violins for each of the cell type classes to check genes used are good markers
curr_plot_path = paste0(plot_path, "cell_type_dotplots/")
dir.create(curr_plot_path)

for(i in names(cell_state_markers)){
  png(paste0(curr_plot_path, i, ".png"), width = (length(cell_state_markers[[i]])+2)*3, height = 15, units = 'cm', res = 200)
  print(DotPlot(seurat_data, features = cell_state_markers[[i]], group.by = "seurat_clusters"))
  graphics.off()
}


# Plot multi-feature plots for each of the cell type classes
curr_plot_path = paste0(plot_path, "cell_type_feature_plots/")
dir.create(curr_plot_path)

for(i in names(cell_state_markers)){
  ncol = ceiling((length(cell_state_markers[[i]])+1)/8)+1
  nrow = ceiling((length(cell_state_markers[[i]])+1)/ncol)

  png(paste0(curr_plot_path, i, '.png'), width = ncol*10, height = nrow*10, units = "cm", res = 200)
  MultiFeaturePlot(seurat_data, plot_stage = TRUE, stage_col = "stage", plot_celltype = TRUE, celltype_col = "scHelper_cell_type",
                   gene_list = cell_state_markers[[i]], n_col = ncol, label = '')
  graphics.off()
}
