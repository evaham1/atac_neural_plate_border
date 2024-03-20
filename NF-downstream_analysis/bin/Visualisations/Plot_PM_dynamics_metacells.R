#!/usr/bin/env Rscript

print("Calculate the dynamics of PM accessiblity across metacells by latent time weighted by lineage probs (GAMs)")

############################## Load libraries #######################################
library(optparse)
library(future)
library(pheatmap)
library(tidyverse)
library(Antler)
library(RColorBrewer)
library(scHelper)
library(ComplexHeatmap) # Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic data. DOI: 10.1093/bioinformatics/btw313
library(data.table)
library(Seurat)

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
    
    # interactive local paths
    #data_path = "./local_test_data/peak_count_matrix_to_cluster/"
    #plot_path = "./local_test_data/clustered_peaks/plots"
    #rds_path = "./clustered_peaks/rds_files"
    
    # # interactive NEMO paths
    # data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/"
    # # output paths:
    # rds_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/rds_files/"
    # plot_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/plots/"
    # PMs_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/2_peak_clustering/PMs/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  dir.create(PMs_path, recursive = T)
}

########################       CELL STATE COLOURS    ########################################
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'Non-neural',
                              'PPR', 'aPPR', 'pPPR', 'Placodal',
                              'eNPB', 'NPB', 'aNPB', 'pNPB',
                              'NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'Neural',
                              'node', 'streak', 'PGC', 'BI', 'meso', 'endo', 'Contam',
                              'MIXED', 'Unmapped')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                "#A5718D", "#3F918C", "#ed5e5f", "#9792A3",
                                "#7C8483", "#EAEAEA")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                       'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo',
                                       'Neural', 'Placodal', 'Non-neural', 'Contam',
                                       'MIXED', 'Unmapped')
########################       STAGE COLOURS     ###########################################
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
names(stage_colours) <- stage_order
############################################################################################


########################################################################################################
#                                 Read in data and clean up                               #
########################################################################################################

########## COMBINED SEACELL METADATA ############# - with average latent time and lineage probabilities 

print("Reading in SEACell metadata...")

metadata <- read.csv(paste0(data_path, "Combined_SEACell_integrated_metadata.csv"), row.names = 'ATAC')

# Add stage to metadata using SEACell IDs
substrRight <- function(x, n){
  sapply(x, function(xx)
    substr(xx, (nchar(xx)-n+1), nchar(xx))
  )
}
metadata <- metadata %>% mutate(stage = substrRight(rownames(metadata), 3))
metadata <- metadata[,-1]

# Change cell names to match matrix
rownames(metadata) <- gsub('-', '_', rownames(metadata))

# Check metadata
print(head(metadata))

print("Metadata read in!")

########## NORMALISED COUNTS MATRIX ############# - are these scaled?

print("Reading in normalised count data...")

# read in SEACells data
SEACells_normalised_summarised <- fread(paste0(data_path, "Filtered_normalised_summarised_counts.csv"), header = TRUE)
print("Normalised data read in!")

# Extract SEACell IDs from first column
SEACells_IDs <- SEACells_normalised_summarised$V1
length(SEACells_IDs)

# Check for duplicates
table(duplicated(SEACells_IDs))

# Clean up df
SEACells_normalised_summarised <- SEACells_normalised_summarised[,-1]

# Turn into numeric matrix for downstream processing
SEACells_normalised_summarised_numeric <- as.matrix(sapply(SEACells_normalised_summarised, as.numeric))

# Add SEACell IDs as rownames
rownames(SEACells_normalised_summarised_numeric) <- SEACells_IDs

# change cell names for Antler
rownames(SEACells_normalised_summarised_numeric) <- gsub('-', '_', rownames(SEACells_normalised_summarised_numeric))

# Check resulting matrix
print(dim(SEACells_normalised_summarised_numeric))
print("Preview of summarised count df:")
print(SEACells_normalised_summarised_numeric[1:4, 1:4])

# Overwrite cleaned data
SEACells_normalised_summarised <- SEACells_normalised_summarised_numeric

print("Normalised data read in!")


########## PEAK MODULES ############# 



########################################################################################################
#                                 Plot GAMs                               #
########################################################################################################

#### !!! edit the below code which was alex's to use for seurat RNA, instead for matrix directly ATAC

# DefaultAssay(seurat_data) <- "RNA"

# # Iteratively get expression data for each gene module and bind to tidy dataframe
# plot_data <- data.frame()
# for(module in names(gms)){
#   temp <- GetAssayData(seurat_data, assay = 'RNA', slot = 'scale.data')
#   temp <- temp[1:10, 1:10]
  
#   seurat_data <- AddMetaData(seurat_data, metadata = runif(nrow(seurat_data@meta.data)), col.name = "latent_time")
#   seurat_data <- AddMetaData(seurat_data, metadata = runif(nrow(seurat_data@meta.data)), col.name = "lineage_NC_probability")
#   seurat_data <- AddMetaData(seurat_data, metadata = runif(nrow(seurat_data@meta.data)), col.name = "lineage_neural_probability")
#   seurat_data <- AddMetaData(seurat_data, metadata = runif(nrow(seurat_data@meta.data)), col.name = "lineage_placodal_probability")
  
#   temp <- merge(t(temp), seurat_data@meta.data[,c('latent_time', 'lineage_NC_probability', 'lineage_neural_probability', 'lineage_placodal_probability'), drop=FALSE], by=0)
#   plot_data <- temp %>%
#     column_to_rownames('Row.names') %>%
#     pivot_longer(!c(latent_time, lineage_NC_probability, lineage_neural_probability, lineage_placodal_probability)) %>%
#     dplyr::rename(scaled_expression = value) %>%
#     dplyr::rename(gene = name) %>%
#     pivot_longer(cols = !c(latent_time, gene, scaled_expression)) %>%
#     # dplyr::mutate(module = module) %>%
#     dplyr::rename(lineage_probability = value) %>%
#     dplyr::rename(lineage = name) %>%
#     dplyr::group_by(lineage) %>%
#     dplyr::mutate(lineage = unlist(strsplit(lineage, '_'))[2]) %>%
#     bind_rows(plot_data) %>%
#     ungroup()
# }

## final plot data should look like this, and be way more rows than cells originally
plot_data
# A tibble: 65,220 × 5
   latent_time gene               scaled_expression lineage  lineage_probability
         <dbl> <chr>                          <dbl> <chr>                  <dbl>
 1       0.368 ENSGALG00000054818             0     NC                     0.166
 2       0.368 ENSGALG00000054818             0     neural                 0.838
 3       0.368 ENSGALG00000054818             0     placodal               0.865
 4       0.368 ENSGALG00000053455            -0.113 NC                     0.166
 5       0.368 ENSGALG00000053455            -0.113 neural                 0.838
 6       0.368 ENSGALG00000053455            -0.113 placodal               0.865
 7       0.368 ENSGALG00000045540            -0.175 NC                     0.166
 8       0.368 ENSGALG00000045540            -0.175 neural                 0.838
 9       0.368 ENSGALG00000045540            -0.175 placodal               0.865
10       0.368 ENSGALG00000051297            -0.294 NC                     0.166
# ℹ 65,210 more rows
# ℹ Use `print(n = ...)` to see more rows

### then plotting is something like this:

  plot = ggplot(filter(plot_data), aes(x = latent_time, y = scaled_expression)) +
    geom_smooth(method="gam", se=FALSE, mapping = aes(weight = lineage_probability, linetype = lineage, group=lineage)) +
    xlab("Latent time") + ylab("Scaled expression") +
    # facet_wrap(~lineage, dir = 'v') +
    theme_classic()
  
  print(plot)