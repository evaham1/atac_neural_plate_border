#!/usr/bin/env Rscript

print("Calculate Peak modules using Antler")

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
    data_path = "./local_test_data/peak_count_matrix_to_cluster/"
    plot_path = "./local_test_data/clustered_peaks/plots"
    rds_path = "./clustered_peaks/rds_files"
    
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
}

########################################################################################################
#                                 Read in data and clean up                               #
########################################################################################################

# read in SEACells data
SEACells_summarised <- fread(paste0(data_path, "Filtered_summarised_counts.csv"), header = TRUE)
print("Data read in!")

# Check input
print("Preview of input df:")
print(SEACells_summarised[1:4, 1:4])
print(dim(SEACells_summarised))

# Extract SEACell IDs from first column
SEACells_IDs <- SEACells_summarised$V1
print(head(SEACells_IDs))
length(SEACells_IDs)

# Clean up df
SEACells_summarised <- SEACells_summarised[,-1]
print("Preview of input df after cleanup:")
print(SEACells_summarised[1:4, 1:4])
dim(SEACells_summarised)

# Turn into numeric matrix for downstream processing
SEACells_summarised_numeric <- as.matrix(sapply(SEACells_summarised, as.numeric))  

# Add SEACell IDs as rownames
rownames(SEACells_summarised_numeric) <- SEACells_IDs

# change cell names for Antler
rownames(SEACells_summarised_numeric) <- gsub('-', '_', rownames(SEACells_summarised_numeric))
rownames(SEACells_summarised_numeric) <- gsub('#', '', rownames(SEACells_summarised_numeric))

# Check resulting matrix
print(dim(SEACells_summarised_numeric))
print("Preview of summarised count df:")
print(SEACells_summarised_numeric[1:4, 1:4])

# Overwrite cleaned data
SEACells_summarised <- SEACells_summarised_numeric

print("Data read in!")

########################################################################################################
#                                 Generate Antler Inputs                               #
########################################################################################################

# generate fake metadata needed for Antler
pheno_data <- data.frame(row.names = rownames(SEACells_summarised),
                          "timepoint" = rep(1, nrow(SEACells_summarised)),
                          "treatment" = rep("null", nrow(SEACells_summarised)),
                          "replicate_id" = rep(1, nrow(SEACells_summarised))
)

# create antler folder
antler_path = "./antler/"
dir.create(antler_path)

# save pheno data
write.table(pheno_data, file = paste0(antler_path, "phenoData.csv"), row.names = T, sep = "\t", col.names = T)

# save count data
write.table(t(SEACells_summarised), file = paste0(antler_path, "assayData.csv"), row.names = T, sep = "\t", col.names = T, quote = F)

########################################################################################################
#                            Load Antler data and generate correlation matrix                          #
########################################################################################################

# Create Antler object
antler_data <- Antler$new(output_folder = plot_path, num_cores = ncores)
antler_data$load_dataset(folder_path = antler_path)

# Normalize data
antler_data$normalize(method = 'CPM')

###########################
# Calculate GMs unbiasedly
antler_data$gene_modules$identify(
  name                  = "unbiasedPMs",
  corr_t                = 0.3,  # the Spearman correlation threshold
  corr_min              = 3,    # min. number of genes a gene must correlate with
  mod_min_cell          = 10,   # min. number of cells expressing the module
  mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
  process_plots         = TRUE)

# rename peak modules
names(antler_data$gene_modules$lists$unbiasedPMs$content) <- paste0("GM", 1:length(antler_data$gene_modules$lists$unbiasedPMs$content))

# how many peak modules were generated
print(paste0("Number of peak modules made: ", length(antler_data$gene_modules$lists$unbiasedPMs$content)))

########## DE PMs ############## <- NEED TO WORK ON THIS ONCE HAVE MY METADATA
# # Subset peak modules that have at least 50% of peaks DE > 0.25 logFC & FDR < 0.001 between scHelperCellTypes
# gms <- DEGeneModules(seurat_data, antler_data$gene_modules$get("unbiasedGMs"), logfc = 0.25, pval = 0.001, selected_gene_proportion = 0.5, active_ident = meta_col,
#                      ident_1 = c('aPPR', 'pPPR', 'PPR'), ident_2 = c('NC', 'dNC'))
# 
# # save unbiasedGMs_DE in antler object
# antler_data$gene_modules$set(name= "unbiasedGMs_DE", content = gms)

########################################################################################################
#                                        Visualisations                                                #
########################################################################################################

### look through Alex's scripts and adapt the heatmap functions so can plot them here



