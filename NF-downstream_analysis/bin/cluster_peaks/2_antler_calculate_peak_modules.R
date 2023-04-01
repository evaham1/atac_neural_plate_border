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
    
    # output from SEACells - summarised by metacells
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/SEACells/4_exported_SEACells_data/rds_files/"
    label = "summarised_counts_1000.csv"
    plot_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/Peak_modules/plots/"
    rds_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/Peak_modules/rds_files/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    label = "AnnData_summarised_by_metacells_peak_counts.csv"
    ncores = opt$cores
    
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

########################################################################################################
#                                 Read in data and create Antler inputs                               #
########################################################################################################

# read in SEACells data
SEACells_summarised <- fread(paste0(data_path, "Filtered_summarised_counts.csv"), header = TRUE)
#SEACells_summarised <- as.matrix(fread(paste0(data_path, "ss8_summarised_by_metacells_counts.csv"), header = TRUE), rownames = 1)
print("Data read in!")

# Check input
print("Preview of input df:")
print(SEACells_summarised[1:4, 1:4])
print(dim(SEACells_summarised))

## make into numeric matrix adjust colnames/rownmames

# change cell names for Antler
rownames(SEACells_summarised) <- gsub('-', '_', rownames(SEACells_summarised))
rownames(SEACells_summarised) <- gsub('#', '', rownames(SEACells_summarised))

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
if (is.null(names(antler_data$gene_modules$lists$unbiasedPMs$content))) {
  names(antler_data$gene_modules$lists$unbiasedPMs$content) <- paste0("GM", 1:length(antler_data$gene_modules$lists$unbiasedPMs$content))
}

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



