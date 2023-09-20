#!/usr/bin/env Rscript

print("Transfer latent time and lineage probabilities from RNA object to ATAC object")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(Seurat)
library(SummarizedExperiment)

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
    addArchRThreads(threads = 1)
    
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Single_cell_integration/" # ATAC
    data_path = "./output/NF-downstream_analysis/rna_objects/" # RNA
    rds_path = "./output/NF-downstream_analysis/Processing/FullData/Transfer_latent_time/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/"
    rds_path = "./rds_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

set.seed(42)

############################## Read in data #######################################

# read in ATAC data ( in input/rds_files)
ArchR <- loadArchRProject(path = paste0(data_path, "rds_files/", "FullData_Save-ArchR"), force = FALSE, showLogo = TRUE)
paste0("Memory Size = ", round(object.size(ArchR) / 10^6, 3), " MB")

# read in RNA data (in input)
RNA_obj <- readRDS(paste0(data_path, "seurat_label_transfer_latent_time.RDS"))


############################## Prepare for data transfer #######################################

# extract metadata from RNA object
RNA_metadata <- FetchData(RNA_obj, vars = c("latent_time", "lineage_NC_probability", "lineage_neural_probability", "lineage_placodal_probability"))
RNA_metadata <- rownames_to_column(RNA_metadata, var = "predictedCell_Un")
head(RNA_metadata)
dim(RNA_metadata) # 17992     5

# extract the matched RNA cell IDs from the ATAC object
ATAC_metadata <- getCellColData(ArchR, select = "predictedCell_Un")
ATAC_metadata <- rownames_to_column(as.data.frame(ATAC_metadata), var = "ATAC_Cell_ID")
head(ATAC_metadata)
dim(ATAC_metadata) # 86217     2

# match the two dataframes using the matched RNA cell IDs
merged_metadata <- merge(ATAC_metadata, RNA_metadata, by = "predictedCell_Un", all.x = TRUE)
head(merged_metadata)
dim(merged_metadata) # 86217     6

# make new column indicating if that cell is not part of the latent time calculation
sum(is.na(merged_metadata$latent_time)) # 4316 cells NA for latent time
merged_metadata <- merged_metadata %>% 
  mutate(in_trajectory = ifelse(is.na(latent_time), FALSE, TRUE)) %>%
  mutate_all(~ifelse(is.na(.), 0, .))

############################## Add data to ArchR object #######################################

# add this metadata to ArchR using matched RNA cell ID
ArchR <- addCellColData(ArchR, data = merged_metadata$latent_time, name = "rna_latent_time", cells = merged_metadata$ATAC_Cell_ID, force = T)
ArchR <- addCellColData(ArchR, data = merged_metadata$lineage_neural_probability, name = "rna_lineage_neural_probability", cells = merged_metadata$ATAC_Cell_ID, force = T)
ArchR <- addCellColData(ArchR, data = merged_metadata$lineage_NC_probability, name = "rna_lineage_NC_probability", cells = merged_metadata$ATAC_Cell_ID, force = T)
ArchR <- addCellColData(ArchR, data = merged_metadata$lineage_placodal_probability, name = "rna_lineage_placodal_probability", cells = merged_metadata$ATAC_Cell_ID, force = T)
ArchR <- addCellColData(ArchR, data = merged_metadata$in_trajectory, name = "in_rna_lineage", cells = merged_metadata$ATAC_Cell_ID, force = T)

############################## Plots #######################################

# plot these transferred labels

png(paste0(plot_path, 'UMAP_latent_time.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "rna_latent_time", plotAs = "points", size = 2, baseSize = 0, 
              labelSize = 0, legendSize = 0, labelAsFactors = FALSE, pal = ArchRPalettes$beach,
              highlightCells = merged_metadata[merged_metadata$in_trajectory == TRUE, "ATAC_Cell_ID"]) + 
  theme_ArchR(legendTextSize = 12, baseSize = 16, plotMarginCm = 0.5)
graphics.off()

png(paste0(plot_path, 'UMAP_lineage_neural_probability.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "rna_lineage_neural_probability", plotAs = "points", size = 2, baseSize = 0, 
              labelSize = 0, legendSize = 0, labelAsFactors = FALSE, pal = ArchRPalettes$beach,
              highlightCells = merged_metadata[merged_metadata$in_trajectory == TRUE, "ATAC_Cell_ID"]) + 
  theme_ArchR(legendTextSize = 12, baseSize = 16, plotMarginCm = 0.5)
graphics.off()

png(paste0(plot_path, 'UMAP_lineage_NC_probability.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "rna_lineage_NC_probability", plotAs = "points", size = 2, baseSize = 0, 
              labelSize = 0, legendSize = 0, labelAsFactors = FALSE, pal = ArchRPalettes$beach,
              highlightCells = merged_metadata[merged_metadata$in_trajectory == TRUE, "ATAC_Cell_ID"]) + 
  theme_ArchR(legendTextSize = 12, baseSize = 16, plotMarginCm = 0.5)
graphics.off()

png(paste0(plot_path, 'UMAP_lineage_placodal_probability.png'), height = 20, width = 20, units = 'cm', res = 400)
plotEmbedding(ArchR, name = "rna_lineage_placodal_probability", plotAs = "points", size = 2, baseSize = 0, 
              labelSize = 0, legendSize = 0, labelAsFactors = FALSE, pal = ArchRPalettes$beach,
              highlightCells = merged_metadata[merged_metadata$in_trajectory == TRUE, "ATAC_Cell_ID"]) + 
  theme_ArchR(legendTextSize = 12, baseSize = 16, plotMarginCm = 0.5)
graphics.off()


############################## Save #######################################

# save ArchR object
saveArchRProject(ArchRProj = ArchR, outputDirectory = paste0(rds_path, "TransferLabel_Save-ArchR"), load = FALSE)

############################## Trajectory exploration plots #######################################

# ArchR <- addImputeWeights(ArchR)
# 
# trajMM  <- getTrajectory(ArchR, name = "rna_lineage_placodal_probability", useMatrix = "GeneScoreMatrix", log2Norm = FALSE)
# p1 <- plotTrajectoryHeatmap(trajMM, varCutOff = 0.1, pal = paletteContinuous(set = "solarExtra"))
# p1[[1]]
# 
# trajMM  <- getTrajectory(ArchR, name = "rna_lineage_neural_probability", useMatrix = "GeneScoreMatrix", log2Norm = FALSE)
# p1 <- plotTrajectoryHeatmap(trajMM, varCutOff = 0.1, pal = paletteContinuous(set = "solarExtra"))
# p1[[1]]
# 
# p1 <- plotTrajectory(ArchR, trajectory = "rna_lineage_placodal_probability", colorBy = "cellColData", name = "rna_lineage_placodal_probability", 
#                      addArrow = F)
# p1[[1]]
# 
# 
# trajMM  <- getTrajectory(ArchR, name = "rna_lineage_placodal_probability", useMatrix = "GeneScoreMatrix", log2Norm = FALSE)
# p2 <- plotTrajectoryHeatmap(trajMM,  pal = paletteContinuous(set = "horizonExtra"))
# 
# p2 <- plotTrajectory(ArchR, trajectory = "rna_lineage_placodal_probability", colorBy = "GeneScoreMatrix", name = "SIX1", continuousSet = "blueYellow", addArrow = F)
# p2[[1]]
# p2[[2]]
# 
# p2 <- plotTrajectory(ArchR, trajectory = "rna_lineage_placodal_probability", colorBy = "cellColData", name = "rna_lineage_placodal_probability", continuousSet = "blueYellow")
# 
# p2 <- plotTrajectory(ArchR, trajectory = "rna_latent_time", colorBy = "cellColData", name = "rna_latent_time", continuousSet = "blueYellow", addArrow = FALSE)
# p2[[1]]
# 
# 
# plotEmbedding(ArchR, colorBy = "GeneScoreMatrix", name = "SIX1", plotAs = "points", size = 2, baseSize = 0,
#               labelSize = 0, legendSize = 0, labelAsFactors = FALSE)
# 
# ArchR <- addImputeWeights(ArchR, seed = 1)