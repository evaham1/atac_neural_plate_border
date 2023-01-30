##### Putting the SEACell assignments onto ArchR and checking their purity
# also adding metadata to the SEACells themselves and exporting this for visualisations of peak modules

# load libraries
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(hexbin)
library(gridExtra)
library(grid)
library(parallel)
library(data.table)
library(doSNOW)

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
    
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/SEACells/4_exported_SEACells_data/rds_files/"
    label = "summarised_counts_1000.csv"
    
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

######    Set colours
###### stage colours
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

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

############################## Read in SEACells assignments and add to ArchR object #######################################
data_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/SEACells/4_exported_SEACells_data/rds_files/"
SEACells_metadata <- read_csv(paste0(data_path, "AnnData_metacells_assigned_cell_metadata.csv"))
SEACells_cell_assignments <- SEACells_metadata %>% select(index, SEACell)
dim(SEACells_cell_assignments)

dim(getCellColData(ArchR))
ArchR <- addCellColData(ArchRProj = ArchR, data = SEACells_cell_assignments$SEACell,
                            cells = SEACells_cell_assignments$index, name = "SEACell")
getCellColData(ArchR)

###### plot seacells
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

plotEmbedding(ArchR, 
              name = "SEACell",
              plotAs = "points", size = ifelse(length(unique(ArchR$stage)) == 1, 1.8, 1),
              baseSize = 0, labelSize = 0, legendSize = 0,
              randomize = TRUE)

### add some plots here to show purity of seacells

## save ArchR object with metacell assignments