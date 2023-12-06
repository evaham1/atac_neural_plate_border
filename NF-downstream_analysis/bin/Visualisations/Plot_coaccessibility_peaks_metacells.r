#!/usr/bin/env Rscript

print("plot coaccessibility ArchR on SEACells matrix as heatmap")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(parallel)
library(scHelper)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-g", "--group_by"), action = "store", type = "character", help = "How to group cells", default = "clusters",),
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
    
    # ArchR object with no coaccessibility data
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Metacell_ID_purity/rds_files/"
    # co-accessibility csv and RDS files
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Single_cell_integration/"
    # # ArchR object with coaccessibility data
    # data_path = "./output/NF-downstream_analysis/Processing/FullData/Single_cell_integration/"
    
    # rds_path = "./output/NF-downstream_analysis/Processing/ss8/Coaccessibility/rds_files/"
    # plot_path = "./output/NF-downstream_analysis/Processing/ss8/Coaccessibility/plots/"
    
    addArchRThreads(threads = 1) 
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    txt_path = "./text_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    #addArchRThreads(threads = ncores)
    addArchRThreads(threads = 1) 
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  dir.create(txt_path, recursive = T)
  
  
}

set.seed(42)

########################       CELL STATE COLOURS    ########################################
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak', 
                              'PGC', 'BI', 'meso', 'endo', 'MIXED', 'Unmapped',
                              'Neural', 'Placodal', 'Non-neural', 'Contam')
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



############################## Read in data #######################################

print("reading in interactions data...")

# read in P2G dataframe
p2g_df <- read.csv(paste0("output/NF-downstream_analysis/Processing/FullData/Single_cell_integration/csv_files/Peak_to_gene_linkage_df_250000_distance.csv"))
head(p2g_df)

# read in SEACells matrix
SEACells_peak_matrix <- fread('./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/Unfiltered_normalised_summarised_counts.csv', sep = ",")
SEACells_IDs <- scan("./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/SEACell_IDs.txt", character(), quote = "")
table(duplicated(SEACells_IDs))
rownames(SEACells_peak_matrix) <- SEACells_IDs
SEACells_peak_matrix[1:2, 1:2]

SEACells_peak_matrix <- as.matrix(sapply(SEACells_peak_matrix, as.numeric))
rownames(SEACells_peak_matrix) <- SEACells_IDs
SEACells_peak_matrix[1:2, 1:2]
dim(SEACells_peak_matrix)

# Read in metadata for all SEACells
SEACells_metadata <- fread('./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/')
SEACells_metadata <- SEACells_metadata %>% mutate(stage = sub(".*-", "", ATAC))
head(SEACells_metadata)


############################## Plot heatmap of peaks that interact with placodal gene modules #######################################

# PPR specific GMs
# GM12; AKR1D1, ATP1A1, ATP1B1, ATP2B1, B3GNT7, CCDC25, CGNL1, DAG1, EMILIN2, ENSGALG00000004814, EPCAM, FAM184B, FAM89A, GATA2, GATA3, IRAK2, IRF6, MAP7, METRNL, MPZL3, NET1, PLEKHA5, POGLUT2, PPFIBP1, RGN, SLC16A10, SLC25A4, TSPAN13, TTC39A, UNC5B, Z-TJP2
# GM14; BRINP1, CITED4, DLX5, DLX6, ENSGALG00000023936, ENSGALG00000040010, FN1, HESX1, HIF1A, KCNAB1, LAMB1, LRP11, NFKB1, PAX6, PITX1, PITX2, SEM1, SFRP1, SH3D19, SHISA2, SIX3, SPON1, SST, WDR1, Z-HAPLN1
# GM13; AKAP12, ASS1, BASP1, CD99, ENSGALG00000011296, ENSGALG00000041054, ENSGALG00000042443, EYA2, FERMT2, LGMN, METTL24, NR2F2, NUCKS1, SIX1

PPR_genes <- c("AKR1D1", "ATP1A1", "ATP1B1", "ATP2B1", "B3GNT7", "CCDC25", "CGNL1", "DAG1", "EMILIN2", "ENSGALG00000004814", "EPCAM", "FAM184B", "FAM89A", "GATA2", "GATA3", "IRAK2", "IRF6", "MAP7", "METRNL", "MPZL3", "NET1", "PLEKHA5", "POGLUT2", "PPFIBP1", "RGN", "SLC16A10", "SLC25A4", "TSPAN13", "TTC39A", "UNC5B", "Z-TJP2",
               "BRINP1", "CITED4", "DLX5", "DLX6", "ENSGALG00000023936", "ENSGALG00000040010", "FN1", "HESX1", "HIF1A", "KCNAB1", "LAMB1", "LRP11", "NFKB1", "PAX6", "PITX1", "PITX2", "SEM1", "SFRP1", "SH3D19", "SHISA2", "SIX3", "SPON1", "SST", "WDR1", "Z-HAPLN1",
               "AKAP12", "ASS1", "BASP1", "CD99", "ENSGALG00000011296", "ENSGALG00000041054", "ENSGALG00000042443", "EYA2", "FERMT2", "LGMN", "METTL24", "NR2F2", "NUCKS1", "SIX1")





############################## Plot heatmap of peaks that interact with NC gene modules #######################################

# # NC specific GMs
# # GM40; AGTRAP, BRINP2, CDH11, CMTM8, ENSGALG00000001136, FRZB, GADD45A, HUNK, LARP7, LMX1B, MRAS, MSX1, SFRP2, SOX11, SPSB4, Z-FST, ZEB2, ZIC1
# # GM42; BMP5, CDH6, CSRNP1, DRAXIN, EN1, ENSGALG00000013505, FOXD3, NKD1, NRP2, OLFM1, PAX7, SNAI2, SOX9, TFAP2B, TMEM132C, TSPAN18, WLS, WNT6, Z-ENC1, ZFHX4, ZNF423
# # GM43; COL9A3, ENSGALG00000037717, ENSGALG00000053185, ERMN, ETS1, LMO4, PPP1R1C, RASL11B, RFTN2, SOX10, SOX8, TNC, Z-MEF2C, Z-PLK2
# # GM44; CXCR4, ENSGALG00000030512, ENSGALG00000031427, FABP7, FKBP11, GLIPR2, ID1, ID2, MEOX1, MYL4, OLFML3, SOX5, WNT1
# 
NC_genes <- c("AGTRAP", "BRINP2", "CDH11", "CMTM8", "ENSGALG00000001136", "FRZB", "GADD45A", "HUNK", "LARP7", "LMX1B", "MRAS", "MSX1", "SFRP2", "SOX11", "SPSB4", "Z-FST", "ZEB2", "ZIC1",
              "BMP5", "CDH6", "CSRNP1", "DRAXIN", "EN1", "ENSGALG00000013505", "FOXD3", "NKD1", "NRP2", "OLFM1", "PAX7", "SNAI2", "SOX9", "TFAP2B", "TMEM132C", "TSPAN18", "WLS", "WNT6", "Z-ENC1", "ZFHX4", "ZNF423",
              "CDON", "COTL1", "ENSGALG00000048488", "MFAP2", "PRTG", "TUBB3",
              "COL9A3", "ENSGALG00000037717", "ENSGALG00000053185", "ERMN", "ETS1", "LMO4", "PPP1R1C", "RASL11B", "RFTN2", "SOX10", "SOX8", "TNC", "Z-MEF2C", "Z-PLK2",
              "CXCR4", "ENSGALG00000030512", "ENSGALG00000031427", "FABP7", "FKBP11", "GLIPR2", "ID1", "ID2", "MEOX1", "MYL4", "OLFML3", "SOX5", "WNT1")


