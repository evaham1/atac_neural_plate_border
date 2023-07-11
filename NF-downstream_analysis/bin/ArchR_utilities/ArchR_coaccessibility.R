#!/usr/bin/env Rscript

print("coaccessibility ArchR")

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
    data_path = "./output/NF-downstream_analysis/Processing/FullData/Clustering/rds_files/"
    # smaller object
    data_path = "./output/NF-downstream_analysis/Processing/ss8/Peak_call//rds_files/"
    
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

############################## Read in ArchR project #######################################

# If files are not in rds_files subdirectory look in input dir 
label <- unique(sub('_.*', '', list.files(data_path)))
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

# see what is in the ArchR object already
print("ArchR object info: ")
print(ArchR)
getPeakSet(ArchR)
getAvailableMatrices(ArchR)


###############################################################################################
############################## CO-ACCESSIBILITY BETWEEN PEAKS #################################

# calculate co-accessibility between all peaks
ArchR <- addCoAccessibility(ArchR)

# extract interactions - returns indexes of queryHits and subjectHits
cA <- getCoAccessibility(ArchR, corCutOff = 0.5, returnLoops = FALSE)
cA
  # DataFrame with 120270 rows and 11 columns
  # queryHits subjectHits seqnames correlation Variability1 Variability2     TStat        Pval         FDR VarQuantile1 VarQuantile2
  # <integer>   <integer>    <Rle>   <numeric>    <numeric>    <numeric> <numeric>   <numeric>   <numeric>    <numeric>    <numeric>
  #   1              3           4     chr1    0.548725   0.00437754   0.00683964   14.5441 4.15759e-40 4.52151e-38     0.911185     0.965430
  # 2              4           3     chr1    0.548725   0.00683964   0.00437754   14.5441 4.15759e-40 4.52151e-38     0.965430     0.911185
  # 3              4           5     chr1    0.517190   0.00683964   0.00356568   13.3901 4.49249e-35 3.64027e-33     0.965430     0.870967
  # 4              5           4     chr1    0.517190   0.00356568   0.00683964   13.3901 4.49249e-35 3.64027e-33     0.870967     0.965430
  # 5             27          40     chr1    0.761607   0.01690577   0.00855042   26.0418 1.47916e-94 2.12498e-91     0.995825     0.978303
coacessibility_df <- as.data.frame(cA)

# Need to use indices from df to extract granges and therefore informative peak IDs
coacessibility_df <- coacessibility_df %>% 
  mutate(query_PeakID = paste0(seqnames(metadata(cA)[[1]][queryHits]), "-", start(metadata(cA)[[1]][queryHits]), "-", end(metadata(cA)[[1]][queryHits]))) %>%
  mutate(subject_PeakID = paste0(seqnames(metadata(cA)[[1]][subjectHits]), "-", start(metadata(cA)[[1]][subjectHits]), "-", end(metadata(cA)[[1]][subjectHits])))

head(coacessibility_df)

# sanity check that all interaction Peak IDs are in the ArchR peakset
table(coacessibility_df$subject_PeakID %in% getPeakSet(ArchR)$name)

# save df
write.csv(coacessibility_df, file = paste0(rds_path, "Peak_coaccessibility_df.csv"))



#### Browser tracks
p <- plotBrowserTrack(
  ArchRProj = ArchR,
  groupBy = "clusters", 
  geneSymbol = "SIX1", 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(ArchR)
)
grid::grid.newpage()
grid::grid.draw(p[[1]])



#########################################################################################################
############################## CO-ACCESSIBILITY BETWEEN PEAKS AND GENES #################################

# calculate gene-to-peak co-accessibility using GeneIntegrationMatrix
ArchR <- addPeak2GeneLinks(ArchR)

# extract resulting interactions
p2g <- getPeak2GeneLinks(ArchR, corCutOff = 0.5, returnLoops = FALSE)
p2g_df <- as.data.frame(p2g)
head(p2g)

# need to add correct Peak IDs and gene names to df
p2g_df <- p2g_df %>% 
  mutate(PeakID = paste0(seqnames(metadata(p2g)$peakSet[idxATAC]), "-", start(metadata(p2g)$peakSet[idxATAC]), "-", end(metadata(p2g)$peakSet[idxATAC]))) %>%
  mutate(gene_name = metadata(p2g)$geneSet[idxRNA]$name)

head(p2g_df)

# sanity check that all interaction Peak IDs are in the ArchR peakset
table(p2g_df$PeakID %in% getPeakSet(ArchR)$name)

# save df
write.csv(p2g_df, file = paste0(rds_path, "Peak_to_gene_linkage_df.csv"))

## Heatmap of linkage across clusters
p <- plotPeak2GeneHeatmap(ArchRProj = ArchR, groupBy = "clusters")
png(paste0(plot_path, 'Peak_to_gene_linkage_clusters_heatmap.png'), height = 80, width = 60, units = 'cm', res = 400)
print(p)
graphics.off()

p <- plotPeak2GeneHeatmap(ArchRProj = ArchR, groupBy = "stage")
png(paste0(plot_path, 'Peak_to_gene_linkage_stage_heatmap.png'), height = 80, width = 60, units = 'cm', res = 400)
print(p)
graphics.off()


###########################################################################################
############################## BROWSER TRACKS P2G LINKAGE #################################

# extract resulting interactions
head(p2g)
head(p2g_df)

# extract interactions from gene of interest using df
gene = "SIX1"
interactions <- p2g_df %>% filter(gene_name %in% gene)
print(interactions)

# extract the peak IDs that interact with that gene
interacting_peaks <- unique(interactions$PeakID)

# make own loops GRanges for plotting these filtered interactions
interactions

# extract gene TSS coordinate from metadata of p2g linkage
gene_metadata <- metadata(p2g)$geneSet
gene_coord <- start(ranges(gene_metadata[which(gene_metadata$name %in% gene)]))
chr <- as.character(seqnames(gene_metadata[which(gene_metadata$name %in% gene)]))

# extract interacting peak coordinates from unique peak IDs
split <- unlist(lapply(interacting_peaks, FUN = function(x) strsplit(x, split = "-")))
npeaks <- length(split) / 3
peak_coord <- as.numeric(split[1:npeaks*3]) - 250

# put this together into a df
df <- data.frame(start=gene_coord, end=peak_coord)
df_ordered <- as.data.frame(t(as.data.frame(apply(df, 1, function(x) x[order(as.vector(x))]))))
colnames(df_ordered) <- c("start", "end")
loops_of_interest_df <- df_ordered %>% mutate(chr = chr) %>% 
  dplyr::select(chr, start, end)
loops_of_interest_df

# filter all interactions from P2GL for these ones
all_interactions <- getPeak2GeneLinks(ArchR)[[1]]
filtered_granges <- subset(all_interactions, start %in% loops_of_interest_df$start & end %in% loops_of_interest_df$end)
length(filtered_granges) == nrow(loops_of_interest_df)

# identify the range you need to plot to see all these interactions
max_coordinate <- max(end(filtered_granges))
min_coordinate <- min(start(filtered_granges))

# identify how far away from the TSS these coordinates are
distance <- max(abs(gene_coord - max_coordinate), abs(gene_coord - min_coordinate)) + 500
print(distance)

# plot all the interactions pertaining to one gene around that gene
p <- plotBrowserTrack(
  ArchRProj = ArchR,
  groupBy = "clusters", 
  geneSymbol = gene, 
  upstream = distance,
  downstream = distance,
  loops = filtered_granges
)
grid::grid.newpage()
grid::grid.draw(p[[1]])