##### Input files: 1) summarised data .csv from SEACells computation (table of summarised peak counts across metacells)
##### and any 2) feature_metadata.csv from SEACells computation (table with all PeakSet information)
##### Normlises peak counts across seacells, then filters based on annotation and variability

# load libraries 
library(getopt)
library(optparse)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(parallel)
library(data.table)

############################## Set up script options #######################################
# Read in command line opts
option_list <- list(
  make_option(c("-r", "--runtype"), action = "store", type = "character", help = "Specify whether running through through 'nextflow' in order to switch paths"),
  make_option(c("-c", "--cores"), action = "store", type = "integer", help = "Number of CPUs"),
  make_option(c("-n", "--nPeaks"), action = "store", type = "integer", help = "Number of top variable peaks to include", default = 100),
  make_option(c("", "--verbose"), action = "store", type = "logical", help = "Verbose", default = TRUE)
)

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(opt$verbose) print(opt)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    ncores = 8
    
    # combined SEACell outputs
    data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/0_combining_outputs/csv_files/"
    plot_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/plots/"
    rds_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/csv_files/"
    
    ncores = opt$cores
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n")) 
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Functions #######################################

# function to calculate variance and pick out top n features
calc_top_variable_features <- function(counts, n = 1000){
  # calculate variance
  variance <- apply(counts, 2, var)
  # rank and pick top features
  features <- names(variance[rank(-variance) < n+1])
  
  return(features)
}

############################## Pass through SEACell metadata #######################################

## Read in SEACell metadata
metadata <- fread(paste0(data_path, "Combined_SEACell_integrated_metadata.csv"), header = TRUE)
metadata <- metadata[,-1]
print("Preview of SEACell metadata df:")
print(head(metadata))

# Write out unaltered SEACell metadata
write.csv(metadata, paste0(rds_path, "Combined_SEACell_integrated_metadata.csv"), col.names = TRUE)

############################## Read in PeakSet data + pass through #######################################

## Read in peak metadata
peak_metadata <- fread(paste0(data_path, "Feature_metadata.csv"), header = TRUE)
peak_metadata <- peak_metadata[,-c(1:2)]
print("Preview of peakset df:")
print(head(peak_metadata))

# Write out unaltered peak metadata
write.csv(peak_metadata, paste0(rds_path, "Feature_metadata.csv"), col.names = TRUE)

############################## Read in SEACells summarised data #######################################

# Read in summarised peak counts across all seacells
SEACells_summarised <- fread(paste0(data_path, "Combined_summarised_by_metacells_counts.csv"), header = TRUE)
SEACells_summarised <- SEACells_summarised[, -1]
print("Data read in!")

# Check input
print(dim(SEACells_summarised))
print("Preview of input df:")
print(SEACells_summarised[1:4, 1:4])

# Extract SEACell IDs from first column and print out as a text file
SEACells_IDs <- dplyr::pull(SEACells_summarised, 1)
print(head(SEACells_IDs))
length(SEACells_IDs)
write(SEACells_IDs, paste0(rds_path, "SEACell_IDs.txt"))

# Clean up df
SEACells_summarised <- SEACells_summarised[,-1]
print("Preview of input df after cleanup:")
print(SEACells_summarised[1:4, 1:4])
dim(SEACells_summarised)

# Turn into numeric matrix for downstream processing
SEACells_summarised_numeric <- as.matrix(sapply(SEACells_summarised, as.numeric))  

# Adjust peak IDs name format
colnames(SEACells_summarised_numeric) <- gsub(':', '-', colnames(SEACells_summarised_numeric))

# Add SEACell IDs as rownames
rownames(SEACells_summarised_numeric) <- SEACells_IDs

# Check resulting matrix
print(dim(SEACells_summarised_numeric))
print("Preview of summarised count df:")
print(SEACells_summarised_numeric[1:4, 1:4])

# Write out unfiltered raw counts matrix
fwrite(SEACells_summarised_numeric, paste0(rds_path, "Unfiltered_raw_summarised_counts.csv"), sep = ",")

############################## 1) Normalise #######################################
## normalise each metacell by the total number of cut sites

print("Normalising counts...")

# Calculate total number of reads per cell
total_counts <- rowSums(SEACells_summarised_numeric)

# Convert to counts per 10,000 (like normalising done in Seurat) 
normalised_counts <- SEACells_summarised_numeric / (total_counts / 1e4)

print("Preview of normalised counts df:")
print(normalised_counts[1:2, 1:2])
print(normalised_counts[2145:2156, 1:2])
dim(normalised_counts)

## calculate new variance and plot 
variance_before <- apply(SEACells_summarised_numeric, 2, var)
print("Before normalising:")
print(summary(variance_before))

png(paste0(plot_path, "hist_variance_before_normalising.png"), width=60, height=40, units = 'cm', res = 200)
hist(variance_before, breaks = 1000, xlim = c(0, 100))
graphics.off()

variance_after <- apply(normalised_counts, 2, var)
print("After normalising:")
print(summary(variance_after))

png(paste0(plot_path, "hist_variance_after_normalising.png"), width=60, height=40, units = 'cm', res = 200)
hist(variance_after, breaks = 1000, xlim = c(0, 0.0005))
graphics.off()

print("Counts normalised!")

# Write out unfiltered normalised counts matrix
fwrite(normalised_counts, paste0(rds_path, "Unfiltered_normalised_summarised_counts.csv"), sep = ",")

################ 2) Filter peaks by annotation #########################
## only include intergenic and intronic peaks

print("Filtering by peak annotation...")

included_peak_set <- peak_metadata[which(peak_metadata$peakType %in% c("Distal", "Intronic")), ]
included_peaks <- included_peak_set$name

# filter normalised summarised counts to only include these peaks
annot_filtered_matrix <- normalised_counts[, which(colnames(normalised_counts) %in% included_peaks)]
dim(annot_filtered_matrix)

# plot how many peaks filtered
peak_counts <- data.frame(unfiltered = dim(normalised_counts)[2], filtered = dim(annot_filtered_matrix)[2])

png(paste0(plot_path, "annotation_filtering.png"), width=60, height=40, units = 'cm', res = 200)
grid.arrange(top=textGrob("Remaining Peak Counts", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
               tableGrob(peak_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

############## 3) Filter peaks by variance #######################################
## pick peaks with most variance across all cells

print("Filtering by variance...")

## pick peaks with top 10,000 variance (as test do top 100) determined by n param
top_features <- calc_top_variable_features(annot_filtered_matrix, n = opt$nPeaks)
print(paste0("Number of top variable features: ", length(top_features)))

## filter normalised matrix based on these features
variable_filtered_matrix <- annot_filtered_matrix[, top_features]

# plot how many peaks filtered
peak_counts <- data.frame(unfiltered = dim(annot_filtered_matrix)[2], filtered = dim(variable_filtered_matrix)[2])

png(paste0(plot_path, "variable_filtering.png"), width=60, height=40, units = 'cm', res = 200)
grid.arrange(top=textGrob("Remaining Peak Counts", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
               tableGrob(peak_counts, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# check new matrix
print("Preview of filtered normalised counts df:")
print(variable_filtered_matrix[1:2, 1:2])
dim(variable_filtered_matrix)

# Check new variance
variance <- apply(variable_filtered_matrix, 2, var)
print(summary(variance))

png(paste0(plot_path, "hist_variance_after_normalising_after_filtering.png"), width=60, height=40, units = 'cm', res = 200)
hist(variance, breaks = 1000)
graphics.off()

print("Peaks filtered!")

############################## Plot how many peaks filtered in each step #######################################

df <- data.frame(Filtering_step = c("Unfiltered", "Annotation_filter", "Variance_filter"),
                 Remaining_peak_count = c(dim(normalised_counts)[2], dim(annot_filtered_matrix)[2], dim(variable_filtered_matrix)[2])
)
png(paste0(plot_path, 'how_many_peaks_filtered_table.png'), height = 5, width = 12, units = 'cm', res = 400)
grid.arrange(top=textGrob(" ", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, theme = ttheme_minimal()))
graphics.off()

############################## Save filtered normalised summarised count data #######################################

write.csv(variable_filtered_matrix, paste0(rds_path, "Filtered_normalised_summarised_counts.csv"), row.names = TRUE, col.names = TRUE)

############################## Filter raw counts matrix and save #######################################

# Use same peaks to filter raw data matrix
filtered_raw_matrix <- SEACells_summarised_numeric[, top_features]

# check new matrix
print("Preview of filtered raw counts df:")
print(filtered_raw_matrix[1:2, 1:2])
dim(filtered_raw_matrix)

# write to file
write.csv(filtered_raw_matrix, paste0(rds_path, "Filtered_raw_summarised_counts.csv"), row.names = TRUE, col.names = TRUE)