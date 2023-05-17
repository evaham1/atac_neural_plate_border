# Prepare skinny data for shiny
library(tidyverse)
library(data.table)

# working dir as base of repository
data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/"

#####  Normalised SEACells peak count matrix

# read in SEACells data
SEACells_normalised_summarised <- fread(paste0(data_path, "Filtered_normalised_summarised_counts.csv"), header = TRUE)

# Extract SEACell IDs from first column
SEACells_IDs <- SEACells_normalised_summarised$V1
length(SEACells_IDs)

# Check for duplicates
table(duplicated(SEACells_IDs))

#Clean up df
SEACells_normalised_summarised <- SEACells_normalised_summarised[,-1]

# Add SEACell IDs as rownames
rownames(SEACells_normalised_summarised) <- SEACells_IDs

# Check resulting matrix
print(dim(SEACells_normalised_summarised))
print("Preview of summarised count df:")
print(SEACells_normalised_summarised[1:4, 1:4])

# Make tiny data with all SEACells and only 100 peaks
tiny_df <- SEACells_normalised_summarised[, 1:100]
rownames(tiny_df) <- rownames(SEACells_normalised_summarised)
print(tiny_df[1:4, 1:4])

# write to file
write.csv(tiny_df, "./output/Rshiny_input_SEACells_matrix.csv", row.names = TRUE, col.names = TRUE)

#####  SEACells metadata
metadata <- read.csv(paste0(data_path, "Combined_SEACell_integrated_metadata.csv"), row.names = 'ATAC')
substrRight <- function(x, n){
  sapply(x, function(xx)
    substr(xx, (nchar(xx)-n+1), nchar(xx))
  )
}
metadata <- metadata %>% mutate(stage = substrRight(rownames(metadata), 3))
metadata <- metadata[,-1]

# Check metadata
print(head(metadata))

# write to file
write.csv(tiny_df, "./output/Rshiny_input_SEACells_metadata.csv", row.names = TRUE, col.names = TRUE)
