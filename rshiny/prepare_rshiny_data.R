# Prepare skinny data for shiny
library(tidyverse)
library(data.table)

# read in normalised SEACells data with all peaks and all SEACells (output from 1_peak_filtering)
data_path = "./output/NF-downstream_analysis/Downstream_processing/Cluster_peaks/1_peak_filtering/rds_files/"
SEACells_normalised_summarised <- fread(paste0(data_path, "Unfiltered_normalised_summarised_counts.csv"), header = TRUE)
print("Normalised data read in!")

# Extract SEACell IDs from first column
SEACells_IDs <- SEACells_normalised_summarised$V1
length(SEACells_IDs)

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

# Make tiny data
tiny_df <- SEACells_normalised_summarised[1:100, ]

# write to file
write.csv(tiny_df, "output/Rshiny_input.csv", row.names = TRUE, col.names = TRUE)