#!/usr/bin/env Rscript

### Script to predict gene activity from ATAC data
print("3_gene_activity")

############################## Load libraries #######################################
library(getopt)
library(Signac)
library(Seurat)
library(future)
library(tidyverse)
library(grid)
library(gridExtra)
library(clustree)
library(ggplot2)
library(dplyr)
library(scHelper)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
library(GenomicRanges)

############################## Set up script options #######################################
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    setwd("~/NF-downstream_analysis")
    ncores = 8

    plot_path = "../output/NF-downstream_analysis/gene_activity/plots/"
    rds_path = "../output/NF-downstream_analysis/gene_activity/rds_files/"
    data_path = "../output/NF-downstream_analysis/test_input/"   
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 305* 1024^3)
    plan()
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

############################## Read in Seurat RDS object and fragment files #######################################

seurat <- readRDS(paste0(data_path, "rds_files/seurat_all_filtered.RDS"))
print(seurat)

# read in fragment files
paths <- list.dirs(paste0(data_path, "cellranger_atac_output/"), recursive = FALSE, full.names = TRUE)
input <- data.frame(sample = sub('.*/', '', paths), 
                    matrix_path = paste0(paths, "/outs/filtered_peak_bc_matrix.h5"),
                    metadata_path = paste0(paths, "/outs/singlecell.csv"),
                    fragments_path = paste0(paths, "/outs/fragments.tsv.gz"))
new.paths <- as.list(input$fragments_path)
frags <- Fragments(seurat)  # get list of fragment objects
Fragments(seurat) <- NULL  # remove fragment information from assay

for (i in seq_along(frags)) {
  frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
}
Fragments(seurat) <- frags # assign updated list back to the object

print("data read in")

######################################## Annotations information ###############################################
annotations_plot_path = paste0(plot_path, "annotations/")
dir.create(annotations_plot_path, recursive = T)

annotation <- Annotation(seurat[["peaks"]])
colnames(mcols(annotation))
#[1] "source"             "type"               "score"              "phase"             
#[5] "gene_id"            "gene_version"       "gene_source"        "gene_biotype"      
#[9] "transcript_id"      "transcript_version" "transcript_source"  "transcript_biotype"
#[13] "exon_number"        "exon_id"            "exon_version"       "protein_id"        
#[17] "protein_version"    "gene_name"          "transcript_name"    "tag"
dim(mcols(annotation)) # 24356, 20
unique(mcols(annotation)$gene_biotype)
#  [1] "protein_coding"       "snRNA"                "lncRNA"               "miRNA"               
# [5] "pseudogene"           "misc_RNA"             "snoRNA"               "rRNA"                
# [9] "scaRNA"               "ribozyme"             "processed_pseudogene" "IG_V_gene"           
# [13] "sRNA"                 "Mt_tRNA"              "Mt_rRNA"      
sum(is.na(mcols(annotation)$gene_id)) #0
sum(is.na(mcols(annotation)$gene_name)) #10501
sum(mcols(annotation)$gene_biotype == "protein_coding") # 16779

datalist <- c()
for (i in unique(mcols(annotation)$gene_biotype)){
  sum <- sum(mcols(annotation)$gene_biotype == i)
  datalist[[i]] <- sum # add it to your list
}
pie_chart = do.call(rbind, datalist)

png(paste0(annotations_plot_path, "Gene_biotypes.png"), width=30, height=20, units = 'cm', res = 200)
pie(pie_chart, labels = pie_chart, main = "Gene biotype annotations from GTF",
    col = rainbow(length(pie_chart)), radius = 1)
legend("bottomleft", rownames(pie_chart), cex = 0.8,
       fill = rainbow(length(pie_chart)))
graphics.off()

protein_coding_annotations <- mcols(annotation)[grep("protein_coding", mcols(annotation)$gene_biotype), ]
dim(protein_coding_annotations) # 16779 20
sum(is.na(protein_coding_annotations$gene_id)) #0
gene_name_NA <- sum(is.na(protein_coding_annotations$gene_name)) #4303
gene_name <- length(protein_coding_annotations[,1]) - sum(is.na(protein_coding_annotations$gene_name))

png(paste0(annotations_plot_path, "protein_coding_gene_names.png"), width=30, height=20, units = 'cm', res = 200)
pie(c(gene_name_NA, gene_name), labels = c(gene_name_NA, gene_name), 
    main = "Proportion of genes which have a gene name",
    col = rainbow(2), radius = 1)
legend("bottomleft", c("NAs", "named genes"), cex = 0.8,
       fill = rainbow(2))
graphics.off()

print("annotations investigated")

######################################## ESTIMATE GEX #####################################################

gex_plot_path = paste0(plot_path, "gene_activity/")
dir.create(gex_plot_path, recursive = T)

gene.activities.names <- GeneActivity(seurat)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
seurat[['RNA']] <- CreateAssayObject(counts = gene.activities.names)
seurat <- NormalizeData(
  object = seurat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat$nCount_RNA)
)

DefaultAssay(seurat) <- 'RNA'

print("gex predicted")

######################################## Gex information ###############################################

length(rownames(gene.activities.names)) # 16731
gene_name_NA <- sum(grepl("NA.", rownames(gene.activities.names), fixed = TRUE)) # 4274
rownames(gene.activities.names)[grep("NA.", rownames(gene.activities.names), fixed = TRUE)]
gene_name <- length(rownames(gene.activities.names)) - gene_name_NA

png(paste0(gex_plot_path, "protein_coding_gene_names.png"), width=30, height=20, units = 'cm', res = 200)
pie(c(gene_name_NA, gene_name), labels = c(gene_name_NA, gene_name), 
    main = "Proportion of genes which have a gene name",
    col = rainbow(2), radius = 1)
legend("bottomleft", c("NAs", "named genes"), cex = 0.8,
       fill = rainbow(2))
graphics.off()

df <- data.frame(total_counts = rowSums(gene.activities.names))
png(paste0(gex_plot_path, "gene_counts_boxplot.png"), width=30, height=20, units = 'cm', res = 200)
ggplot(df, aes(x = total_counts)) + geom_boxplot()
graphics.off()

png(paste0(gex_plot_path, "gene_counts_hist_xmax20.png"), width=30, height=20, units = 'cm', res = 200)
ggplot(df, aes(x = total_counts)) + geom_histogram(binwidth = 1) + xlim(-1, 20)
graphics.off()

print("gex plots made")

######################################## Top variable genes ###############################################

seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 10000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png(paste0(plot_path, 'MostVariable.png'), height = 15, width = 20, units = 'cm', res = 400)
plot2
graphics.off()

print("variable genes found")

######################################## Scaling ###############################################

# save seurat object with predicted RNA assay
saveRDS(seurat, paste0(rds_path, "seurat_GeneActivity.RDS"), compress = FALSE)