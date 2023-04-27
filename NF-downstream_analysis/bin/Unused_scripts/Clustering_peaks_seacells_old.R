##### Filtering peaks before clustering them

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


############################## Read in SEACells summarised data #######################################
data_path = "./output/NF-downstream_analysis/Downstream_processing/Peak_clustering/SEACells/4_exported_SEACells_data/rds_files/"
SEACells_summarised <- read_csv(paste0(data_path, "AnnData_summarised_by_metacells_peak_counts.csv"))
dim(SEACells_summarised)

SEACells_summarised <- as.matrix(fread(paste0(data_path, "summarised_counts_1000.csv"), header = TRUE), rownames = 1)

hist(colSums(SEACells_summarised), breaks = 1000)
summary(colSums(SEACells_summarised))
dim(SEACells_summarised)

### 1) Normalising
### 2) Filtering - variance + annotations
### 3) Calculate correlation matrix
### 4) Filter features based on correlation (corr threshold and min number of higihly correlated features)
### 5) Identify how many peak modules explain the overall variance
### 6) Cluster corr matrix into that number of modules

## normalise each metacell by the total number of cut sites
normalised_counts <- t(apply(SEACells_summarised, 1, function(x) x/sum(x))) * 1000
normalised_counts[1:2, 1:2]
dim(normalised_counts)


## could add a minimum number of total cut sites here


# variance <- apply(SEACells_summarised, 2, var)
# hist(variance, breaks = 1000)
# summary(variance)
# 


## pick peaks with top 10,000 variance
# as test do top 100

# function to calculate variance and pick out top n features
calc_top_variable_features <- function(counts, n = 1000){
  # calculate variance
  variance <- apply(counts, 2, var)
  # rank and pick top features
  features <- names(variance[rank(-variance) < n+1])
  
  return(features)
}

# variance <- apply(SEACells_summarised, 2, var)
# hist(variance, breaks = 1000)

top_features <- calc_top_variable_features(normalised_counts, n = 100)
length(top_features)
filtered_normalised_counts <- normalised_counts[, top_features]
dim(filtered_normalised_counts)
# variance <- apply(filtered_counts, 2, var)
# hist(variance, breaks = 1000)
# summary(variance)



### correlation matrix
corr_mat <- cor(filtered_normalised_counts, method = "spearman")
dim(corr_mat)
#hist(correlation)

# function to filter features based on correlation threshold
calc_top_variable_features <- function(counts, corr_t = 0.7, corr_min = 20){
  # calculate correlation
  corr_matrix <- cor(counts)
  # filter
  filtered_features <- colSums(corr_matrix > corr_t) > corr_min
  filtered_counts <- corr_matrix[filtered_features, filtered_features]
  
  return(filtered_counts)
}

filtered_corr_mat <- calc_top_variable_features(corr_mat, corr_t = 0.5, corr_min = 8)
dim(filtered_corr_mat)
#dim(filtered_counts)

filtered_normalised_counts <- filtered_normalised_counts[, colnames(filtered_corr_mat)]

dim(filtered_corr_mat)
dim(filtered_normalised_counts)

# cluster the correlation matrix
# to do this we want to iteratively cluster and then calculate the variance
# we can then find the 'elbow' in which adding more clusters doesnt reduce intracluster variance anymore
# this heuristically allows us to determine total number of clusters

diss_matrix <- as.dist(1 - filtered_corr_mat)
tree <- hclust(diss_matrix, method="complete")
plot(tree)

dataset = filtered_normalised_counts
krange = seq(1, 30, 5)
k = krange[2]

hc.mod_ids = cutree(tree, k = k)
distortion = 0
for (ik in 1:k) {
  clust = dataset[which(hc.mod_ids == ik), , 
                  drop = F]
  clust_mean = colMeans(clust)
  distortion <- distortion + sum(apply(clust, 
                                       1, function(r) {
                                         rc = r - clust_mean
                                         return(sum(rc^2))
                                       }))
}


test(filtered_normalised_counts, hc = tree, num_mod_max = 20, numcores = 5, display = TRUE, pdf_plot_path = "test.pdf")

test <- function (dataset, hc, num_mod_max, numcores, display = FALSE, 
          pdf_plot_path = NULL, method = c("min_area", "piecewise")) 
{
  method <- match.arg(method)
  # krange = c(seq(1, 99, 10), seq(100, 999, 100), seq(1000, 
  #                                                    1e+05, 1000))
  # krange <- krange[which(krange < nrow(dataset))]
  # krange <- krange[which(krange <= num_mod_max)]
  # krange <- c(krange, nrow(dataset)
  krange = seq(1, 30, 5)
  if (numcores == 1) {
    distortion.vals <- unlist(lapply(krange, function(k) {
      getDistortion(k, hc, dataset)
    }))
  }
  else {
    cl <- parallel::makeCluster(min(numcores, length(krange)), 
                                outfile = "/dev/null", methods = FALSE)
    parallel::clusterExport(cl, c("hc", "dataset", "krange"), 
                            envir = environment())
    doSNOW::registerDoSNOW(cl)
    cat("\nRun parallel distortion jobs\n")
    pb <- txtProgressBar(max = length(krange), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    distortion.vals = unlist(foreach(k = krange, .options.snow = opts, 
                                     .noexport = ls(), .packages = c("Matrix")) %dopar% 
                               {
                                 if (any(is.na(dataset))) {
                                   print(paste(k, " contains NA values"))
                                 }
                                 hc.mod_ids = cutree(hc, k = k)
                                 distortion = 0
                                 for (ik in 1:k) {
                                   clust = dataset[which(hc.mod_ids == ik), , 
                                                   drop = F]
                                   clust_mean = colMeans(clust)
                                   distortion <- distortion + sum(apply(clust, 
                                                                        1, function(r) {
                                                                          rc = r - clust_mean
                                                                          return(sum(rc^2))
                                                                        }))
                                 }
                                 print(paste(k, "before"))
                                 print(paste(k, ":", distortion))
                                 return(distortion)
                               })
    parallel::stopCluster(cl)
  }
  distortion.df = data.frame(x = krange, y = distortion.vals)
  num_mod_max_idx <- length(which(distortion.df$x <= num_mod_max))
  if (method == "piecewise") {
    num_modules <- as.integer(elbow_pos_piecewise_heuristic(distortion.df$x[1:num_mod_max_idx], 
                                                            distortion.df$y[1:num_mod_max_idx], plot = display, 
                                                            save.pdf.path = pdf_plot_path)[1])
  }
  else if (method == "min_area") {
    num_modules <- as.integer(elbow_pos_area_heuristic(distortion.df$x[1:num_mod_max_idx], 
                                                       distortion.df$y[1:num_mod_max_idx], plot = display, 
                                                       save.pdf.path = pdf_plot_path))
  }
  return(num_modules)
}

elbow_pos_area_heuristic <- function(
    x,
    y,
    plot=FALSE,
    save.pdf.path=NULL) {
  
  x.interp = seq(min(x), max(x), 1)
  
  # normalize y to match x range
  norm_coeff = (max(x)-1) / max(y)
  y2 <- norm_coeff * y
  y.interp = approx(x, y2, x.interp)$y
  num_points = length(x.interp)
  
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  
  segments_length = c(0, unlist(lapply(seq(num_points-1), function(i){euc.dist(c(x.interp[i], y.interp[i]), c(x.interp[i+1], y.interp[i+1]))})))
  
  cumulative_segments_length = cumsum(segments_length)
  
  num_points2 = 100
  
  x.interp2 = approx(cumulative_segments_length, 
                     seq(length(cumulative_segments_length)),
                     seq(0, max(cumulative_segments_length), length.out=num_points2)
  )$y
  
  y.interp2 = approx(x,y2,x.interp2)$y
  
  # Point M (xm, ym) and line passing through A (xa, ya) and B (xb, yb)
  point_line_distance <- function(xa, ya, xb, yb, xm, ym){
    
    nx = -(yb-ya)
    ny = (xb-xa)
    
    n_norm = sqrt(nx^2+ny^2)
    
    d = abs(nx * (xa-xm) + ny *(ya-ym)) / n_norm
    return(d)
  }
}


