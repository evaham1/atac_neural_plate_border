#!/usr/bin/env Rscript

print("Downstream analysis on GRN")

############################## Load libraries #######################################
library(getopt)
library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
library(Seurat)
library(Signac)
library(scMEGA)
library(SummarizedExperiment)
library(igraph)
library(ggraph)
library(clusterProfiler)
library(org.Gg.eg.db)
library(pheatmap)
library(ComplexHeatmap)
library(VennDiagram)
library(RColorBrewer)

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
    
    data_path = "./output/NF-downstream_analysis/Processing/FullData/scMEGA/MEGA_GRNi/csv_files/"
    
    plot_path = "./output/NF-downstream_analysis/Processing/FullData/scMEGA/MEGA_GRN_vis/plots/"
    csv_path = "./output/NF-downstream_analysis/Processing/FullData/scMEGA/MEGA_GRN_vis/csv_files/"
    
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/rds_files/"
    rds_path = "./rds_files/"
    csv_path = "./csv_files/"
    ncores = opt$cores
    
    addArchRThreads(threads = ncores)
    
  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  dir.create(csv_path, recursive = T)
}

set.seed(42)

############################## EDITED FUNCTIONS #######################################

GetTrajectory_updated <- function (object = NULL, trajectory.name = "Trajectory", assay = NULL, 
                                   slot = "counts", groupEvery = 1, log2Norm = TRUE, scaleTo = 10000, 
                                   smoothWindow = 11) 
{
  if (is.null(assay) | !assay %in% Assays(object)) {
    stop("Please provide an available assay!")
  }
  if (!(trajectory.name %in% colnames(object@meta.data))) {
    stop(glue::glue("Cannot find trajecotry {trajectory.name}!"))
  }
  trajectory <- object@meta.data[trajectory.name]
  trajectory <- trajectory[!is.na(trajectory[, 1]), , drop = FALSE]
  breaks <- seq(0, 100, groupEvery)
  if (!all(is.numeric(trajectory[, 1]))) {
    stop("Trajectory must be a numeric. Did you add the trajectory with addTrajectory?")
  }
  if (!all(trajectory[, 1] >= 0 & trajectory[, 1] <= 100)) {
    stop("Trajectory values must be between 0 and 100. Did you add the trajectory with addTrajectory?")
  }
  groupList <- lapply(seq_along(breaks), function(x) {
    if (x == 1) {
      NULL
    }
    else {
      rownames(trajectory)[which(trajectory[, 1] > breaks[x - 
                                                            1] & trajectory[, 1] <= breaks[x])]
    }
  })[-1]
  names(groupList) <- paste0("T.", breaks[-length(breaks)], 
                             "_", breaks[-1])
  message("Creating Trajectory Group Matrix..")
  
  ### LayerData fails due to seurat package version issue so have hacked this section myself to extract and deal with
  # both ATAC/RNA count data but also motif data (NB motif data must be in a slot called 'chromvar')
  #data.use <- GetAssayData(object, assay = assay, slot = slot)
  #data.use <- LayerData(object, assay = assay)
  if (assay %in% c("chromvar", "target")){
    data.use <- object[[assay]]@data
  } else {
    data.use <- object[[assay]]$counts
  }
  
  # Adapted due to errors when there is only one cell in a group
  groupMat <- lapply(1:length(groupList), function(x) {
    cell_names <- groupList[[x]]
    if (length(cell_names) == 1) {
      mat <- data.use[, cell_names]
    } else {
      mat <- Matrix::rowMeans(data.use[, cell_names]) # rowMeans doesnt work if there is only one cell in group
    }
  }) %>% Reduce(cbind, .)
  
  # Added this to avoid the NAs that arise when there are no cells in a group
  groupMat[is.na(groupMat)] <- 0
  
  colnames(groupMat) <- names(groupList)
  if (!is.null(scaleTo)) {
    if (any(groupMat < 0)) {
      message("Some values are below 0, this could be the Motif activity matrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!")
    }
    else {
      groupMat <- t(t(groupMat)/colSums(groupMat)) * scaleTo
    }
  }
  
  # Added this because NaNs arrive now
  groupMat[is.nan(groupMat)] <- 0
  
  if (log2Norm) {
    if (any(groupMat < 0)) {
      message("Some values are below 0, this could be a Motif activity matrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!")
    }
    else {
      groupMat <- log2(groupMat + 1)
    }
  }
  if (!is.null(smoothWindow)) {
    message("Smoothing...")
    smoothGroupMat <- as.matrix(t(apply(groupMat, 1, function(x) centerRollMean(x, 
                                                                                k = smoothWindow))))
    colnames(smoothGroupMat) <- paste0(colnames(groupMat))
    colnames(groupMat) <- paste0(colnames(groupMat))
    seTrajectory <- SummarizedExperiment(assays = SimpleList(smoothMat = as.matrix(smoothGroupMat), 
                                                             mat = as.matrix(groupMat)))
  }
  else {
    colnames(groupMat) <- paste0(colnames(groupMat))
    seTrajectory <- SummarizedExperiment(assays = SimpleList(mat = as.matrix(groupMat)))
  }
  return(seTrajectory)
}


centerRollMean <- function(v = NULL, k = NULL){
  o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
  if(k%%2==0){
    o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else if(k%%2==1){
    o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else{
    stop("Error!")
  }
  o2
}

PseudotimePlot_updated <- function (object, tf.use, tf.assay = "chromvar", rna.assay = "RNA", 
                                    atac.assay = "ATAC", target.assay = "target", trajectory.name = "Trajectory", 
                                    groupEvery = 1) 
{
  trajMM <- suppressMessages(GetTrajectory_updated(object, assay = tf.assay, 
                                                   trajectory.name = trajectory.name, groupEvery = groupEvery, 
                                                   slot = "data", smoothWindow = 7, log2Norm = FALSE))
  rownames(trajMM) <- object@assays[[atac.assay]]@motifs@motif.names
  df.tf.activity <- assay(trajMM)
  df.tf.activity <- t(scale(t(df.tf.activity)))
  df.tf.activity <- as.data.frame(df.tf.activity)
  colnames(df.tf.activity) <- seq(0, 100, groupEvery)[-1]
  df.tf.activity$tf <- toupper(rownames(df.tf.activity))
  df.tf.activity <- tidyr::pivot_longer(df.tf.activity, -tf, 
                                        names_to = "pseudotime", values_to = "value")
  df.tf.activity$pseudotime <- as.numeric(df.tf.activity$pseudotime)
  trajGEX <- suppressMessages(GetTrajectory_updated(object, assay = rna.assay, 
                                                    trajectory.name = trajectory.name, groupEvery = groupEvery, 
                                                    slot = "data", smoothWindow = 7, log2Norm = TRUE))
  df.tf.expression <- assay(trajGEX)
  df.tf.expression <- t(scale(t(df.tf.expression)))
  df.tf.expression <- as.data.frame(df.tf.expression)
  colnames(df.tf.expression) <- seq(0, 100, groupEvery)[-1]
  df.tf.expression$tf <- toupper(rownames(df.tf.expression))
  df.tf.expression <- tidyr::pivot_longer(df.tf.expression, 
                                          -tf, names_to = "pseudotime", values_to = "value")
  df.tf.expression$pseudotime <- as.numeric(df.tf.expression$pseudotime)
  traj.target <- suppressMessages(GetTrajectory_updated(object, assay = target.assay, 
                                                        trajectory.name = trajectory.name, groupEvery = groupEvery, 
                                                        slot = "data", smoothWindow = 7, log2Norm = FALSE))
  df.target <- assay(traj.target)
  df.target <- t(scale(t(df.target)))
  df.target <- as.data.frame(df.target)
  colnames(df.target) <- seq(0, 100, groupEvery)[-1]
  df.target$tf <- toupper(rownames(df.target))
  df.target <- tidyr::pivot_longer(df.target, -tf, names_to = "pseudotime", 
                                   values_to = "value")
  df.target$pseudotime <- as.numeric(df.target$pseudotime)
  df.tf.activity$data <- "TF activity"
  df.tf.expression$data <- "TF expression"
  df.target$data <- "Targets expression"
  df.tf <- rbind(df.tf.activity, df.tf.expression, df.target)
  df.plot <- subset(df.tf, tf == tf.use)
  p <- ggplot(df.plot, aes(x = pseudotime, y = value, color = data)) + 
    geom_smooth(method = "loess", se = FALSE) + ggtitle(tf.use) + 
    cowplot::theme_cowplot() + ylab("") + theme(legend.title = element_blank())
  return(p)
}

AddModuleScore_updated <- function (object, features, pool = NULL, nbin = 24, ctrl = 100, 
                                    k = FALSE, assay = NULL, name = "Cluster", seed = 1, search = FALSE, 
                                    ...) 
{
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  # GetAssayData doesnt work package issues
  # assay.data <- GetAssayData(object = object)
  assay.data <- object[[assay]]$counts
  
  features.old <- features
  if (k) {
    .NotYetUsed(arg = "k")
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == 
                                         i))
    }
    cluster.length <- length(x = features)
  }
  else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {
        warning("The following features are not present in the object: ", 
                paste(missing.features, collapse = ", "), ifelse(test = search, 
                                                                 yes = ", attempting to find updated synonyms", 
                                                                 no = ", not searching for symbol synonyms"), 
                call. = FALSE, immediate. = TRUE)
        if (search) {
          tryCatch(expr = {
            updated.features <- UpdateSymbolList(symbols = missing.features, 
                                                 ...)
            names(x = updated.features) <- missing.features
            for (miss in names(x = updated.features)) {
              index <- which(x == miss)
              x[index] <- updated.features[miss]
            }
          }, error = function(...) {
            warning("Could not reach HGNC's gene names database", 
                    call. = FALSE, immediate. = TRUE)
          })
          missing.features <- setdiff(x = x, y = rownames(x = object))
          if (length(x = missing.features) > 0) {
            warning("The following features are still not present in the object: ", 
                    paste(missing.features, collapse = ", "), 
                    call. = FALSE, immediate. = TRUE)
          }
        }
      }
      return(intersect(x = x, y = rownames(x = object)))
    })
    cluster.length <- length(x = features)
  }
  if (!all(LengthCheck(values = features))) {
    warning(paste("Could not find enough features in the object from the following feature lists:", 
                  paste(names(x = which(x = !LengthCheck(values = features)))), 
                  "Attempting to match case..."))
    features <- lapply(X = features.old, FUN = CaseMatch, 
                       match = rownames(x = object))
  }
  if (!all(LengthCheck(values = features))) {
    stop(paste("The following feature lists do not have enough features present in the object:", 
               paste(names(x = which(x = !LengthCheck(values = features)))), 
               "exiting..."))
  }
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                         n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                                                                              data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                        ncol = ncol(x = object))
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, 
    ])
  }
  features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
                            ncol = ncol(x = object))
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}

AddTargetAssay_updated <- function (object, target.assay = "target", rna.assay = "RNA", 
                                    df.grn = NULL) 
{
  if (is.na(df.grn)) {
    stop("Cannot find the gene regulatory network!")
  }
  df.genes <- split(df.grn$gene, df.grn$tf)
  object <- AddModuleScore_updated(object, features = df.genes, assay = rna.assay, 
                                   name = "tf_target_", ctrl = 80)
  target_gex <- object@meta.data %>% as.data.frame() %>% dplyr::select(contains("tf_target_"))
  colnames(target_gex) <- names(df.genes)
  object[["target"]] <- CreateAssayObject(data = t(target_gex))
  return(object)
}

LengthCheck <- function(values, cutoff = 0) {
  return(vapply(
    X = values,
    FUN = function(x) {
      return(length(x = x) > cutoff)
    },
    FUN.VALUE = logical(1)
  ))
}

## function extract df showing shared and unique target genes for a list of TFs
extract_target_genes_df <- function(TF_list, grn_df){
  # extract target genes into list called targets_list
  targets_list <- list()
  for (TF in TF_list){
    df.temp <- grn_df %>% 
      filter(tf == TF)
    genes <- unique(df.temp$gene)
    targets_list[[TF]] <- genes
  }
  
  # combine target lists into df called all_targets
  all_targets <- data.frame(
    ID = unique(unlist(targets_list))
  )
  for (i in 1:length(targets_list)){
    print(i)
    targets <- targets_list[[i]]
    df2 <- data.frame(
      ID = targets,
      TF = rep(1, length(targets))
    )
    colnames(df2)[2] <- TF_list[i]
    all_targets <- merge(all_targets, df2, by = "ID", all = TRUE)
  }
  all_targets <- all_targets %>%
    mutate_all(~ ifelse(is.na(.), 0, .))
  all_targets <- column_to_rownames(all_targets, var = "ID")
}


############################## Read in data #######################################

TF_names <- read_tsv("./output/NF-downstream_analysis/Processing/FullData/scMEGA/MEGA_GRNi/Known_TF_names.txt", col_names = FALSE)
TF_names <- TF_names$X1

temp_data_path = paste0(data_path, "./placodal_lineage/")

# read in csvs
df.tfs <- read.csv(paste0(temp_data_path, "TF_correlations_all.csv"))
df.p2g <- read.csv(paste0(temp_data_path, "Variable_genes_with_matched_enhancers.csv"))
df.tf.gene <- read.csv(paste0(temp_data_path, "TF_to_gene_correlations.csv"))

# read in networks
df.grn.un <- read.csv(paste0(temp_data_path, "GRN_unfiltered.csv"))
df.grn <- read_tsv(paste0(temp_data_path, "GRN_filtered.txt"))
df.grn.pos <- read_tsv(paste0(temp_data_path, "GRN_filtered_pos_corr.txt"))

######################################################################################
##############################    FILTERED GRN     ###################################
######################################################################################

temp_plot_path = paste0(plot_path, "filtered_network/")
dir.create(temp_plot_path, recursive = T)

############################## Plot filtered GRN #######################################

# define colors for nodes representing TFs (i.e., regulators)
df.tfs <- df.tfs[order(df.tfs$time_point), ]
tfs.timepoint <- df.tfs$time_point
names(tfs.timepoint) <- df.tfs$tfs

# plot whole GRN
p <- GRNPlot(df.grn,
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 42, 
             plot.importance = TRUE,
             min.importance = 2,
             remove.isolated = FALSE)

png(paste0(temp_plot_path, 'Network_filtered.png'), height = 30, width = 45, units = 'cm', res = 400)
print(p)
graphics.off()

# only plot TFs
p <- GRNPlot(df.grn,
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             plot.importance = TRUE,
             genes.use = df.tfs$tfs,
             remove.isolated = TRUE)

png(paste0(temp_plot_path, 'Network_filtered_TFs.png'), height = 30, width = 45, units = 'cm', res = 400)
print(p)
graphics.off()

############################## Top factors of filtered network #######################################

# make network
dgg <- graph.edgelist(as.matrix(df.grn[,1:2]), directed = T)

########## EDGES
# extract number of edges
edges <- as.data.frame(igraph::degree(dgg))
edges <- rownames_to_column(edges, var = "node")
colnames(edges)[2] <- "nEdges"
edges <- edges %>% arrange(desc(nEdges))
print(edges[1:20,])

png(paste0(temp_plot_path, 'Edges_hist.png'), height = 10, width = 20, units = 'cm', res = 400)
hist(edges$nEdges, breaks = 100)
graphics.off()

factors <- edges[1:10, 1]

# plot TF-target corr
df.tf.gene.subset <- df.tf.gene %>%
  dplyr::filter(tf %in% factors)
df.tfs.subset <- df.tfs %>%
  dplyr::filter(tfs %in% factors)
ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = 1)

png(paste0(temp_plot_path, 'Top_edges_TF_gene_corr_heatmap.png'), height = 10, width = 20, units = 'cm', res = 400)
ht
graphics.off()

########## scMEGA IMPORTANCE
# have to extract manually top 40 factors:
factors <- c("ZNF384", "TCF3", "HMBOX1", "TEAD4", "TFAP2A",
             "ETV1", "MEIS1", "GATA3", "KLF6", "E2F8",
             "FOXK1", "RFX2", "MTF1", "ETV6", "REL",
             "IRF7", "ZBTB7A", "KLF11", "E2F6", "MYCN", 
             "THRB", "OTX1", "REB1", "MNT", "MSX1", 
             "PPARD", "TGIF1", "ZNF143", "POU2E1","NFYC", 
             "SP8", "MEF2D", "JUN", "NR2C2", "KLF3", 
             "ZBTB14", "MAFK", "HES5", "TFAP2E", "LIN54")

# plot TF-target corr of top 40 factors
df.tf.gene.subset <- df.tf.gene %>%
  dplyr::filter(tf %in% factors)
df.tfs.subset <- df.tfs %>%
  dplyr::filter(tfs %in% factors)
ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = 1)

png(paste0(temp_plot_path, 'Top_importance_TF_gene_corr_heatmap.png'), height = 10, width = 20, units = 'cm', res = 400)
ht
graphics.off()

######################################################################################
##############################    POS CORR GRN     ###################################
######################################################################################

temp_plot_path = paste0(plot_path, "filtered__pos_corr_network/")
dir.create(temp_plot_path, recursive = T)

temp_csv_path = paste0(csv_path, "filtered__pos_corr_network/")
dir.create(temp_csv_path, recursive = T)

############################## Plot filtered pos corr GRN #######################################

# plot whole GRN
p <- GRNPlot(df.grn.pos,
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 42, 
             plot.importance = TRUE,
             min.importance = 2,
             remove.isolated = FALSE)

png(paste0(temp_plot_path, 'Network_filtered_pos_corr.png'), height = 30, width = 45, units = 'cm', res = 400)
print(p)
graphics.off()

# only plot TFs
p <- GRNPlot(df.grn.pos,
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             plot.importance = TRUE,
             genes.use = df.tfs$tfs,
             remove.isolated = TRUE)

png(paste0(temp_plot_path, 'Network_filtered_pos_corrTFs.png'), height = 30, width = 45, units = 'cm', res = 400)
print(p)
graphics.off()

############################## Top factors of filtered network #######################################

# make network
dgg <- graph.edgelist(as.matrix(df.grn.pos[,1:2]), directed = T)

########## EDGES
temp_plot_path_subset = paste0(temp_plot_path, "top_edges/")
dir.create(temp_plot_path_subset, recursive = T)
k = 6

edges <- as.data.frame(igraph::degree(dgg))
edges <- rownames_to_column(edges, var = "node")
colnames(edges)[2] <- "nEdges"
edges <- edges %>% arrange(desc(nEdges))
print(edges[1:20,])

png(paste0(temp_plot_path_subset, 'Edges_hist.png'), height = 10, width = 20, units = 'cm', res = 400)
hist(edges$nEdges, breaks = 100)
graphics.off()

factors <- edges[1:10, 1]

# plot corr heatmap
df.tf.gene.subset <- df.tf.gene %>%
  dplyr::filter(tf %in% factors)
df.tfs.subset <- df.tfs %>%
  dplyr::filter(tfs %in% factors)
ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = 1)

png(paste0(temp_plot_path_subset, 'TF_gene_corr_heatmap.png'), height = 10, width = 20, units = 'cm', res = 400)
ht
graphics.off()

# Create target gene heatmap
target_genes_df <- extract_target_genes_df(factors, df.grn.pos)
colnames(target_genes_df) <- paste0(colnames(target_genes_df), " - ", colSums((target_genes_df)))
hm <- pheatmap::pheatmap(target_genes_df,
                         color = c("grey", "purple"),  # Color scheme
                         cluster_rows = TRUE,  # Do not cluster rows
                         cluster_cols = TRUE,  # Do not cluster columns
                         fontsize_row = 0.1,  # Font size for row labels
                         fontsize_col = 10,   # Font size for column labels
                         cutree_rows = k)
png(paste0(temp_plot_path_subset, 'Targets_heatmap.png'), height = 18, width = 10, units = 'cm', res = 400)
hm
graphics.off()

# GO analysis on each cluster of targets
df_row_cluster = data.frame(cluster = cutree(hm$tree_row, k = k))
for (i in 1:k){
  print(i)
  targets <- rownames(df_row_cluster %>% dplyr::filter(cluster == i))
  go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
  if (nrow(as.data.frame(go_output)) > 0){
    png(paste0(temp_plot_path_subset, 'Target_genes_cluster_', i, '_GO_plot.png'), height = 10, width = 20, units = 'cm', res = 400)
    print(plot(barplot(go_output, showCategory = 20)))
    graphics.off()
  }
}

# GO analysis of each TF's target genes
for (i in length(factors)){
  TF <- factors[i]
  targets <- rownames(target_genes_df)[as.logical(target_genes_df[,i])]
  go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
  if (nrow(as.data.frame(go_output)) > 0){
    png(paste0(temp_plot_path_subset, 'Target_genes_TF_', TF, '_GO_plot.png'), height = 10, width = 20, units = 'cm', res = 400)
    print(plot(barplot(go_output, showCategory = 20)))
    graphics.off()
  }
}

# subset network
df.grn.pos.selected <- df.grn.pos %>%
  dplyr::filter(tf %in% factors) %>%
  dplyr::filter(gene %in% factors)
nrow(df.grn.pos.selected) # 15,345 interactions

# TF network numbers
df <- data.frame(
  nSource = length(unique(df.grn.pos.selected$tf)),
  nTarget = length(unique(df.grn.pos.selected$gene)),
  nBoth = sum(unique(df.grn.pos.selected$tf) %in% unique(df.grn.pos.selected$gene)),
  nInteractions = nrow(df.grn.pos.selected),
  nPositiveInteractions = length(which(df.grn.pos.selected$correlation > 0)),
  nNegativeInteractions = length(which(df.grn.pos.selected$correlation < 0))
)
png(paste0(temp_plot_path_subset, 'Subset_network_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# save network
write_tsv(df.grn.pos.selected, file = paste0(temp_plot_path_subset, "GRN_subset.txt"))

########## scMEGA IMPORTANCE
# have to extract manually top 15 factors:
factors <- c("TCF3", "HMBOX1", "TEAD4", "FOXK1", "GATA3",
             "ETV1", "REL", "TFAP2A", "KLF6", "THRB",
             "PPARD", "TGIF1", "MSX1", "HES5", "MEF2D")

temp_plot_path_subset = paste0(temp_plot_path, "top_importance/")
dir.create(temp_plot_path_subset, recursive = T)
k = 8

# plot corr heatmap
df.tf.gene.subset <- df.tf.gene %>%
  dplyr::filter(tf %in% factors)
df.tfs.subset <- df.tfs %>%
  dplyr::filter(tfs %in% factors)
ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = 1)

png(paste0(temp_plot_path_subset, 'TF_gene_corr_heatmap.png'), height = 10, width = 20, units = 'cm', res = 400)
ht
graphics.off()

# Create target gene heatmap
target_genes_df <- extract_target_genes_df(factors, df.grn.pos)
colnames(target_genes_df) <- paste0(colnames(target_genes_df), " - ", colSums((target_genes_df)))
hm <- pheatmap::pheatmap(target_genes_df,
                         color = c("grey", "purple"),  # Color scheme
                         cluster_rows = TRUE,  # Do not cluster rows
                         cluster_cols = TRUE,  # Do not cluster columns
                         fontsize_row = 0.1,  # Font size for row labels
                         fontsize_col = 10,   # Font size for column labels
                         cutree_rows = k)
png(paste0(temp_plot_path_subset, 'Targets_heatmap.png'), height = 18, width = 10, units = 'cm', res = 400)
hm
graphics.off()

# GO analysis on each cluster of targets
df_row_cluster = data.frame(cluster = cutree(hm$tree_row, k = k))
for (i in 1:k){
  print(i)
  targets <- rownames(df_row_cluster %>% dplyr::filter(cluster == i))
  go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
  if (nrow(as.data.frame(go_output)) > 0){
    png(paste0(temp_plot_path_subset, 'Target_genes_cluster_', i, '_GO_plot.png'), height = 10, width = 20, units = 'cm', res = 400)
    print(plot(barplot(go_output, showCategory = 20)))
    graphics.off()
  }
}

# GO analysis of each TF's target genes
for (i in length(factors)){
  TF <- factors[i]
  targets <- rownames(target_genes_df)[as.logical(target_genes_df[,i])]
  go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
  if (nrow(as.data.frame(go_output)) > 0){
    png(paste0(temp_plot_path_subset, 'Target_genes_TF_', TF, '_GO_plot.png'), height = 10, width = 20, units = 'cm', res = 400)
    print(plot(barplot(go_output, showCategory = 20)))
    graphics.off()
  }
}

# subset network
df.grn.pos.selected <- df.grn.pos %>%
  dplyr::filter(tf %in% factors) %>%
  dplyr::filter(gene %in% factors)
nrow(df.grn.pos.selected) # 15,345 interactions

# TF network numbers
df <- data.frame(
  nSource = length(unique(df.grn.pos.selected$tf)),
  nTarget = length(unique(df.grn.pos.selected$gene)),
  nBoth = sum(unique(df.grn.pos.selected$tf) %in% unique(df.grn.pos.selected$gene)),
  nInteractions = nrow(df.grn.pos.selected),
  nPositiveInteractions = length(which(df.grn.pos.selected$correlation > 0)),
  nNegativeInteractions = length(which(df.grn.pos.selected$correlation < 0))
)
png(paste0(temp_plot_path_subset, 'Subset_network_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# save network
write_tsv(df.grn.pos.selected, file = paste0(temp_plot_path_subset, "GRN_subset.txt"))
