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
    
    data_path = "./output/NF-downstream_analysis/Processing/FullData/scMEGA/MEGA_GRNi_motif_cutoff/"
    
    plot_path = "./output/NF-downstream_analysis/Processing/FullData/scMEGA/MEGA_GRN_vis/plots/"
    csv_path = "./output/NF-downstream_analysis/Processing/FullData/scMEGA/MEGA_GRN_vis/csv_files/"
    
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    data_path = "./input/"
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

GRNPlot_updated <- function (df.grn, tfs.use = NULL, show.tf.labels = TRUE, tfs.timepoint = NULL, 
          genes.cluster = NULL, genes.use = NULL, genes.highlight = NULL, 
          cols.highlight = "#984ea3", seed = 42, plot.importance = TRUE, return.importance = TRUE,
          min.importance = 2, remove.isolated = FALSE) 
{
  if (is.null(tfs.timepoint)) {
    stop("Need time point for each TF!")
  }
  if (!is.null(tfs.use)) {
    df.grn <- subset(df.grn, tf %in% tfs.use)
  }
  if (!is.null(genes.use)) {
    df.grn <- subset(df.grn, gene %in% genes.use)
  }
  tf.list <- unique(df.grn$tf)
  gene.list <- setdiff(unique(df.grn$gene), tf.list)
  g <- igraph::graph_from_data_frame(df.grn, directed = TRUE)
  if (remove.isolated) {
    isolated <- which(degree(g) == 0)
    g <- igraph::delete.vertices(g, isolated)
  }
  pagerank <- page_rank(g, weights = E(g)$weights)
  bet <- betweenness(g, weights = E(g)$weights, normalized = TRUE)
  df_measure <- data.frame(tf = V(g)$name, pagerank = pagerank$vector, 
                           betweenness = bet) %>% subset(tf %in% df.grn$tf) %>% 
    mutate(pagerank = scale(pagerank)[, 1]) %>% mutate(betweenness = scale(betweenness)[, 
                                                                                        1])
  min.page <- min(df_measure$pagerank)
  min.bet <- min(df_measure$betweenness)
  df_measure$importance <- sqrt((df_measure$pagerank - min.page)^2 + 
                                  (df_measure$betweenness - min.bet)^2)
  if (plot.importance) {
    p <- ggplot(data = df_measure) + aes(x = reorder(tf, 
                                                     -importance), y = importance) + geom_point() + xlab("TFs") + 
      ylab("Importance") + cowplot::theme_cowplot() + theme(axis.text.x = element_text(angle = 60, 
                                                                                       hjust = 1))
    print(p)
  }
  if (return.importance) {
    return(df_measure)
  }
  df_measure_sub <- subset(df_measure, importance > 2)
  tf_size <- df_measure$importance
  names(tf_size) <- df_measure$tf
  gene_size <- rep(min(df_measure$importance), length(unique(df.grn$gene)))
  names(gene_size) <- gene.list
  v_size <- c(tf_size, gene_size)
  V(g)$size <- v_size[V(g)$name]
  cols.tf <- ArchR::paletteContinuous(set = "blueYellow", n = length(tfs.timepoint))
  names(cols.tf) <- names(tfs.timepoint)
  if (is.null(genes.cluster)) {
    cols.gene <- rep("gray", length(gene.list))
    names(cols.gene) <- gene.list
  }
  else {
    genes.cluster <- genes.cluster %>% subset(gene %in% gene.list)
    cols <- ArchR::paletteDiscrete(values = as.character(genes.cluster$cluster))
    df.gene <- lapply(1:length(cols), function(x) {
      df <- subset(genes.cluster, cluster == x)
      df$color <- rep(cols[[x]], nrow(df))
      return(df)
    }) %>% Reduce(rbind, .)
    cols.gene <- df.gene$color
    names(cols.gene) <- df.gene$gene
  }
  v_color <- c(cols.tf, cols.gene)
  v_color <- v_color[V(g)$name]
  tf_alpha <- rep(1, length(tf.list))
  gene_alpha <- rep(0.5, length(gene.list))
  names(tf_alpha) <- tf.list
  names(gene_alpha) <- gene.list
  v_alpha <- c(tf.list, gene.list)
  V(g)$alpha <- v_alpha[V(g)$name]
  set.seed(seed)
  layout <- layout_with_fr(g, weights = E(g)$weights, dim = 2, 
                           niter = 1000)
  p <- ggraph(g, layout = layout) + geom_edge_link(edge_colour = "gray", 
                                                   edge_alpha = 0.25) + geom_node_point(aes(size = V(g)$size, 
                                                                                            color = as.factor(name), alpha = V(g)$alpha), show.legend = FALSE) + 
    scale_size(range = c(1, 10)) + scale_color_manual(values = v_color)
  if (show.tf.labels) {
    p <- p + geom_node_label(aes(filter = V(g)$name %in% 
                                   tf.list, label = V(g)$name), repel = TRUE, hjust = "inward", 
                             color = "#ff7f00", size = 5, show.legend = FALSE, 
                             max.overlaps = Inf)
  }
  if (!is.null(genes.highlight)) {
    p <- p + geom_node_label(aes(filter = V(g)$name %in% 
                                   genes.highlight, label = V(g)$name), repel = TRUE, 
                             hjust = "inward", size = 5, color = cols.highlight, 
                             show.legend = FALSE)
  }
  p <- p + theme_void()
  return(p)
}

## function to write out gene lists in nested list
export_gene_list <- function(gene_list, publish_dir){
  for(name in names(gene_list)){
    print(name)
    write(paste0(name, ": ", paste0(gene_list[[name]], collapse = ", ")), file = paste0(publish_dir, '.txt'), append = TRUE)
  }
}

############################## Read in data #######################################

print("reading in data...")

TF_names <- read_tsv(paste0(data_path, "csv_files/Known_TF_names.txt"), col_names = FALSE)
TF_names <- TF_names$X1

temp_data_path = paste0(data_path, "./csv_files/placodal_lineage/")

# read in csvs
df.tfs <- read.csv(paste0(temp_data_path, "TF_correlations_all.csv"))
df.p2g <- read.csv(paste0(temp_data_path, "Target_genes_with_matched_enhancers.csv"))
df.tf.gene <- read.csv(paste0(temp_data_path, "TF_to_gene_correlations.csv"))

# read in networks
df.grn.un <- read.csv(paste0(temp_data_path, "GRN_initial.txt"))
df.grn <- read_tsv(paste0(temp_data_path, "GRN_filtered.txt"))
df.grn.pos <- read_tsv(paste0(temp_data_path, "GRN_filtered_pos_corr.txt"))

# read in data object
obj.traj <- readRDS(paste0(data_path, "rds_files/Placodal_traj_obj.RDS"))

# read in RNA data object
seurat <- readRDS(paste0(data_path, "seurat_label_transfer_minus_HH4.RDS"))

######################################################################################
##############################    FILTERED GRN     ###################################
######################################################################################

print("Filtered network plotting...")

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

print("Filtered network analysis...")

########## EDGES

print("Edges analysis...")

temp_plot_path_subset = paste0(temp_plot_path, "top_edges/")
dir.create(temp_plot_path_subset, recursive = T)
k = 1

# analyse edges
dgg <- graph.edgelist(as.matrix(df.grn[,1:2]), directed = T)
edges <- as.data.frame(igraph::degree(dgg))
edges <- rownames_to_column(edges, var = "node")
colnames(edges)[2] <- "nEdges"
edges <- edges %>% arrange(desc(nEdges))
print(edges[1:20,])
write_tsv(edges[1:100,], file = paste0(temp_plot_path_subset, "Edges_top_100_df.txt"))

png(paste0(temp_plot_path_subset, 'Edges_hist.png'), height = 10, width = 20, units = 'cm', res = 400)
hist(edges$nEdges, breaks = 100)
graphics.off()

factors <- edges[1:10, 1]
print("Top 10 edges factors:")
print(factors)

# # Pseudotime plots
# obj.traj <- AddTargetAssay_updated(object = obj.traj, df.grn = df.grn)
# for (TF in factors){
#   print(TF)
#   p <- PseudotimePlot_updated(object = obj.traj, tf.use = TF, trajectory.name = "lineage_placodal_probability")
#   png(paste0(temp_plot_path_subset, 'Pseudotime_plot_', TF, '.png'), height = 8, width = 18, units = 'cm', res = 400)
#   print(p)
#   graphics.off()
# }

# plot corr heatmap
df.tf.gene.subset <- df.tf.gene %>%
  dplyr::filter(tf %in% factors)
df.tfs.subset <- df.tfs %>%
  dplyr::filter(tfs %in% factors)
ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = 1)

png(paste0(temp_plot_path_subset, 'TF_gene_corr_heatmap.png'), height = 10, width = 20, units = 'cm', res = 400)
ht
graphics.off()

# # Create target gene heatmap
# target_genes_df <- extract_target_genes_df(factors, df.grn)
# colnames(target_genes_df) <- paste0(colnames(target_genes_df), " - ", colSums((target_genes_df)))
# hm <- pheatmap::pheatmap(t(target_genes_df),
#                          color = c("grey", "purple"),  # Color scheme
#                          cluster_rows = TRUE,  # Do not cluster rows
#                          cluster_cols = TRUE,  # Do not cluster columns
#                          fontsize_row = 10,  # Font size for row labels
#                          fontsize_col = 0.0001,   # Font size for column labels
#                          cutree_cols = k)
# png(paste0(temp_plot_path_subset, 'Targets_heatmap.png'), height = 10, width = 18, units = 'cm', res = 400)
# hm
# graphics.off()

# # GO analysis on each cluster of targets
# df_row_cluster = data.frame(cluster = cutree(hm$tree_col, k = k))
# for (i in 1:k){
#   print(i)
#   targets <- rownames(df_row_cluster %>% dplyr::filter(cluster == i))
#   go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
#   if (nrow(as.data.frame(go_output)) > 0){
#     png(paste0(temp_plot_path_subset, 'Target_genes_cluster_', i, '_GO_plot.png'), height = 30, width = 20, units = 'cm', res = 400)
#     print(plot(barplot(go_output, showCategory = 20)))
#     graphics.off()
#   }
# }

# # GO analysis of each TF's target genes
# for (i in length(factors)){
#   TF <- factors[i]
#   targets <- rownames(target_genes_df)[as.logical(target_genes_df[,i])]
#   go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
#   if (nrow(as.data.frame(go_output)) > 0){
#     png(paste0(temp_plot_path_subset, 'Target_genes_TF_', TF, '_GO_plot.png'), height = 30, width = 20, units = 'cm', res = 400)
#     print(plot(barplot(go_output, showCategory = 20)))
#     graphics.off()
#   }
# }

# # subset network
# df.grn.selected <- df.grn %>%
#   dplyr::filter(tf %in% factors) %>%
#   dplyr::filter(gene %in% factors)
# nrow(df.grn.selected)

# # TF network numbers
# df <- data.frame(
#   nSource = length(unique(df.grn.selected$tf)),
#   nTarget = length(unique(df.grn.selected$gene)),
#   nBoth = sum(unique(df.grn.selected$tf) %in% unique(df.grn.selected$gene)),
#   nInteractions = nrow(df.grn.selected),
#   nPositiveInteractions = length(which(df.grn.selected$correlation > 0)),
#   nNegativeInteractions = length(which(df.grn.selected$correlation < 0))
# )
# png(paste0(temp_plot_path_subset, 'Subset_network_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(df, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# # save network
# write_tsv(df.grn.selected, file = paste0(temp_plot_path_subset, "GRN_subset.txt"))

########## scMEGA IMPORTANCE

print("Importance analysis...")

temp_plot_path_subset = paste0(temp_plot_path, "top_importance/")
dir.create(temp_plot_path_subset, recursive = T)
k = 1

# plot importance ranking + save the df
png(paste0(temp_plot_path_subset, 'Top_importance_plot.png'), height = 10, width = 120, units = 'cm', res = 400)
importance_df <- GRNPlot_updated(df.grn,
                                 tfs.timepoint = tfs.timepoint,
                                 show.tf.labels = TRUE,
                                 seed = 42, 
                                 plot.importance = TRUE,
                                 min.importance = 2,
                                 remove.isolated = FALSE,
                                 return.importance = TRUE)
graphics.off()
importance_df <- arrange(importance_df, by = desc(importance))
write_tsv(importance_df, file = paste0(temp_plot_path_subset, "Importance_df.txt"))

# extract top 20 factors:
factors <- importance_df[1:20, 1]
print("Top 20 importance factors:")
print(factors)

# # Pseudotime plots
# for (TF in factors){
#   print(TF)
#   p <- PseudotimePlot_updated(object = obj.traj, tf.use = TF, trajectory.name = "lineage_placodal_probability")
#   png(paste0(temp_plot_path_subset, 'Pseudotime_plot_', TF, '.png'), height = 8, width = 18, units = 'cm', res = 400)
#   print(p)
#   graphics.off()
# }

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
target_genes_df <- extract_target_genes_df(factors, df.grn)
colnames(target_genes_df) <- paste0(colnames(target_genes_df), " - ", colSums((target_genes_df)))
hm <- pheatmap::pheatmap(t(target_genes_df),
                         color = c("grey", "purple"),  # Color scheme
                         cluster_rows = TRUE,  # Do not cluster rows
                         cluster_cols = TRUE,  # Do not cluster columns
                         fontsize_row = 10,  # Font size for row labels
                         fontsize_col = 0.0001,   # Font size for column labels
                         cutree_cols = k)
png(paste0(temp_plot_path_subset, 'Targets_heatmap.png'), height = 10, width = 18, units = 'cm', res = 400)
hm
graphics.off()

# GO analysis on each cluster of targets
df_row_cluster = data.frame(cluster = cutree(hm$tree_col, k = k))
for (i in 1:k){
  print(i)
  targets <- rownames(df_row_cluster %>% dplyr::filter(cluster == i))
  go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
  if (nrow(as.data.frame(go_output)) > 0){
    png(paste0(temp_plot_path_subset, 'Target_genes_cluster_', i, '_GO_plot.png'), height = 30, width = 20, units = 'cm', res = 400)
    print(plot(barplot(go_output, showCategory = 20)))
    graphics.off()
  }
}

# GO analysis of each TF's target genes
for (i in length(factors)){
  TF <- factors[i]
  print(TF)
  targets <- rownames(target_genes_df)[as.logical(target_genes_df[,i])]
  go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
  if (nrow(as.data.frame(go_output)) > 0){
    png(paste0(temp_plot_path_subset, 'Target_genes_TF_', TF, '_GO_plot.png'), height = 30, width = 20, units = 'cm', res = 400)
    print(plot(barplot(go_output, showCategory = 20)))
    graphics.off()
  }
}

# subset network
df.grn.selected <- df.grn %>%
  dplyr::filter(tf %in% factors) %>%
  dplyr::filter(gene %in% factors)
nrow(df.grn.selected) # 15,345 interactions

# TF network numbers
df <- data.frame(
  nSource = length(unique(df.grn.selected$tf)),
  nTarget = length(unique(df.grn.selected$gene)),
  nBoth = sum(unique(df.grn.selected$tf) %in% unique(df.grn.selected$gene)),
  nInteractions = nrow(df.grn.selected),
  nPositiveInteractions = length(which(df.grn.selected$correlation > 0)),
  nNegativeInteractions = length(which(df.grn.selected$correlation < 0))
)
png(paste0(temp_plot_path_subset, 'Subset_network_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
             tableGrob(df, rows=NULL, theme = ttheme_minimal()))
graphics.off()

# save network
write_tsv(df.grn.selected, file = paste0(temp_plot_path_subset, "GRN_subset.txt"))


######################################################################################
##############################    POS CORR GRN     ###################################
######################################################################################

temp_plot_path = paste0(plot_path, "filtered_pos_corr_network/")
dir.create(temp_plot_path, recursive = T)

temp_csv_path = paste0(csv_path, "filtered_pos_corr_network/")
dir.create(temp_csv_path, recursive = T)

############################## Plot filtered pos corr GRN #######################################

print("Pos corr network plotting...")

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

############################## Top factors of filtered pos corr network #######################################

print("Pos corr network analysis...")

########## EDGES

print("Edges analysis...")

temp_plot_path_subset = paste0(temp_plot_path, "top_edges/")
dir.create(temp_plot_path_subset, recursive = T)
k = 3

# analyse edges
dgg <- graph.edgelist(as.matrix(df.grn.pos[,1:2]), directed = T)
edges <- as.data.frame(igraph::degree(dgg))
edges <- rownames_to_column(edges, var = "node")
colnames(edges)[2] <- "nEdges"
edges <- edges %>% arrange(desc(nEdges))
print(edges[1:20,])
write_tsv(edges[1:100,], file = paste0(temp_plot_path_subset, "Edges_top_100_df.txt"))

png(paste0(temp_plot_path_subset, 'Edges_hist.png'), height = 10, width = 20, units = 'cm', res = 400)
hist(edges$nEdges, breaks = 100)
graphics.off()

factors <- edges[1:10, 1]
print("Top 10 edges factors:")
print(factors)

# # Pseudotime plots
# obj.traj <- AddTargetAssay_updated(object = obj.traj, df.grn = df.grn.pos)
# for (TF in factors){
#   print(TF)
#   p <- PseudotimePlot_updated(object = obj.traj, tf.use = TF, trajectory.name = "lineage_placodal_probability")
#   png(paste0(temp_plot_path_subset, 'Pseudotime_plot_', TF, '.png'), height = 8, width = 18, units = 'cm', res = 400)
#   print(p)
#   graphics.off()
# }

# plot corr heatmap
df.tf.gene.subset <- df.tf.gene %>%
  dplyr::filter(tf %in% factors)
df.tfs.subset <- df.tfs %>%
  dplyr::filter(tfs %in% factors)
ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = 1)

png(paste0(temp_plot_path_subset, 'TF_gene_corr_heatmap.png'), height = 10, width = 20, units = 'cm', res = 400)
ht
graphics.off()

# # Create target gene heatmap
# target_genes_df <- extract_target_genes_df(factors, df.grn.pos)
# colnames(target_genes_df) <- paste0(colnames(target_genes_df), " - ", colSums((target_genes_df)))
# hm <- pheatmap::pheatmap(t(target_genes_df),
#                          color = c("grey", "purple"),  # Color scheme
#                          cluster_rows = TRUE,  # Do not cluster rows
#                          cluster_cols = TRUE,  # Do not cluster columns
#                          fontsize_row = 10,  # Font size for row labels
#                          fontsize_col = 0.0001,   # Font size for column labels
#                          cutree_cols = k)
# png(paste0(temp_plot_path_subset, 'Targets_heatmap.png'), height = 10, width = 18, units = 'cm', res = 400)
# hm
# graphics.off()

# # GO analysis on each cluster of targets
# df_row_cluster = data.frame(cluster = cutree(hm$tree_col, k = k))
# for (i in 1:k){
#   print(i)
#   targets <- rownames(df_row_cluster %>% dplyr::filter(cluster == i))
#   go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
#   if (nrow(as.data.frame(go_output)) > 0){
#     png(paste0(temp_plot_path_subset, 'Target_genes_cluster_', i, '_GO_plot.png'), height = 30, width = 20, units = 'cm', res = 400)
#     print(plot(barplot(go_output, showCategory = 20)))
#     graphics.off()
#   }
# }

# # GO analysis of each TF's target genes
# for (i in length(factors)){
#   TF <- factors[i]
#   targets <- rownames(target_genes_df)[as.logical(target_genes_df[,i])]
#   go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
#   if (nrow(as.data.frame(go_output)) > 0){
#     png(paste0(temp_plot_path_subset, 'Target_genes_TF_', TF, '_GO_plot.png'), height = 30, width = 20, units = 'cm', res = 400)
#     print(plot(barplot(go_output, showCategory = 20)))
#     graphics.off()
#   }
# }

# # subset network
# df.grn.pos.selected <- df.grn.pos %>%
#   dplyr::filter(tf %in% factors) %>%
#   dplyr::filter(gene %in% factors)
# nrow(df.grn.pos.selected)

# # TF network numbers
# df <- data.frame(
#   nSource = length(unique(df.grn.pos.selected$tf)),
#   nTarget = length(unique(df.grn.pos.selected$gene)),
#   nBoth = sum(unique(df.grn.selected$tf) %in% unique(df.grn.pos.selected$gene)),
#   nInteractions = nrow(df.grn.selected),
#   nPositiveInteractions = length(which(df.grn.pos.selected$correlation > 0)),
#   nNegativeInteractions = length(which(df.grn.pos.selected$correlation < 0))
# )
# png(paste0(temp_plot_path_subset, 'Subset_network_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
# grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
#              tableGrob(df, rows=NULL, theme = ttheme_minimal()))
# graphics.off()

# # save network
# write_tsv(df.grn.pos.selected, file = paste0(temp_plot_path_subset, "GRN_subset.txt"))

########## scMEGA IMPORTANCE

print("Importance analysis...")

temp_plot_path_subset = paste0(temp_plot_path, "top_importance/")
dir.create(temp_plot_path_subset, recursive = T)
k = 10

# plot importance ranking + save the df
png(paste0(temp_plot_path_subset, 'Importance_plot.png'), height = 10, width = 120, units = 'cm', res = 400)
importance_df <- GRNPlot_updated(df.grn.pos,
                                 tfs.timepoint = tfs.timepoint,
                                 show.tf.labels = TRUE,
                                 seed = 42, 
                                 plot.importance = TRUE,
                                 min.importance = 2,
                                 remove.isolated = FALSE,
                                 return.importance = TRUE)
graphics.off()
importance_df <- arrange(importance_df, by = desc(importance))
write_tsv(importance_df, file = paste0(temp_plot_path_subset, "Analysis_metrics_df.txt"))

# plot the individual metrics - betweeness
betweeness_df <- importance_df %>% 
  dplyr::select(c("tf", "betweenness")) %>%
  arrange(by = desc(betweenness))
png(paste0(temp_plot_path_subset, 'Betweeness_plot.png'), height = 10, width = 40, units = 'cm', res = 400)
barplot(betweeness_df$betweenness, names.arg = betweeness_df$tf, las = 2, cex.names=.5)
graphics.off()

# plot the individual metrics - pagerank
pagerank_df <- importance_df %>% 
  dplyr::select(c("tf", "pagerank")) %>%
  arrange(by = desc(pagerank))
png(paste0(temp_plot_path_subset, 'Pagerank_plot.png'), height = 10, width = 40, units = 'cm', res = 400)
barplot(pagerank_df$pagerank, names.arg = pagerank_df$tf, las = 2, cex.names=.5)
graphics.off()

# extract top 20 factors:
factors <- importance_df[1:20, 1]
print("Top 20 importance factors:")
print(factors)

# # Pseudotime plots
# for (TF in factors){
#   print(TF)
#   p <- PseudotimePlot_updated(object = obj.traj, tf.use = TF, trajectory.name = "lineage_placodal_probability")
#   png(paste0(temp_plot_path_subset, 'Pseudotime_plot_', TF, '.png'), height = 8, width = 18, units = 'cm', res = 400)
#   print(p)
#   graphics.off()
# }

# plot corr heatmap
df.tf.gene.subset <- df.tf.gene %>%
  dplyr::filter(tf %in% factors)
df.tfs.subset <- df.tfs %>%
  dplyr::filter(tfs %in% factors)
ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = k)
ht <- draw(ht)

png(paste0(temp_plot_path_subset, 'TF_gene_corr_heatmap.png'), height = 10, width = 20, units = 'cm', res = 400)
ht
graphics.off()

# extract target gene clusters
mat.cor <- df.tf.gene.subset %>% as.data.frame() %>% 
  select(c(tf, gene, correlation)) %>% 
  tidyr::pivot_wider(names_from = tf, values_from = correlation) %>% 
  textshape::column_to_rownames("gene")
target_gene_clusters <- list()
for (cluster_name in names(row_order(ht))){
  print(paste0("cluster: ", cluster_name))
  indices <- row_order(ht)[[cluster_name]]
  target_gene_clusters[[cluster_name]] <- rownames(mat.cor)[indices]
  print(length(target_gene_clusters[[cluster_name]]))
  
  # print expression of these genes
  seurat <- AddModuleScore(object = seurat, features = list(target_gene_clusters[[cluster_name]]), name = "cluster")
  png(paste0(temp_plot_path_subset, 'FeaturePlot_of_gene_cluster_from_corr_', cluster_name, '.png'), height = 10, width = 20, units = 'cm', res = 400)
  FeaturePlot(seurat, features = "cluster1", pt.size = 1.5)
  graphics.off()
  
}
export_gene_list(target_gene_clusters, publish_dir = paste0(temp_plot_path_subset, "target_gene_clusters_from_corr"))

# Create target gene heatmap
target_genes_df <- extract_target_genes_df(factors, df.grn.pos)
colnames(target_genes_df) <- paste0(colnames(target_genes_df), " - ", colSums((target_genes_df)))
hm <- pheatmap::pheatmap(t(target_genes_df),
                         color = c("grey", "purple"),  # Color scheme
                         cluster_rows = TRUE,  # Do not cluster rows
                         cluster_cols = TRUE,  # Do not cluster columns
                         fontsize_row = 10,  # Font size for row labels
                         fontsize_col = 0.0001,   # Font size for column labels
                         cutree_cols = k)
png(paste0(temp_plot_path_subset, 'Targets_heatmap.png'), height = 10, width = 18, units = 'cm', res = 400)
hm
graphics.off()

# Extract each cluster of targets
df_row_cluster = data.frame(cluster = cutree(hm$tree_col, k = k))
target_gene_direct_clusters <- list()
for (i in 1:k){
  print(i)
  targets <- rownames(df_row_cluster %>% dplyr::filter(cluster == i))
  target_gene_direct_clusters[[i]] <- targets
  # go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
  # if (nrow(as.data.frame(go_output)) > 0){
  #   png(paste0(temp_plot_path_subset, 'Target_genes_cluster_', i, '_GO_plot.png'), height = 30, width = 20, units = 'cm', res = 400)
  #   print(plot(barplot(go_output, showCategory = 20)))
  #   graphics.off()
  # }
}
export_gene_list(target_gene_direct_clusters, publish_dir = paste0(temp_plot_path_subset, "target_gene_clusters_from_direct_interactions"))

# Extract each of the TF's target genes
target_gene_direct <- list()
for (i in 1:length(factors)){
  print(i)
  TF <- factors[i]
  print(TF)
  targets <- rownames(target_genes_df)[as.logical(target_genes_df[,i])]
  target_gene_direct[[TF]] <- targets
  # go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
  # if (nrow(as.data.frame(go_output)) > 0){
  #   png(paste0(temp_plot_path_subset, 'Target_genes_TF_', TF, '_GO_plot.png'), height = 30, width = 20, units = 'cm', res = 400)
  #   print(plot(barplot(go_output, showCategory = 20)))
  #   graphics.off()
  # }
}
export_gene_list(target_gene_direct, publish_dir = paste0(temp_plot_path_subset, "target_genes_from_direct_interactions"))

# subset network
df.grn.pos.selected <- df.grn %>%
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
