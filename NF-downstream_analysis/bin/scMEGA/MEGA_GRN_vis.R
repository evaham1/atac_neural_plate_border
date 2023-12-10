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
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Ggallus.UCSC.galGal6)
library(pheatmap)

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
    
    # addArchRThreads(threads = ncores)
    addArchRThreads(threads = 1)
    
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

TenxPheatmap <- function(data, metadata, col_order = metadata, custom_order = NULL, custom_order_column = NULL, 
                         assay = "RNA", slot = "scale.data", selected_genes,
                         main = '', hide_annotation = NULL, show_rownames = TRUE, hclust_rows = FALSE, hclust_cols = FALSE, gaps_col = NULL,
                         cell_subset = NULL, treeheight_row = 0, use_seurat_colours = TRUE,  annotation_colours = NULL,
                         col_ann_order = rev(metadata)){
  
  if(!is.null(cell_subset)){
    data <- subset(data, cells = cell_subset)
  } else {}
  
  # if there are any character cols in metadata convert to factor
  data@meta.data[sapply(data@meta.data, is.character)] <- lapply(data@meta.data[sapply(data@meta.data, is.character)], as.factor)
  
  # reset levels in seurat_clusters metadata to default numerical order as default
  if("seurat_clusters" %in% metadata){
    data@meta.data[,"seurat_clusters"] <-  factor(data@meta.data[,"seurat_clusters"], levels = sort(unique(as.numeric(as.character(data@meta.data[,"seurat_clusters"])))))
  } else {}
  
  # initialise heatmap metadata
  HM.col <- droplevels(data@meta.data[, metadata, drop=FALSE])
  
  # order HM metadata based on col_order variable
  HM.col <- HM.col[do.call('order', c(HM.col[col_order], list(decreasing=FALSE))), , drop = FALSE]
  
  # order HM metadata based on custom order
  if(!is.null(custom_order)){
    if(is.null(custom_order_column)){
      "custom_order column must be specified \n"
    } else {}
    if(!setequal(custom_order, unique(HM.col[[custom_order_column]]))){
      stop("custom_order factors missing from custom_order_column \n\n")
    } else {}
    HM.col[[custom_order_column]] <-  factor(HM.col[[custom_order_column]], levels = custom_order)
    HM.col <- HM.col[order(HM.col[[custom_order_column]]),,drop = FALSE]  
  }
  
  # gaps_col specifies a metadata column which column gaps are calculated from
  if(!is.null(gaps_col)) {
    if(class(gaps_col) != "character"){
      stop("gaps_col must be a metadata column name")
    } else {
      gaps_col = cumsum(rle(as.vector(HM.col[[gaps_col]]))[["lengths"]])
    }
  } else {
  }
  
  # hide as many annotations in metadata as desired with hide_annotation
  if(!is.null(hide_annotation)){
    HM.col[,hide_annotation] <- NULL
  } else{}
  
  if(use_seurat_colours == FALSE){
    # must have given annotation colours
    if(is.null(annotation_colours)){
      "annotation_colours column must be specified \n"
    } else {}
    ### use custom colours
    cols <- list()
    cols[[metadata]] <- annotation_colours
    cols[[metadata]] <- cols[[metadata]][match(levels(HM.col[[metadata]]), names(cols[[metadata]]))]
  } else {
    # set colours ggplot default colours, as in Seurat::DimPlot
    ann_colours <- list()
    for(column in metadata){
      ann_colours[[column]] <- setNames(ggPlotColours(n = length(levels(data@meta.data[,column]))), levels(data@meta.data[,column]))
      
      # change levels of HM col so that heatmap annotations are in the same order as plotted
      ann_colours[[column]] <- ann_colours[[column]][match(levels(HM.col[[column]]), names(ann_colours[[column]]))]
    }
  }
  
  # extract data to plot from seurat object
  new.dat <- t(as.matrix(x = GetAssayData(object = data, assay = assay, slot = slot)[selected_genes, rownames(HM.col), drop = FALSE]))
  if(!is.null(cell_subset)){
    cat("rescaling data as cells have been subset \n")
    new.dat <- t(scale(t(new.dat)))
  } else {}
  new.dat <- replace(new.dat, new.dat >= 2, 2)
  new.dat <- replace(new.dat, new.dat <= -2, -2)
  
  print(pheatmap(t(new.dat), color = PurpleAndYellow(),
                 cluster_rows = hclust_rows, cluster_cols = hclust_cols, show_colnames = F,
                 annotation_col = HM.col[, rev(col_ann_order), drop = FALSE], fontsize = 22, fontsize_row = 12, gaps_col = gaps_col,
                 main = main, show_rownames = show_rownames, annotation_colors = cols, treeheight_row = treeheight_row))
}


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

############################## Set colours and annotation order #######################################

scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak', 
                              'PGC', 'BI', 'meso', 'endo',
                              'Neural', 'Placodal', 'Non-neural', 'Contam')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", 
                                         "#53A651", "#6D8470", "#87638F", "#A5548D", "#C96555", "#ED761C", 
                                         "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30", "#CC9F2C", "#AD6428", 
                                         "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                         "#786D73", "#581845", "#9792A3", "#BBB3CB",
                                         "#A5718D", "#3F918C", "#ed5e5f", "#9792A3")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                                                                'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                                                                'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 
                                                                                'MB','vFB', 'aNP', 'node', 'FB', 'pEpi',
                                                                                'PGC', 'BI', 'meso', 'endo',
                                                                                'Neural', 'Placodal', 'Non-neural', 'Contam')


############################## Read in data #######################################

print("reading in data...")

TF_names <- read_tsv(paste0(data_path, "csv_files/Known_TF_names.txt"), col_names = FALSE)
TF_names <- TF_names$X1

temp_data_path = paste0(data_path, "./csv_files/placodal_lineage/")

# read in csvs
df.tfs <- read.csv(paste0(temp_data_path, "TF_correlations_all.csv"))
df.p2g <- read.csv(paste0(temp_data_path, "Target_nodes_with_matched_enhancers.csv"))
df.tf.gene <- read.csv(paste0(temp_data_path, "TF_to_gene_correlations.csv"))

# read in networks
df.grn.un <- read.csv(paste0(temp_data_path, "GRN_initial.txt"))
df.grn <- read_tsv(paste0(temp_data_path, "GRN_filtered.txt"))
df.grn.pos <- read_tsv(paste0(temp_data_path, "GRN_filtered_pos_corr.txt"))

# read in data object
obj.traj <- readRDS(paste0(data_path, "rds_files/Placodal_traj_obj.RDS"))

# read in RNA data object
seurat <- readRDS(paste0(data_path, "seurat_label_transfer_minus_HH4.RDS"))

# read in ATAC data object
# ArchR <- loadArchRProject(path = paste0(data_path, "ss8_Save-ArchR"), force = FALSE, showLogo = TRUE)
ArchR <- loadArchRProject(path = paste0(data_path, "FullData_Save-ArchR"), force = FALSE, showLogo = TRUE)
getCellColData(ArchR)

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

# ########## EDGES

# print("Edges analysis...")

# temp_plot_path_subset = paste0(temp_plot_path, "top_edges/")
# dir.create(temp_plot_path_subset, recursive = T)
# k = 1

# # analyse edges
# dgg <- graph.edgelist(as.matrix(df.grn[,1:2]), directed = T)
# edges <- as.data.frame(igraph::degree(dgg))
# edges <- rownames_to_column(edges, var = "node")
# colnames(edges)[2] <- "nEdges"
# edges <- edges %>% arrange(desc(nEdges))
# print(edges[1:20,])
# write_tsv(edges[1:100,], file = paste0(temp_plot_path_subset, "Edges_top_100_df.txt"))

# png(paste0(temp_plot_path_subset, 'Edges_hist.png'), height = 10, width = 20, units = 'cm', res = 400)
# hist(edges$nEdges, breaks = 100)
# graphics.off()

# factors <- edges[1:10, 1]
# print("Top 10 edges factors:")
# print(factors)

# # # Pseudotime plots
# # obj.traj <- AddTargetAssay_updated(object = obj.traj, df.grn = df.grn)
# # for (TF in factors){
# #   print(TF)
# #   p <- PseudotimePlot_updated(object = obj.traj, tf.use = TF, trajectory.name = "lineage_placodal_probability")
# #   png(paste0(temp_plot_path_subset, 'Pseudotime_plot_', TF, '.png'), height = 8, width = 18, units = 'cm', res = 400)
# #   print(p)
# #   graphics.off()
# # }

# # plot corr heatmap
# df.tf.gene.subset <- df.tf.gene %>%
#   dplyr::filter(tf %in% factors)
# df.tfs.subset <- df.tfs %>%
#   dplyr::filter(tfs %in% factors)
# ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = 1)

# png(paste0(temp_plot_path_subset, 'TF_gene_corr_heatmap.png'), height = 10, width = 20, units = 'cm', res = 400)
# ht
# graphics.off()

# # # Create target gene heatmap
# # target_genes_df <- extract_target_genes_df(factors, df.grn)
# # colnames(target_genes_df) <- paste0(colnames(target_genes_df), " - ", colSums((target_genes_df)))
# # hm <- pheatmap::pheatmap(t(target_genes_df),
# #                          color = c("grey", "purple"),  # Color scheme
# #                          cluster_rows = TRUE,  # Do not cluster rows
# #                          cluster_cols = TRUE,  # Do not cluster columns
# #                          fontsize_row = 10,  # Font size for row labels
# #                          fontsize_col = 0.0001,   # Font size for column labels
# #                          cutree_cols = k)
# # png(paste0(temp_plot_path_subset, 'Targets_heatmap.png'), height = 10, width = 18, units = 'cm', res = 400)
# # hm
# # graphics.off()

# # # GO analysis on each cluster of targets
# # df_row_cluster = data.frame(cluster = cutree(hm$tree_col, k = k))
# # for (i in 1:k){
# #   print(i)
# #   targets <- rownames(df_row_cluster %>% dplyr::filter(cluster == i))
# #   go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
# #   if (nrow(as.data.frame(go_output)) > 0){
# #     png(paste0(temp_plot_path_subset, 'Target_genes_cluster_', i, '_GO_plot.png'), height = 30, width = 20, units = 'cm', res = 400)
# #     print(plot(barplot(go_output, showCategory = 20)))
# #     graphics.off()
# #   }
# # }

# # # GO analysis of each TF's target genes
# # for (i in length(factors)){
# #   TF <- factors[i]
# #   targets <- rownames(target_genes_df)[as.logical(target_genes_df[,i])]
# #   go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
# #   if (nrow(as.data.frame(go_output)) > 0){
# #     png(paste0(temp_plot_path_subset, 'Target_genes_TF_', TF, '_GO_plot.png'), height = 30, width = 20, units = 'cm', res = 400)
# #     print(plot(barplot(go_output, showCategory = 20)))
# #     graphics.off()
# #   }
# # }

# # # subset network
# # df.grn.selected <- df.grn %>%
# #   dplyr::filter(tf %in% factors) %>%
# #   dplyr::filter(gene %in% factors)
# # nrow(df.grn.selected)

# # # TF network numbers
# # df <- data.frame(
# #   nSource = length(unique(df.grn.selected$tf)),
# #   nTarget = length(unique(df.grn.selected$gene)),
# #   nBoth = sum(unique(df.grn.selected$tf) %in% unique(df.grn.selected$gene)),
# #   nInteractions = nrow(df.grn.selected),
# #   nPositiveInteractions = length(which(df.grn.selected$correlation > 0)),
# #   nNegativeInteractions = length(which(df.grn.selected$correlation < 0))
# # )
# # png(paste0(temp_plot_path_subset, 'Subset_network_numbers.png'), height = 8, width = 18, units = 'cm', res = 400)
# # grid.arrange(top=textGrob("Network numbers", gp=gpar(fontsize=12, fontface = "bold"), hjust = 0.5, vjust = 3),
# #              tableGrob(df, rows=NULL, theme = ttheme_minimal()))
# # graphics.off()

# # # save network
# # write_tsv(df.grn.selected, file = paste0(temp_plot_path_subset, "GRN_subset.txt"))

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

# # # Pseudotime plots
# # for (TF in factors){
# #   print(TF)
# #   p <- PseudotimePlot_updated(object = obj.traj, tf.use = TF, trajectory.name = "lineage_placodal_probability")
# #   png(paste0(temp_plot_path_subset, 'Pseudotime_plot_', TF, '.png'), height = 8, width = 18, units = 'cm', res = 400)
# #   print(p)
# #   graphics.off()
# # }

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
#   print(TF)
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
# nrow(df.grn.selected) # 15,345 interactions

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

# ### plot the distribution of stages across the placodal trajectory
# trajectory_df <- data.frame(
#   traj = obj.traj$lineage_placodal_probability,
#   stage = obj.traj$stage
# )

# table(trajectory_df$stage)
# summary(trajectory_df$traj)

# df_HH5 <- trajectory_df %>% dplyr::filter(stage == "HH5")
# hist(df_HH5$traj)

# df_HH6 <- trajectory_df %>% dplyr::filter(stage == "HH6")
# hist(df_HH6$traj)

# df_HH7 <- trajectory_df %>% dplyr::filter(stage == "HH7")
# hist(df_HH7$traj)

# df_ss4 <- trajectory_df %>% dplyr::filter(stage == "ss4")
# hist(df_ss4$traj)

# df_ss8 <- trajectory_df %>% dplyr::filter(stage == "ss8")
# hist(df_ss8$traj, breaks = 100)

# stage_cols = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
# ggplot(trajectory_df, aes(x = traj, fill = stage)) + 
#   geom_histogram(binwidth=5) + 
#   facet_grid(stage ~ .) +
#   theme_minimal() +
#   scale_fill_manual(values = stage_cols) +
#   guides(fill = guide_legend(title = "Annotation")) +
#   theme(text = element_text(size = 25))


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

# # plot corr heatmap
# df.tf.gene.subset <- df.tf.gene %>%
#   dplyr::filter(tf %in% factors)
# df.tfs.subset <- df.tfs %>%
#   dplyr::filter(tfs %in% factors)
# ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = 1)

# png(paste0(temp_plot_path_subset, 'TF_gene_corr_heatmap.png'), height = 10, width = 20, units = 'cm', res = 400)
# ht
# graphics.off()

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
k = 12

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
factors <- importance_df[1:25, 1]
print("Top 20 importance factors:")
print(factors)

# plot expression of these factors
png(paste0(temp_plot_path_subset, 'Factors_expression_heatmap.png'), height = 10, width = 40, units = 'cm', res = 400)
DoHeatmap(object = seurat, features = factors, group.by = "scHelper_cell_type")
graphics.off()

# # Pseudotime plots of these factors
# obj.traj <- AddTargetAssay_updated(object = obj.traj, df.grn = df.grn.pos)
# for (TF in factors){
#   print(TF)
#   p <- PseudotimePlot_updated(object = obj.traj, tf.use = TF, trajectory.name = "lineage_placodal_probability")
#   png(paste0(temp_plot_path_subset, 'Factors_pseudotime_plot_', TF, '.png'), height = 8, width = 18, units = 'cm', res = 400)
#   print(p)
#   graphics.off()
# }

# ######## TEMP: Chromvar and footprinting of these factors
# # download motif database
# motifList <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))

# # rename each motif to have TF name
# name_vector <- c()
# for (i in 1:length(motifList)){
#   name <- name(motifList[[i]])
#   name_vector <- c(name_vector, name)
# }
# names(motifList) <- name_vector

# # annotate peaks in ArchR object with these motifs
# ArchR <- addMotifAnnotations(ArchR, name = "Motif", motifPWMs = motifList, cutOff = 1e-05, force = T)
# print("Motifs matrix added to ArchR object!")

# # chromvar plots for these factors
# ArchR <- addBgdPeaks(ArchR)
# ArchR <- addDeviationsMatrix(ArchR, peakAnnotation = "Motif", force = TRUE)
# ArchR <- addImputeWeights(ArchR)
# print("Making chromvar plots...")
# for (TF in factors){
#   print(TF)
#   markerMotif <- getFeatures(ArchR, select = TF, useMatrix = "MotifMatrix")
#   if(length(markerMotif) == 0){stop("Motif of that TF not found!")}
  
#   # Plot chromvar scores on UMAP
#   p <- plotEmbedding(ArchR, colorBy = "MotifMatrix", name = markerMotif, embedding = "UMAP", 
#                      imputeWeights = getImputeWeights(ArchR), plotAs = "points", size = 1.8,)
#   png(paste0(temp_plot_path_subset, 'Factors_chromvar_UMAP_Full_', TF, '.png'), height = 12, width = 10, units = 'cm', res = 400)
#   print(p)
#   graphics.off()
  
# }

# # set up for footprinting
# print("setting up for footprinting...")
# motifPositions <- getPositions(ArchR)
# ArchR <- addGroupCoverages(ArchRProj = ArchR, groupBy = "scHelper_cell_type_broad")

# # Footprinting of these factors
# print("Footprinting...")
# for (TF in factors){
#   print(TF)
#   seFoot <- getFootprints(ArchR, positions = motifPositions[TF], groupBy = "scHelper_cell_type_broad")
#   p <- plotFootprints(seFoot, names = TF, normMethod = "Subtract", plotName = "Footprints-Subtract-Bias",
#                       smoothWindow = 10, baseSize = 16, plot = FALSE)
#   png(paste0(temp_plot_path_subset, 'Factors_footprint_Full_', TF, '.png'), height = 20, width = 20, units = 'cm', res = 400)
#   grid::grid.newpage()
#   grid::grid.draw(p[[1]])
#   graphics.off()
# }
# ########

# plot corr heatmap of factors against all target nodes in the network
dir.create(paste0(temp_plot_path_subset, 'Target_gene_corr/'), recursive = T)

df.tf.gene.subset <- df.tf.gene %>%
  dplyr::filter(tf %in% factors)
df.tfs.subset <- df.tfs %>%
  dplyr::filter(tfs %in% factors)
ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = 1)

png(paste0(temp_plot_path_subset, 'Target_gene_corr/TF_gene_corr_heatmap.png'), height = 20, width = 20, units = 'cm', res = 400)
ht <- draw(ht)
graphics.off()

# # extract target gene clusters
# mat.cor <- df.tf.gene.subset %>% as.data.frame() %>% 
#   select(c(tf, gene, correlation)) %>% 
#   tidyr::pivot_wider(names_from = tf, values_from = correlation) %>% 
#   textshape::column_to_rownames("gene")
# # loop through each of these clusters
# target_gene_clusters <- list()
# for (cluster_name in names(row_order(ht))){
#   print(paste0("cluster: ", cluster_name))
#   indices <- row_order(ht)[[cluster_name]]
#   target_gene_clusters[[cluster_name]] <- rownames(mat.cor)[indices]
#   print(length(target_gene_clusters[[cluster_name]]))
#   #print expression of cluster of genes on featureplot
#   seurat <- AddModuleScore(object = seurat, features = list(target_gene_clusters[[cluster_name]]), name = paste0("cluster", cluster_name))
#   png(paste0(temp_plot_path_subset, 'Target_gene_corr/FeaturePlot_of_gene_cluster_from_corr_', cluster_name, '.png'), height = 10, width = 12, units = 'cm', res = 400)
#   print(FeaturePlot(seurat, features = paste0("cluster", cluster_name, "1"), pt.size = 1.5))
#   graphics.off()
# }
# export_gene_list(target_gene_clusters, publish_dir = paste0(temp_plot_path_subset, "Target_gene_corr/target_gene_clusters_from_corr"))

# Create target gene heatmap
dir.create(paste0(temp_plot_path_subset, 'Target_gene_direct_targets/'), recursive = T)
target_genes_df <- extract_target_genes_df(factors, df.grn.pos)
colnames(target_genes_df) <- paste0(colnames(target_genes_df), " - ", colSums((target_genes_df)))
hm <- pheatmap::pheatmap(t(target_genes_df),
                         color = c("grey", "purple"),  # Color scheme
                         cluster_rows = TRUE,  # Do not cluster rows
                         cluster_cols = TRUE,  # Do not cluster columns
                         fontsize_row = 10,  # Font size for row labels
                         fontsize_col = 0.0001,   # Font size for column labels
                         cutree_cols = k)
png(paste0(temp_plot_path_subset, 'Target_gene_direct_targets/Targets_heatmap.png'), height = 10, width = 18, units = 'cm', res = 400)
print(hm)
graphics.off()

# Extract each of the TF's target genes and plot their expression on a heatmap
order <- scHelper_cell_type_order[scHelper_cell_type_order %in% unique(seurat@meta.data$scHelper_cell_type)]
cols <- scHelper_cell_type_colours[names(scHelper_cell_type_colours) %in% order]

target_gene_direct <- list()
for (i in 1:length(factors)){
  print(i)
  TF <- factors[i]
  print(TF)
  targets <- rownames(target_genes_df)[as.logical(target_genes_df[,i])]
  target_gene_direct[[TF]] <- targets

  # GO analysis for each set of target genes
  go_output <- enrichGO(targets, OrgDb = org.Gg.eg.db, keyType = "SYMBOL", ont = "BP")
  if (nrow(as.data.frame(go_output)) > 0){
    png(paste0(temp_plot_path_subset, 'Target_gene_direct_targets/Target_genes_TF_', TF, '_GO_plot.png'), height = 30, width = 20, units = 'cm', res = 400)
    print(plot(barplot(go_output, showCategory = 20)))
    graphics.off()
  }

  # heatmap of gene expression for each set of target genes
  # height = length(targets)/7
  png(paste0(temp_plot_path_subset, 'Target_gene_direct_targets/Target_genes_heatmap_', TF, '.png'), height = 40, width = 40, units = 'cm', res = 400)
  print(TenxPheatmap(seurat, metadata = "scHelper_cell_type", selected_genes = targets,
             custom_order = order, custom_order_column = "scHelper_cell_type",
             gaps_col = "scHelper_cell_type",
             use_seurat_colours = FALSE, annotation_colours = cols,
             hclust_rows = TRUE, treeheight_row = 5, show_rownames = FALSE))
  graphics.off()
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

# plot corr heatmap of factors to test and known placodal factors
dir.create(paste0(temp_plot_path_subset, 'Selected_factors/'), recursive = T)

factors <- c("TFAP2E", "ZEB1", "TEAD3", "TBXT", "NKX2-3", "FOXK2", "EOMES", "TCF3",
             "TFAP2A", "TFAP2B", "TFAP2C", "GATA2", "GATA3", "DLX5", "DLX6", "SIX1")

k = 8

df.tf.gene.subset <- df.tf.gene %>%
  dplyr::filter(tf %in% factors)
df.tfs.subset <- df.tfs %>%
  dplyr::filter(tfs %in% factors)
ht <- GRNHeatmap(df.tf.gene.subset, tf.timepoint = df.tfs.subset$time_point, km = k)

png(paste0(temp_plot_path_subset, 'Selected_factors/TF_gene_corr_heatmap.png'), height = 20, width = 20, units = 'cm', res = 400)
ht <- draw(ht)
graphics.off()

# extract target gene clusters
mat.cor <- df.tf.gene.subset %>% as.data.frame() %>% 
  select(c(tf, gene, correlation)) %>% 
  tidyr::pivot_wider(names_from = tf, values_from = correlation) %>% 
  textshape::column_to_rownames("gene")
# loop through each of these clusters
target_gene_clusters <- list()
for (cluster_name in names(row_order(ht))){
  print(paste0("cluster: ", cluster_name))
  indices <- row_order(ht)[[cluster_name]]
  target_gene_clusters[[cluster_name]] <- rownames(mat.cor)[indices]
  print(length(target_gene_clusters[[cluster_name]]))
  #print expression of cluster of genes on featureplot
  seurat <- AddModuleScore(object = seurat, features = list(target_gene_clusters[[cluster_name]]), name = paste0("cluster", cluster_name))
  png(paste0(temp_plot_path_subset, 'Selected_factors/FeaturePlot_of_gene_cluster_from_corr_', cluster_name, '.png'), height = 10, width = 12, units = 'cm', res = 400)
  print(FeaturePlot(seurat, features = paste0("cluster", cluster_name, "1"), pt.size = 1.5))
  graphics.off()
}
export_gene_list(target_gene_clusters, publish_dir = paste0(temp_plot_path_subset, "Selected_factors/target_gene_clusters_from_corr"))

# plot target heatmap of selected factors + known factors
target_genes_df <- extract_target_genes_df(factors, df.grn.pos)
colnames(target_genes_df) <- paste0(colnames(target_genes_df), " - ", colSums((target_genes_df)))
hm <- pheatmap::pheatmap(t(target_genes_df),
                         color = c("grey", "purple"),  # Color scheme
                         cluster_rows = TRUE,  # Do not cluster rows
                         cluster_cols = TRUE,  # Do not cluster columns
                         fontsize_row = 10,  # Font size for row labels
                         fontsize_col = 0.0001,   # Font size for column labels
                         cutree_cols = k)
png(paste0(temp_plot_path_subset, 'Selected_factors/Targets_heatmap_plus_known_factors.png'), height = 10, width = 18, units = 'cm', res = 400)
print(hm)
graphics.off()

# plot target heatmap of just selected factors
factors <- c("TFAP2E", "ZEB1", "TEAD3", "TBXT", "NKX2-3", "FOXK2", "EOMES", "TCF3")
target_genes_df <- extract_target_genes_df(factors, df.grn.pos)
colnames(target_genes_df) <- paste0(colnames(target_genes_df), " - ", colSums((target_genes_df)))
hm <- pheatmap::pheatmap(t(target_genes_df),
                         color = c("grey", "purple"),  # Color scheme
                         cluster_rows = TRUE,  # Do not cluster rows
                         cluster_cols = TRUE,  # Do not cluster columns
                         fontsize_row = 10,  # Font size for row labels
                         fontsize_col = 0.0001,   # Font size for column labels
                         cutree_cols = k)
png(paste0(temp_plot_path_subset, 'Selected_factors/Targets_heatmap.png'), height = 10, width = 18, units = 'cm', res = 400)
print(hm)
graphics.off()