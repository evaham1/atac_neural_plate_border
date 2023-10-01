#!/usr/bin/env Rscript

print("Run scMEGA on paired RNA/ATAC object to infer GRN")

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
    
    # full data
    data_path = "./output/NF-downstream_analysis/Processing/FullData/scMEGA/MEGA_cell_pairing_and_chromvar/rds_files/" # paired object with chromvar run
    
    
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


############################## Read in seurat object #######################################

print("reading in data...")

## read in paired seurat object
obj.pair <- readRDS(paste0(data_path, "paired_object_chromvar.RDS"))
obj.pair

#### TEMP: Check lineage trajectories are there

head(obj.pair@meta.data)

# # and then for now subsample data so can run quicker
# test <- subset(obj.pair, downsample = 100)

######################################################################################
##############################    PLACODAL     #######################################
######################################################################################

############################## Read in csv files #######################################

temp_plot_path = "./plots/placodal_lineage/"
dir.create(temp_plot_path, recursive = T)

temp_csv_path = "./csv_files/placodal_lineage/"

df.p2g <- read.csv(paste0(temp_csv_path, "Variable_genes_and_matched_enhancers.csv"))
df.tfs <- read.csv(paste0(temp_csv_path, "TF_correlations.csv"))
df.tf.gene <- read.csv(paste0(temp_csv_path, "TF_to_gene_correlations.csv"))
df.grn <- read.csv(paste0(temp_csv_path, "GRN_data.csv"))

############################## Plot GRN maps using scMEGA #######################################

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

png(paste0(temp_plot_path, 'Network_all.png'), height = 30, width = 45, units = 'cm', res = 400)
print(p)
graphics.off()


# only plot TFs
p <- GRNPlot(df.grn,
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             plot.importance = TRUE,
             genes.use = df.tfs$tfs,
             remove.isolated = TRUE)

png(paste0(temp_plot_path, 'Network_TFs.png'), height = 30, width = 45, units = 'cm', res = 400)
print(p)
graphics.off()

############################## Plot whole GRN map using igraphs #######################################

# filter grn to only include TFs
df.grn.tfs <- df.grn[df.grn$gene %in% df.grn$tf, ]
dim(df.grn.tfs)

# filter grn to only include positive interactions
df.grn.tfs <- df.grn.tfs %>% dplyr::filter(correlation > 0.1)
dim(df.grn.tfs)

# extract data to plot
links <- df.grn.tfs %>% dplyr::select(c("tf", "gene", "correlation"))

# timepoint expression information of TFs
nodes <- rownames_to_column(as.data.frame(tfs.timepoint), var = "tf")
length(unique(nodes$tf))
nodes <- nodes %>% 
  dplyr::filter(tf %in% unique(links$tf)) %>% # remove TFs with no links
  mutate(tfs.timepoint = round(tfs.timepoint, digits = 0))

# colours for latent time
cols <- data.frame(col = paletteContinuous(set = "beach", n = 100),
                   val = seq(1:100))
cols_matched <- merge(cols, nodes, by.x = "val", by.y = "tfs.timepoint")

# Turn it into igraph object
network <- graph_from_data_frame(d = links, vertices = nodes, directed = T) 

# Plot network - its different everytime, cant seem to fix randomness
png(paste0(temp_plot_path, 'Network_TFs_igraph.png'), height = 30, width = 45, units = 'cm', res = 400)
plot(network, edge.arrow.size=0.5, edge.width = E(network)$correlation*5,
          vertex.color = cols_matched$col,
          vertex.label.color = "black", vertex.label.family = "Helvetica", vertex.label.cex = 0.8,
     layout = layout_nicely(network))
graphics.off()

## annoyingly there doesnt seem to be an easy way to prevent node overlap
# plot(network, edge.arrow.size=0.2, edge.width = E(network)$correlation*5,
#      vertex.color = cols_matched$col,
#      vertex.label.color = "black", vertex.label.family = "Helvetica", vertex.label.cex = 0.8,
#      layout = layout_as_tree(network))

############################## Plot known GRN map using igraphs #######################################

# Plot network of just known TFs
df.grn.known.tfs <- df.grn.tfs %>% 
  filter(tf %in% c("ZNF462", "DLX5", "DLX6", "GATA2", "GATA3", "SIX1", "EYA2", "SIX4",
                   "FOXI1", "FOXI3", "TFAP2A", "ZIC1", "ZIC5", "PRDM1", "PAX3", "PAX7",
                   "MSX1", "AXUD1", "FOXD3", "SNAI2"))

# extract data to plot
links <- df.grn.known.tfs %>% dplyr::select(c("tf", "gene", "correlation"))

# timepoint expression information of TFs
nodes <- rownames_to_column(as.data.frame(tfs.timepoint), var = "tf")
length(unique(nodes$tf))
nodes <- nodes %>% 
  dplyr::filter(tf %in% c(links$tf, links$gene)) %>% # remove TFs with no links
  mutate(tfs.timepoint = round(tfs.timepoint, digits = 0))

# colours for latent time
cols <- data.frame(col = paletteContinuous(set = "beach", n = 100),
                   val = seq(1:100))
cols_matched <- merge(cols, nodes, by.x = "val", by.y = "tfs.timepoint")

# Turn it into igraph object
network <- graph_from_data_frame(d = links, vertices = nodes, directed = T) 

# Plot network - its different everytime, cant seem to fix randomness
png(paste0(temp_plot_path, 'Network_known_TFs_igraph.png'), height = 30, width = 45, units = 'cm', res = 400)
plot(network, edge.arrow.size=1, edge.width = E(network)$correlation*5,
     vertex.color = cols_matched$col,
     vertex.label.color = "black", vertex.label.family = "Helvetica", vertex.label.cex = 0.8,
     layout = layout_nicely(network))
graphics.off()

############################## Plot GRN map with high centrality using igraphs #######################################

# TEAD3     TEAD3   6.083333
# JUN         JUN   7.666667
# ATF6       ATF6  25.000000
# GATA3     GATA3  28.416667
# REL         REL  41.333333
# KLF6       KLF6  41.416667
# ZBTB7A   ZBTB7A  53.750000
# KLF11     KLF11  95.833333
# RFX2       RFX2  99.500000
# PBX1       PBX1 106.000000

# Plot network of just known TFs
df.grn.known.tfs <- df.grn.tfs %>% 
  filter(tf %in% c("TEAD3", "JUN", "ATF6", "GATA3", "REL", "KLF6", "ZBTB7A", "KLF11",
                   "RFX2", "PBX1"))

# extract data to plot
links <- df.grn.known.tfs %>% dplyr::select(c("tf", "gene", "correlation"))

# timepoint expression information of TFs
nodes <- rownames_to_column(as.data.frame(tfs.timepoint), var = "tf")
length(unique(nodes$tf))
nodes <- nodes %>% 
  dplyr::filter(tf %in% c(links$tf, links$gene)) %>% # remove TFs with no links
  mutate(tfs.timepoint = round(tfs.timepoint, digits = 0))

# colours for latent time
cols <- data.frame(col = paletteContinuous(set = "beach", n = 100),
                   val = seq(1:100))
cols_matched <- merge(cols, nodes, by.x = "val", by.y = "tfs.timepoint")

# Turn it into igraph object
network <- graph_from_data_frame(d = links, vertices = nodes, directed = T) 

# Plot network - its different everytime, cant seem to fix randomness
png(paste0(temp_plot_path, 'Network_high_betweeness_TFs_igraph.png'), height = 30, width = 45, units = 'cm', res = 400)
plot(network, edge.arrow.size=1, edge.width = E(network)$correlation*5,
     vertex.color = cols_matched$col,
     vertex.label.color = "black", vertex.label.family = "Helvetica", vertex.label.cex = 0.8,
     layout = layout_nicely(network))
graphics.off()


############################## Network analysis #######################################

# degree: how many connections
# eigenvector centrality: vertices that are in turn connected to other vertices that are highly interconnected in the network
# betweenness: how frequently a vertex lies on the shortest path(s) between any two vertices in the network.

# local degree centrality, ie the number of links held by a node
degree.cent <- centr_degree(network, mode = "all")
results <- data.frame(tf = nodes$tf,
                      res = degree.cent$res)
results <- results[order(results$res), ]
print(results)

# closeness centrality, closeness to all other nodes in the network
closeness.cent <- closeness(network, mode="all")
results <- data.frame(tf = nodes$tf,
                      res = closeness.cent)
results <- results[order(results$res), ]
print(results)

# betweeness, how frequently the vertex lies on shortest path between ay two vertices in network
# ie how critical this node is to flow of information through the network
betweeness <- betweenness(network, directed = "TRUE")
results <- data.frame(tf = nodes$tf,
                      res = betweeness)
results <- results[order(results$res), ]
print(results)

# grouping factors
x <- fastgreedy.community(network)


############################## Plot lineage dynamics #######################################

obj.pair <- AddTargetAssay_updated(test, df.grn = df.grn)

# plot dynamics of each TF in network
for (TF in df.tfs$tfs){
  print(TF)
  
  png(paste0(temp_plot_path, TF, '_dynamics_plot.png'), height = 15, width = 20, units = 'cm', res = 400)
  PseudotimePlot_updated(obj.pair, trajectory.name = "rna_lineage_placodal_probability",
                         tf.use = TF)
  graphics.off()
}



############################## Interactive visualisation #######################################

# Libraries
library(tidyverse)
library(viridis)
library(patchwork)
library(ggraph)
library(igraph)
library(networkD3)

# Plot
simpleNetwork(links,     
              Source = 1,                 # column number of source
              Target = 2,                 # column number of target
              height = 880,               # height of frame area in pixels
              width = 1980,
              linkDistance = 50,         # distance between node. Increase this value to have more space between nodes
              charge = -4,              # numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value)
              fontSize = 5,              # size of the node names
              fontFamily = "serif",       # font og node names
              linkColour = "#666",        # colour of edges, MUST be a common colour for the whole graph
              nodeColour = "#69b3a2",     # colour of nodes, MUST be a common colour for the whole graph
              opacity = 0.9,              # opacity of nodes. 0=transparent. 1=no transparency
              zoom = T                    # Can you zoom on the figure?
)


############################## Playing around with ggraph #######################################

package <- data.frame(
  name = c("igraph", "ggraph", "dplyr", "ggplot", "tidygraph")
)

tie <- data.frame(
  from = c("igraph", "igraph", "ggplot", "igraph", "dplyr", "ggraph"),
  to =   c("ggraph", "tidygraph", "ggraph", "tidygraph", "tidygraph", "tidygraph")
)

g = graph_from_data_frame(links, directed = TRUE, vertices = nodes)

# use arrows for directions
ggraph(g, layout = 'graphopt') + 
  geom_edge_link(aes(start_cap = label_rect(node1.name), end_cap = label_rect(node2.name)), 
                 arrow = arrow(type = "closed", length = unit(3, 'mm'))) + 
  geom_node_text(aes(label = name)) +
  geom_node_point(aes(fill = tfs.timepoint, size = 20),shape = 21) +
  theme_graph()

 
    
ggraph(g, layout = 'graphopt') + 
  geom_edge_link(aes(edge_width = correlation),edge_colour = "grey66") +
    geom_node_text(aes(label = name)) +
    geom_node_point(aes(colour = tfs.timepoint), size = 4)
  theme_graph()
  
  got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", 
                   "#50594B", "#8968CD", "#9ACD32", "#feb24c")
  size = 2
  ggraph(g, layout = "stress")+
    geom_edge_link(aes(edge_width = correlation),edge_colour = "grey66",
                    arrow = arrow(type = "open", length = unit(4, 'mm')),
                    start_cap = label_rect(node1.name), end_cap = label_rect(node2.name))+
    geom_node_text(aes(label = name)) +
    geom_node_point(aes(fill = tfs.timepoint, size = 2),shape = 21)+
    geom_node_text(aes(filter = size >= 26, label = name),family="serif")+
    scale_fill_continuous()+
    scale_edge_width(range = c(0.2,3))+
    scale_size(range = c(1,8))+
    theme_graph()+
    theme(legend.position = "right")
  