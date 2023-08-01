# library(getopt)
# library(optparse)
library(ArchR)
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
library(GenomicFeatures)
# library(parallel)
# library(scHelper)
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(scMEGA))
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Ggallus.UCSC.galGal6)
library(SummarizedExperiment)
suppressMessages(library(igraph))
suppressMessages(library(ggraph))

## set cols
###### schelper cell type colours
scHelper_cell_type_order <- c('EE', 'NNE', 'pEpi', 'PPR', 'aPPR', 'pPPR',
                              'eNPB', 'NPB', 'aNPB', 'pNPB','NC', 'dNC',
                              'eN', 'eCN', 'NP', 'pNP', 'HB', 'iNP', 'MB', 
                              'aNP', 'FB', 'vFB', 'node', 'streak', 
                              'PGC', 'BI', 'meso', 'endo')
scHelper_cell_type_colours <- c("#ed5e5f", "#A73C52", "#6B5F88", "#3780B3", "#3F918C", "#47A266", "#53A651", "#6D8470",
                                "#87638F", "#A5548D", "#C96555", "#ED761C", "#FF9508", "#FFC11A", "#FFEE2C", "#EBDA30",
                                "#CC9F2C", "#AD6428", "#BB614F", "#D77083", "#F37FB8", "#DA88B3", "#B990A6", "#b3b3b3",
                                "#786D73", "#581845", "#9792A3", "#BBB3CB")
names(scHelper_cell_type_colours) <- c('NNE', 'HB', 'eNPB', 'PPR', 'aPPR', 'streak',
                                       'pPPR', 'NPB', 'aNPB', 'pNPB','eCN', 'dNC',
                                       'eN', 'NC', 'NP', 'pNP', 'EE', 'iNP', 'MB', 
                                       'vFB', 'aNP', 'node', 'FB', 'pEpi',
                                       'PGC', 'BI', 'meso', 'endo')



# set paths
rds_path = "./output/NF-downstream_analysis/Processing/ss8/scMEGA/rds_files/"

## read in paired seurat object
obj.pair <- readRDS(paste0(rds_path, "paired_object.RDS"))
cols <- scHelper_cell_type_colours[as.character(unique(obj.pair$scHelper_cell_type))]

# UMAP
p1 <- DimPlot(obj.pair, group.by = "scHelper_cell_type", shuffle = TRUE, label = TRUE, reduction = "umap_harmony", cols = cols)
p1

## create trajectory - need to specify order of cell groupings
obj.pair <- AddTrajectory(object = obj.pair, 
                          trajectory = c("FB", "vFB", "aPPR", "pPPR"),
                          group.by = "scHelper_cell_type", 
                          reduction = "umap_harmony",
                          dims = 1:2, 
                          use.all = TRUE)

# we only incluce the cells that are in this trajectory
obj.pair <- obj.pair[, !is.na(obj.pair$Trajectory)]

p1 <- DimPlot(obj.pair, reduction = "umap_harmony", 
              group.by = "scHelper_cell_type",
              cols = cols)
p2 <- TrajectoryPlot(object = obj.pair, 
                     reduction = "umap_harmony",
                     continuousSet = "blueYellow",
                     size = 1,
                     addArrow = FALSE)
p1 + p2

## add motif information
# download motif database
motifList <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))

# rename each motif to have TF name
### THIS MAKES CHROMVAR FAIL!
# name_vector <- c()
# for (i in 1:length(motifList)){
#   name <- name(motifList[[i]])
#   name_vector <- c(name_vector, name)
# }
# names(motifList) <- name_vector
# head(names(motifList))

# add motif information to ATAC data
obj.pair <- AddMotifs(
  object = obj.pair,
  genome = BSgenome.Ggallus.UCSC.galGal6,
  pfm = motifList,
  assay = "ATAC"
)

# run chromvar
obj_chromvar <- RunChromVAR(
  object = obj.pair,
  genome = BSgenome.Ggallus.UCSC.galGal6,
  assay = 'ATAC'
)

saveRDS(obj_chromvar, paste0(rds_path, "paired_object_chromvar.RDS"), compress = FALSE)
obj_chromvar <- readRDS(paste0(rds_path, "paired_object_chromvar.RDS"))

# inspect chromvar results
DefaultAssay(obj_chromvar) <- 'chromvar'
p1 <- FeaturePlot(
  object = obj,
  features = "SIX1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1

# select the TFs that correlate in 'activity' and expression
res <- SelectTFs(object = obj_chromvar, trajectory.name = "Trajectory", return.heatmap = TRUE,
                 cor.cutoff = 0.1)

# plot TF activity dynamics across trajectory
df.cor <- res$tfs
ht <- res$heatmap
draw(ht)

## select genes that are associated with peaks
res <- SelectGenes(object = obj_chromvar,
                   labelTop1 = 0,
                   labelTop2 = 0)

# plot the dynamics of these genes across trajectory
df.p2g <- res$p2g
ht <- res$heatmap
draw(ht)

## GRN inference
tf.gene.cor <- GetTFGeneCorrelation(object = obj_chromvar, 
                                    tf.use = df.cor$tfs, 
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar", 
                                    gene.assay = "RNA",
                                    trajectory.name = "Trajectory")
## plot regulons? heatmap
ht <- GRNHeatmap(tf.gene.cor, 
                 tf.timepoint = df.cor$time_point)
ht

# associate genes with TFs to build network
motif.matching <- obj_chromvar@assays$ATAC@motifs@data
colnames(motif.matching) <- obj_chromvar@assays$ATAC@motifs@motif.names
motif.matching <- motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]

df.grn <- GetGRN(motif.matching = motif.matching, 
                 df.cor = tf.gene.cor, 
                 df.p2g = df.p2g)

## visualise network
# define colors for nodes representing TFs (i.e., regulators)
df.cor <- df.cor[order(df.cor$time_point), ]
tfs.timepoint <- df.cor$time_point
names(tfs.timepoint) <- df.cor$tfs

# plot the graph, here we can highlight some genes
df.grn2 <- df.grn %>%
  subset(correlation > 0.8) %>%
  dplyr::select(c(tf, gene, correlation)) %>%
  rename(weights = correlation)

unique(df.grn2$tf)

p <- GRNPlot(df.grn2, 
             tfs.use = c("SIX1", "DLX5", "DLX6", "GATA2", "GATA3", "SNAI2", "SOX13", "SOX2", "SOX21",
                         "TFAP2A", "TFAP2B", "TFAP2C", "ZIC1", "ZIC3"),
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 42, 
             plot.importance = TRUE,
             min.importance = 2,
             remove.isolated = FALSE)

options(repr.plot.height = 20, repr.plot.width = 20)

print(p)

