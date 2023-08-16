library(ArchR)

data_path = "./output/NF-downstream_analysis/Processing/FullData/remove_contam/rds_files/TransferLabels_Save-ArchR"

ArchR <- loadArchRProject(path = paste0(data_path), force = FALSE, showLogo = TRUE)

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



stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order



cols <- scHelper_cell_type_colours[unique(ArchR$scHelper_cell_type_broad)]
p1 <- plotEmbedding(ArchR, name = "scHelper_cell_type_broad", plotAs = "points", size = 2, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = cols, labelAsFactors = FALSE)

p2 <- plotEmbedding(ArchR, name = "stage", plotAs = "points", size = 2, baseSize = 0, 
              labelSize = 0, legendSize = 0, pal = stage_colours, labelAsFactors = FALSE)


p1
p2

ArchR <- addTrajectory(
  ArchRProj = ArchR, 
  name = "Trajectory", 
  groupBy = "stage",
  trajectory = stage_order, 
  embedding = "UMAP", 
  force = TRUE
)

p <- plotTrajectory(ArchR, trajectory = "Trajectory", colorBy = "cellColData", name = "Trajectory")
p[[1]]

############################## Create cell groupings for trajectory #######################################

###### add trajectory groupings
## HH5 and HH6, then split HH7, ss4 and ss8 by neural vs non-neural
stages <- as.data.frame(ArchR$stage)
cell_types <- as.data.frame(ArchR$scHelper_cell_type_broad)
metadata <- stages %>% mutate(cell_types = cell_types$`ArchR$scHelper_cell_type_broad`)
colnames(metadata) <- c("stage", "cell_type")
head(metadata)

## first merge placodal and non-neural
unique(metadata$cell_type)
metadata <- metadata %>% mutate(broad = cell_type)
metadata <- metadata %>% mutate(
  broad = case_when(
    broad == "Neural" ~ cell_type,
    broad == "Non-neural" ~ cell_type,
    broad == "NC" ~ cell_type,
    broad == "Placodal" ~ "Non-neural",
  )
)
head(metadata)
unique(metadata$broad)

## add new metadata called 'Order'
# all HH5 -> HH5
# all HH6 -> HH6
# HH7 split: non-neural and placodal cells = HH7_NN, neural and NC cells = HH7_Neural
# ss4 split: non-neural and placodal cells = ss4_NN, neural and NC cells = ss4_Neural
# ss8 split: non-neural and placodal cells = ss8_NN, neural and NC cells = ss8_Neural
metadata <- metadata %>% mutate(
  order = case_when(
    stage == "HH5" ~ stage,
    stage == "HH6" ~ stage,
    stage == "HH7" ~ paste0(stage, "_", broad),
    stage == "ss4" ~ paste0(stage, "_", broad),
    stage == "ss8" ~ paste0(stage, "_", broad),
  )
)

## add new cell groupings to seurat object
ArchR$order <- metadata$order
ArchR$broad <- metadata$broad

# plot order groupings
plotEmbedding(ArchR, name = "order", plotAs = "points", size = 2, baseSize = 0, 
                    labelSize = 0, legendSize = 0, labelAsFactors = FALSE)

plotEmbedding(ArchR, name = "broad", plotAs = "points", size = 2, baseSize = 0, 
                    labelSize = 0, legendSize = 0, labelAsFactors = FALSE)


unique(ArchR$order)


## add trajectories
# non-neural/placodal
ArchR <- addTrajectory(
  ArchRProj = ArchR, 
  name = "Trajectory", 
  groupBy = "order",
  trajectory = c("HH5", "HH6", "HH7_Non-neural", "ss4_Non-neural", "ss8_Non-neural"), 
  embedding = "UMAP", 
  force = TRUE
)

p <- plotTrajectory(ArchR)
p[[1]]

# neural
ArchR <- addTrajectory(
  ArchRProj = ArchR, 
  name = "Trajectory", 
  groupBy = "order",
  trajectory = c("HH5", "HH6", "HH7_Neural", "ss4_Neural", "ss8_Neural"), 
  embedding = "UMAP", 
  force = TRUE
)


p <- plotTrajectory(ArchR)
p[[1]]

# NC
ArchR <- addTrajectory(
  ArchRProj = ArchR, 
  name = "Trajectory", 
  groupBy = "order",
  trajectory = c("HH5", "HH6", "HH7_NC", "ss4_NC", "ss8_NC"), 
  embedding = "UMAP", 
  force = TRUE
)

p <- plotTrajectory(ArchR)
p[[1]]

ArchR$

