# load and process tutorial data
inputFiles <- getTutorialData(tutorial = "hematopoiesis")
addArchRGenome("hg19")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, 
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  threads = 1
)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE,
  threads = 1
)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- saveArchRProject(ArchRProj = proj)

# add gene score matrix and calculate differential expression
proj <- addGeneScoreMatrix(proj, threads = 1, force = TRUE)
se <- getMarkerFeatures(proj,
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters", threads = 1)

#######   1: inputting the whole SE object and using 'cutoff' to filter   #########

plotMarkerHeatmap(
  seMarker = se, 
  cutOff = "FDR <= 0.01 & Log2FC >= 5",
  nLabel = 3)

 #######   2: inputting manually subsetted se object using same cutoff   #########

# extract df from se object
markerList <- getMarkers(se, cutOff = "FDR <= 0.01 & Log2FC >= 5")
df <- data.frame()
for (i in 1:length(names(markerList))) {
  print(i)
  df_i <- as.data.frame(markerList[i])
  df <- rbind(df, df_i)
}

# subsetting the se object
rowData(se)
subsetSE <- se[which(rowData(se)$name %in% df$name)]

# plot
plotMarkerHeatmap(
  seMarker = subsetSE, 
  cutOff = "FDR = FDR",
  nLabel = 3)