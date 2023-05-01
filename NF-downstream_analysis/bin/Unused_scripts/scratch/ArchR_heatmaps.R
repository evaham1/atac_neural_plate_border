ArchR <- loadArchRProject(path = "./ArchRSubset", force = TRUE, showLogo = TRUE)
markersPeaks <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = opt$group_by)
markersPeaks_untouched <- markersPeaks

# creating the unique ids in slot called 'name'
tmp_peaks = data.frame(ArchR@peakSet)
tmp_diff_peaks = data.frame(rowData(markersPeaks))
diff_peaks_join_peakset = left_join(tmp_diff_peaks, tmp_peaks, 
                                    by = c("seqnames" = "seqnames", "start" = "start", "end" = "end"))
diff_peaks_join_peakset$gene_name = paste(diff_peaks_join_peakset$nearestGene, diff_peaks_join_peakset$distToTSS,sep="_")
diff_peaks_join_peakset$unique_id = paste0(diff_peaks_join_peakset$seqnames, ":", diff_peaks_join_peakset$start, "-", diff_peaks_join_peakset$end)

# adding unique ids to the se object
rowData(markersPeaks) = diff_peaks_join_peakset


#############   1   ###############
### create matrix using cut offs in plotheatmap function
plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 5",
  nLabel = 3)

matrix_1 <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 5",
  nLabel = 3,
  returnMatrix = TRUE)

dim(matrix_1) # 26 x 3
head(matrix_1)
rownames(matrix_1)


#############   2   ###############
### create matrix by subsetting markersPeaks and plotting whole se object
# subsetting cannot be done with idx as they are not unique, need to create unique ids first

# extract df from markersPeaks and use to subset se object
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 5")
df <- data.frame()
for (i in 1:length(names(markerList))) {
  print(i)
  df_i <- as.data.frame(markerList[i])
  df <- rbind(df, df_i)
}
dim(df) # 26 x 26
length(unique(df$unique_id)) # 26
df$unique_id %in% rownames(matrix_1)

# subsetting the se object using unique ids
marker_subset = markersPeaks[rowData(markersPeaks)$unique_id %in% df$unique_id, ]
dim(marker_subset) # 26 x 3

# plot
plotMarkerHeatmap(
  seMarker = marker_subset, 
  cutOff = "FDR == FDR",
  nLabel = 3,
  log2Norm = FALSE)

matrix_2 <- plotMarkerHeatmap(
  seMarker = marker_subset, 
  cutOff = "FDR == FDR",
  nLabel = 3,
  returnMatrix = TRUE,
  log2Norm = FALSE)
dim(matrix_2)

########### comparing the two
sum(rownames(matrix_1) %in% rownames(matrix_2)) # 8475
sum(matrix_1[,1] %in% matrix_2[,1])
matrix_1 == matrix_2

head(matrix_1)
head(matrix_2)


matrix_1 <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 5",
  nLabel = 3,
  returnMatrix = TRUE,
  log2Norm = FALSE)

matrix_2 <- plotMarkerHeatmap(
  seMarker = marker_subset, 
  cutOff = "FDR == FDR",
  nLabel = 3,
  returnMatrix = TRUE,
  log2Norm = FALSE)



