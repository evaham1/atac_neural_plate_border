#### EXTRACT MATRIX FROM ARCHR OBJECT (use getMarkerFeatures for now)

# # extract peak matrix, colnames = cell ids, rownames = peak idx
# peak_se <- getMatrixFromProject(ArchR, useMatrix = "PeakMatrix")
# peak_matrix <- as.matrix(assays(peak_se)$PeakMatrix)
# rownames(peak_matrix) <- rowData(peak_se)$idx
# 
# peak_matrix[1:2, 1:2]
# 
# # extract peaks of interest (work on)
# idxs = c(1:100)
# matrix <- peak_matrix[idxs, ]
# 
# matrix[1:2, 1:2]
# 
# # group cells and calculate average peak score per group
# matrix_t <- as.data.frame(t(matrix)) %>% mutate(group = ArchR$clusters)
# means <- aggregate(. ~ group, matrix_t, mean)
# means_matrix <- t(means) %>% set_colnames(means$group)
# means_matrix <- means_matrix[-1, ]
# means_matrix <- data.frame(apply(means_matrix, 2, function(x) as.numeric(as.character(x))))

# run differential test to get se object
seMarker <- getMarkerFeatures(
  ArchRProj = ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = opt$group_by)

### ADD UNIQUE IDs TO SE OBJECT SO CAN SUBSET IT (FOR PEAKS ONLY?)

add_unique_ids_to_se <- function(seMarker, ArchR) {
  tmp_peaks = data.frame(ArchR@peakSet)
  tmp_diff_peaks = data.frame(rowData(seMarker))
  diff_peaks_join_peakset = left_join(tmp_diff_peaks, tmp_peaks, 
                                      by = c("seqnames" = "seqnames", "start" = "start", "end" = "end"))
  diff_peaks_join_peakset$gene_name = paste(diff_peaks_join_peakset$nearestGene, diff_peaks_join_peakset$distToTSS,sep="_")
  diff_peaks_join_peakset$unique_id = paste0(diff_peaks_join_peakset$seqnames, ":", diff_peaks_join_peakset$start, "-", diff_peaks_join_peakset$end)
  
  rowData(seMarker) = diff_peaks_join_peakset
  return(seMarker)
}

seMarker_named <- add_unique_ids_to_se(seMarker, ArchR)

seMarker <- seMarker_named

### EXTRACT MEANS FROM SE OBJECT, rownames = unique_id

extract_means_from_se <- function(seMarker) {
  mat <- as.data.frame(SummarizedExperiment::assays(seMarker)[["Mean"]])
  rownames(mat) <- rowData(seMarker)$unique_id
  
  return(mat)
}

matrix <- extract_means_from_se(seMarker_named)

#### SCALING FUNCTION - feature dependent so must be done before subsetting

Log2norm <- function(mat, scaleTo = 10^4) {
  # normalising means for depth of cluster
  #scaleTo is just to increase values x1000
  mat <- log2(t(t(mat)/colSums(mat)) * scaleTo + 1)
  return(mat)
}

normalised_matrix <- Log2norm(matrix)

### GET LIST OF IDXS TO PLOT - either using cut off or top n per cell group
extract_ids <- function(seMarker, cutOff = "FDR <= 1 & Log2FC >= 0", top_n = TRUE, n = 10, group_name = "clusters") {
  
  markerList <- getMarkers(seMarker, cutOff = cutOff) # extract features that pass threshold
  
  df <- data.frame() # merged all features into a df
  for (i in 1:length(names(markerList))) {
    print(i)
    df_i <- as.data.frame(markerList[i])
    df <- rbind(df, df_i)
  }
  
  if (top_n == FALSE){
    ids <- df$unique_id
  } else {
    df <- df %>%
      group_by(group_name) %>%
      top_n(n, Log2FC) %>%
      dplyr::arrange(Log2FC, .by_group = TRUE)
    ids <- unique(df$unique_id)
  }
  
  return(ids)
}

ids <- extract_ids(seMarker, cutOff = "FDR <= 0.01 & Log2FC >= 5", top_n = FALSE)
length(ids)


#### USE IDXS TO SUBSET SCALED MATRIX

subset_matrix <- function(mat, ids) {
  subsetted_matrix <- mat[ids, ]
  return(subsetted_matrix)
}

mat <- subset_matrix(normalised_matrix, ids)

### HEATMAP
pal <- paletteContinuous(set = "viridis", n = 100)
pal <- viridis::magma(100)
pal <- paletteContinuous(set = "solarExtra", n = 100)
pal <- viridisLite::mako(256)
stage_order <- c("HH5", "HH6", "HH7", "ss4", "ss8")
stage_colours = c("#8DA0CB", "#66C2A5", "#A6D854", "#FFD92F", "#FC8D62")
names(stage_colours) <- stage_order

marker_heatmap <- function(mat, pal = NULL, 
                           labelRows = TRUE, clusterRows = TRUE, showRowDendrogram = TRUE, fontSizeRows = 12,
                           labelCols = TRUE, clusterCols = TRUE, showColDendrogram = TRUE, fontSizeCols = 12,
                           ) {
  
  # scale each feature independently and add min/max limits
  limits <- c(-2, 2) # could make this user-defined
  mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), 
               `/`)
  mat[mat > max(limits)] <- max(limits)
  mat[mat < min(limits)] <- min(limits)
  
  # colours - set default if NULL
  if (is.null(pal) == TRUE) {
    pal <- paletteContinuous(set = "solarExtra", n = 100)
  }

  # legend
  legend <- list(at = c(0, 1),
    labels = c(round(min(limits),2), round(max(limits),2)),
    color_bar = "continuous",
    legend_direction = "horizontal",
    legend_width = unit(3, "cm"),
    title = "Z-scores"
    )
  
  # # set top annotation (for now stage, can add scHelper later?)
  # colData <- data.frame(clusters = c("C7", "C5", "C6"), stage = c("ss8", "ss8", "ss4"), stage_cols = c("red", "red", "blue"))
  # topAnno <- HeatmapAnnotation(stage = anno_block(gp = gpar(fill = stage_colours),
  #                                                 labels = unique(colData$stage))
  # )
  
  
  # topAnno <- HeatmapAnnotation(
  #   df = colData,
  #   col = colorMap, 
  #   show_legend = showLegend,
  #   show_annotation_name = TRUE,
  #   gp = gpar(col = "NA"),
  #   annotation_legend_param =
  #     list(
  #       nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
  #     ),
  #   foo = anno_check_version_cols(at = customColLabel, labels = customColLabelIDs, labels_gp = gpar(fontsize = fontSizeLabels))
  # )
  # 
  # topAnno <- HeatmapAnnotation(stage = anno_block(gp = gpar(fill = plot_data$ann_colours$stage),
  #                                      labels = levels(plot_data$col_ann$stage),
  #                                      labels_gp = gpar(col = "white", fontsize = 50, fontface='bold'))),
  
  Heatmap(
    matrix = mat,
    col = pal,
    heatmap_legend_param = legend,
    top_annotation = topAnno, 
    # add raster stuff?
    
    #Column Options
    show_column_names = labelCols,
    cluster_columns = clusterCols,
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    column_names_gp = gpar(fontsize = fontSizeCols),
    column_names_max_height = unit(100, "mm"),
    column_split = colData$stage,
    
    #Row Options
    show_row_names = labelRows,
    cluster_rows = clusterRows,
    show_row_dend = showRowDendrogram,
    clustering_method_rows = "ward.D2",
    row_names_gp = gpar(fontsize = fontSizeRows)
    #row_split = row_split_params
  )
  
  return(Heatmap)
  
}
  

    
    
    



