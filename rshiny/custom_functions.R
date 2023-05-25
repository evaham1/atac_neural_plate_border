
## Function to make shiny heatmap

plot_shiny_heatmap <- function(stage, cell_types, peaks){
  
  #peaks <- c(PPR_hichip_peaks, NC_hichip_peaks)
  #peaks <- c(PPR_peaks. shared_peaks)
  
  # init data from global data
  matrix <- SEACells_peak_matrix
  metadata <- as.data.frame(SEACells_metadata)
  metadata <- column_to_rownames(metadata, var = "ATAC")
  
  # subset metadata to only include SEACells that are in that chosen cell type
  if (!stage == "Full Data"){metadata <- metadata %>% filter(stage == !!stage)}
  metadata <- metadata %>% filter(scHelper_cell_type %in% !!cell_types)
  
  # subset matrix to only include cells in the metadata
  matrix <- matrix[which(rownames(matrix) %in% rownames(metadata)), ]
  
  # extract cols and order based on seacells subset
  order <- scHelper_cell_type_order[scHelper_cell_type_order %in% metadata$scHelper_cell_type]
  scHelper_cell_type_colours <- scHelper_cell_type_colours[order]
  
  # filter matrix by selected peaks
  peaks <- peaks[peaks %in% colnames(matrix)]
  matrix <- matrix[, which(colnames(matrix) %in% peaks)]
  
  # make heatmap
  plot_data <- PrepPeakModuleHeatmap(matrix, metadata, 
                                     col_order = c('stage', 'scHelper_cell_type'), custom_order_column = "scHelper_cell_type", custom_order = order, 
                                     hclust_SEACells = TRUE, hclust_SEACells_within_groups = TRUE,
                                     peak_modules = peaks, peak_row_annotation = FALSE,
                                     log_path = NULL)
  plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = TRUE,
                  show_column_names = FALSE, column_title = NULL, show_row_names = TRUE, row_title_gp = gpar(fontsize = 10), row_title_rot = 90,
                  column_split = plot_data$col_ann$stage,
                  bottom_annotation = CreateCellTypeAnnotation(plot_data, scHelper_cell_type_colours),
                  top_annotation = CreateStageAnnotation(plot_data, stage_colours),
                  col = PurpleAndYellow())
  
  # return the plot
  return(plot)
  
}

## Function to extract cell types avaliable for a given stage

check_cell_types <- function(stage){
  metadata <- SEACells_metadata
  if (!stage == "Full Data"){metadata <- metadata %>% filter(stage == !!stage)}
  cell_types <- unique(metadata$scHelper_cell_type)
  return(cell_types)
}




###################################################################################################################
########### functions to go into schelper package


PrepPeakModuleHeatmap <- function (peak_normalised_matrix, cell_metadata, 
                                   col_order, custom_order_column = NULL, custom_order = NULL, 
                                   hclust_SEACells = FALSE, hclust_SEACells_within_groups = TRUE,
                                   peak_modules, peak_row_annotation = TRUE,
                                   scale_data = TRUE,
                                   log_path = paste0(plot_path, "PrepPeakModuleHeatmap_logs/")) 
{
  
  ### Cell-level ordering and annotations ###
  
  # Initiate column anndata
  col_ann <- cell_metadata %>% mutate_if(is.character, as.factor)
  
  # If 'custom_order' is set use this to reorder cells
  if (!is.null(custom_order)) {
    if (!setequal(custom_order, unique(col_ann[[custom_order_column]]))) {
      stop("custom_order factors missing from custom_order_column \n\n")
    }
    col_ann[[custom_order_column]] <- factor(col_ann[[custom_order_column]], levels = custom_order)
    col_ann <- col_ann[order(col_ann[[custom_order_column]]),]
  }
  
  # If 'col_order' is use these columns to order cells
  if (!is.null(col_order)) {
    col_ann <- col_ann[do.call("order", c(col_ann[col_order], list(decreasing = FALSE))), , drop = FALSE]
  }
  
  # Optionally hclust SEACells
  if (hclust_SEACells == TRUE) {
    
    # Hclust SEACells within SEACell groups eg within scHelper_cell_type groups, which is specified by last col_order value
    if (hclust_SEACells_within_groups == TRUE) {
      cell_groups <- split(col_ann, col_ann[[tail(col_order, n=1)]])
      CellGroups_ordered_SEACells <- c()
      for (i in names(cell_groups)) {
        mat <- peak_normalised_matrix[rownames(cell_groups[[i]]), ]
        dist_mat <- dist(mat, method = "euclidean")
        hclust_avg <- hclust(dist_mat, method = "average")
        if (!is.null(log_path)){
          dir.create(log_path, recursive = T)
          png(paste0(log_path, i, '_SEACells_dendogram.png'),  width = 60, height = 40, units = 'cm', res = 400)
          plot(hclust_avg, main = paste0("SEACells dendogram for ", i))
          graphics.off()
        }
        ordered_SEACells <- hclust_avg$labels[c(hclust_avg$order)]
        CellGroups_ordered_SEACells[[i]] <- ordered_SEACells
      }
      col_ann <- col_ann[order(match(rownames(col_ann), unlist(CellGroups_ordered_SEACells))), , drop = FALSE]
      
      # Hclust SEACells across all SEACells, then if there is a col_order secondarily order by the last col_order value
    } else {
      dist_mat <- dist(peak_normalised_matrix, method = "euclidean")
      hclust_avg <- hclust(dist_mat, method = "average")
      if (!is.null(col_order)){ hclust_avg <- with(col_ann, reorder(hclust_avg, as.numeric(col_ann[[tail(col_order, n=1)]])))}
      if (!is.null(log_path)){
        dir.create(log_path, recursive = T)
        png(paste0(log_path, 'SEACells_dendogram.png'), width = 120, height = 40, units = 'cm', res = 400)
        plot(hclust_avg, main = "SEACells dendogram")
        graphics.off()
        }
    ordered_SEACells <- hclust_avg$labels[c(hclust_avg$order)]
    col_ann <- col_ann[order(match(rownames(col_ann), ordered_SEACells)), , drop = FALSE]
  }
  
  ### Peak-level ordering and annotations ###
  
  # Optionally annotate peaks by their modules
  if (peak_row_annotation == TRUE) {
    row_ann <- stack(peak_modules) %>% dplyr::rename(`Peak Modules` = ind) %>%
      column_to_rownames("values")
  } else {
    row_ann <- NA
  }
  
  ### Prepare data for plotting ###
  
  # Order matrix by row and column annotation orders
  plot_data <- t(peak_normalised_matrix)[unlist(peak_modules), rownames(col_ann)]
  
  # Optionally scale
  if (scale_data) {
    cat("Scaling data \n")
    plot_data <- t(scale(t(plot_data)))
    plot_data <- replace(plot_data, plot_data >= 2, 2)
    plot_data <- replace(plot_data, plot_data <= -2, -2)
  }
  
  ### Output plotting data and annotations ###
  
  output <- list(plot_data = plot_data,
                 row_ann = row_ann,
                 col_ann = col_ann)
  return(output)
  
  }
}

## Function to create bottom annotation for Complex Heatmap
# plot_data has to be generated by `PrepPeakModuleHeatmap` function and include $col_ann, scHelper_cell_type_colors should be named and ordered vector of colours
CreateCellTypeAnnotation <- function(plot_data, scHelper_cell_type_colors){
  return(
    HeatmapAnnotation(scHelper_cell_type = anno_simple(x = as.character(plot_data$col_ann$scHelper_cell_type),
                                                       col = scHelper_cell_type_colors, height = unit(0.5, "cm")), show_annotation_name = FALSE,
                      labels = anno_mark(at = cumsum(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths) - floor(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths/2),
                                         labels = rle(as.character(plot_data$col_ann$scHelper_cell_type))$values,
                                         which = "column", side = 'bottom',
                                         labels_gp = gpar(fontsize = 10), lines_gp = gpar(lwd=2)))
  )
}

## Function to create stage top annotation from ComplexHeatmap
# Needs plot_data generated from PrepPeakModuleHeatmap and named ordered vector of stage colours
CreateStageAnnotation <- function(plot_data, stage_colours){
  return(
    HeatmapAnnotation(stage = anno_block(gp = gpar(fill = stage_colours),
                                         labels = levels(plot_data$col_ann$stage),
                                         labels_gp = gpar(col = "white", fontsize = 20, fontface='bold')),
                      simple_anno_size = unit(1, "cm"),
                      annotation_label = "stage", gp = gpar(fontsize = 20))
  )
}
