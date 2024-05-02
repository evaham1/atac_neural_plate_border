
## Function to make shiny heatmap

plot_shiny_heatmap <- function(stage, cell_types, input_peaks){
  
  # init data from global data
  matrix <- SEACells_peak_matrix
  metadata <- as.data.frame(SEACells_metadata)
  metadata <- column_to_rownames(metadata, var = "ATAC")
  
  ######  SEACELLS ########
  # subset metadata to only include SEACells that are in that chosen cell type
  if (!stage == "Full Data"){metadata <- metadata %>% filter(stage == !!stage)}
  metadata <- metadata %>% filter(scHelper_cell_type %in% !!cell_types)
  
  # subset matrix to only include cells in the metadata
  matrix <- matrix[which(rownames(matrix) %in% rownames(metadata)), ]
  
  # extract cols and order based on seacells subset
  order <- scHelper_cell_type_order[scHelper_cell_type_order %in% metadata$scHelper_cell_type]
  scHelper_cell_type_colours <- scHelper_cell_type_colours[order]
  
  ######  Peaks ########
  # if input is PMs, extract peaks
  if (length(grep("PM", input_peaks)) == length(input_peaks)){
    peaks <- ExtractModuleContents(PMs_df, module_names = input_peaks)
    peak_row_annotation <- TRUE
    peak_modules <- list()
    for (PM in input_peaks){
      peak_modules[[PM]] <- ExtractModuleContents(PMs_df, module_names = PM)
    }
  } else{
    if (length(grep("chr", input_peaks)) == length(input_peaks)){
      peaks <- input_peaks
      peak_row_annotation <- FALSE
      peak_modules <- peaks
    } else {
      peaks <- hichip_peaks_list[[input_peaks]]
      peak_modules <- peaks
    }
    
  }
  
  # filter matrix by selected peaks
  peaks <- peaks[peaks %in% colnames(matrix)]
  matrix <- matrix[, which(colnames(matrix) %in% peaks)]
  
  ######  Heatmap ########
  # prep heatmap data
  plot_data <- PrepPeakModuleHeatmap(matrix, metadata, 
                                     col_order = c('stage', 'scHelper_cell_type'), custom_order_column = "scHelper_cell_type", custom_order = order, 
                                     hclust_SEACells = TRUE, hclust_SEACells_within_groups = TRUE,
                                     # peak annotation depends on if plotting individual peaks or PMs
                                     peak_modules = peak_modules, peak_row_annotation = peak_row_annotation,
                                     log_path = NULL)
  
  # if input is PMs, different plotting params
  if (length(grep("chr", input_peaks)) != length(input_peaks)){
    cluster_rows <- FALSE
    show_row_names <- FALSE
    row_split <- plot_data$row_ann$`Peak Modules`
  } else{
    cluster_rows <- TRUE
    show_row_names <- TRUE
    row_split <- NULL
  }
  
  # plot heatmap
  plot <- Heatmap(plot_data$plot_data, cluster_columns = FALSE, cluster_rows = cluster_rows,
                  show_column_names = FALSE, column_title = NULL, 
                  show_row_names = show_row_names, row_title_gp = gpar(fontsize = 10), row_title_rot = 0,
                  column_split = plot_data$col_ann$stage, row_split = row_split,
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