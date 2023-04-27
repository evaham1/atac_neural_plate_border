##################### Distribution of labels across stages ##################################

if (length(unique(ArchR$stage)) > 1){

  plot_path = "./plots/labels_by_stage_distribution/"
  dir.create(plot_path, recursive = T)
  
  png(paste0(plot_path, 'counts_by_stage_table.png'), height = 25, width = 40, units = 'cm', res = 400)
  cell_counting(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "stage", scHelper_cell_type_order = scHelper_cell_type_order)
  graphics.off()
  
  # visualise distribution across stages: confusion matrix
  png(paste0(plot_path, "stage_distribution.png"), width=25, height=20, units = 'cm', res = 200)
  cell_counts_heatmap(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "stage")
  graphics.off()
  
  # visualise distribution across stages: piecharts
  counts <- cell_counting(ArchR = ArchR, group1 = "scHelper_cell_type_old", group2 = "stage", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
  png(paste0(plot_path, "label_by_stage_piecharts_unscaled.png"), width=50, height=40, units = 'cm', res = 200)
  cell_counts_piecharts(counts, col = scHelper_cell_type_colours)
  graphics.off()
  
  png(paste0(plot_path, "label_by_stage_piecharts_scaled.png"), width=50, height=40, units = 'cm', res = 200)
  cell_counts_piecharts(counts, col = scHelper_cell_type_colours, scale = TRUE)
  graphics.off()
  
##################### Distribution of stages across labels ##################################
  
  # visualise distribution across stages: piecharts
  counts <- cell_counting(ArchR = ArchR, group1 = "stage", group2 = "scHelper_cell_type_old", print_table = FALSE, scHelper_cell_type_order = scHelper_cell_type_order)
  png(paste0(plot_path, "stage_by_label_piecharts_unscaled.png"), width=50, height=40, units = 'cm', res = 200)
  cell_counts_piecharts(counts, col = stage_colours)
  graphics.off()
  
  png(paste0(plot_path, "stage_by_label_piecharts_scaled.png"), width=50, height=40, units = 'cm', res = 200)
  cell_counts_piecharts(counts, col = stage_colours, scale = TRUE)
  graphics.off()
  
##################### Distribution of rna stages across atac stages ##################################

  rna_stages <- plotEmbedding(ArchR, name = "rna_stage", plotAs = "points", size = 1.8, baseSize = 0, 
                            labelSize = 8, legendSize = 0, labelAsFactors = FALSE, pal = stage_colours)
  atac_stages <- plotEmbedding(ArchR, name = "stage", plotAs = "points", size = 1.8, baseSize = 0, 
                             labelSize = 8, legendSize = 0, labelAsFactors = FALSE, pal = stage_colours)
  png(paste0(plot_path, 'UMAPs_rna_stages_VS_atac_stages.png'), height = 20, width = 40, units = 'cm', res = 400)
  print(rna_stages + atac_stages)
  graphics.off()

  png(paste0(plot_path, "rna_atac_stage_distribution.png"), width=25, height=20, units = 'cm', res = 200)
  cell_counts_heatmap(ArchR = ArchR, group1 = "rna_stage", group2 = "stage")
  graphics.off()

}