### will need to read in transfer_labels object

sc_Helper_cell_states_df <- data.frame(cell_ids = colnames(seurat_data),
                                       scHelper_cell_type = seurat_data@meta.data$scHelper_cell_type)


write.csv(sc_Helper_cell_states_df, "./TEST.csv")