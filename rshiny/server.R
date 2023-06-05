
server <- function(input, output, session){
  
  ####################################################################
  # Generate test UMAP
  output$dimplot_test <- renderPlot(DimPlot(SEACells_seurat, group.by = req(input$dimplot_groupby)) + my_theme, height = function() {session$clientData$output_dimplot_test_width * 0.8})
  
  
  # ####################################################################
  # Generate heatmap
  output$heatmap <- renderPlot(
    plot_shiny_heatmap(input$heatmap_stage, input$heatmap_celltype, input$heatmap_peaks), height = function() {session$clientData$output_heatmap_width * 0.8}
    )

  ## Observe stage input and use this to update which cell types can be selected to visualise
  observeEvent(input$heatmap_stage, {
    updateSelectInput(session, "heatmap_celltype", choices = check_cell_types(input$heatmap_stage), selected = "")
  })
  
  ## Observe what type of peak input and use this to update which peak options (IDs, PMs or genes) you can input
  observeEvent(input$how_to_choose_peaks, {
    if (input$how_to_choose_peaks == "Individual peaks") {
      updateSelectInput(session, "heatmap_peaks", choices = colnames(SEACells_peak_matrix), selected = "")
      output$peak_module_selection <- renderUI({NULL})
    }
    if (input$how_to_choose_peaks == "Peak modules") {
      output$peak_module_selection <- renderUI({
        radioButtons("peak_module_selection", "Select which set of peak modules:",
                     choices = c("FullData", "HH5", "HH6", "HH7", "ss4", "ss8"),
                     selected = "")
      })
    }
    if (input$how_to_choose_peaks == "From HiChip") {
      updateSelectInput(session, "heatmap_peaks", choices = c("PPR", "NC"), selected = "")
    }
  })
  
  observeEvent(input$peak_module_selection, {
    updateSelectInput(session, "heatmap_peaks", choices = rownames(PMs_df)[grep(input$peak_module_selection, rownames(PMs_df))], selected = "")
  })
  
  
  
}

