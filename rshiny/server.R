
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
  
}

