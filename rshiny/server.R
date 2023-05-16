
server <- function(input, output, session){
  
  ####################################################################
  # Generate heatmap
  output$heatmap <- renderPlot(DimPlot(dat_list[[input$subset_featureplots]], group.by = input$group) +
                                 my_theme,
                               height = function() {session$clientData$output_dimplot_width * 0.8})
}