# library(shinythemes)

tab_umap <- tabItem(tabName = "SEACell_UMAP",
                    fluidRow(
                      column(12,
                             radioButtons("dimplot_groupby", "How to colour UMAP", dimplot_groupby_options, inline = TRUE, selected = 'seurat_clusters', width = '800')
                      )
                    ),
                    fluidRow(
                      column(5,
                             box(
                               plotOutput("dimplot_test"),
                               width = 12,
                               style='height:35vw'
                             )
                      )
                    )
)

tab_heatmap <- tabItem(tabName = "SEACell_heatmaps",
                       fluidRow(
                         column(5, radioButtons("heatmap_stage", "Select stage to visualise", data_subsets, inline = TRUE, selected = 'Full Data', width = '800'))
                         ),
                       fluidRow(
                         column(12, selectInput("heatmap_celltype", "Select cell type to visualise", choices = NULL, multiple = TRUE, width = "250"))
                       ),
                       fluidRow(
                         column(5, radioButtons("how_to_choose_peaks", "How do you want to choose your peaks?", c("Individual peaks", "Peak modules", "From HiChip"), inline = TRUE, selected = 'Individual peaks', width = '800'))
                       ),
                       fluidRow(
                         column(5, uiOutput("peak_module_selection"))
                       ),
                       fluidRow(
                         column(12, selectInput("heatmap_peaks", "Select peaks/PMs to visualise", choices = NULL, multiple = TRUE, width = "250"))
                       ),
                       fluidRow(
                         column(12, plotOutput("heatmap", width = "1000", height = "3000"))
                         )
)



ui <- dashboardPage(
  #fullscreen = TRUE,
  header = dashboardHeader(
    title = dashboardBrand(
      title = "10x ATAC Neural Plate Border",
      href = "https://github.com/evaham1/atac_neural_plate_border"
    )

  ),


  dashboardSidebar(
    # tags$style("@import url(https://use.fontawesome.com/releases/v5.7.2/css/all.css);"),

    sidebarMenu(
      # menuItem("Home", tabName = "dashboard", icon = icon('home')),
      menuItem("Heatmaps", tabName = "SEACell_heatmaps", icon = icon("border-none")),
      # menuItem("Genome tracks", tabName = "SEACell_genome_tracks", icon = icon('chart-line')),
      # menuItem("UMAP", tabName = "SEACell_UMAP", icon = icon("braille")),
      menuItem("Test", tabName = "SEACell_UMAP", icon = icon("arrows-alt"))
    )
  ),

  dashboardBody(
    # tags$head(
    #   tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    # ),
    ## do these need to be in the same order as above??
    tabItems(
      tab_umap,
      # tab_test1,
      tab_heatmap
      # tab_umap,
      # tab_test
  )
)
)




