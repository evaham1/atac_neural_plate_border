library(shinythemes)
tab_home          <- tabItem(tabName = "home",
                             fluidRow(
                               column(8, offset = 2,
                                      includeMarkdown("home.md")
                               )
                             )
)

tab_heatmaps <- tabItem(tabName = "heatmaps",
                             fluidRow(
                               column(12,
                                      radioButtons("subset_heatmap", "Select dataset to visualise", data_subsets, inline = TRUE, selected = 'Full data', width = '800')
                               )
                             ),
                               column(5,
                                      box(
                                        selectizeInput("peak_id", "Select Peak", choices = NULL, width = "250"),
                                        plotOutput("heatmap"),
                                        width = 12,
                                        style='height:35vw'
                                      )
                               )
                             )

ui <- dashboardPage(
  fullscreen = TRUE,
  header = dashboardHeader(
    title = dashboardBrand(
      title = "10x scATAC Neural Plate Border",
      href = "https://github.com/evaham1/atac_neural_plate_border"
    )
    # skin = "light",
    # status = "white",
    # border = TRUE,
    # sidebarIcon = icon("bars"),
    
  ),
  
  
  
  dashboardSidebar(
    # tags$style("@import url(https://use.fontawesome.com/releases/v5.7.2/css/all.css);"),
    
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon('home')),
      menuItem("Heatmaps", tabName = "heatmaps", icon = icon("border-none")),
      menuItem("Lineage Dynamics", tabName = "lineage_dynamics", icon = icon('chart-line')),
      menuItem("UMAP co-expression", tabName = "coexpression_umaps", icon = icon("braille")),
      menuItem("Differential expression", tabName = "dea", icon = icon("arrows-alt"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    tabItems(
      tab_home,
      tab_heatmaps,
      tab_lineage_dynamics,
      tab_coexpression_umaps,
      tab_dea
    )
  )
)
