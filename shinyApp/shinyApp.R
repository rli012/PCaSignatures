
library(shiny)
library(shinydashboard)
#shinyUI(ui)



##### ui

header=dashboardHeader(title = 'PCa Transcriptomes')

sidebar=dashboardSidebar(
  sidebarMenu(
    tags$div(tags$h4("Gene-Level"),
             style = 'color: white'),
    menuItem("Differential Expression", tabName = 'degene', icon = icon("dna"),
             menuSubItem("Boxplot" , tabName = "boxplot", icon = icon("bar-chart")),
             menuSubItem("Heatmap" , tabName = "heatmap", icon = icon("scatter-chart"))
    ),
    menuItem("Survival Analysis" , tabName = "kmplot", icon = icon("pencil")),
    
    tags$div(tags$h4("Dataset-Level")),
    menuItem("Differential Expression", tabName = 'dedataset', icon = icon("database")),
    menuItem("Pathway Analysis", tabName = 'pathway', icon = icon("database")),
    menuItem("Co-expression Analysis", tabName = 'coexpression', icon = icon("chart-network"))
  )
)


body=dashboardBody()


ui <- dashboardPage(title='PCa Transcriptomes', skin = 'green', header, sidebar, body)


##### Server

server <- function(input, output) { }

##### help

gene.default <- 'INFG'



shinyApp(
  ui = ui,
  server = server
)

