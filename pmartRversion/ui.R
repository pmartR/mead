kovera_k <- 0
kovera_A <- 0
source("./functions/helper_functions.R")

shinyUI(fluidPage(
  title = (windowTitle = "mead"),
  titlePanel(div(img(src = "Honey_Jar.png", height = 33, width = 22), "mead")),
  
  fluidRow(
    column(width = 10, offset = 1, 
           tags$hr(),
           h3("Load Data"),
           fluidRow(
             column(width = 4, fileInput('biom', 'Choose BIOM File')),
             column(width = 4, fileInput('fasta', 'Choose FASTA file')),
             column(width = 4, fileInput('qiime', 'Choose Sample Metadata QIIME'))
           )
    )
  ),
  h3("Full Metadata"),
  tags$table(
    DT::dataTableOutput("sample_metadata"),
    tags$head(
      tags$title(h4("Uploaded Metadata View"))
    )),
  h3("Sample Filtering"),
  h4("Keep samples with specific metadata"),
  uiOutput("plots"),
  br(),
  h4("Keep samples above a minimum number of reads"),
  fluidRow(
    column(width = 5, h5("Keep samples with at least")),
    column(width = 2, numericInputRow(inputId = "n", label = "", value = 0, min = 0, max = , width = 150)),
    column(width = 2, h5("reads"))
  ),
  plotOutput("sample_counts_plot"),
  hr(),
  actionButton("reset_button", label = "Reset Filter", icon = icon("trash")),
  h3("OTU Filtering"),
  fluidRow(
    column(width = 3, h5("Keep OTUs with more than"), align = 'center'),
    column(width = 2,  numericInputRow("filter_kOverA_count_threshold", "",
                                       value = kovera_A, min = 0, step = 1)),
    column(width = 3, h5("counts in at least")),
    column(width = 2, uiOutput("filter_ui_kOverA_k")),
    column(width = 2, h5("samples"))
  ),
  plotOutput("read_counts_plot")
)) #end page