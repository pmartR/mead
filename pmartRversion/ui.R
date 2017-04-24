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
  h3("Metadata Filtering"),
  uiOutput("plots"),
  hr(),
  actionButton("reset_button", label = "Reset Filter", icon = icon("trash"))
  
)) #end page