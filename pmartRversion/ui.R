kovera_k <- 0
kovera_A <- 0
source("./functions/helper_functions.R")

shinyUI(navbarPage(
  title = (windowTitle = "mead"),
 # titlePanel(div(img(src = "Honey_Jar.png", height = 33, width = 22), "mead")),
  
  tabPanel("Data and Filtering", 
  fluidRow(
    column(width = 12, 
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
  fluidRow(
    actionButton("metadata_reset_button", label = "Reset Meta Filter", icon = icon("trash")),
    actionButton("metadata_filter_go", label = "Apply Meta Filter", icon = icon("bar-chart"))
  ),
  br(),
  h4("Keep samples above a minimum number of reads"),
  fluidRow(
    column(width = 5, h5("Keep samples with at least")),
    column(width = 2, numericInput(inputId = "n", label = "", value = 0, min = 0, width = 150)),
    column(width = 2, h5("reads"))
  ),
  plotOutput("sample_counts_plot"),
  fluidRow(
    actionButton("sample_reset_button", label = "Reset Sample Filter", icon = icon("trash")),
    actionButton("sample_filter_go", label = "Apply Sample Filter", icon = icon("bar-chart"))
  ),
  hr(),
  h3("OTU Filtering"),
  fluidRow(
    column(width = 3, h5("Keep OTUs with more than"), align = 'center'),
    column(width = 2,  numericInputRow("filter_kOverA_count_threshold", "",
                                       value = kovera_A, min = 0, step = 1)),
    column(width = 3, h5("counts in at least")),
    column(width = 2, uiOutput("filter_ui_kOverA_k")),
    column(width = 2, h5("samples"))
  ),
  plotOutput("read_counts_plot"),
  fluidRow(
    actionButton("otu_reset_button", label = "Reset OTU Filter", icon = icon("trash")),
    actionButton("otu_filter_go", label = "Apply OTU Filter", icon = icon("bar-chart"))
  )
  ),
  
  # tabPanel("Group Designation",
  #          #sidebarLayout(
  #            sidebarPanel(
  #              h3("Groups"),
  #              uiOutput("mainEffects"),
  #              fluidRow(
  #                actionButton("groupDF_reset_button", label = "Reset Groupings", icon = icon("trash")),
  #                actionButton("groupDF_go", label = "Apply Groupings", icon = icon("check"))
  #              )
  #           # )
  # 
  #          ),
  #          column(width=8, DT::dataTableOutput("group_DF"))
  # ),
  # 
  # tabPanel("Normalization",
  #          br()),
  
 tabPanel("Community Metrics",
          h3("Groupings"),
          uiOutput("group1"),
          uiOutput("group2"),
          # fluidRow(
          #   actionButton("groupDF_reset_button", label = "Reset Groupings", icon = icon("trash")),
          #   actionButton("groupDF_go", label = "Apply Groupings", icon = icon("check"))
          # ),
          DT::dataTableOutput("group_DF"),
          h3("Plot Parameters"),
          uiOutput("xaxis"),
          uiOutput("color"),
          h3("Alpha Diversity"),
          uiOutput("adiv_index"),
          plotOutput("adiv_plot"),
          br(),
          h3("Richness"),
          uiOutput("rich_index"),
          plotOutput("rich_plot"))
  
)) #end page