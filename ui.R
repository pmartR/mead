# KoverA Default Filtering Settings
kovera_A <- 0
kovera_k <- 0
source("./functions/helper_functions.R")

shinyUI(navbarPage(title = div(img(src = "Honey_Jar.png", height = 33, width = 22),
                               "mead"),
                   theme = "spacelab.css",
#---------------------------------------- Load Data Tab ----------------------------------------#                   
                   tabPanel("Load Data",
                            sidebarLayout(
                              sidebarPanel(
                                h3("Metadata"),
                                fileInput('qiime', 'Choose qiime sample data'),
                                tags$hr(),
                                h3("Data"),
                                fileInput('biom', 'Choose BIOM File'),
                                fileInput('fasta', 'Choose FASTA file'),
                                fileInput('tree', 'Choose TREE file'),
                                actionButton('go_button','Begin Analysis', icon = icon("flash")),
                                br(),
                                br(),
                                downloadButton('downloadOTUtable', "Download OTU table")
                              ),
                            mainPanel(
                                 plotOutput("library_sizes"),
                                 uiOutput("plots"),
                                 h4("Uploaded Metadata View"),
                                 dataTableOutput("sample_metadata")
                            )
                          )
                    ),
#---------------------------------------- kOverA FIltering Tab ----------------------------------------# 
                   tabPanel( "kOverA Filtering",
                             fluidPage(
                               fluidRow(
                              # Sidebar 
                                sidebarPanel(
                                  actionButton("actionb_filter", "Execute Filter", icon("filter")),
                                  h4('kOverA OTU Filtering'),
                                  fluidRow(
                                    column(width = 12,
                                               div(class = "col-md-6",
                                                   numericInputRow("filter_kOverA_count_threshold", "A: Count value minimum",
                                                                   value = kovera_A, min = 0, step = 1, class = "col-md-12")), 
                                               div(class = "col-md-6", uiOutput("filter_ui_kOverA_k"))
                                    )
                                  )
                                ),
                               # Main Panel
                                column(
                                 width = 8, offset = 0, 
                                 h4("Histograms Before and After Filtering"),
                                 plotOutput("filter_summary_plot"),
                                 h4("Data Summaries"),
                                 fluidRow(
                                   column(width = 6,
                                          p("Original"),
                                          htmlOutput('contents')
                                   ),
                                   column(width = 6,
                                          p("Filtered Data:"),
                                          htmlOutput('filtered_contents')
                                   )
                                  )
                                )
                               )
                             )
                   ),
#---------------------------------------- Alpha Diversity Tab ----------------------------------------# 
                   tabPanel("Alpha Diversity",
                                sidebarLayout(
                                  sidebarPanel(
                                    h3("Diversity"),
                                    uiOutput("rich_uix_color"),
                                    selectInput("dist_measures", label = "Measures", choices = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), multiple = TRUE),
                                    downloadButton('downloadDiversityImage', "Download plot"),
                                    h3("Richness"),
                                    checkboxGroupInput("split", "Split?", choices = c("TRUE", "FALSE"), selected = "TRUE"),
                                    downloadButton('downloadRichnessEstimates', "Download csv of richness estimates")
                                    
                                  ),
                                  mainPanel(plotOutput("richness_plot", width = "100%"))
                                )
                    ),
#---------------------------------------- Ordination Tab ----------------------------------------# 
                    tabPanel("Ordination",
                      sidebarLayout(
                        sidebarPanel(
                          downloadButton('downloadOrdinationImage', "Download plot")
                        ),
                      mainPanel(plotOutput("ordination_plot",  height = 500, width = 700))
                      )
                    )

))