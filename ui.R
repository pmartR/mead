

shinyUI(navbarPage("preDE Analysis Tool",
                   tabPanel("Load Data",
                            sidebarLayout(
                              sidebarPanel(
                                fileInput('qiime', 'Choose qiime sample data'),
                                tags$hr(),
                                h3("Related Files"),
                                fileInput('biom', 'Choose BIOM File'),
                                fileInput('fasta', 'Choose FASTA file'),
                                fileInput('tree', 'Choose TREE file')
                              ),
                            mainPanel(
                              flowLayout(
#                                 plotOutput("lib_size_hist"),
#                                 br(),
                                downloadButton('downloadOTUtable', "Download OTU table")
                              )
                            )
                          )
                    ),
                   tabPanel("Alpha Diversity",
                                sidebarLayout(
                                  sidebarPanel(
                                    h3("Diversity"),
                                    selectInput("dist_measures", label = "Measures", choices = c("Observed", "Simpson", "Shannon", "InvSimpson"), multiple = TRUE),
                                    downloadButton('downloadDiversityImage', "Download plot"),
                                    h3("Richness"),
                                    checkboxGroupInput("split", "Split?", choices = c("TRUE", "FALSE"), selected = "TRUE"),
                                    downloadButton('downloadRichnessEstimates', "Download csv of richness estimates")
                                    
                                  ),
                                  mainPanel(plotOutput("richness_plot",  height = 500, width = 700))
                                )
                              ),
                    tabPanel("Ordination",
                      sidebarLayout(
                        sidebarPanel(
                          downloadButton('downloadOrdinationImage', "Download plot")
                        ),
                      mainPanel(plotOutput("ordination_plot",  height = 500, width = 700))
                      )
                    )

))