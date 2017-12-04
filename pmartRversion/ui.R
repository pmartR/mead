# Get dependencies imported if any are missing
#list of all packages
packs <- installed.packages()[,"Package"]
#list of required packages
dependencies <- c("shiny", "shinyjs", "lazyeval", "V8", "dplyr", "DT", "ggplot2", "vegan", "goeveg", "plotly")

#check for missing packages
missing <- dependencies[!(dependencies %in% packs)]

#install missing ones
if (length(missing) > 0) install.packages(missing, dependencies = TRUE)

# add in custom packages
if (!("filterWidget" %in% packs)) devtools::install("filterWidget")
if (!("pmartRseq" %in% packs)) devtools::install("../../pmartRseq/")

# source packages and functions
library(shiny)
library(shinyjs) 
library(V8)
library(lazyeval)
library(dplyr)
library(DT)
library(filterWidget)
library(ggplot2)
library(pmartRseq)
#library(phyloseq)
library(vegan)
library(goeveg)
library(plotly)
source("./functions/helper_functions.R")
source("./functions/test_functions.R")

kovera_k <- 0
kovera_A <- 0
source("./functions/helper_functions.R")
# Get dependencies imported if any are missing
#list of all packages
packs <- installed.packages()[,"Package"]
#list of required packages
dependencies <- c("shiny", "shinyjs", "lazyeval", "V8", "dplyr", "DT", "ggplot2", "vegan", "goeveg", "plotly", "DESeq2")

#check for missing packages
missing <- dependencies[!(dependencies %in% packs)]

#install missing ones
if (length(missing) > 0) install.packages(missing, dependencies = TRUE)

# add in custom packages and installs
if (!("filterWidget" %in% packs)) devtools::install("filterWidget")
if (!("pmartRseq" %in% packs)) devtools::install("../../pmartRseq/")
if (!("DESeq2" %in% packs)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
} 

# source packages and functions
library(shiny)
library(shinyjs) 
library(V8)
library(lazyeval)
library(dplyr)
library(DT)
library(filterWidget)
library(ggplot2)
library(pmartRseq)
#library(phyloseq)
library(vegan)
library(goeveg)
library(plotly)
library(DESeq2)
source("./functions/helper_functions.R")
source("./functions/test_functions.R")
# create a reset session button using a js method
jsResetCode <- "shinyjs.reset = function() {history.go(0)}"
shinyUI(navbarPage(title = (windowTitle = "mead"),




  
 # titlePanel(div(img(src = "Honey_Jar.png", height = 33, width = 22), "mead")),
 

 tabPanel("Data and Filtering", 
          tags$head(
          ),
          fluidRow(column(width=12,
                          useShinyjs(),                                           
                          extendShinyjs(text = jsResetCode),                     
                          actionButton("reset_button", "Reset All") 
                          )
            
          ),

          fluidRow(
            column(width = 12,
                   tags$hr(),
                   h3("Load Data"),
                   fluidRow(
                     column(width = 6, fileInput('biom', 'Choose BIOM File')),
                     # column(width = 4, fileInput('fasta', 'Choose FASTA file')),
                     column(width = 6, fileInput('qiime', 'Choose Sample Metadata QIIME'))
                   )
            )
          ),
          #actionButton("Upload", "Upload ze data!"),
  h3("Metadata View"),
  fluidRow(
    column(width = 12, tags$table(
      DT::dataTableOutput("sample_metadata"),
      tags$head(
        tags$title(h4("Uploaded Metadata View"))
      )))
  ),
  hr(),
  h3("Filtering"),
  h4("Sample Metadata"),
  p("Filter samples on numeric ranges or groupable characteristics"),
  br(),
  fluidRow(
    column(width = 5,  uiOutput("boxes")),
    column(width = 6, uiOutput("plots"))
    ),
  fluidRow(
    actionButton("metadata_reset_button", label = "Reset Meta Filter", icon = icon("trash")),
    actionButton("metadata_filter_go", label = "Apply Meta Filter", icon = icon("bar-chart"))
  ),
  br(),
  hr(),
  # h3("Taxonomic Rollup"),
  # p("Roll up from OTU-level data to a specified taxonomic level"),
  # uiOutput("taxa_level"),
  # br(),
  # hr(),
  fluidRow(
    column(width = 4, 
           h4("Sample Reads"),
           p("Keep samples above a minimum number of reads"),
           fluidRow(
             column(width = 5, h5("Keep samples with at least")),
             column(width = 2, numericInput(inputId = "n", label = "", value = 0, min = 0, width = 150)),
             column(width = 2, h5("reads"))
           ),
           plotOutput("sample_counts_plot"),
           fluidRow(
             actionButton("sample_reset_button", label = "Reset Sample Filter", icon = icon("trash")),
             actionButton("sample_filter_go", label = "Apply Sample Filter", icon = icon("bar-chart"))
           )
           ),
    column(width = 4,
           verbatimTextOutput("summ_filt")
    ),
    column(width = 4, 
           h4("OTU Counts Per Sample"),
           p("Keep OTUs with a minimum number of counts per sample"),
           fluidRow(
             column(width = 3, h5("Keep OTUs with more than"), align = 'center'),
             column(width = 2,  numericInputRow("filter_kOverA_count_threshold", "",
                                                value = kovera_A, min = 0, step = 1)),
             column(width = 3, h5("counts in at least")),
             column(width = 2, uiOutput("filter_ui_kOverA_k")),
             column(width = 2, h5("samples"))
           ),
           br(),
           plotOutput("read_counts_plot"),
           fluidRow(
             actionButton("otu_reset_button", label = "Reset OTU Filter", icon = icon("trash")),
             actionButton("otu_filter_go", label = "Apply OTU Filter", icon = icon("bar-chart"))
           )
           )
           )
  ),

  
  tabPanel("Group Designation",
           #verbatimTextOutput("summ_filt"),
           #verbatimTextOutput("nrow_edata"),
           h2("Groupings"),
           h3("Main Effects"),
           # fluidRow(
           #    column(width=3, uiOutput("group1")),
           #    column(width=3, uiOutput("group2"))
           #    ),
           uiOutput("gdfMainEffect"),
           #actionButton("covs","Add Covariates"),
           # h3("Covariates"),
           # fluidRow(
           #   column(width=3, uiOutput("cov1")),
           #   column(width=3, uiOutput("cov2"))
           # ),
           #uiOutput("cov1"),
           #uiOutput("cov2"),
           # fluidRow(
           #   actionButton("groupDF_reset_button", label = "Reset Groupings", icon = icon("trash")),
           #   actionButton("groupDF_go", label = "Apply Groupings", icon = icon("check"))
           # ),
           p("The following table shows which group each sample belongs to."),
           DT::dataTableOutput("group_DF"),
           br(),
           p("The following table shows the number of samples in each group."),
           DT::dataTableOutput("group_tab"),
           br()
   ),
   
 tabPanel("Outliers",
          p("Use the Jaccard Index to look for other outliers in the dataset."),
          br(),
          plotlyOutput("jac_plot")
   ),
 
  tabPanel("Normalization",
           h2("Which normalization function to use?"),
           uiOutput("normFunc"),
           #actionButton("normGo","Normalize Data"),
           DT::dataTableOutput("normData"),
           br(),
           h4("Richness vs Abundance"),
           p("The normalization should reduce the correlation between richness
             and abundance. If it doesn't appear to, a different normalization
             function might be preferable."),
           fluidRow(
             splitLayout(cellWidths = c("50%","50%"), 
                         plotlyOutput("ra_raw"),
                         plotlyOutput("ra_norm"))
           ),
           uiOutput("norm_class"),
           plotlyOutput("norm_plot")
  ),
  
  tabPanel("Community Metrics",
           h3("Plot Parameters"),
           uiOutput("xaxis"),
           uiOutput("color"),
           h3("Alpha Diversity"),
           uiOutput("adiv_index"),
           plotOutput("adiv_plot"),
           verbatimTextOutput("adiv_summary"),
           br(),
           h3("Richness"),
           uiOutput("rich_index"),
           plotOutput("rich_plot"),
           verbatimTextOutput("rich_summary"),
           br(),
           h3("Abundance"),
           plotOutput("abun_plot"),
           verbatimTextOutput("abun_summary"),
           br(),
           h3("Evenness"),
           uiOutput("even_index"),
           plotOutput("even_plot"),
           verbatimTextOutput("even_summary"),
           br(),
           h3("Effective Species"),
           plotOutput("effsp_plot"),
           verbatimTextOutput("effsp_summary")
  ),
 
  tabPanel("Ordination",
           h3("Parameters"),
           uiOutput("beta_index"),
           actionButton(
             inputId = "submit_goe",
             label = "Submit"
           ),
           plotOutput("dimcheck"),
           #uiOutput("ord_method"),
           fluidRow(
             column(width=3, uiOutput("k")),
             column(width=3, uiOutput("ord_x"))
           ),
           fluidRow(
             column(width=3, uiOutput("ord_colors")),
             column(width=3, uiOutput("ord_y"))
           ),
           actionButton(
             inputId = "submit_ord",
             label = "Submit"
           ),
           uiOutput("ellipses"),
           #DT::dataTableOutput("beta"),
           #DT::dataTableOutput("mydist"),
           plotOutput("ord_plot")
  ),
 
 tabPanel("Differential Abundance",
          h3("Parameters"),
          uiOutput("da_index"),
          uiOutput("pval_adjust"),
          uiOutput("pval_thresh"),
          actionLink("selectall", "Select All"),
          uiOutput("comparisons"),
          actionButton(
            inputId = "submit_da",
            label = "Submit"
              ),
          DT::dataTableOutput("da_res"),
          verbatimTextOutput("da_summary"),
          plotOutput("flag_plot"),
          plotOutput("logfc_plot"),
          plotOutput("plot_all_da")
  ),
 
 tabPanel("Indicator Species",
          h3("Parameters"),
          uiOutput("within"),
          uiOutput("is_pval_thresh"),
          actionButton(
            inputId = "submit_is",
            label = "Submit"
          ),
          DT::dataTableOutput("indsp_results"),
          verbatimTextOutput("indsp_summary"),
          uiOutput("indsp_xaxis"),
          uiOutput("indsp_group"),
          plotOutput("indsp_plot")
  ),
 
 tabPanel("Statistics Results",
          conditionalPanel(
            condition = "input.submit_is == TRUE && input.submit_da == TRUE",
            DT::dataTableOutput("stats_res")
          )
          #plotOutput("newisplot")
          #if(exists(indsp_res()) & exists(diffabun_res())){
          #  DT::dataTableOutput("stats_res")
          # }else{
          #   p("This page is for combining the results of differential abundance analysis and indicator species analysis.")
          # }
  ),
 
 tabPanel("Differential Abundance - ALDEx2",
          h3("Parameters"),
          uiOutput("pa_mainEffects"),
          uiOutput("pa_randomEffect"),
          uiOutput("pa_Interactions"),
          uiOutput("mcsamples"),
          actionButton(
            inputId = "submit_pa",
            label = "Submit"
          ),
          DT::dataTableOutput("pa_res"),
          verbatimTextOutput("pa_summary"),
          plotOutput("pa_pval_plot"),
          plotOutput("pa_flag_plot")
 ),
 
 tabPanel("Download",
          uiOutput("files_to_download"),
          downloadButton("downloadData","Download")
  ),
  
 tabPanel("Meg's Tab",
          uiOutput("megs_output")
  )
  
)) #end page