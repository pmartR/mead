kovera_k <- 0
kovera_A <- 0
source("./functions/helper_functions.R")

shinyUI(navbarPage(title = (windowTitle = "mead"),

  
 # titlePanel(div(img(src = "Honey_Jar.png", height = 33, width = 22), "mead")),
 

 tabPanel("Data and Filtering", 
          tags$head(
            # tags$style(HTML("
            #                  .shiny-output-error-validation {
            #                  color: #D8000C;
            #                  }
            #                  "))
          ),
          
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
  # tags$table(
  #   DT::dataTableOutput("new_samples"),
  #   tags$head(
  #     tags$title(h4("Uploaded Metadata View"))
  #   )),
  h3("Sample Filtering"),
  h4("Keep samples with specific metadata"),
  # uiOutput("plots"),
  br(),
  fluidRow(
    column(width = 6, uiOutput("boxes")),
    column(width = 6, uiOutput("plots"))
    ),
  #uiOutput("boxes"),
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
  
  tabPanel("Group Designation",
           h2("Groupings"),
           h3("Main Effects"),
           fluidRow(
              column(width=3, uiOutput("group1")),
              column(width=3, uiOutput("group2"))
              ),
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
          plotOutput("jac_plot")
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
                         plotOutput("ra_raw"),
                         plotOutput("ra_norm"))
           ),
           plotOutput("norm_plot")
  ),
  
  tabPanel("Community Metrics",
           h3("Plot Parameters"),
           uiOutput("xaxis"),
           uiOutput("color"),
           h3("Alpha Diversity"),
           uiOutput("adiv_index"),
           plotOutput("adiv_plot"),
           br(),
           h3("Richness"),
           uiOutput("rich_index"),
           plotOutput("rich_plot")),
 
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
           plotOutput("ord_plot")),
 
 tabPanel("Differential Abundance",
          h3("Parameters"),
          uiOutput("da_index"),
          uiOutput("pval_adjust"),
          uiOutput("pval_thresh"),
          uiOutput("comparisons"),
          actionButton(
            inputId = "submit_da",
            label = "Submit"
              ),
          DT::dataTableOutput("da_res"),
          plotOutput("flag_plot"),
          plotOutput("logfc_plot"),
          plotOutput("plot_all_da")),
 
 tabPanel("Indicator Species",
          h3("Parameters"),
          uiOutput("within"),
          uiOutput("is_pval_thresh"),
          actionButton(
            inputId = "submit_is",
            label = "Submit"
          ),
          DT::dataTableOutput("indsp_results"),
          uiOutput("indsp_xaxis"),
          uiOutput("indsp_group"),
          plotOutput("indsp_plot")),
 
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
  
 tabPanel("Meg's Tab",
          uiOutput("megs_output"))
  
)) #end page