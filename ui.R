kovera_k <- 0
kovera_A <- 0
source("./functions/helper_functions.R")
# Get dependencies imported if any are missing
#list of all packages
packs <- installed.packages()[,"Package"]
#list of required packages
dependencies <- c("shiny", "shinyjs", "shinycssloaders", "lazyeval","lubridate", "V8", "dplyr", "DT", "ggplot2", "vegan", "goeveg", "plotly", "DESeq2", "edgeR", "indicspecies", "ALDEx2", "fdrtool")

# add in custom packages and installs
if (!("filterWidget" %in% packs)) devtools::install("filterWidget")
if (!("pmartRseq" %in% packs)) devtools::install("../../pmartRseq/")
if (!("DESeq2" %in% packs)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
}
if (!("edgeR" %in% packs)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
}
if (!("ALDEx2" %in% packs)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("ALDEx2")
}

#check for missing packages
missing <- dependencies[!(dependencies %in% packs)]

#install missing ones
if (length(missing) > 0) install.packages(missing, dependencies = TRUE)

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
library(indicspecies)
library(plotly)
library(DESeq2)
library(fdrtool)
library(shinycssloaders)
source("./functions/helper_functions.R")
source("./functions/test_functions.R")
# create a reset session button using a js method
#jsResetCode <- "shinyjs.reset = function() {history.go(0)}"
tagList(
  navbarPage(title = div(img(src = "Honey_Jar.png", height = 30, width = 15), "mead"),
        windowTitle = "mead",
        #theme = "yeti.css",

# titlePanel(div(img(src = "Honey_Jar.png", height = 33, width = 22), "mead")),
 
#---------------- Load Data Page ----------------#
 tabPanel("Data", 
          tags$head(
          ),
          # fluidRow(column(width = 12,
          #                 useShinyjs(),                                           
          #                 extendShinyjs(text = jsResetCode),                     
          #                 actionButton("reset_button", "Reset All") 
          #                 )
          #   
          # ),

          fluidRow(
            column(width = 12,
                   tags$hr(),
                   h3("Load Data"),
                   fluidRow(
                     column(width = 4, fileInput('e_data', 'Choose Data File')),
                     # column(width = 4, fileInput('fasta', 'Choose FASTA file')),
                     column(width = 4, fileInput('f_data', 'Choose Sample Metadata File')),
                     column(width = 4, fileInput('e_meta', 'Optional, Choose Feature Metadata File'))
                   )
            )
          ),

          # fluidRow(
          #   column(width = 12,
          #          tags$hr(),
          #          h3("Specify Column Names"),
          #          fluidRow(
          #            column(width = 4, uiOutput("edata_cname")),
          #            column(width = 4, uiOutput("fdata_cname")),
          #            column(width = 4, uiOutput("taxa_cname"))
          #          )
          #   )
          # ),
          # actionButton("Upload", "Upload ze data!"),
          #-------- Guess which colnames are the proper identifiers --------#
          h3("Selected Identifiers"),
          fluidRow(
            column(width = 6,
                    h4("The selected expression data identifier is:"),
                    h3(textOutput("e_data_cname"))),
            column(width = 6,
                   h4("The selected metadata identifier is:"),
                   h3(textOutput("f_data_cname")))
          ),
          h4("If the selected identifiers are not correct, use the dropdowns below to change them"),
          fluidRow(
            column(width = 6,
            uiOutput("new_edata_cname")
            ),
            column(width = 6,
                   uiOutput("new_fdata_cname")
            )
          ),
          hr(),
          h3("Data View"),
          fluidRow(
            column(width = 12, tags$table(
              DT::dataTableOutput("sample_data"),
              tags$head(
                tags$title(h4("Uploaded Data View"))
              )))
          ),
          br(),
          hr(),
          h2("Taxonomic Level"),
          uiOutput("rollup"),
          br(),
          hr(),
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
  ),#end data upload

#---------------- All Filtering ----------------#
  tabPanel("Filtering",
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
                    radioButtons(inputId = "count_filter_fun", label = "Choose filter type:",
                                 choices = list("mean", "percent", "max", "sum", "nonmiss", "ka"),
                                 selected = "mean", inline = TRUE),
                    h4("OTU Counts Per Sample"),
                    p("Keep OTUs with some minimum number of counts"),
                    fluidRow(
                      column(width = 3, h5("Keep OTU Threshold"), align = 'center'),
                      column(width = 2,  numericInputRow("filter_count_threshold", "",
                                                         value = kovera_A, min = 0, step = 1)),
                      column(width = 2, uiOutput("filter_ui_kOverA_k"))
                      ),
                    br(),
                    plotOutput("read_counts_plot"),
                    fluidRow(
                      actionButton("otu_reset_button", label = "Reset OTU Filter", icon = icon("trash")),
                      actionButton("otu_filter_go", label = "Apply OTU Filter", icon = icon("bar-chart"))
                    )
             )
           ),
           br(),
           hr(),
           h4("Taxonomic Filtering"),
           uiOutput("criteria"),
           uiOutput("keep_taxa"),
           verbatimTextOutput("taxa_counts"),
           fluidRow(
             actionButton("taxa_reset_button", label = "Reset Taxa Filter", icon = icon("trash")),
             actionButton("taxa_filter_go", label = "Apply Taxa Filter", icon = icon("bar-chart"))
           ),
           br()
    ),

  # tabPanel("Group Designation",
  #          #verbatimTextOutput("summ_filt"),
  #          #verbatimTextOutput("nrow_edata"),
  #
  #  ),

 tabPanel("Outliers",
          p("Use the Jaccard Index to look for other outliers in the dataset."),
          br(),
          #plotOutput("jac_plot"),
          fluidRow(
            # column(width = 6,
            #        plotOutput("outlier_abundance_plot", height=300,  brush = brushOpts(id = "plot_brush"))
            #        ),
            # column(width = 6,
            #        h4("Brushed points"),
            #        tableOutput("plot_brushedpoints")
            #        )
            column(width = 4,
                   plotlyOutput("outlier_abundance_plot", height = 300)
            ),
            column(width = 4,
                   plotlyOutput("outlier_jaccard_plot", height = 300)
            ),
            column(width = 4,
                   plotlyOutput("outlier_richness_plot", height = 300)
            )

            ),
          br(),
          fluidRow(
            column(width = 6, tableOutput("audies")),
            column(width = 6, actionButton("remove_outliers", label = "Remove Outliers"))
          )
   ),#end outliers

#---------------- Normalization ----------------#
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
  ),#end normalization

#---------------- Community Metrics ----------------#  
  tabPanel("Community Metrics",
           h3("Plot Parameters"),
           uiOutput("xaxis"),
           uiOutput("color"),
           h3("Alpha Diversity"),
           uiOutput("adiv_index"),
           plotlyOutput("adiv_plot"),
           verbatimTextOutput("adiv_summary"),
           br(),
           h3("Richness"),
           uiOutput("rich_index"),
           plotlyOutput("rich_plot"),
           verbatimTextOutput("rich_summary"),
           br(),
           h3("Abundance"),
           plotlyOutput("abun_plot"),
           verbatimTextOutput("abun_summary"),
           br(),
           h3("Evenness"),
           uiOutput("even_index"),
           plotlyOutput("even_plot"),
           verbatimTextOutput("even_summary"),
           br(),
           h3("Effective Species"),
           plotlyOutput("effsp_plot"),
           verbatimTextOutput("effsp_summary")
  ),# end community metrics

#---------------- Ordination ----------------# 
  tabPanel("Ordination",
           h3("Parameters"),
           uiOutput("beta_index"),
           actionButton(
             inputId = "submit_goe",
             label = "Submit"
           ),
           withSpinner(plotOutput("dimcheck")),
           #uiOutput("ord_method"),
           fluidRow(
             column(width = 3, uiOutput("k")),
             column(width = 3, uiOutput("ord_x"))
           ),
           fluidRow(
             column(width = 3, uiOutput("ord_colors")),
             column(width = 3, uiOutput("ord_y"))
           ),
           actionButton(
             inputId = "submit_ord",
             label = "Submit"
           ),
           uiOutput("ellipses"),
           #DT::dataTableOutput("beta"),
           #DT::dataTableOutput("mydist"),
           withSpinner(plotlyOutput("ord_plot"))
  ),# end ordination

#---------------- Defferential Abundance ----------------#
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
          withSpinner(dataTableOutput("da_res")),
          verbatimTextOutput("da_summary"),
          plotOutput("flag_plot"),
          plotOutput("logfc_plot"),
          plotOutput("plot_all_da")
  ),#end differential abundance

#---------------- Indicator Species ----------------#
 tabPanel("Indicator Species",
          h3("Parameters"),
          uiOutput("within"),
          uiOutput("max_indsp_grp"),
          uiOutput("is_pval_thresh"),
          actionButton(
            inputId = "submit_is",
            label = "Submit"
          ),
          withSpinner(dataTableOutput("indsp_results")),
          verbatimTextOutput("indsp_summary"),
          uiOutput("indsp_xaxis"),
          uiOutput("indsp_group"),
          plotlyOutput("indsp_plot", width = "100%")
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

#---------------- ALDEx2 Differential Abundance ----------------#
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
          withSpinner(dataTableOutput("pa_res")),
          verbatimTextOutput("pa_summary"),
          plotOutput("pa_pval_plot"),
          plotOutput("pa_flag_plot")
 ), #end ALDEx2 diff abundance

#---------------- Network Analysis -------------# 
 tabPanel("Network Analysis",
          h3("Parameters"),
          fluidRow(
            column(width = 4, uiOutput("na_group")),
            column(width = 4, uiOutput("na_size")),
            column(width = 4, uiOutput("na_colour"))
          ),
          fluidRow(
            column(width = 4, uiOutput("na_group_var")),
            column(width = 4, uiOutput("na_coeff"))
          ),
          fluidRow(
            column(width = 4, uiOutput("na_missingval")),
            column(width = 4, uiOutput("na_qval"))
          ),
          # uiOutput("na_group"),
          # uiOutput("na_group_var"),
          # uiOutput("na_missingval"),
          # uiOutput("na_coeff"),
          # uiOutput("na_qval"),
          #uiOutput("na_colour"),
          # uiOutput("na_size"),
          actionButton(
            inputId = "submit_na",
            label = "Submit"
          ),
          withSpinner(plotOutput("na_network_plot")),
          br(),
          hr(),
          h3("Modules"),
          uiOutput("na_cluster"),
          uiOutput("na_mod_size"),
          actionButton(
            inputId = "submit_modules",
            label = "Submit Module Detection"
          ),
          withSpinner(plotOutput("na_mod_plot")),
          br(),
          hr(),
          h3("Environmental Variables"),
          uiOutput("na_envvars"),
          uiOutput("env_pval"),
          actionButton(
            inputId = "submit_envvars",
            label = "Submit Environmental Variables"
          ),
          withSpinner(plotOutput("na_envvars_plot"))
          ),

 tabPanel("Download",
          uiOutput("files_to_download"),
          downloadButton("downloadData","Download")
  )#,
#---------- Placeholder for Meg's Awesome Viz----------#
 #  
 # tabPanel("Meg's Tab",
 #          uiOutput("megs_output")
 #  )
 #  
) #end page
)