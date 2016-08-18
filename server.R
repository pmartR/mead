
library(shiny)
library(phyloseq)
library(ggplot2)
theme_set(theme_bw())
source("https://bioconductor.org/biocLite.R")
source("./functions/helper_functions.R")
shinyServer(function(input, output) {

  #---------- BIOM -------------#
  BIOM <- reactive({
    if(is.null(input$biom$datapath)){
      return(NULL)
    }
      BIOM_file <- import_biom(BIOMfilename = input$biom$datapath, 
                               refseqfilename = input$fasta$datapath, 
                               treefilename = input$tree$datapath, 
                               parseFunction = parse_taxonomy_greengenes)
    return(BIOM_file)
  })

  #--------- Metadata -----------#
  Scoping_only_meta <- reactive({
    if(is.null(input$qiime)) {
      return(NULL)
    }
    if(!is.null(input$qiime)) {
      return(import_qiime_sample_data(input$qiime$datapath))
    }
  })
  #--------- Merge BIOM and Metadata -----------#
  full_data <- reactive({
    if (is.null(c(BIOM(),Scoping_only_meta()))) {
      return(NULL)
    }
    if (!is.null(c(BIOM(),Scoping_only_meta()))) {
      return(merge_phyloseq(BIOM(), Scoping_only_meta()))
    }
  })

  #--------- Prune for alpha diversity -----------#
  pruned_data <- reactive({
    if (!is.null(full_data())) {
      dat <- full_data()
      temp <- prune_taxa(taxa_sums(dat) > 0, dat)
    }
    if (is.null(full_data())) {
      temp <- NULL
    }
    return(temp)
  })
  
  # borrowed historgram function from phyloseq
  lib_size_hist = reactive({
    if (is.null(full_data())) {
      return(NULL)
    }
    xlab = "Number of Reads (Counts)"
    ylab = "Number of Libraries"
    return(sums_hist(sample_sums(full_data()), xlab, ylab))
  })
  otu_sum_hist = reactive({
    if (is.null(full_data())) {
      return(NULL)
    }
    xlab = "Number of Reads (Counts)"
    ylab = "Number of OTUs"
    return(sums_hist(taxa_sums(full_data()), xlab, ylab))    
  })
  
  #--------- Data-reactive variable lists --------------#
  
  rankNames = reactive({
    rankNames = as.list(rank_names(full_data(), errorIfNULL=FALSE))
    names(rankNames) <- rankNames
    return(rankNames)
  })
  
  variNames = reactive({
    variNames = as.list(sample_variables(full_data(), errorIfNULL=FALSE))
    names(variNames) <- variNames
    return(variNames)
  })
  
  vars = function(type="both", withnull=TRUE, singles=FALSE){
    if (!type %in% c("both", "taxa", "samples")) {
      stop("incorrect `type` specification when accessing variables for UI.")
    }
    returnvars = NULL
    if (type == "samples") {
      if (singles) {
        returnvars <- c(list(Sample = "Sample"), variNames())
      } else {
        returnvars <- variNames()
      }
    }
    if (type == "taxa") {
      if (singles) {
        returnvars <- c(rankNames(), list(OTU = "OTU"))
      } else {
        returnvars <- rankNames()
      }
    } 
    if ( type == "both") {
      # Include all variables
      if (singles) {
        returnvars <- c(rankNames(), variNames(), list(OTU = "OTU", Sample = "Sample"))
      } else {
        returnvars <- c(rankNames(), variNames())
      }
    }
    if (withnull) {
      # Put NULL first so that it is default when `select` not specified
      returnvars <- c(list("NULL" = "NULL"), returnvars)
    }
    return(returnvars)
  }
  
  # A generic selectInput UI. Plan is to pass a reactive argument to `choices`.
  uivar = function(id, label="Variable:", choices, selected="NULL"){
    selectInput(inputId = id, label = label, choices = choices, selected = selected)
  }
  
  #---------- Outputs -------------------------#

  output$library_sizes <- renderPlot({
    if (is.null(full_data())) {
      return(NULL)
    } 
    if (!is.null(full_data())) {
      p = lib_size_hist() + ggtitle("Library Sizes")
      q = otu_sum_hist() + ggtitle("OTU Totals")
      gridExtra::grid.arrange(p, q, ncol = 2)
    }
  })
  
  output$sample_metadata <- renderDataTable({
    if (is.null(Scoping_only_meta())) {
      return(NULL)
    }
    Scoping_only_meta()
  })
  
  output$downloadOTUtable <- downloadHandler(
    filename = "OTU_Sample_Table.csv",
    content = function(file) {
      write.csv(otu_table(full_data()), file)
    }
  )
  
  #---------------- Alpha Diversity ------------------#
  output$rich_uix_color <- renderUI({
    selectInput("color_rich",
                label = "Color",
                choices = c(list("samples"), vars("samples")),
                selected = "NULL")
  })
  
  output$downloadRichnessEstimates <- downloadHandler(
    filename = "alpha_diversity_estimates.csv",
    content = function(file) {
    write.csv(estimate_richness(pruned_data(), split = as.logical(input$split), measures =input$dist_measures), file)
  })
  
  output$richness_plot <- renderPlot({
    if (!is.null(pruned_data())) {
      plot_data <- pruned_data()
      p <- plot_richness(plot_data, measures = input$dist_measures)
    if (!is.null(input$color_rich)) {
      p$mapping$colour <- as.symbol(input$color_rich)
      p <- update_labels(p, list(colour = input$color_rich))
      return(p)
    }

    return()
    }
    if (is.null(pruned_data)) {
      return(NULL)
    }
  })
  output$downloadDiversityImage <- 
    downloadHandler(
      filename = "alpha_diversity_plot.png",
      content = function(file) {
        plot_data <- pruned_data()
        plot <- estimate_richness(plot_richness(plot_data, measures = input$dist_measures))
        ggsave(plot, filename = file, dpi = 400)
      })
  
  #---------------- Ordination ------------------#
  
  output$ordination_plot <- renderPlot({
    if (!is.null(pruned_data())) {
      plot_data <- pruned_data()
      ord_plot_data <- ordinate(plot_data, "NMDS")
      return(plot_ordination(plot_data, ord_plot_data))
    }
    if (is.null(pruned_data)) {
      return(NULL)
    }
  })
  output$downloadOrdinationImage <- 
    downloadHandler(
      filename = "ordination_plot.png",
      content = function(file) {
        plot_data <- pruned_data()
        plot <- estimate_richness(plot_richness(plot_data, measures = input$dist_measures))
        ggsave(plot, filename = file, dpi = 400)
      })
})

 