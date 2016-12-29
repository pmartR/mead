
library(shiny)
library(phyloseq)
library(ggplot2)
library(filterWidget)
theme_set(theme_bw())
source("https://bioconductor.org/biocLite.R")
source("./functions/helper_functions.R")
max_plots <- 100

get_plot_output_list <- function(meta_data) {
  # Insert plot output objects the list
  input_n <- ncol(meta_data)
  plot_output_list <- lapply(1:input_n, function(i) {
    plotname <- paste("plot", i, sep="")
    plot_output_object <- filterWidget::filterWidgetOutput(plotname, height = 100, width = 300)
    plot_output_object <- renderfilterWidget({
      filterWidget(as.character(names(meta_data)[i]), meta_data)
    })
  })
  do.call(tagList, plot_output_list) # needed to display properly.

  return(plot_output_list)
}

shinyServer(function(input, output) {

  #---------- BIOM -------------#
  BIOM <- eventReactive( input$go_button, {
    if (is.null(input$biom$datapath)) {
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
    if (is.null(input$qiime)) {
      return(NULL)
    }
    if (!is.null(input$qiime)) {
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

  #--------- Data-reactive variable lists also borrowed from phyloseq --------------#

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

  vars = function(type = "both", withnull = TRUE, singles = FALSE) {
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
    if (type == "both") {
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

  #---------------------------------------- Load Data Tab ----------------------------------------#

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
  # Insert the right number of plot output objects into the web page

  observeEvent(Scoping_only_meta(), {
    output$plots <- renderUI({ get_plot_output_list(Scoping_only_meta()) })
  })
# Get subset based on selection

 #co2_subset <- prune_samples(sample_names(full_data()) %in% event.data$X.SampleID, full_data())

 output$metatable <- renderTable({

   return(event.data())
 })

  # If NULL dont do anything
 # if(is.null(event.data) == T) return(NULL)
  #---------------------------------------- kOverA Filtering Tab ----------------------------------------#


  kovera_k <- 0

  maxSamples = reactive({
  # Create logical indicated the samples to keep, or dummy logical if nonsense input
  if (inherits(full_data(), "phyloseq")) {
    return(nsamples(full_data()))
  } else {
    return(NULL)
  }
  })

  output$filter_ui_kOverA_k <- renderUI({
    numericInputRow("filter_kOverA_sample_threshold", "k: number of samples in which a taxa exceeded A",
                 min = 0, max = maxSamples(), value = kovera_k, step = 1, class = "col-md-12")
  })

  filtered_data = reactive({
    #' On click of filter-button, filters the imported data with genefilter to k
    #' elements that exceed A.
    ps0 = full_data()
    if (input$actionb_filter == 0) {
      # Don't execute filter if filter-button has never been clicked.
      if (inherits(ps0, "phyloseq")) {
        return(ps0)
      } else {
        return(NULL)
      }
    }
    # Isolate all filter code so that button click is required for update
    isolate({
      if (inherits(ps0, "phyloseq")) {
          if (input$filter_kOverA_sample_threshold > 1) {
            # kOverA OTU Filtering
            flist = genefilter::filterfun(
              genefilter::kOverA(input$filter_kOverA_sample_threshold,
                                 input$filter_kOverA_count_threshold, na.rm=TRUE)
            )
            koatry = try(ps0 <- filter_taxa(ps0, flist, prune=TRUE))
            if (inherits(koatry, "try-error")) {
              warning("kOverA parameters resulted in an error, kOverA filtering skipped.")
            }
          }
        return(ps0)
      } else {
        return(NULL)
      }
    })
  })

  output$contents <- renderUI({
  output_phyloseq_print_html(full_data())
  })

  output$filtered_contents <- renderUI({
    output_phyloseq_print_html(filtered_data())
  })

  output$sample_variables <- renderText({return(
    paste0(sample_variables(full_data(), errorIfNULL=FALSE), collapse=", ")
  )})

  output$rank_names <- renderText({return(
    paste0(rank_names(full_data(), errorIfNULL=FALSE), collapse=", ")
  )})

  output$filter_summary_plot <- renderPlot({
    plib0 = lib_size_hist() + ggtitle("Original Data")
    potu0 = otu_sum_hist() + ggtitle("Original Data")

    if (inherits(filtered_data(), "phyloseq")) {
      potu1 = sums_hist(taxa_sums(filtered_data()), xlab = "Number of Reads (Counts)",
                        ylab = "Number of OTUs"
      ) +
        ggtitle("Filtered Data")
      plib1 = sums_hist(sample_sums(filtered_data()), xlab = "Number of Reads (Counts)",
                        ylab = "Number of Libraries"
      ) +
        ggtitle("Filtered Data")
    } else {
      potu1 = plib1 = fail_gen()
    }
    gridExtra::grid.arrange(plib0, potu0, plib1, potu1, ncol=2)
  })


  #---------------------------------------- Alpha Diversity Tab ----------------------------------------#
  output$rich_uix_color <- renderUI({
    selectInput("color_rich",
                label = "Color",
                choices = c(list("samples"), vars("samples")),
                selected = "NULL")
  })

  output$downloadRichnessEstimates <- downloadHandler(
    filename = "alpha_diversity_estimates.csv",
    content = function(file) {
    write.csv(estimate_richness(pruned_data(), split = as.logical(input$split), measures = input$dist_measures), file)
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

  #---------------------------------------- Ordination Tab ----------------------------------------#

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

