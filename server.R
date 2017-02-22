library(shiny)
library(phyloseq)
library(ggplot2)
library(filterWidget)
library(DT)
theme_set(theme_bw())
source("https://bioconductor.org/biocLite.R")
source("./functions/helper_functions.R")
source("./functions/test_functions.R")

shinyServer(function(input, output) {
  
  #---------- BIOM -------------#
  biom_obj <- reactive({
    validate(
      need(input$biom$datapath != "", "Please select a biom file")
    )
    if (is.null(input$biom$datapath)) {
      return(NULL)
    }else{
      BIOM_file <- import_biom(BIOMfilename = input$biom$datapath, 
                               refseqfilename = input$fasta$datapath, 
                               treefilename = input$tree$datapath, 
                               parseFunction = parse_taxonomy_greengenes)
      return(BIOM_file)
    }
  })
  
  #--------- Metadata -----------#
  metadata_obj <- reactive({
    validate(
      need(input$qiime$datapath != "", "Please select a QIIME file")
    )
    if (is.null(input$qiime)) {
      return(NULL)
    }else{
      temp <- import_qiime_sample_data(input$qiime$datapath)
      inds <- lapply(temp, function(x) sum(is.na(x)) == length(x) )
      results <- temp[,!(unlist(inds))]
      return(results)
    }
  })
  
  
  #--------- Merge BIOM and Metadata -----------#
  full_data <- reactive({
    #--------- Warnings and errors -----------#
    # complete mismatch in biom and QIIME
    validate(
      try(biom_not_matching_metadata(biom_input = sample_names(biom_obj()), metadata_input = sample_names(metadata_obj()))),
      try(metadata_not_matching_biom(biom_input = sample_names(biom_obj()), metadata_input = sample_names(metadata_obj())))
    )
    # partial mismatch in biom and QIIME
    # using shiny notifications
    metadata_warning <- NULL
    # observeEvent(input$biom, {
    #   if (!is.null(metadata_warning))
    #     return()
    #   id <<- showNotification("Warning: Samples in biom not in QIIME")
    # })
    if (try(metadata_mismatching_biom(biom_input = sample_names(biom_obj()), metadata_input = sample_names(metadata_obj())))) {
      metadata_warning <<- showNotification("Warning: Samples in QIIME not in biom", duration = NA, type = "error")
    }
    if (try(biom_mismatching_metadata(biom_input = sample_names(biom_obj()), metadata_input = sample_names(metadata_obj())))) {
      metadata_warning <<- showNotification("Warning: Samples in biom not in QIIME", duration = NA, type = "error")
    }
    if (any(is.null(c(biom_obj(),metadata_obj())))) {
      return(NULL)
    }else{return(merge_phyloseq(biom_obj(), metadata_obj()))}
  })
  
  #--------- Prune for alpha diversity -----------#
  pruned_data <- reactive({
    if (is.null( full_data())) {
      temp <- NULL
    }else{
      dat <-  full_data()
      temp <- prune_taxa(taxa_sums(dat) > 0, dat)
    }
    return(temp)
  })
  
  # borrowed historgram function from phyloseq
  lib_size_hist = reactive({
    if (is.null( full_data())) {
      return(NULL)
    }else{
      xlab = "Number of Reads (Counts)"
      ylab = "Number of Libraries"
      return(sums_hist(sample_sums( full_data()), xlab, ylab))
    }
  })
  otu_sum_hist = reactive({
    if (is.null( full_data())) {
      return(NULL)
    }else{
      xlab = "Number of Reads (Counts)"
      ylab = "Number of OTUs"
      return(sums_hist(taxa_sums( full_data()), xlab, ylab))    
    }
  })
  
  #--------- Data-reactive variable lists also borrowed from phyloseq --------------#
  
  rankNames = reactive({
    rankNames = as.list(rank_names( full_data(), errorIfNULL=FALSE))
    names(rankNames) <- rankNames
    return(rankNames)
  })
  
  variNames = reactive({
    variNames = as.list(sample_variables( full_data(), errorIfNULL=FALSE))
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
    if (is.null( full_data())) {
      return(NULL)
    }else{
      p = lib_size_hist() + ggtitle("Library Sizes")
      q = otu_sum_hist() + ggtitle("OTU Totals")
      gridExtra::grid.arrange(p, q, ncol = 2)
    } 
  })
  outputOptions(output, "library_sizes", suspendWhenHidden = FALSE)
  
  observeEvent(input$qiime, 
               output$sample_metadata <- DT::renderDataTable(
                 data.frame(metadata_obj()), rownames = FALSE, class = 'cell-border stripe compact hover',
                 options = list(columnDefs = list(list(
                   targets = c(1:(ncol(metadata_obj()) - 1)),
                   render = JS(
                     "function(data, type, row, meta) {",
                     "return type === 'display' && data.length > 10 ?",
                     "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
                     "}")
                 ))), callback = JS('table.page(3).draw(false);'))
  )
  
  output$downloadOTUtable <- downloadHandler(
    filename = "OTU_Sample_Table.csv",
    content = function(file) {
      write.csv(otu_table(full_data()), file)
    }
  )
  
  #---- Charts for filtering -------#
  new_metadata_obj <- eventReactive(input$selected_indices, {
    return(metadata_obj()[input$selected_indices + 1, ]) # initialize object
  })
  observeEvent(metadata_obj(), {
    # If the metadata has been loaded, create the charts
    output$plots <- renderUI({get_plot_output_list(metadata_obj())})
    # Return a new metadata object that will be used for filtering
    new_metadata_obj <- reactive({
      return(metadata_obj()[input$selected_indices + 1, ]) # add one because javascript is zero-indexed
    })
    
    # If the user brushes a chart, input$selected_indices will change
    # If that change is observed, subset new_metadata_obj and 
    # create subsetted charts and table
    observeEvent(input$selected_indices, {
      new_metadata_obj <- reactive({
        return(metadata_obj()[input$selected_indices + 1, ]) # add one because javascript is zero-indexed
      })
      output$plots <- renderUI({get_plot_output_list(new_metadata_obj())})
      output$new_samples <- DT::renderDataTable(
        data.frame(new_metadata_obj()), rownames = FALSE, class = 'cell-border stripe compact hover',
        options = list(columnDefs = list(list(
          targets = c(1:(ncol(metadata_obj()) - 1)),
          render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 10 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
            "}")
        ))), callback = JS('table.page(3).draw(false);'))
    })
    # render the first charts as soon as the data becomes available
    outputOptions(output, "plots", suspendWhenHidden = FALSE)
    outputOptions(output, "library_sizes", suspendWhenHidden = FALSE)
    outputOptions(output, "sample_metadata", suspendWhenHidden = FALSE)
  })
  #------- Reset everything or keep the filtered subset -------#
  observeEvent(input$reset_button,{
    output$plots <- renderUI({get_plot_output_list(metadata_obj())})
    output$new_samples <- DT::renderDataTable(
      data.frame(metadata_obj()), rownames = FALSE, class = 'cell-border stripe compact hover',
      options = list(columnDefs = list(list(
        targets = c(1:(ncol(metadata_obj()) - 1)),
        render = JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 10 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
          "}")
      ))), callback = JS('table.page(3).draw(false);'))
  })
  
  #--------- Merge BIOM and Metadata -----------#
  meta_filtered_data <- reactive({
    if(is.null(input$selected_indices)){
      return(full_data())
    }else{
      return(prune_samples(as.character(sample_names(new_metadata_obj())), full_data()))
    }
  })
  
  #---------------------------------------- kOverA Filtering Tab ----------------------------------------# 
  kovera_k <- 0
  observeEvent(meta_filtered_data(), {
    maxSamples = reactive({
      # Create logical indicating the samples to keep, or dummy logical if nonsense input
      if (inherits( meta_filtered_data(), "phyloseq")) {
        return(nsamples( meta_filtered_data()))
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
      ps0 =  meta_filtered_data()
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
      output_phyloseq_print_html( meta_filtered_data())
    })
    
    output$filtered_contents <- renderUI({
      output_phyloseq_print_html(filtered_data())
    })
    
    output$sample_variables <- renderText({return(
      paste0(sample_variables( meta_filtered_data(), errorIfNULL=FALSE), collapse=", ")
    )})
    
    output$rank_names <- renderText({return(
      paste0(rank_names( meta_filtered_data(), errorIfNULL=FALSE), collapse=", ")
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
