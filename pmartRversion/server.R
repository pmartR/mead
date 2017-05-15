
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(DT)
library(dplyr)
library(filterWidget)
library(ggplot2)
library(pmartRseq)
source("./functions/helper_functions.R")
source("./functions/test_functions.R")

shinyServer(function(input, output) {

  #---------- rRNA object import -------------#
  rRNAobj <- reactive({
  validate(
    need(input$biom != "", "Please select a biom file")
  )
    validate(
      need(input$qiime != "", "Please select a qiime file")
    )
    return(pmartRseq::as.rRNAdata(e_data = as.character(input$biom$datapath), f_data = as.character(input$qiime$datapath), edata_cname = "OTU"))
  }) #end rRNAobj
  
  #--------- Metadata Object for filtering -----------#
  metadata_obj <- reactive({
  #------- TODO: need to display an error if the rRNA object isn't created! ---------#
    validate(
      need(!(is.null(rRNAobj()$f_data)), message = "rRNA object fail. Check to make sure file paths are correct")
    )
    # validate(
    #   need(input$qiime$datapath != "", "Please select a QIIME file")
    # )
    if (is.null(input$qiime)) {
      return(NULL)
    }else{
      temp <- rRNAobj()$f_data
      inds <- lapply(temp, function(x) sum(is.na(x)) == length(x) )
      results <- temp[,!(unlist(inds))]
      return(results)
    }
  })
  
  observeEvent(input$qiime, 
               output$sample_metadata <- DT::renderDataTable(
                 data.frame(  metadata_obj()), rownames = FALSE, class = 'cell-border stripe compact hover',
                 options = list(columnDefs = list(list(
                   targets = c(1:(ncol((metadata_obj())) - 1)),
                   render = JS(
                     "function(data, type, row, meta) {",
                     "return type === 'display' && data.length > 10 ?",
                     "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
                     "}")
                 ))), callback = JS('table.page(3).draw(false);'))
  )
  
  #---- Charts for f_meta filtering -------#

  #initialize new_metadata_obj object (REMOVE)
  # new_metadata_obj <- eventReactive(metadata_obj(), {
  #   return(metadata_obj()[input$selected_indices + 1, ])
  # })
  # 
  observeEvent(metadata_obj(), {
    # If the metadata has been loaded, create the charts
    output$plots <- renderUI({get_plot_output_list(metadata_obj())})
    
    # Initialize new_metadata_obj object. If no metadata filtering then it's a copy of metadata_obj.
    # If metadata filtering, will be overwritten with filters 
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
    #outputOptions(output, "library_sizes", suspendWhenHidden = FALSE)
    #outputOptions(output, "sample_metadata", suspendWhenHidden = FALSE)
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
  #filtered_data_obj <- pmartRseq::as.rRNAdata()
  # end sample metadata filtering
  
  #--------- Object for filtering -----------#
  filtered_rRna_obj <- reactive({
    #------- TODO: need to display an error if the rRNA object isn't created! ---------#
    validate(
      need(!(is.null(rRNAobj()$e_data)), message = "rRNA object fail. Check to make sure file paths are correct")
    )
    if (is.null(input$biom)) {
      return(NULL)
    }else{
      results <- rRNAobj()
      return(results)
    }
  })
  
  # filtered_rRna_obj <- reactive({
  #   validate(need(!(is.null(biom_obj())), message = "biom object needed"))
  #   return(pmartRseq::as.rRNAdata(e_data = biom_obj(), f_data = metadata_obj(), edata_cname = "OTU", fdata_cname = names(metadata_obj()[1]))) 
  #   # If the metadata has been loaded, create the charts
  #   })

  
  #--------------- k over a  and library read count filtering  ---------#
  kovera_k <- 0
  observeEvent(rRNAobj(), {
    maxSamples = reactive({
      # Create logical indicating the samples to keep, or dummy logical if nonsense input
      validate(
        need(!(is.null(metadata_obj())), message = "please import sample metadata")
      )
      if (!is.null(metadata_obj())) {
        return(nrow(metadata_obj()))
      } else {
        return(NULL)
      }
    })
    
    output$filter_ui_kOverA_k <- renderUI({
      numericInputRow("filter_kOverA_sample_threshold", "",
                      min = 0, max = maxSamples(), value = kovera_k, step = 1)
    })
  })
  
  #-------------- Library read filtering -----------#
  sample_filter_obj <- reactive({
    pmartRseq::sample_based_filter(rRNAobj(), fn = "sum")
  })
  output$sample_counts_plot <- renderPlot({
    validate(
      need(input$n >= 0, message = "Enter a count minimum >= 0")
    )
    plot(sample_filter_obj(), min_num = input$n)
  })
  

  
  # output$sample_counts_plot <- renderPlot({
  #   plot_data <- pmartRseq::count_based_filter(omicsData = filtered_rRna_obj(), fn = "persample")
  #   p <- ggplot(plot_data, aes(x = sumSamps))+
  #     geom_histogram(color = "black", fill = "black", bins =  nrow(plot_data))+
  #     geom_vline(xintercept = as.numeric(input$n), color = "red")+
  #     ylab("Samples")+
  #     xlab("OTU Reads")+
  #     theme_bw()
  #   return(p)
  # })

  #-------------- OTU read filtering -----------#  
  otu_filter_obj <- reactive({
    pmartRseq::count_based_filter(rRNAobj(), fn = "sum")
  })
  
  output$read_counts_plot <- renderPlot({
    validate(
      need( input$filter_kOverA_count_threshold >= 0, message = "Enter a count minimum >= 0")
    )
    plot(otu_filter_obj(), min_num = input$filter_kOverA_count_threshold)
  })
  
  filtered_rRna_obj <- eventReactive(input$otu_filter_go, {
    pmartRseq::applyFilt(otu_filter_obj(), rRNAobj(), min_num = input$filter_kOverA_count_threshold)
  })
  

  
  # output$read_counts_plot <- renderPlot({
  #   plot_data <- pmartRseq::count_based_filter(omicsData = filtered_rRna_obj(), fn = "sum")
  #   p <- ggplot(plot_data, aes(x = sumOTUs))+
  #     geom_histogram(color = "black", fill = "black")+
  #     scale_x_log10()+
  #     geom_vline(xintercept = log10(as.numeric(input$filter_kOverA_count_threshold)), color = "red")+
  #     ylab("OTUs")+
  #     xlab("OTU Reads")+
  #     theme_bw()
  #   return(p)
  #   
  #   # Initialize new_metadata_obj object. If no metadata filtering then it's a copy of metadata_obj.
  #   # If metadata filtering, will be overwritten with filters 
  #   # new_biom_obj <- reactive({
  #   #   return(biom_obj()[input$selected_indices + 1, ]) # add one because javascript is zero-indexed
  #   # })
  # })
  
  #-------- Charts for e_meta filtering -----------#
  #new_biom_obj <- reactive()
  #-------- Create filtered object -----------#
  #filtered_rRna_obj <- as.rRNAdata(e_data = new_biom_obj(), f_data = new_metadata_obj())  
 
  # 
  # output$sample_counts_plot <- renderPlot({
  #   plot_data <- sample_counts()
  #   ggplot(plot_data, aes(x = sample_count))+
  #     geom_histogram()+
  #     theme(axis.text.x = element_text(angle = 90))+
  #     ylab("Libraries")+
  #     xlab("OTU Reads")
  # })
  # 
  # output$otu_counts_plot <- renderPlot({
  #   plot_data <- sample_counts()
  #   ggplot(plot_data, aes(x = sample_count))+
  #     geom_histogram()+
  #     theme(axis.text.x = element_text(angle = 90))+
  #     ylab("OTU Reads")+
  #     xlab("OTUs")
  # })
}) #end server
