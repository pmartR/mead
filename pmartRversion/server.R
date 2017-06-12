
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

filtered_rRNA_obj <- list()

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

  #--------- somewhere to store the applied filters ---------#
  filters <- reactiveValues(otu = list(), sample = list())

  observeEvent(input$otu_filter_go, {
    filters$otu[[input$otu_filter_go]] <- otu_filter_obj()
      # isolate(filtered_rRNA_obj <- pmartRseq::applyFilt(filter_object = filters$otu[[input$otu_filter_go]],
      #                                             omicsData = rRNAobj(), 
      #                                             num_samps = input$filter_kOverA_sample_threshold, 
      #                                             upper_lim = input$filter_kOverA_count_threshold))
      # print(str(filtered_rRNA_obj))
    })
  
  observeEvent(input$sample_filter_go, {
    filters$sample[[input$sample_filter_go]] <- sample_filter_obj()
    # isolate(filtered_rRNA_obj <- pmartRseq::applyFilt(filter_object = filters$otu[[input$otu_filter_go]],
    #                                             omicsData = rRNAobj(), 
    #                                             num_samps = input$filter_kOverA_sample_threshold, 
    #                                             upper_lim = input$filter_kOverA_count_threshold))
    # print(str(filtered_rRNA_obj))
  })
    
   # second_filt <-  pmartRseq::applyFilt(filter_object = filters$sample[[input$sample_filter_go]], omicsData = first_filt, upper_lim = input$n)
                                       
  #})
  
  
  # filtered_rRNA_obj <- reactive({
  #   if (values$default == 0) {
  #     return(rRNAobj())
  #            }
  #   else{
  #     # p <<- pmartRseq::applyFilt(filter_object = otu_filter_obj(),
  #     #                                omicsData = filtered_rRNA_obj(),
  #     #                                num_samps = input$filter_kOverA_sample_threshold, upper_lim = input$filter_kOverA_count_threshold)  
  #       return(p)
  #     }
  #   })
  #-------- filter history support -----------#
  #filters <- reactiveValues(otu_filt = list(), sample_filt = list())
  
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
  # end sample metadata filtering
  
 
  #--------- Object for filtering -----------#

  # filtered_rRNA_obj <- eventReactive(input$sample_filter_go,{
  #   current <<-  pmartRseq::applyFilt(filter_object = sample_filter_obj(),
  #                                     omicsData = current,
  #                                     upper_lim = input$n)
  #   print(str(current))
  #   current
  # })
  
  # observeEvent(input$otu_filter_go,{
  #   current <- isolate(filtered_rRNA_obj)
  #   filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = otu_filter_obj(),
  #                                              omicsData = current,
  #                                              num_samps = input$filter_kOverA_sample_threshold, upper_lim = input$filter_kOverA_count_threshold)
  #   print(str(filtered_rRNA_obj))
  #   return(filtered_rRNA_obj)
  # })
  # 
  # observeEvent(input$sample_filter_go,{
  #   current <- isolate(filtered_rRNA_obj)
  #   filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = sample_filter_obj(),
  #                                              omicsData =  current,
  #                                              upper_lim = input$n)
  # })
  # 
  # filtered_rRNA_obj <- reactive({
  #     current <<- pmartRseq::applyFilt(filter_object = otu_filter_obj(),
  #                                                 omicsData = current,
  #                                                 num_samps = input$filter_kOverA_sample_threshold, upper_lim = input$filter_kOverA_count_threshold)
  #     current})
  # 
  # filtered_rRNA_obj <- reactive({
  #   current <<- pmartRseq::applyFilt(filter_object = sample_filter_obj(),
  #                                                           omicsData = current, upper_lim = input$n)
  #     current})
  
  # otu_counter <- eventReactive(input$otu_filter_go,{
  #   input$otu_filter_go
  # })
  # 
  # sample_counter <- eventReactive(input$sample_filter_go,{
  #   input$sample_filter_go
  # })
  # 
  # filters$otu_filt[[otu_counter()]] <- otu_filter_obj()
  # filters$sample_filt[[sample_counter()]] <- sample_filter_obj()
  # 
  

  # eventReactive(input$otu_filter_go,{
  #   print("Happening")
  #   temp1 <- pmartRseq::applyFilt(filter_object = otu_filter_obj(),
  #                                omicsData = filtered_rRNA_obj(),
  #                                num_samps = input$filter_kOverA_sample_threshold, upper_lim = input$filter_kOverA_count_threshold)
  # }, ignoreNULL = TRUE)
  


  # observeEvent(input$sample_filter_go,{
  #   temp2 <- pmartRseq::applyFilt(filter_object = sample_filter_obj(), omicsData = filtered_rRNA_obj(), upper_lim = input$n)
  #   filtered_rRNA_obj <- reactive({
  #     return(temp2)
  #   })
  # }, ignoreNULL = TRUE)
  
  
  #------------ Apply filters on action button -----------#
  # fix this to include A
  # filtered_otu_obj <- isolate({
  #   eventReactive(input$otu_filter_go, {
  #   validate(
  #     need(input$filter_kOverA_count_threshold > 0, message = "Please use a filter threshold greater than 0")
  #   )
  #   pmartRseq::applyFilt(filter_object = otu_filter_obj(), omicsData = filtered_rRNA_obj(), upper_lim = input$filter_kOverA_count_threshold)
  # })
  # })
  # 
  # filtered_sample_obj <- eventReactive(input$sample_filter_go, {
  #     validate(
  #       need(input$filter_kOverA_count_threshold > 0, message = "Please use a filter threshold greater than 0")
  #     )
  #     pmartRseq::applyFilt(filter_object = sample_filter_obj(), omicsData = filtered_rRNA_obj(), min_num = input$n)
  #   })
  
  #if(input$otu_filter_go){
 #   filtered_rRna_obj <- reactiveValues(filtered_otu_obj())
 # }
  # filtered_rRna_obj <- eventReactive(input$sample_filter_go, {
  #   validate(
  #     need(input$n > 0, message = "Please use a filter threshold greater than 0")
  #   )
  #   pmartRseq::applyFilt(filter_object = sample_filter_obj(), omicsData = rRNAobj(), upper_lim = as.numeric(input$n))
  # })
  # filtered_rRna_obj <- reactive({
  #   #------- TODO: need to display an error if the rRNA object isn't created! ---------#
  #   validate(
  #     need(!(is.null(rRNAobj()$e_data)), message = "rRNA object fail. Check to make sure file paths are correct")
  #   )
  #   if (is.null(input$biom)) {
  #     return(NULL)
  #   }else{
  #     results <- rRNAobj()
  #     return(results)
  #   }
  # })
  
  # filtered_rRna_obj <- reactive({
  #   validate(need(!(is.null(biom_obj())), message = "biom object needed"))
  #   return(pmartRseq::as.rRNAdata(e_data = biom_obj(), f_data = metadata_obj(), edata_cname = "OTU", fdata_cname = names(metadata_obj()[1]))) 
  #   # If the metadata has been loaded, create the charts
  #   })

  
  #--------------- k over a  and library read count filtering  ---------#
  kovera_k <- 1
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
                      min = 1, max = maxSamples(), value = kovera_k, step = 1)
    })
  })
  
  otu_filter_obj <- reactive({
      return(pmartRseq::count_based_filter(rRNAobj(), fn = "ka"))
  })
  
  output$read_counts_plot <- renderPlot({
    validate(
      need( input$filter_kOverA_count_threshold >= 0, message = "Enter a count minimum >= 0"),
      need( input$filter_kOverA_sample_threshold >= 0, message = "Enter a sample minimum >= 0")
    )
    plot(otu_filter_obj(), min_num = input$filter_kOverA_count_threshold, min_samp = input$filter_kOverA_sample_threshold)
  })
  
  
  #-------------- Library read filtering -----------#
  sample_filter_obj <- reactive({
      return(pmartRseq::sample_based_filter(omicsData = rRNAobj(), fn = "sum"))
  })
  
  output$sample_counts_plot <- renderPlot({
    validate(
      need(input$n >= 0, message = "Enter a count minimum >= 0")
    )
    plot(sample_filter_obj(), min_num = input$n)
  })
  

}) #end server
