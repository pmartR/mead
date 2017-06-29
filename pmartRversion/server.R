
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

shinyServer(function(input, output, session) {

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

  #-------- filter history support -----------#
  filters <- reactiveValues(otu = list(), sample = list())
  
  #--------- kovera observer ---------#
  observeEvent(input$otu_filter_go,{
    filters$otu[[input$otu_filter_go]] <- otu_filter_obj()
    })
  
  #--------- sample observer ---------#
  observeEvent(input$sample_filter_go, {
    filters$sample[[input$sample_filter_go]] <- sample_filter_obj()
  })
  
  #--------- filter application observer ---------#
  filtered_rRNA_obj <- NULL
  observe({
    # no kovera filter yet
    if (input$otu_filter_go == 0) {
      filt1 <- rRNAobj()
    } 
    if (input$otu_filter_go != 0) {
      # apply k over a filter
      isolate({
        filt1 <- pmartRseq::applyFilt(filter_object = filters$otu[[input$otu_filter_go]],
                                      omicsData = rRNAobj(),
                                      num_samps = input$filter_kOverA_sample_threshold,
                                      upper_lim = input$filter_kOverA_count_threshold)  
      })
    }
    
    # no sample filter yet
    if (input$sample_filter_go == 0) {
      filtered_rRNA_obj <<- filt1
      return(filtered_rRNA_obj)
    } else {
      # apply sample filter
      isolate({
        filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$sample[[input$sample_filter_go]],
                                                  omicsData = filt1,
                                                  upper_lim = input$n)
        return(filtered_rRNA_obj)
      })
    }
    
 
    
  })

  #----------------- observe resets ----------#
  observeEvent(input$otu_reset_button, {
    kovera_k <- 1
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
    filt1 <- rRNAobj()
    # no sample filter yet
    if (input$sample_filter_go == 0) {
      filtered_rRNA_obj <<- filt1
      return(filtered_rRNA_obj)
    } else {
      # apply sample filter
      isolate({
        filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$sample[[input$sample_filter_go]],
                                                   omicsData = filt1,
                                                   upper_lim = input$n)
        return(filtered_rRNA_obj)
      })
    }
  })
  
  observeEvent(input$sample_reset_button, {
    updateNumericInput(session, "n",
                       value = 0)
    if (input$otu_filter_go == 0) {
      filt1 <- rRNAobj()
    } 
    if (input$otu_filter_go != 0) {
      # apply k over a filter
      isolate({
        filt1 <- pmartRseq::applyFilt(filter_object = filters$otu[[input$otu_filter_go]],
                                      omicsData = rRNAobj(),
                                      num_samps = input$filter_kOverA_sample_threshold,
                                      upper_lim = input$filter_kOverA_count_threshold)  
      })
    }
    
    # no sample filter yet
      filtered_rRNA_obj <<- filt1
      return(filtered_rRNA_obj)
  })
  
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
  # 
  observeEvent(metadata_obj(), {
    # If the metadata has been loaded, create the charts
    output$plots <- renderUI({get_plot_output_list(metadata_obj())})
    #output$boxes <- renderUI({get_checkbox_output_list(metadata_obj())})
    
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
  
  #--------------- k over a filtering  -----------------#
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
    validate(
      need(length(filtered_rRNA_obj) > 0 , message = "Upload data first")
    )
    return(pmartRseq::count_based_filter(rRNAobj(), fn = "ka"))
  })
  
  output$read_counts_plot <- renderPlot({
    validate(
      need( input$filter_kOverA_count_threshold >= 0, message = "Enter a count minimum >= 0"),
      need( input$filter_kOverA_sample_threshold >= 0, message = "Enter a sample minimum >= 0")
    )
    if (input$sample_filter_go == 0 & input$otu_filter_go == 0) {
      plot(otu_filter_obj(), min_num = input$filter_kOverA_count_threshold, min_samp = input$filter_kOverA_sample_threshold)
    } else{
      otu_filt_obj <- pmartRseq::count_based_filter(filtered_rRNA_obj, fn = "ka")
      plot(otu_filt_obj, min_num = input$filter_kOverA_count_threshold, min_samp = input$filter_kOverA_sample_threshold)
      
    }
  })
  
  
  #-------------- Library read filtering -----------#

  sample_filter_obj <- reactive({
    return(pmartRseq::sample_based_filter(omicsData = rRNAobj(), fn = "sum"))
  })
  

  output$sample_counts_plot <- renderPlot({
    validate(
      need(input$n >= 0, message = "Enter a count minimum >= 0")
    )
    if (input$sample_filter_go == 0 & input$otu_filter_go == 0) {
      plot(sample_filter_obj(), min_num = input$n)
    } else {
      sample_filt_obj <- pmartRseq::sample_based_filter(omicsData = filtered_rRNA_obj, fn = "sum")
      plot(sample_filt_obj, min_num = input$n)
    }
  })
  
  #------------ reactive filtered data for downstream processing --------------#

  filtered_data <- reactive({
    observe({
      filtered_rRNA_obj
      rRNA_filtered <- jsonlite::toJSON(list(filtered_rRNA_obj$e_data, filtered_rRNA_obj$e_meta, filtered_rRNA_obj$f_data))
      write(rRNA_filtered, file = "rRNA_filtered.json")
    })
    return(isolate(filtered_rRNA_obj))
  })
  

  
  ################ Group Designation Tab #################
  output$mainEffects <- renderUI({
    checkboxGroupInput("mainEffects",
                label = "Grouping Main Effects",
                choices = as.list(colnames(filtered_data()$f_data)),
                selected = "NULL")
  })
  
  groupDF <- reactive({
    validate(
      need(length(input$mainEffects) <= 2, "There can only be a maximum of 2 grouping variables")
    )
    
    validate(
      need(length(input$mainEffects) > 0, "There needs to be at least one grouping variable")
    )
    
    # if (input$groupDF_reset_button != 0) {
    #   input$groupDF_go = 0
    #   attr(rRNAobj(), "group_DF") <- NULL
    # }
    
    if (input$groupDF_go == 0) {
      return()
    }else{
      return(pmartRseq::group_designation(filtered_data(), main_effects = input$mainEffects))
    }
  })
  
  # output$group_DF <- renderTable({
  #   attr(groupDF(), "group_DF")
  # })
  
  observeEvent(input$groupDF_go, 
               output$group_DF <- DT::renderDataTable(
                 attr(groupDF(), "group_DF")
                 ))
  
  ################ Community Metrics Tab #################
  
  output$xaxis <- renderUI({
    selectInput("xaxis",
                label = "x-axis",
                choices = c(colnames(attr(groupDF(),"group_DF"))),
                selected = "Group")
  })
  
  output$color <- renderUI({
    selectInput("color",
                label = "color",
                choices = c(colnames(attr(groupDF(),"group_DF"))),
                selected = "Group")
  })
  
  #----------- alpha diversity example ----------#
  
  output$adiv_index <- renderUI({
    checkboxGroupInput("adiv_index",
                       label = "Alpha Diversity Index",
                       choices = list("Shannon"="shannon","Simpson"="simpson","InverseSimpson"="invsimpson"),
                       selected = c("shannon","simpson","invsimpson"))
  })

  
  
  a_div <- reactive({
    validate(
      need(length(input$adiv_index) > 0, "There needs to be at least one alpha diversity index")
    )
    
    return(pmartRseq::alphaDiv_calc(groupDF(), index=input$adiv_index))
  })
  
  output$adiv_plot <- renderPlot({
    plot(a_div(), x_axis=input$xaxis, color=input$color)
  })

               
  #----------- richness example ----------#
  output$rich_index <- renderUI({
    checkboxGroupInput("rich_index",
                        label = "Richness Index",
                        choices = list("Observed"="observed","Chao1"="chao1","ACE"="ace","Breakaway"="break"),
                        selected = c("observed","chao1","ace"))
   })
   
  rich <- reactive({
    validate(
      need(length(input$rich_index) > 0, "There needs to be at least one richness index")
    )
     
    return(pmartRseq::richness_calc(groupDF(), index=input$rich_index))
  })
   
  output$rich_plot <- renderPlot({
    plot(rich(), x_axis=input$xaxis, color=input$color)
  })
  
}) #end server
