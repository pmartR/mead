
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
      need(grepl(".biom", as.character(input$biom$name)), "Biom file must be in .biom format"), errorClass = "red"
    )

    validate(
      need(input$qiime != "", "Please select a qiime file")
    )
    return(pmartRseq::as.rRNAdata(e_data = as.character(input$biom$datapath), f_data = as.character(input$qiime$datapath), edata_cname = "OTU"))
  }) #end rRNAobj

  #-------- filter history support -----------#
  filters <- reactiveValues(otu = list(), sample = list(), metadata = list())
  
  #--------- kovera observer ---------#
  observeEvent(input$otu_filter_go,{
    filters$otu[[input$otu_filter_go]] <- otu_filter_obj()
    })
  
  #--------- sample observer ---------#
  observeEvent(input$sample_filter_go, {
    filters$sample[[input$sample_filter_go]] <- sample_filter_obj()
  })
  
  #--------- metadata observer ---------#
  observeEvent(input$metadata_filter_go, {
    filters$metadata[[input$metadata_filter_go]] <- sample_metadata_filter_obj()
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
    
    # no sample metadata filter yet
    if (input$metadata_filter_go == 0) {
      filtered_rRNA_obj <<- filt1
      return(filtered_rRNA_obj)
    } 
    if (input$metadata_filter_go > 0) {
      # apply a metadata filter
     isolate({
        # get the intersection of the samples left after checkbox groups are removed
        temp <- filtered_rRNA_obj$f_data
        column_class <- sapply(temp, class)
        categorical <- temp[, which(column_class %in% c("character", "factor", "logical"))]
        # If categorical check for non-uniqueness
        non_unique_columns <- which(unlist(lapply(categorical, function(x) length(unique(x)) != nrow(categorical) & length(unique(x)) != 1)))
        if(length(non_unique_columns > 0)){
          check_boxes <- names(categorical)[non_unique_columns]
        } else {
          check_boxes <- names(categorical)
        }
        sample_names <- list()
        for (i in 1:length(check_boxes)) {
          print(check_boxes[i])
          #   temp <- pmartRseq::applyFilt(sample_metadata_filter_obj(), rRNAobj(), samps_to_remove = )
          sample_names[[i]] <- dplyr::group_by_(temp, check_boxes[i]) %>%
          dplyr::filter_(check_boxes[i] %in% input[[check_boxes[i]]]) %>%
          dplyr::select_(attr(filtered_rRNA_obj, "cnames")$fdata_cname)#find the column name and associated check box
        }
        sample_names <- Reduce(intersect, sample_names)
        
        
        filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$metadata[[input$metadata_filter_go]],
                                                   omicsData = filt1,
                                                   samps_to_remove = sample_names)
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
      temp <- filtered_rRNA_obj$f_data
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
  sample_metadata_filter_obj <- reactive({
    return(pmartRseq::sample_based_filter(omicsData = rRNAobj(), fn = "criteria"))
  })
  
  observeEvent(metadata_obj(), {
    # If the metadata has been loaded, create the charts
    output$plots <- renderUI({get_plot_output_list(metadata_obj())})
    output$boxes <- renderUI({get_checkbox_output_list(metadata_obj())})
    # Initialize new_metadata_obj object. If no metadata filtering then it's a copy of metadata_obj.
    # If metadata filtering, will be overwritten with filters 

    # If the user brushes a chart, input$selected_indices will change
    # If that change is observed, subset new_metadata_obj and 
    # create subsetted charts and table
    observeEvent(input$selected_indices, {

      new_metadata_obj <- reactive({
        return(metadata_obj()[input$selected_indices + 1, ]) # add one because javascript is zero-indexed
      })
      output$plots <- renderUI({get_plot_output_list(new_metadata_obj())})
      output$new_samples <- DT::renderDataTable(
        data.frame( new_metadata_obj()), rownames = FALSE, class = 'cell-border stripe compact hover',
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
  observeEvent(input$metadata_reset_button,{
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
    # observe({
    #   filtered_rRNA_obj
    ################# JSON object for Meg #################
       rRNA_filtered <- jsonlite::toJSON(list(filtered_rRNA_obj$e_data, filtered_rRNA_obj$e_meta, filtered_rRNA_obj$f_data))
    #   #write(rRNA_filtered, file = "rRNA_filtered.json")
    # })
    return(isolate(filtered_rRNA_obj))
  })
  

  
  # ################ Group Designation Tab #################

    output$group1 <- renderUI({
      selectInput("group1",
                  label = "Main Effect 1",
                  choices = colnames(filtered_data()$f_data))
    })
    
    
    output$group2 <- renderUI({
      selectInput("group2",
                  label = "Main Effect 2",
                  choices = c("NA",colnames(filtered_data()$f_data)),
                  selected = NULL)
      #groupDesignation
    })
    
    # output$covs <- renderUI({
    #   checkboxInput("covs",
    #                 label = "Any covariates?")
    # })
    # input$cov1 <- NA
    # input$cov2 <- NA
    #observeEvent(input$covs, {
      output$cov1 <- renderUI({
        selectInput("cov1",
                    label = "Covariate 1",
                    choices = c("NA",colnames(filtered_data()$f_data)),
                    selected = NULL)
      })
      
      output$cov2 <- renderUI({
        selectInput("cov2",
                    label = "Covariate 2",
                    choices = c("NA",colnames(filtered_data()$f_data)),
                    selected = NULL)
        #groupDesignation
      })
    #})
    

    
    groupDF <- reactive({
      if(input$group2 == "NA"){
        mainEffects <- input$group1
      }else{
        mainEffects <- c(input$group1, input$group2)
      }
      if(!is.na(input$cov1) | !is.na(input$cov2)){
        covariates <- c(input$cov1, input$cov2)
        if(any(is.na(covariates))){
          covariates <- covariates[!is.na(covariates)]
        }
      }
      validate(
        need(length(mainEffects) > 0, "There needs to be at least one grouping variable")
      )
      
      return(pmartRseq::group_designation(filtered_data(), main_effects = mainEffects, covariates = covariates))
      
    })
    
    #observeEvent(input$groupDF_go,
    output$group_DF <- DT::renderDataTable(attr(groupDF(), "group_DF"))
    #)
    
    
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
