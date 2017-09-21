
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(lazyeval)
library(dplyr)
library(DT)
library(filterWidget)
library(ggplot2)
library(pmartRseq)
#library(phyloseq)
library(vegan)
library(goeveg)
source("./functions/helper_functions.R")
source("./functions/test_functions.R")
#source("./Report/R/report.R")

#Sys.setenv(R_ZIPCMD="/usr/bin/zip")

filtered_rRNA_obj <- list()

shinyServer(function(input, output, session) {
  
  #---------- rRNA object import -------------#
  rRNAobj <- reactive({
    validate(
      need(input$biom != "", "Please select a biom file")
    )
    validate(
      need(grepl(".biom", as.character(input$biom$name)), "Biom file must be in .biom format")
    )
    
    validate(
      need(input$qiime != "", "Please select a qiime file")
    )
    return(pmartRseq::as.seqData(e_data = as.character(input$biom$datapath), f_data = as.character(input$qiime$datapath), edata_cname = "OTU", data_type = 'rRNA'))
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
      #return(filtered_rRNA_obj)
    } 
    if (input$sample_filter_go != 0) {
      # apply sample filter
       #browser()
      #isolate({
      filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$sample[[input$sample_filter_go]],
                                                 omicsData = filt1,
                                                 upper_lim = input$n)
      filt1 <- filtered_rRNA_obj
      # return(filtered_rRNA_obj)
      # })
    }
    
    # no sample metadata filter yet
    if (input$metadata_filter_go == 0) {
      filtered_rRNA_obj <<- filt1
      # return(filtered_rRNA_obj)
    } 
    
    if (input$metadata_filter_go != 0) {
      
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
          # build in a check for no samples to remove
          # if (length(check_boxes[i]) == 0){
          #   return(filtered_rRNA_obj)
          # }
          sample_names[[i]] <- dplyr::filter_(temp, interp(~v%in%input[[check_boxes[i]]], v=as.name(check_boxes[i]))) %>%
            dplyr::select_(attr(filtered_rRNA_obj, "cnames")$fdata_cname)#find the column name and associated check box
        }
        sample_names <- Reduce(intersect, sample_names)
        to_remove <- temp[ which(!(as.character(temp[, attr(filtered_rRNA_obj, "cnames")$fdata_cname]) %in% as.character(unlist(sample_names)))),
                           attr(filtered_rRNA_obj, "cnames")$fdata_cname]
        
        
        filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$metadata[[input$metadata_filter_go]],
                                                   omicsData = filt1,
                                                   samps_to_remove = to_remove)
        
        # return(filtered_rRNA_obj)
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
    }
    if (input$sample_filter_go != 0) {
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
    input$metadata_filter_go
    input$metadata_reset_button
    if (is.null(input$qiime)) {
      return(NULL)
    }else{
      temp <- filtered_rRNA_obj$f_data
      inds <- lapply(temp, function(x) sum(is.na(x)) == length(x) )
      results <- temp[,!(unlist(inds))]
      return(results)
    }
  })
  
  
  output$sample_metadata <- DT::renderDataTable(expr = 
    data.frame(metadata_obj()), rownames = FALSE, class = 'cell-border stripe compact hover',
    options = list(columnDefs = list(list(
      targets = c(1:(ncol((metadata_obj())) - 1)),
      render = JS(
        "function(data, type, row, meta) {",
        "return type === 'display' && data.length > 10 ?",
        "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
        "}")
    ))), callback = JS('table.page(3).draw(false);'))
  
  
  #---- Charts for f_meta filtering -------#
  # 
  sample_metadata_filter_obj <- reactive({
    return(pmartRseq::sample_based_filter(omicsData = rRNAobj(), fn = "criteria"))
  })
  observeEvent(input$metadata_reset_button, {
    output$sample_metadata <- DT::renderDataTable(
      rRNAobj()$f_data, rownames = FALSE, class = 'cell-border stripe compact hover',
      options = list(columnDefs = list(list(
        targets = c(1:(ncol((metadata_obj())) - 1)),
        render = JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 10 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
          "}")
      ))), callback = JS('table.page(3).draw(false);'))  })
  
  observeEvent(rRNAobj(), {
    # If the metadata has been loaded, create the charts
    output$plots <- renderUI({get_plot_output_list(metadata_obj())})
    output$boxes <- renderUI({get_checkbox_output_list(rRNAobj()$f_data)})
    # Initialize new_metadata_obj object. If no metadata filtering then it's a copy of metadata_obj.
    # If metadata filtering, will be overwritten with filters 
    
    # If the user brushes a chart, input$selected_indices will change
    # If that change is observed, subset new_metadata_obj and 
    # create subsetted charts and table
    observeEvent(input$selected_indices, {
      new_metadata_obj <- reactive({
        print(input$selected_indices)
        
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
  observeEvent(input$metadata_reset_button, {
    output$boxes <- renderUI({get_checkbox_output_list(rRNAobj()$f_data)})
  })
  
  # end sample metadata filtering
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
    input$metadata_filter_go
    # observe({
    #   filtered_rRNA_obj
    ################# JSON object for Meg #################
    rRNA_filtered <- jsonlite::toJSON(list(filtered_rRNA_obj$e_data, filtered_rRNA_obj$e_meta, filtered_rRNA_obj$f_data))
    #   #write(rRNA_filtered, file = "rRNA_filtered.json")
    # })
    return(isolate(filtered_rRNA_obj))
  })
  
  
  # ################ Group Designation Tab #################
  
  output$summ_filt <- renderPrint({
    #browser()
    summary(filtered_data())
  })
  
  output$nrow_edata <- renderPrint({
    nrow(filtered_data()$e_data)
  })
  
  
  group_vars <- reactive({
    intersect(which(lapply(apply(filtered_data()$f_data, 2, function(z) table(z))[unlist(lapply(apply(filtered_data()$f_data, 2, function(x) table(x)), function(y) any(is.finite(y))))], function(w) max(w)) > 2), which(apply(filtered_data()$f_data, 2, function(v) length(unique(v))) > 2))
    
  })
  
  # Input main effects used for groupings
  output$group1 <- renderUI({
    selectInput("group1",
                label = "Main Effect 1",
                choices = colnames(filtered_data()$f_data)[group_vars()])
  })
  
  # Can have up to 2 main effects
  output$group2 <- renderUI({
    selectInput("group2",
                label = "Main Effect 2",
                choices = c("NA",colnames(filtered_data()$f_data)[group_vars()]),
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
  #   output$cov1 <- renderUI({
  #     selectInput("cov1",
  #                 label = "Covariate 1",
  #                 choices = c("NA",colnames(filtered_data()$f_data)),
  #                 selected = NULL)
  #   })
  #   
  #   output$cov2 <- renderUI({
  #     selectInput("cov2",
  #                 label = "Covariate 2",
  #                 choices = c("NA",colnames(filtered_data()$f_data)),
  #                 selected = NULL)
  #     #groupDesignation
  #   })
  # #})
  
  
  # Create groups with main effects
  groupDF <- reactive({
    if(input$group2 == "NA"){
      mainEffects <- input$group1
    }else{
      mainEffects <- c(input$group1, input$group2)
    }
    # if(!is.na(input$cov1) | !is.na(input$cov2)){
    #   covariates <- c(input$cov1, input$cov2)
    #   if(any(is.na(covariates))){
    #     covariates <- covariates[!is.na(covariates)]
    #   }
    # }
    validate(
      need(length(mainEffects) > 0, "There needs to be at least one grouping variable")
    )
    
    return(pmartRseq::group_designation(filtered_data(), main_effects = mainEffects))
    
  })
  
  # Show the groupings data frame
  #observeEvent(input$groupDF_go,
  group_df_tab <- reactive({
    attr(groupDF(), "group_DF")
  })
  output$group_DF <- DT::renderDataTable(group_df_tab())
  #)
  
  # Also output a table showing the number of reps in each group
  group_freq_tab <- reactive({
    as.data.frame(table(attr(groupDF(),"group_DF")$Group))
  })
  output$group_tab <- DT::renderDataTable(group_freq_tab())
  
  
  # ################ Outliers Tab #################
  
  outlier_jaccard <- reactive({
    pmartRseq::jaccard_calc(omicsData = groupDF())
  })
  
  jac_plot_obj <- reactive({
    plot(outlier_jaccard())
  })
  
  output$jac_plot <- renderPlot({
    #plot(outlier_jaccard())
    print(jac_plot_obj())
  })
  
  
  # ################ Normalization Tab #################
  
  # Select which normalization function to use
  output$normFunc <- renderUI({
    selectInput("normFunc",
                label = "Normalization Function",
                choices = c("percentile","tss","rarefy","poisson","deseq","tmm","css","none"),
                selected = "css")
  })
  
  # Create normalized data
  normalized_data <- reactive({
    validate(need(length(input$normFunc) == 1, "Need to specify a normalization function."))
    validate(need(input$normFunc %in% c("percentile","tss","rarefy","poisson","deseq","tmm","css","none"), "Normalization function must be one of the options specified."))
    
    return(pmartRseq::split_emeta(pmartRseq::normalize_data(omicsData=groupDF(), norm_fn=input$normFunc, normalize=TRUE), cname="OTU", split1=NULL, numcol=7, split2="__", num=2, newnames=NULL))
  })
  
  # Look at normalized results
  output$normData <- DT::renderDataTable(normalized_data()$e_data, rownames = FALSE)
  
  output$norm_class <- renderUI({
    selectInput("norm_class",
                label = "Taxonomic level to plot",
                choices = colnames(normalized_data()$e_meta)[-which(colnames(normalized_data()$e_meta)==attr(normalized_data(),"cnames")$edata_cname)])
  })
  
  norm_plot_obj <- reactive({
    plot(normalized_data(), class=input$norm_class)
  })
  # Try to make a stacked bar plot - not working right now
  output$norm_plot <- renderPlot({
    #browser()
    #plot(normalized_data(), class="Phylum")
    print(norm_plot_obj())
  })
  
  # Calculate abundance on normalized data
  abun_norm <- reactive({
    return(suppressWarnings(pmartRseq::abundance_calc(normalized_data())))
  })
  
  # Calculate abundance on raw data
  abun_raw <- reactive({
    return(suppressWarnings(pmartRseq::abundance_calc(groupDF())))
  })
  
  # Calculate richness on normalized data
  rich_norm <- reactive({
    return(suppressWarnings(pmartRseq::richness_calc(normalized_data(), index="observed")))
  })
  
  # Calculate richness on raw data
  rich_raw <- reactive({
    return(suppressWarnings(pmartRseq::richness_calc(groupDF(), index="observed")))
  })
  
  ra_raw_plot <- reactive({
    plot(abun_raw(), rich_raw(), plot_title="Raw Data")
  })
  # Create a plot of raw abundance vs raw richness
  output$ra_raw <- renderPlot({
    #plot(abun_raw(), rich_raw(), plot_title="Raw Data")
    print(ra_raw_plot())
  })
  
  ra_norm_plot <- reactive({
    plot(abun_norm(), rich_norm(), plot_title="Normalized Data")
  })
  # Create a plot of normalized abundance vs normalized richness to see if there is a reduction in correlation
  output$ra_norm <- renderPlot({
    #plot(abun_norm(), rich_norm(), plot_title="Normalized Data")
    print(ra_norm_plot())
  })
  
  ################ Community Metrics Tab #################
  
  # Select what variable to put on x-axis in community metrics plots
  output$xaxis <- renderUI({
    selectInput("xaxis",
                label = "x-axis",
                choices = c(colnames(attr(groupDF(),"group_DF"))),
                selected = "Group")
  })
  
  # Select what variable to color by in community metrics plots
  output$color <- renderUI({
    selectInput("color",
                label = "color",
                choices = c(colnames(attr(groupDF(),"group_DF"))),
                selected = "Group")
  })
  
  #----------- alpha diversity example ----------#
  
  # Which alpha diversity indices to calculate
  output$adiv_index <- renderUI({
    checkboxGroupInput("adiv_index",
                       label = "Alpha Diversity Index",
                       choices = list("Shannon"="shannon","Simpson"="simpson","InverseSimpson"="invsimpson"),
                       selected = c("shannon","simpson","invsimpson"))
  })
  
  # Calculate alpha diversity
  a_div <- reactive({
    validate(
      need(length(input$adiv_index) > 0, "There needs to be at least one alpha diversity index")
    )
    
    return(pmartRseq::alphaDiv_calc(groupDF(), index=input$adiv_index))
  })
  
  adiv_plot_obj <- reactive({
    plot(a_div(), x_axis=input$xaxis, color=input$color)
  })
  # Show alpha diversity plot
  output$adiv_plot <- renderPlot({
    #plot(a_div(), x_axis=input$xaxis, color=input$color)
    print(adiv_plot_obj())
  })
  
  output$adiv_summary <- renderPrint({
    summary(a_div())
  })
  
  
  #----------- richness example ----------#
  
  # Which richness indices to calculate
  output$rich_index <- renderUI({
    checkboxGroupInput("rich_index",
                       label = "Richness Index",
                       choices = list("Observed"="observed","Chao1"="chao1","ACE"="ace","Breakaway"="break"),
                       selected = c("observed","chao1","ace"))
  })
  
  # Calculate richness
  rich <- reactive({
    validate(
      need(length(input$rich_index) > 0, "There needs to be at least one richness index")
    )
    
    return(pmartRseq::richness_calc(groupDF(), index=input$rich_index))
  })
  
  rich_plot_obj <- reactive({
    plot(rich(), x_axis=input$xaxis, color=input$color)
  })
  # Show richness plot
  output$rich_plot <- renderPlot({
    #plot(rich(), x_axis=input$xaxis, color=input$color)
    print(rich_plot_obj())
  })
  
  output$rich_summary <- renderPrint({
    summary(rich())
  })
  
  
  ################ Beta Diversity Tab #################
  #----------- phyloseq beta diversity ----------#
  # output$beta_index <- renderUI({
  #   selectInput("beta_index",
  #                      label = "Beta Diversity Index",
  #                      choices = unname(unlist(phyloseq::distanceMethodList)),
  #                      selected = "bray")
  # })
  # 
  # output$ord_method <- renderUI({
  #   selectInput("ord_method",
  #                 label = "Ordination Method",
  #                 choices = list("DCA"="DCA","CCA"="CCA","RDA"="RDA","CAP"="CAP","DPCoA"="DPCoA","NMDS"="NMDS","MDS"="MDS","PCoA"="PCoA"),
  #                 selected = "DCA")
  # })
  # 
  # output$ord_color <- renderUI({
  #   selectInput("ord_color",
  #               label = "Color variable for ordination plot",
  #               choices = c("NA",colnames(attr(groupDF(),"group_DF"))[-c(which(colnames(attr(groupDF(),"group_DF"))=="Group"), which(colnames(attr(groupDF(),"group_DF"))==attr(groupDF(),"cnames")$fdata_cname))]),
  #               selected = NULL)
  # })
  # 
  # phylo <- reactive({
  #   return(mintR_to_phyloseq(normalized_data()))
  # })
  # 
  # beta <- reactive({
  #   validate(
  #     need(length(input$beta_index) > 0, "There needs to be at least one beta diversity index")
  #   )
  #   return(phyloseq::distance(physeq = phylo(), method = input$beta_index))
  # })
  # 
  # #output$beta <- DT::renderDataTable(beta())
  # 
  # mydist <- reactive({
  #   validate(
  #     need(length(input$beta_index) > 0, "There needs to be at least one beta diversity index")
  #   )
  #   validate(
  #     need(length(input$ord_method) > 0, "There needs to be one ordination method")
  #   )
  #   
  #   return(phyloseq::ordinate(physeq = phylo(), method = input$ord_method, distance = beta()))
  # })
  # 
  # #output$mydist <- DT::renderDataTable(mydist())
  # 
  # output$ord_plot <- renderPlot({
  #   validate(
  #     need(length(input$ord_color) > 0, "There needs to be one factor to color by in the ordination plot")
  #   )
  #   phyloseq::plot_ordination(physeq = phylo(), ordination = mydist(), color = input$ord_color)
  # })
  
  #----------- vegan beta diversity ----------#
  
  # Select which beta diversity index to calculate
  output$beta_index <- renderUI({
    selectInput("beta_index",
                label = "Beta Diversity Index",
                choices = list("manhattan","euclidean","canberra","bray","kulczynski","jaccard","gower","altGower","morisita","horn","mountford","raup","binomial","chao","cao","mahalanobis"),
                selected = "bray")
  })
  
  # Translate seqData to something that vegan can use
  vegdata <- reactive({
    return(pmartRseq::pmartRseq_to_vegan(normalized_data()))
  })
  
  observeEvent(input$submit_goe, {
    dimcheck_obj <<- reactive({
      goeveg::dimcheckMDS(vegdata(), distance = input$beta_index, autotransform = FALSE)
    })
    output$dimcheck <- renderPlot({
      #goeveg::dimcheckMDS(vegdata(), distance = input$beta_index, autotransform = FALSE)
      print(dimcheck_obj())
    })
  })
  
  # output$ord_method <- renderUI({
  #   selectInput("ord_method",
  #                 label = "Ordination Method",
  #                 choices = list("NMDS","PCA"),
  #                 selected = "NMDS")
  # })
  
  # Select the number of dimensions to use 
  output$k <- renderUI({
    numericInput("k",
                 label = "Number of dimensions",
                 value = 4)
  })
  
  # Select what variable to color by
  output$ord_colors <- renderUI({
    selectInput("ord_colors",
                label = "Ordination Group Colors",
                choices = colnames(attr(normalized_data(),"group_DF")),
                selected = "Group")
  })
  
  observeEvent(input$submit_ord, {
    
    output$ellipses <- renderUI({
      checkboxInput("ellipses",
                    label = "NMDS Ellipses",
                    value = TRUE)
    })
    
    output$ord_x <- renderUI({
      selectInput("ord_x",
                  label = "NMDS x-axis",
                  choices = paste("NMDS",seq(1,input$k,1),sep=""),
                  selected = "NMDS1")
    })
    
    output$ord_y <- renderUI({
      selectInput("ord_y",
                  label = "NMDS y-axis",
                  choices = paste("NMDS",seq(1,input$k,1),sep=""),
                  selected = "NMDS2")
    })
    
    
    # Use vegan to calculate scores for beta diversity index
    vegmds <- reactive({
      validate(
        need(length(input$beta_index) == 1, "There needs to be one beta diversity index.")
      )
      validate(
        need(input$k >= 1, "The dimension values needs to be greater than 0.")
      )
      
      return(vegan::metaMDS(vegdata(), distance = input$beta_index, k = input$k, autotransform = FALSE))
    })
    
    ord_plot_obj <<- reactive({
      pmartRseq::pmartRseq_NMDS(res = vegmds(), omicsData = normalized_data(), grp = input$ord_colors, k = input$k, 
                                x_axis = input$ord_x, y_axis = input$ord_y, ellipses=input$ellipses)
    })
    
    # Plot showing beta diversity
    output$ord_plot <- renderPlot({
      #if(input$ord_method == "NMDS"){
      # pmartRseq::pmartRseq_NMDS(res = vegmds(), 
      #           grp = as.factor(attr(normalized_data(),"group_DF")[match(rownames(vegdata()), attr(normalized_data(),"group_DF")[,attr(normalized_data(),"cnames")$fdata_cname]),input$ord_colors]),ellipses=input$ellipses)
      #pmartRseq::pmartRseq_NMDS(res = vegmds(), omicsData = normalized_data(), grp = input$ord_colors, k = input$k, 
      #                         x_axis = input$ord_x, y_axis = input$ord_y, ellipses=input$ellipses)
      print(ord_plot_obj())
      # }else if(input$ord_method == "PCA"){
      #   mead_PCA(XX = vegmds(),
      #            ZZ = as.factor(attr(normalized_data(),"group_DF")[match(rownames(vegdata()), attr(normalized_data(),"group_DF")[,attr(normalized_data(),"cnames")$fdata_cname]),input$ord_colors]))
      # }
    })
  })
  
  
  ################ Differential Abundance Tab #################
  #----------- differential abundance ----------#
  
  # Select which differential abundance test to use
  output$da_index <- renderUI({
    selectInput("da_index",
                label = "Differential Abundance Test",
                choices = list("DESeq2 Wald Test"="dw",
                               "DESeq2 Likelihood Ratio"="dl",
                               "EdgeR Likelihood Ratio Test"="el",
                               "EdgeR with QCML Test"="eq",
                               "EdgeR with QL F-Test"="ef"),
                selected = "dw")
  })
  
  # Select which p-value adjustment method to use
  output$pval_adjust <- renderUI({
    selectInput("pval_adjust",
                label = "P-Value Adjustment Method",
                choices = p.adjust.methods,
                selected = "none")
  })
  
  output$pval_thresh <- renderUI({
    sliderInput("pval_thresh",
                label = "P-value significance threshold",
                min = 0,
                max = 1,
                step = 0.01,
                value = 0.05)
  })
  
  output$comparisons <- renderUI({
    checkboxGroupInput("comparisons",
                       label = "Differential abundance pairwise comparisons",
                       choices = lapply(c(1:(factorial(length(unique(attr(groupDF(),"group_DF")$Group)))/(factorial(2)*factorial(length(unique(attr(groupDF(),"group_DF")$Group))-2)))), function(x) paste(combn(unique(attr(groupDF(),"group_DF")$Group),2)[,x], collapse="  VS  ")),
                       selected = sapply(c(1:(factorial(length(unique(attr(groupDF(),"group_DF")$Group)))/(factorial(2)*factorial(length(unique(attr(groupDF(),"group_DF")$Group))-2)))), function(x) paste(combn(unique(attr(groupDF(),"group_DF")$Group),2)[,x], collapse="  VS  ")))
  })
  
  observeEvent(input$submit_da, {
    # Calculate normalization factors to use in differential abundance test - will use the same that was used on normalization tab
    norm_factors <- reactive({
      validate(need(length(input$normFunc) == 1, "Need to specify a normalization function."))
      validate(need(input$normFunc %in% c("percentile","tss","rarefy","poisson","deseq","tmm","css","none"), "Normalization function must be one of the options specified."))
      
      return(pmartRseq::normalize_data(omicsData=groupDF(), norm_fn=input$normFunc, normalize=FALSE))
    })
    
    comps <- reactive({
      tmp1 <- lapply(input$comparisons, function(x) strsplit(x, "  VS  ")[[1]][1])
      tmp1 <- do.call(cbind, tmp1)
      tmp2 <- lapply(input$comparisons, function(x) strsplit(x, "  VS  ")[[1]][2])
      tmp2 <- do.call(cbind, tmp2)
      tmp <- lapply(c(1:length(tmp1)), function(x) c(tmp1[x],tmp2[x]))
      return(tmp)
    })
    
    # Perform differential abundance analysis
    diffabun_res <<- reactive({
      validate(need(length(input$da_index) == 1, "Need to specify a differential abundance test"))
      validate(need(length(input$pval_adjust) == 1, "Need to specify a p-value adjustment method"))
      
      return(pmartRseq::countSTAT(omicsData = groupDF(), norm_factors = norm_factors()$scale_param, comparisons = comps(), control = NULL, test = input$da_index, pval_adjust = input$pval_adjust, pval_thresh = 0.05))
      
    })
    
    # Look at the results - this is hard to look at, maybe remove?
    output$da_res <- DT::renderDataTable(diffabun_res()$allResults)
    
    output$da_summary <- renderPrint({
      summary(diffabun_res())
    })
    
    da_flag_plot_obj <<- reactive({
      plot(diffabun_res(), type = "flag")
    })
    # Plot showing number differentially abundant in each comparison and direction of change
    output$flag_plot <- renderPlot({
      #plot(diffabun_res(), type = "flag")
      print(da_flag_plot_obj())
    })
    
    da_logfc_plot_obj <<- reactive({
      plot(diffabun_res(), type = "logfc")
    })
    # Heatmap showing the log2foldchanges of differentially abundant features
    output$logfc_plot <- renderPlot({
      #plot(diffabun_res(), type = "logfc")
      print(da_logfc_plot_obj())
    })
    
    plot_all_da_obj <<- reactive({
      pmartRseq::plot_all_diffabun(countSTAT_results = diffabun_res(), omicsData = normalized_data(), x_axis = "taxonomy2", x_lab = "Phylum")
    })
    # Plot showing log fold changes and p-values of all features, grouped by taxonomy
    output$plot_all_da <- renderPlot({
      #pmartRseq::plot_all_diffabun(countSTAT_results = diffabun_res(), omicsData = normalized_data(), x_axis = "taxonomy2", x_lab = "Phylum")
      print(plot_all_da_obj())
    })
  }, autoDestroy = FALSE)
  
  ################ Indicator Species Tab #################
  #----------- indicator species ----------#
  output$within <- renderUI({
    selectInput("within",
                label = "Perform indicator species analysis between groups within a variable",
                choices = c("NA",colnames(attr(normalized_data(),"group_DF"))[-which(colnames(attr(normalized_data(),"group_DF")) %in% c("Group",attr(normalized_data(),"cnames")$fdata_cname))]),
                selected = "NA")
  })
  
  output$is_pval_thresh <- renderUI({
    sliderInput("is_pval_thresh",
                label = "P-value significance threshold",
                min = 0,
                max = 1,
                step = 0.01,
                value = 0.05)
  })
  
  observeEvent(input$submit_is, {
    indsp_res <<- reactive({
      if(input$within == "NA"){
        return(pmartRseq::indsp_calc(omicsData = normalized_data(), within = NULL, pval_thresh = input$is_pval_thresh))
      }else{
        return(pmartRseq::indsp_calc(omicsData = normalized_data(), within = input$within, pval_thresh = input$is_pval_thresh))
      }
    })
    
    output$indsp_results <- DT::renderDataTable(indsp_res())
    
    output$indsp_summary <- renderPrint({
      summary(indsp_res())
    })
    
    
    output$indsp_xaxis <- renderUI({
      selectInput("indsp_xaxis",
                  label = "X-Axis for Indicator Species Plot",
                  choices = colnames(attr(normalized_data(),"group_DF"))[-which(colnames(attr(normalized_data(),"group_DF"))==attr(normalized_data(),"cnames")$fdata_cname)],
                  selected = "Group")
    })
    
    output$indsp_group <- renderUI({
      selectInput("indsp_group",
                  label = "Fill Variable for Indicator Species Plot",
                  choices = colnames(normalized_data()$e_meta),
                  selected = colnames(normalized_data()$e_meta)[3])
    })
    
    indsp_plot_obj <<- reactive({
      pmartRseq::plot_indsp(indsp = indsp_res(), omicsData = normalized_data(), x_axis = input$indsp_xaxis, group = input$indsp_group)
    })
    output$indsp_plot <- renderPlot({
      #pmartRseq::plot_indsp(indsp = indsp_res(), omicsData = normalized_data(), x_axis = input$indsp_xaxis, group = input$indsp_group)
      print(indsp_plot_obj())
    })
  }, autoDestroy = FALSE)
  
  
  ################ Stats Results Tab #################
  #if(exists(indsp_res()) & exists(diffabun_res())){
  #----------- differential abundance ----------#
  diffres <- reactive({
    t1 <- diffabun_res()$allResults[,grep("Flag",colnames(diffabun_res()$allResults))]
    if(is.null(ncol(t1))){
      t1 <- data.frame(OTU=rownames(diffabun_res()$allResults), t1)
      rownames(t1) <- t1$OTU
      colnames(t1)[1] <- attr(normalized_data(), "cnames")$edata_cname
      colnames(t1)[2] <- colnames(diffabun_res()$allResults)[grep("Flag",colnames(diffabun_res()$allResults))]
    }else{
      t1 <- t1[-which(rowSums(abs(t1), na.rm=TRUE) == 0),]
    }
    return(t1)
  })
  
  #----------- indicator species ----------#
  isres <- reactive({
    t1 <- indsp_res()[,grep("Flag", colnames(indsp_res()))]
    if(is.null(ncol(t1))){
      t1 <- data.frame(OTU=rownames(indsp_res()), t1)
      rownames(t1) <- t1$OTU
      colnames(t1)[1] <- attr(normalized_data(), "cnames")$edata_cname
      colnames(t1)[2] <- colnames(indsp_res())[grep("Flag",colnames(indsp_res()))]
    }else{
      t1 <- t1[-which(rowSums(abs(t1), na.rm=TRUE) == 0),]
    }
    return(t1)
  })
  
  #----------- combine results ----------#
  statsres <- reactive({
    if(is.null(attr(indsp_res(), "within"))){
      res <- intersect(rownames(diffres()), rownames(isres()))
      res <- data.frame(res)
      colnames(res)[1] <- attr(normalized_data(), "cnames")$edata_cname
      return(res)
    }else{
      res <- lapply(unique(normalized_data()$f_data[,attr(indsp_res(), "within")]), function(x){
        da <- diffres()[,grep(x, colnames(diffres()))]
        da <- da[-which(rowSums(abs(da)) == 0),]
        is <- isres()[,grep(x, colnames(isres()))]
        is <- data.frame(OTU=rownames(isres()),is)
        rownames(is) <- rownames(isres())
        is <- is[which(abs(is[,-1]) > 0),]
        res <- intersect(rownames(da), rownames(is))
        res <- data.frame(Within=x, OTU=res)
        return(res)
      })
      res <- do.call(rbind, res)
      colnames(res)[which(colnames(res) == "Within")] <- attr(indsp_res(), "within")
      colnames(res)[which(colnames(res) == "OTU")] <- attr(normalized_data(), "cnames")$edata_cname
      return(res)
    }
  })
  
  taxares <- reactive({
    res <- merge(statsres(), normalized_data()$e_meta, by=attr(normalized_data(), "cnames")$edata_cname)
    return(res)
  })
  
  output$stats_res <- DT::renderDataTable(taxares())
  output$newisplot <- renderPlot({
    pmartRseq::plot_indsp(indsp = indsp_res(), omicsData = normalized_data(), x_axis = input$indsp_xaxis, group = input$indsp_group)
  })
  output$newdaplot <- renderPlot({
    pmartRseq::plot_all_diffabun(countSTAT_results = diffabun_res(), omicsData = normalized_data(), x_axis = "taxonomy2", x_lab = "Phylum")
  })
  #}
  
  ################ Download Tab #################
  output$files_to_download <- renderUI({
    checkboxGroupInput("files_to_download",
                       label = "Select which datasets to download",
                       choices = list("raw"="raw","filtered"="filtered","normalized"="normalized",
                                      "diffabun"="diffabun"),
                       selected = c("raw","filtered","normalized","diffabun"))
  })
  
  output$files_to_download <- renderUI({
    checkboxGroupInput("files_to_download",
                       label = "Select which datasets to download",
                       choices = list("raw"="raw","filtered"="filtered","groupings"="groupings",
                                      "outliers"="outliers","normalized"="normalized",
                                      "alphadiv"="alphadiv","richness"="richness","ordination"="ordination",
                                      "diffabun"="diffabun","indicspec"="indicspec",
                                      "combined"="combined","REPORT"="report"),
                       selected = c("raw","filtered","groupings",
                                    "outliers","normalized",
                                    "alphadiv","richness","ordination",
                                    "diffabun","indicspec","combined","report"))
  })
  
  output$downloadData <- downloadHandler(
    filename = "mead_data_analysis.zip",
    content = function(fname){
      if("report" %in% input$files_to_download){
        tempReport <- file.path(tempdir(), "seqData_Report.Rmd")
        file.copy("seqData_Report.Rmd", tempReport, overwrite = TRUE)
      }
      tmpdir <- tempdir()
      setwd(tempdir())
      print(tempdir())
      
      
      
      fs <- vector()
      rep <- list()
      if("raw" %in% input$files_to_download){
        fs <- c(fs, "raw.csv")
        rep$data <- rRNAobj()
        write.csv(rRNAobj()$e_data, file="raw.csv")
      }
      if("filtered" %in% input$files_to_download){
        fs <- c(fs, "filtered.csv")
        rep$data <- filtered_data()
        write.csv(filtered_data()$e_data, file="filtered.csv")
      }
      if("groupings" %in% input$files_to_download){
        fs <- c(fs, "groupings.csv")
        write.csv(group_df_tab(), file="groupings.csv")
      }
      if("normalized" %in% input$files_to_download){
        rep$data <- normalized_data()
        fs <- c(fs, "normalized.csv", "normalized.png", "abun_rich_raw.png", "abun_rich_norm.png")
        write.csv(normalized_data()$e_data, file="normalized.csv")
        ggsave(norm_plot_obj(), filename = "normalized.png", device="png")
        ggsave(ra_raw_plot(), filename="abun_rich_raw.png", device="png")
        ggsave(ra_norm_plot(), filename="abun_rich_norm.png", device="png")
      }
      if("outliers" %in% input$files_to_download){
        rep$jaccard <- outlier_jaccard()
        fs <- c(fs, "outliers.png")
        ggsave(jac_plot_obj(), filename="outliers.png", device="png")
      }
      if("alphadiv" %in% input$files_to_download){
        fs <- c(fs, "alphadiv.csv", "alphadiv.png")
        rep$adiv <- a_div()
        write.csv(a_div(), file="alphadiv.csv")
        ggsave(adiv_plot_obj(), filename="alphadiv.png", device="png")
      }
      if("richness" %in% input$files_to_download){
        fs <- c(fs, "rich.csv", "rich.png")
        rep$rich <- rich()
        write.csv(rich(), file="rich.csv")
        ggsave(rich_plot_obj(), filename="rich.png", device="png")
      }
      if("ordination" %in% input$files_to_download){
        fs <- c(fs, "ordplot.png")
        #ggsave(dimcheck_obj(), filename="dimcheck.png", device="png")
        ggsave(ord_plot_obj(), filename="ordplot.png", device="png")
      }
      if("diffabun" %in% input$files_to_download){
        fs <- c(fs, "diffabun.csv", "daflag.png", "dalogfc.png", "allda.png")
        rep$diffabun <- diffabun_res()
        write.csv(diffabun_res()$allResults, file="diffabun.csv")
        ggsave(da_flag_plot_obj(), filename="daflag.png", device="png")
        ggsave(da_logfc_plot_obj(), filename="dalogfc.png", device="png")
        ggsave(plot_all_da_obj()[[1]], filename="allda.png", device="png")
      }
      if("indicspec" %in% input$files_to_download){
        fs <- c(fs, "indicspec.csv", "indsp.png")
        rep$indsp <- indsp_res()
        write.csv(indsp_res(), file="indicspec.csv")
        ggsave(indsp_plot_obj(), filename="indsp.png", device="png")
      }
      if("combined" %in% input$files_to_download){
        fs <- c(fs, "combined.csv")
        rep$combined <- taxares()
        write.csv(taxares(), file="combined.csv")
      }
      if("report" %in% input$files_to_download){
        fs <- c(fs, "report.docx")
        data <- rep
        classes <- unlist(lapply(data, class))
        print(classes)
        params <- list(data=data, classes=classes)
        rmarkdown::render(tempReport, output_file="report.docx", params=params, envir = new.env(parent = globalenv()))
      }
      
      print(fs)
      
      zip(zipfile=fname, files=fs)
      if(file.exists(paste0(fname,".zip"))){file.rename(paste0(fname,".zip"),fname)}
    },
    contentType = "application/zip"
  )
  
  
}) #end server
