    
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
#library(phyloseq)
library(vegan)
library(goeveg)
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
    return(pmartRseq::as.seqData(e_data = as.character(input$biom$datapath), f_data = as.character(input$qiime$datapath), edata_cname = "OTU", data_type = "rRNA"))
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
    output$boxes <- renderUI({get_checkbox_output_list(metadata_obj())})
    
    # Initialize new_metadata_obj object. If no metadata filtering then it's a copy of metadata_obj.
    # If metadata filtering, will be overwritten with filters 
    new_metadata_obj <- reactive({
      return(metadata_obj()[input$selected_indices + 1, ]) # add one because javascript is zero-indexed
    })
    
    # If the user brushes a chart, input$selected_indices will change
    # If that change is observed, subset new_metadata_obj and 
    # create subsetted charts and table
    observeEvent(input$selected_indices, {
      print(input$box1)
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
    # observe({
    #   filtered_rRNA_obj
    ################# JSON object for Meg #################
       rRNA_filtered <- jsonlite::toJSON(list(filtered_rRNA_obj$e_data, filtered_rRNA_obj$e_meta, filtered_rRNA_obj$f_data))
    #   #write(rRNA_filtered, file = "rRNA_filtered.json")
    # })
    return(isolate(filtered_rRNA_obj))
  })
  
  
  # ################ Group Designation Tab #################

  # Input main effects used for groupings
    output$group1 <- renderUI({
      selectInput("group1",
                  label = "Main Effect 1",
                  choices = colnames(filtered_data()$f_data))
    })
    
  # Can have up to 2 main effects
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
    output$group_DF <- DT::renderDataTable(attr(groupDF(), "group_DF"))
    #)
  
  # Also output a table showing the number of reps in each group
    output$group_tab <- DT::renderDataTable(as.data.frame(table(attr(groupDF(), "group_DF")$Group)))
    
    # ################ Normalization Tab #################
    
  # Select which normalization function to use
    output$normFunc <- renderUI({
      selectInput("normFunc",
                  label = "Normalization Function",
                  choices = c("percentile","tss","rarefy","poisson","deseq","tmm","css"),
                  selected = "css")
    })
    
  # Create normalized data
    normalized_data <- reactive({
      validate(need(length(input$normFunc) == 1, "Need to specify a normalization function."))
      validate(need(input$normFunc %in% c("percentile","tss","rarefy","poisson","deseq","tmm","css"), "Normalization function must be one of the options specified."))
      return(pmartRseq::normalize_data(omicsData=groupDF(), norm_fn=input$normFunc, normalize=TRUE))
    })
    
  # Look at normalized results
    output$normData <- DT::renderDataTable(normalized_data()$e_data, rownames = FALSE)

  # Try to make a stacked bar plot - not working right now
    output$norm_plot <- renderPlot({
      #browser()
      plot(normalized_data(), class="taxonomy2")
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
    
  # Create a plot of raw abundance vs raw richness
    output$ra_raw <- renderPlot({
      plot(abun_raw(), rich_raw(), plot_title="Raw Data")
    })
    
  # Create a plot of normalized abundance vs normalized richness to see if there is a reduction in correlation
    output$ra_norm <- renderPlot({
      plot(abun_norm(), rich_norm(), plot_title="Normalized Data")
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
  
  # Show alpha diversity plot
    output$adiv_plot <- renderPlot({
      plot(a_div(), x_axis=input$xaxis, color=input$color)
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
     
  # Show richness plot
    output$rich_plot <- renderPlot({
      plot(rich(), x_axis=input$xaxis, color=input$color)
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
    
    output$dimcheck <- renderPlot({
      goeveg::dimcheckMDS(vegdata(), distance = input$beta_index)
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
    
    
  # Use vegan to calculate scores for beta diversity index
    vegmds <- reactive({
      validate(
        need(length(input$beta_index) == 1, "There needs to be one beta diversity index.")
      )
      
      return(vegan::metaMDS(vegdata(), distance = input$beta_index, k = input$k, autotransform = FALSE, na.rm = TRUE))
    })
    
  # Plot showing beta diversity
    output$ord_plot <- renderPlot({
      #if(input$ord_method == "NMDS"){
        pmartRseq::pmartRseq_NMDS(res = vegmds(), 
                  grp = as.factor(attr(normalized_data(),"group_DF")[match(rownames(vegdata()), attr(normalized_data(),"group_DF")[,attr(normalized_data(),"cnames")$fdata_cname]),input$ord_colors]))
      # }else if(input$ord_method == "PCA"){
      #   mead_PCA(XX = vegmds(),
      #            ZZ = as.factor(attr(normalized_data(),"group_DF")[match(rownames(vegdata()), attr(normalized_data(),"group_DF")[,attr(normalized_data(),"cnames")$fdata_cname]),input$ord_colors]))
      # }
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
    
  # Calculate normalization factors to use in differential abundance test - will use the same that was used on normalization tab
    norm_factors <- reactive({
      validate(need(length(input$normFunc) == 1, "Need to specify a normalization function."))
      validate(need(input$normFunc %in% c("percentile","tss","rarefy","poisson","deseq","tmm","css"), "Normalization function must be one of the options specified."))
      return(pmartRseq::normalize_data(omicsData=groupDF(), norm_fn=input$normFunc, normalize=FALSE))
    })
    
  # Perform differential abundance analysis
    diffabun_res <- reactive({
      validate(need(length(input$da_index) == 1, "Need to specify a differential abundance test"))
      validate(need(length(input$pval_adjust) == 1, "Need to specify a p-value adjustment method"))
      
      return(pmartRseq::countSTAT(omicsData = groupDF(), norm_factors = norm_factors()$scale_param, comparisons = "all", control = NULL, test = input$da_index, pval_adjust = input$pval_adjust, pval_thresh = 0.05))
    })
    
  # Look at the results - this is hard to look at, maybe remove?
    output$da_res <- DT::renderDataTable(diffabun_res()$allResults)
    
  # Plot showing number differentially abundant in each comparison and direction of change
    output$flag_plot <- renderPlot({
      plot(diffabun_res(), type = "flag")
    })
    
  # Heatmap showing the log2foldchanges of differentially abundant features
    output$logfc_plot <- renderPlot({
      plot(diffabun_res(), type = "logfc")
    })
    
  # Plot showing log fold changes and p-values of all features, grouped by taxonomy
    output$plot_all_da <- renderPlot({
      pmartRseq::plot_all_diffabun(countSTAT_results = diffabun_res(), omicsData = normalized_data(), x_axis = "taxonomy2", x_lab = "Phylum")
    })
}) #end server
