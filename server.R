
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#


#source("./Report/R/report.R")

#Sys.setenv(R_ZIPCMD="/usr/bin/zip")

options(shiny.maxRequestSize=30*1024^2)

filtered_rRNA_obj <- list()

shinyServer(function(input, output, session) {
  #---------- listen for session reset -------#
  observeEvent(input$reset_button, {js$reset()})
  #---------- rRNA object import -------------#
  importObj <- reactive({
    # validate(
    #   need(input$biom != "", "Please select a biom file")
    # )
    validate(
      need(grepl(".biom", as.character(input$e_data$name)) | grepl(".csv", as.character(input$e_data$name)), "Data file must be in .biom or .csv format")
    )

    # validate(
    #   need(input$qiime != "", "Please select a qiime file")
    # )

    validate(
      need(input$e_data != "", "Please select a data file")
    )

    validate(
      need(input$f_data != "", "Please select a sample metadata file")
    )
    if (!grepl(pattern = "\\.biom$", x = input$e_data$datapath)) {
      validate(need(input$e_meta != "", message = "Please upload feature metadata file"))
    }

    # check for e_meta file
    if (is.null(input$e_meta$name) & grepl(".biom", as.character(input$e_data$name))) {
      return(pmartRseq::import_seqData(e_data_filepath = as.character(input$e_data$datapath),
                                       f_data_filepath = as.character(input$f_data$datapath),
             e_meta_filepath = NULL))
    }
    # import e_meta if there is one
    if (!is.null(input$e_meta$name) & !grepl(".biom", as.character(input$e_data$name))) {
      return(pmartRseq::import_seqData(e_data_filepath = as.character(input$e_data$datapath),
                                       f_data_filepath = as.character(input$f_data$datapath),
                                       e_meta_filepath = as.character(input$e_meta$datapath)))
    }
  }) #end rRNAobj import

  # Guess at column identifiers
  output$e_data_cname <- reactive(importObj()$guessed_edata_cname)
  output$f_data_cname <- reactive(importObj()$guessed_fdata_cname)

  output$new_edata_cname <- renderUI({
    selectInput(inputId = "new_e_data_cname", label = "New expression data identifier is:",
                choices = colnames(importObj()$e_data),
                selected = importObj()$guessed_edata_cname,
                               multiple = FALSE)
  })

  output$new_fdata_cname <- renderUI({
    selectInput(inputId = "new_f_data_cname", label = "New metadata identifier is:",
                choices = colnames(importObj()$f_data),
                selected = importObj()$guessed_fdata_cname,
                multiple = FALSE)
  })
  # output$edata_cname <- renderUI({
  #   selectInput("edata_cname",
  #               label = "Identifier column name in data",
  #               choices = colnames(importObj()$e_data),
  #               selected = importObj()$guessed_edata_cname,
  #               multiple = FALSE)
  # })
  #
  # output$fdata_cname <- renderUI({
  #   selectInput("fdata_cname",
  #               label = "Sample column name in sample metadata",
  #               choices = colnames(importObj()$f_data),
  #               selected = importObj()$guessed_fdata_cname,
  #               multiple = FALSE)
  # })
  #
  # output$taxa_cname <- renderUI({
  #   selectInput("taxa_cname",
  #               label = "Taxonomic column of interest in feature metadata",
  #               choices = colnames(importObj()$e_meta),
  #               selected = importObj()$guessed_taxa_cname,
  #               multiple = FALSE)
  # })

  #rRNAobj <- reactive({ return(importObj()) })

  #observeEvent(input$Upload,{

    rRNAobj <- reactive({
      validate(
        need(importObj, message = "Please upload data files to begin analysis"),
        need(!is.null(importObj()$guessed_fdata_cname), message = "Please upload data files to begin analysis"),
        need(!is.null(importObj()$guessed_edata_cname), message = "Please upload data files to begin analysis"),
        need(!is.null(importObj()$guessed_taxa_cname), message = "Please upload data files to begin analysis"),
        need(!is.null(input$new_f_data_cname), message = "Please upload data files to begin analysis")
      )
      tmp <- pmartRseq::as.seqData(e_data = importObj()$e_data,
                                   f_data = importObj()$f_data,
                                   e_meta = importObj()$e_meta,
                                   fdata_cname = input$new_f_data_cname,
                                   edata_cname = input$new_e_data_cname,
                                   taxa_cname = importObj()$guessed_taxa_cname,
                                   data_type = "rRNA")

      if(ncol(tmp$e_meta) <= 2){
        tmp <- pmartRseq::split_emeta(tmp, cname=attr(tmp,"cnames")$taxa_cname, split1=",", numcol=7, split2="__", num=2, newnames=NULL)
      }else{
        tmp <- pmartRseq::split_emeta(tmp, cname=attr(tmp,"cnames")$edata_cname, split1=NULL, numcol=7, split2="__", num=2, newnames=NULL)
      }

      return(tmp)

    })


    output$rollup <- renderUI({
      selectInput("rollup",
                  label = "Which taxonomic level should analysis be performed at?",
                  choices = c("Kingdom","Phylum","Class","Order","Family","Genus","Species",attr(rRNAobj(), "cnames")$edata_cname),
                  selected = attr(rRNAobj(), "cnames")$edata_cname,
                  multiple = FALSE)
    })

    rRNA_agg <- reactive({
      validate(
        need(length(input$rollup) == 1, "Need to specify a taxonomic level")
      )

      return(pmartRseq::taxa_rollup(omicsData = rRNAobj(), level = input$rollup, taxa_levels = NULL))
    })

    output$sample_data <- DT::renderDataTable(expr =
                                                data.frame(rRNA_agg()$e_data), rownames = FALSE, class = 'cell-border stripe compact hover',
                                              options = list(columnDefs = list(list(
                                                targets = c(1:(ncol((rRNA_agg()$e_data)) - 1)),
                                                render = JS(
                                                  "function(data, type, row, meta) {",
                                                  "return type === 'display' && data.length > 10 ?",
                                                  "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
                                                  "}")
                                              ))), callback = JS('table.page(3).draw(false);'))

    output$datsum <- renderPrint({
      validate(
        need(!(is.null(rRNA_agg())), message = "please import sample metadata")
      )
      validate(
        need(!(is.null(rRNA_agg())) , message = "Upload data first")
      )
      summary(rRNA_agg())
    })

    fursumm <- reactive({
      validate(
        need(!(is.null(rRNA_agg())), message = "please import sample metadata")
      )
      validate(
        need(!(is.null(rRNA_agg())), message = "Upload data first")
      )

      temp <- rRNA_agg()$e_meta
      names <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
      temp$TAXONOMY <- apply(temp[,which(colnames(temp) %in% names)], 1, function(x) paste(x, collapse="; "))
      temp <- merge(melt(rRNA_agg()$e_data), temp, by=attr(rRNA_agg(), "cnames")$edata_cname)

      vars <- c(attr(rRNA_agg(), "cnames")$edata_cname, attr(rRNA_agg(), "cnames")$taxa_cname, "TAXONOMY")
      vars <- lapply(vars, as.symbol)

      temp <- temp %>% dplyr::mutate(Tot=sum(value,na.rm=TRUE)) %>%
                       dplyr::group_by_(.dots=vars) %>%
                       dplyr::summarise("Mean Across Samples"=mean(value, na.rm=TRUE),
                                        "Median Across Samples"=median(value, na.rm=TRUE),
                                        "Max Across Samples"=max(value, na.rm=TRUE),
                                        "Min Across Samples"=min(value, na.rm=TRUE),
                                        "Sum Across Samples"=sum(value, na.rm=TRUE),
                                        "Percentage of Total"=sum(value,na.rm=TRUE)/unique(Tot)*100)
      return(temp)

    })

    output$further_summary <- DT::renderDataTable(expr =
                                                data.frame(fursumm()), rownames = FALSE, class = 'cell-border stripe compact hover',
                                              options = list(columnDefs = list(list(
                                                render = JS(
                                                  "function(data, type, row, meta) {",
                                                  "return type === 'display' && data.length > 10 ?",
                                                  "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
                                                  "}")
                                              )), pageLength=5, order=list(3, 'desc')), callback = JS('table.page(3).draw(false);'))


    # ################ Group Designation Tab #################

    group_vars <- reactive({
      intersect(which(lapply(apply(rRNA_agg()$f_data, 2, function(z) table(z))[unlist(lapply(apply(rRNA_agg()$f_data, 2, function(x) table(x)), function(y) any(is.finite(y))))], function(w) max(w)) > 2), which(apply(rRNA_agg()$f_data, 2, function(v) length(unique(v))) >= 2))

    })

    output$gdfMainEffect <- renderUI({
      selectInput("gdfMainEffect",
                  label = "Main Effect(s) to use for Groupings",
                  choices = colnames(rRNA_agg()$f_data)[group_vars()],
                  multiple = TRUE)
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
    #                 choices = c("NA",colnames(rRNA_agg()$f_data)),
    #                 selected = NULL)
    #   })
    #
    #   output$cov2 <- renderUI({
    #     selectInput("cov2",
    #                 label = "Covariate 2",
    #                 choices = c("NA",colnames(rRNA_agg()$f_data)),
    #                 selected = NULL)
    #     #groupDesignation
    #   })
    # #})


    # Create groups with main effects
    groupDF <<- reactive({
      mainEffects <- input$gdfMainEffect
      # if(!is.na(input$cov1) | !is.na(input$cov2)){
      #   covariates <- c(input$cov1, input$cov2)
      #   if(any(is.na(covariates))){
      #     covariates <- covariates[!is.na(covariates)]
      #   }
      # }
      validate(
        need(length(mainEffects) > 0, "There needs to be at least one grouping variable")
      )

      return(pmartRseq::group_designation(rRNA_agg(), main_effects = mainEffects))

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

  #})

  # ################ Filtering Tab #################

  #-------- filter history support -----------#
  filters <- reactiveValues(otu = list(), sample = list(), taxa = list(), metadata = list(), outlier = list())

  #--------- kovera observer ---------#
  observeEvent(input$otu_filter_go,{
    filters$otu[[input$otu_filter_go]] <- otu_filter_obj()
  })

  #--------- sample observer ---------#
  observeEvent(input$sample_filter_go, {
    filters$sample[[input$sample_filter_go]] <- sample_filter_obj()
  })

  #--------- taxa observer ---------#
  observeEvent(input$taxa_filter_go, {
    filters$taxa[[input$taxa_filter_go]] <- taxa_filter_obj()
  })

  #--------- metadata observer ---------#
  observeEvent(input$metadata_filter_go, {
    filters$metadata[[input$metadata_filter_go]] <- sample_metadata_filter_obj()
  })

  #--------- filter application observer ---------#
  filtered_rRNA_obj <- NULL
  #filt1 <- NULL
  observeEvent(groupDF(), {
      filtered_rRNA_obj <<- groupDF()
  }, priority = 10)

  #--------- apply filter applications on click ---------#
  observeEvent(input$otu_filter_go, {
    filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$otu[[input$otu_filter_go]],
                                   omicsData = filtered_rRNA_obj,
                                   num_samps = input$filter_kOverA_sample_threshold,
                                   upper_lim = input$filter_count_threshold)
  })
  observeEvent(input$sample_filter_go, {
    filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$sample[[input$sample_filter_go]],
                                               omicsData = filtered_rRNA_obj,
                                               upper_lim = input$n)
  })
  observeEvent(input$taxa_filter_go, {
    filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$taxa[[input$taxa_filter_go]],
                                               omicsData = filtered_rRNA_obj,
                                               keep_taxa = input$keep_taxa)
  })
  observeEvent(input$metadata_filter_go, {
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
      sample_names[[i]] <- dplyr::filter_(temp, interp(~v%in%input[[check_boxes[i]]], v=as.name(check_boxes[i]))) %>%
        dplyr::select(eval(quote(attr(filtered_rRNA_obj, "cnames")$fdata_cname)))#find the column name and associated check box
    }
    sample_names <- Reduce(dplyr::intersect, sample_names)
    #check if there are no samples to remove
    if (nrow(sample_names) == length(temp[, attr(filtered_rRNA_obj, "cnames")$fdata_cname])) {
      #if no samples to remove return unmodified object
      return(filtered_rRNA_obj)
    } else {
      # remove the samples
      to_remove <- temp[ which(!(as.character(temp[, attr(filtered_rRNA_obj, "cnames")$fdata_cname]) %in% as.character(unlist(sample_names)))),
                         attr(filtered_rRNA_obj, "cnames")$fdata_cname]
      filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$metadata[[input$metadata_filter_go]],
                                                 omicsData = filtered_rRNA_obj,
                                                 samps_to_remove = to_remove)
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
    filt1 <- groupDF()
    # no sample filter yet
    if (input$sample_filter_go == 0) {
      filtered_rRNA_obj <<- filt1
      #return(filtered_rRNA_obj)
    }else{
      # apply sample filter
      isolate({
        filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$sample[[input$sample_filter_go]],
                                                   omicsData = filt1,
                                                   upper_lim = input$n)
        #return(filtered_rRNA_obj)
      })
    }

    if(input$taxa_filter_go == 0){
      filtered_rRNA_obj <<- filtered_rRNA_obj
    }else{
      isolate({
        filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filter$taxa[[input$taxa_filter_go]],
                                      omicsData = filtered_rRNA_obj,
                                      keep_taxa = input$keep_taxa)
      })
    }
    return(filtered_rRNA_obj)
  })

  observeEvent(input$sample_reset_button, {

    updateNumericInput(session, "n",
                       value = 0)
    if (input$otu_filter_go == 0) {
      filt1 <<- groupDF()
    }else{
      # apply k over a filter
      isolate({
        filt1 <<- pmartRseq::applyFilt(filter_object = filters$otu[[input$otu_filter_go]],
                                      omicsData = groupDF(),
                                      num_samps = input$filter_kOverA_sample_threshold,
                                      upper_lim = input$filter_count_threshold)
      })
    }

    if(input$taxa_filter_go == 0){
      filt1 <<- filt1
    }else{
      isolate({
        filt1 <<- pmartRseq::applyFilt(filter_object = filter$taxa[[input$taxa_filter_go]],
                                      omicsData = filt1,
                                      keep_taxa = input$keep_taxa)
      })
    }

    # no sample filter yet
    filtered_rRNA_obj <<- filt1
    return(filtered_rRNA_obj)
  })

  observeEvent(input$taxa_reset_button, {

    # output$criteria <- renderUI({
    #   selectInput("criteria",
    #               label = "Which taxonomic level to use for filtering",
    #               choices = colnames(groupDF()$e_meta),
    #               multiple = FALSE)
    # })
    #
    # keep_taxa = reactive({
    #   # Create logical indicating the samples to keep, or dummy logical if nonsense input
    #   validate(
    #     need(!(is.null(groupDF()$e_meta)), message = "please import feature metadata")
    #   )
    #   if (!is.null(groupDF()$e_meta)) {
    #     return(unique(groupDF()$e_meta[,input$criteria]))
    #   } else {
    #     return(NULL)
    #   }
    #
    # })
    #
    # output$keep_taxa <- renderUI({
    #     selectInput("keep_taxa",
    #                 label = "Which taxa to keep in the analysis",
    #                 choices = c(keep_taxa),
    #                 multiple = TRUE)
    #   })

    output$criteria <- renderUI({
      selectInput("criteria",
                  label = "Which taxonomic level to use for filtering",
                  choices = colnames(groupDF()$e_meta),
                  selected = colnames(groupDF()$e_meta)[2],
                  multiple = FALSE)
    })

    output$keep_taxa <- renderUI({
      selectInput("keep_taxa",
                  label = "Which taxonomies to keep",
                  choices = unique(groupDF()$e_meta[,2]),
                  multiple = TRUE)
    })

    filt1 <- groupDF()
    # no sample filter yet
    if (input$sample_filter_go == 0) {
      filtered_rRNA_obj <<- filt1
      #return(filtered_rRNA_obj)
    }else{
      # apply sample filter
      isolate({
        filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$sample[[input$sample_filter_go]],
                                                   omicsData = filt1,
                                                   upper_lim = input$n)
        #return(filtered_rRNA_obj)
      })
    }

    # no otu filter yet
    if (input$otu_filter_go == 0) {
      filtered_rRNA_obj <<- filtered_rRNA_obj
      #return(filtered_rRNA_obj)
    }else{
      # apply sample filter
      isolate({
        filtered_rRNA_obj <<- pmartRseq::applyFilt(filter_object = filters$otu[[input$otu_filter_go]],
                                      omicsData = filtered_rRNA_obj,
                                      num_samps = input$filter_kOverA_sample_threshold,
                                      upper_lim = input$filter_count_threshold)

      })
    }
    return(filtered_rRNA_obj)
    })

  #})

  #--------- Metadata Object for filtering -----------#
  metadata_obj <- reactive({
    #------- TODO: need to display an error if the rRNA object isn't created! ---------#
    validate(
      need(!(is.null(groupDF()$f_data)), message = "rRNA object fail. Check to make sure file paths are correct")
    )
    input$metadata_filter_go
    input$metadata_reset_button
    input$sample_filter_go
    input$sample_reset_button
    input$otu_filter_go
    input$otu_reset_button
    if (is.null(input$f_data)) {
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
    return(pmartRseq::sample_based_filter(omicsData = groupDF(), fn = "criteria"))
  })
  observeEvent(input$metadata_reset_button, {
    output$sample_metadata <- DT::renderDataTable(
      groupDF()$f_data, rownames = FALSE, class = 'cell-border stripe compact hover',
      options = list(columnDefs = list(list(
        targets = c(1:(ncol((metadata_obj())) - 1)),
        render = JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 10 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
          "}")
      ))), callback = JS('table.page(3).draw(false);'))  })

  observeEvent(groupDF(), {
    # If the metadata has been loaded, create the charts
    output$plots <- renderUI({get_plot_output_list(metadata_obj())})
    output$boxes <- renderUI({get_checkbox_output_list(groupDF()$f_data)})
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
    output$boxes <- renderUI({get_checkbox_output_list(rRNA_agg()$f_data)})
  })

  # end sample metadata filtering
  # end sample metadata filtering
  #--------------- k over a filtering  -----------------#
  kovera_k <- 1
  observeEvent(groupDF(), {
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
      if (input$count_filter_fun == "ka") {
        return(numericInputRow("filter_kOverA_sample_threshold", "in",
                        min = 1, max = maxSamples(), value = kovera_k, step = 1))
      } else {
        return(NULL)
      }

    })
  })

  otu_filter_obj <- reactive({
    validate(
      need(length(filtered_rRNA_obj) > 0 , message = "Upload data first")
    )
    return(pmartRseq::count_based_filter(groupDF(), fn = input$count_filter_fun))
  })


  output$read_counts_plot <- renderPlot({
    validate(
      need( input$filter_count_threshold >= 0, message = "Enter a count minimum >= 0")
      #need( input$filter_kOverA_sample_threshold >= 0, message = "Enter a sample minimum >= 0")
    )
    if (input$sample_filter_go == 0 & input$otu_filter_go == 0 & input$taxa_filter_go == 0) {
      plot(otu_filter_obj(), min_num = input$filter_count_threshold, min_samp = input$filter_kOverA_sample_threshold)
    } else{
      otu_filt_obj <- pmartRseq::count_based_filter(filtered_rRNA_obj, fn = input$count_filter_fun)
      if (input$count_filter_fun == "ka") {
        plot(otu_filt_obj, min_num = input$filter_count_threshold, min_samp = input$filter_kOverA_sample_threshold)
      }
      if (input$count_filter_fun != "ka") {
        plot(otu_filt_obj, min_num = input$filter_count_threshold)
      }


    }
  })


  #-------------- Library read filtering -----------#

  sample_filter_obj <- reactive({
    return(pmartRseq::sample_based_filter(omicsData = groupDF(), fn = "sum"))
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

  #-------------- taxa filtering -----------#

  output$criteria <- renderUI({
      selectInput("criteria",
                  label = "Which taxonomic level to use for filtering",
                  choices = colnames(groupDF()$e_meta),
                  selected = colnames(groupDF()$e_meta)[2],
                  multiple = FALSE)
  })

  output$keep_taxa <- renderUI({
      selectInput("keep_taxa",
                  label = "Which taxonomies to keep",
                  choices = unique(groupDF()$e_meta[,input$criteria]),
                  multiple = TRUE)
  })

  taxa_filter_obj <- reactive({
    if(length(input$criteria) == 0){
      criteria <- colnames(groupDF()$e_meta)[2]
    }else{
      criteria <- input$criteria
    }
    return(pmartRseq::metadata_based_filter(omicsData = groupDF(), criteria = criteria))
  })

  output$taxa_counts <- renderPrint({
    validate(
      need(length(input$keep_taxa) > 0, message = "Need taxa criteria")
    )
    if (input$taxa_filter_go == 0 & input$otu_filter_go == 0 & input$sample_filter_go == 0) {
      #table(taxa_filter_obj()[which(taxa_filter_obj()[,input$criteria] %in% input$keep_taxa),input$criteria])
      #which(taxa_filter_obj()[,input$criteria] %in% input$keep_taxa)
      cat("This keeps ",length(which(taxa_filter_obj()[,input$criteria] %in% input$keep_taxa))," out of a possible ",nrow(taxa_filter_obj()), " features (roughly ", length(which(taxa_filter_obj()[,input$criteria] %in% input$keep_taxa))/nrow(taxa_filter_obj())*100,"%). This correlates to a total number of ",sum(taxa_filter_obj()$Sum[which(taxa_filter_obj()[,input$criteria] %in% input$keep_taxa)], na.rm=TRUE)," sequences kept out of a possible ",sum(taxa_filter_obj()$Sum, na.rm=TRUE)," sequences (roughly ",sum(taxa_filter_obj()$Sum[which(taxa_filter_obj()[,input$criteria] %in% input$keep_taxa)], na.rm=TRUE)/sum(taxa_filter_obj()$Sum, na.rm=TRUE)*100,"%).")
    } else {
      taxa_filter_obj <- pmartRseq::metadata_based_filter(omicsData = filtered_rRNA_obj, criteria = input$criteria)
      cat("This keeps ",length(which(taxa_filter_obj[,input$criteria] %in% input$keep_taxa))," out of a possible ",nrow(taxa_filter_obj), " features (roughly ", length(which(taxa_filter_obj[,input$criteria] %in% input$keep_taxa))/nrow(taxa_filter_obj)*100,"%). This correlates to a total number of ",sum(taxa_filter_obj$Sum[which(taxa_filter_obj[,input$criteria] %in% input$keep_taxa)], na.rm=TRUE)," sequences kept out of a possible ",sum(taxa_filter_obj$Sum, na.rm=TRUE)," sequences (roughly ",sum(taxa_filter_obj$Sum[which(taxa_filter_obj[,input$criteria] %in% input$keep_taxa)], na.rm=TRUE)/sum(taxa_filter_obj$Sum, na.rm=TRUE)*100,"%).")
    }
  })

  #------------ reactive filtered data for downstream processing --------------#
  #--------- outlier observer ---------#
  outlier_filter_obj <- reactive({
    validate(
      need(length(filtered_rRNA_obj) > 0 , message = "Upload data first"),
      need(!is.null(groupDF()), message = "Need group designation first")
    )
    return(pmartRseq::sample_based_filter(omicsData = filtered_rRNA_obj, fn = "criteria"))
  })

  observeEvent(input$remove_outliers, priority = 2, {
    filters$outliers[[input$remove_outliers]] <- outlier_filter_obj()
    filtO <- pmartRseq::applyFilt(omicsData = filtered_rRNA_obj, outlier_filter_obj(), samps_to_remove =  outlier_jaccard()[outlier_jaccard()$jaqsID %in% selected_outliers()[["key"]], "jaqsID"])
    filtered_rRNA_obj <<- filtO
  })

  filtered_data <- reactive({
    event_data("plotly_selected")
    input$gdfMainEffect
    input$metadata_filter_go
    input$metadata_reset_button
    input$sample_filter_go
    input$sample_reset_button
    input$otu_filter_go
    input$otu_reset_button
    input$taxa_filter_go
    input$taxa_reset_button
    input$remove_outliers
    return(filtered_rRNA_obj)
  })



  output$summ_filt <- renderPrint({
    validate(
      need(!(is.null(metadata_obj())), message = "please import sample metadata")
    )
    validate(
      need(!(is.null(filtered_data())) , message = "Upload data first")
    )
    summary(filtered_data())
  })

  output$nrow_edata <- renderPrint({
    nrow(filtered_data()$e_data)
  })



  ################# Outliers Tab #################



  outlier_jaccard <- reactive({
    event_data("plotly_selected")
    jaqs <- pmartRseq::jaccard_calc(omicsData = filtered_data())
    jaqs$jaqsID <- jaqs[, attr(jaqs,"cname")$fdata_cname]
    return(jaqs)
  })

  jac_plot_obj <- reactive({
    plot(outlier_jaccard())
  })

  #define global plotly style
  f <<- list(
    family = "Courier New, monospace",
    size = 18,
    color = "#7f7f7f"
  )

  output$outlier_jaccard_plot <- renderPlotly({
    d <- event_data("plotly_selected")
    p <- plotly::plot_ly(data = outlier_jaccard(),
                         x = ~jaqsID,
                         y = ~Average) %>%
      add_markers(key = ~jaqsID, color = I("black"))
    if (!is.null(d)) {
      m <- outlier_jaccard()[outlier_jaccard()$jaqsID %in% d[["key"]], ]
      p <- add_markers(p, data = m, color = I("red"))
    }


    x <- list(
      title = "",
      titlefont = f
    )
    y <- list(
      title = "Jaccard Dis/similarity",
      titlefont = f
    )
    p$elementId <- NULL
    layout(p, dragmode = "lasso", showlegend = FALSE, xaxis = x, yaxis = y)

  })

  selected_outliers <- reactive({
    input$remove_outliers
    event_data("plotly_selected")
  })

  output$audies <- renderTable({
    outlier_jaccard()[outlier_jaccard()$jaqsID %in% selected_outliers()[["key"]], ]
  })



  abundance <- reactive({
    event_data("plotly_selected")
    input$remove_outliers
    abundances <- abundance_calc(omicsData = filtered_data())
    abundances$abundID <- rownames(abundances)
    return(abundances)
  })

  output$outlier_abundance_plot <- renderPlotly({
    d <- event_data("plotly_selected")
    p <- plotly::plot_ly(data = data.frame(abundance()),
                         x = ~abundID,
                         y = ~abundance) %>%
      add_markers(key = ~abundID, color = I("black"))
    if (!is.null(d)) {
      m <- abundance()[abundance()$abundID %in% d[["key"]], ]
      p <- add_markers(p, data = m, color = I("red"))
    }
    x <- list(
      title = "",
      titlefont = f
    )
    y <- list(
      title = "Abundance",
      titlefont = f
    )
    p$elementId <- NULL
    layout(p, dragmode = "lasso", showlegend = FALSE, xaxis = x, yaxis = y)  })


  output$outlier_richness_plot <- renderPlotly({
    d <- event_data("plotly_selected")
    long_richness <- reshape2::melt(rich_raw())
    p <- plotly::plot_ly(data = long_richness,
                         x = ~variable,
                         y = ~value) %>%
      add_markers(key = ~variable, color = I("black"))
    if (!is.null(d)) {
      m <- long_richness[long_richness$variable %in% d[["key"]], ]
      p <- add_markers(p, data = m, color = I("red"))
    }
    x <- list(
      title = "",
      titlefont = f
    )
    y <- list(
      title = "Richness",
      titlefont = f
    )
    p$elementId <- NULL
    layout(p, dragmode = "lasso", showlegend = FALSE, xaxis = x, yaxis = y)  })




  # ################ Normalization Tab #################

  # Select which normalization function to use
  output$normFunc <- renderUI({
    selectInput("normFunc",
                label = "Normalization Function",
                choices = c("percentile","tss","rarefy","poisson","deseq","tmm","css","log","clr","none"),
                selected = "css")
  })

  # Create normalized data
  normalized_data <- reactive({
    validate(need(length(input$normFunc) == 1, "Need to specify a normalization function."))
    validate(need(input$normFunc %in% c("percentile","tss","rarefy","poisson","deseq","tmm","css","log","clr","none"), "Normalization function must be one of the options specified."))

    #return(pmartRseq::split_emeta(pmartRseq::normalize_data(omicsData=filtered_data(), norm_fn=input$normFunc, normalize=TRUE), cname="OTU", split1=NULL, numcol=7, split2="__", num=2, newnames=NULL))
    return(pmartRseq::normalize_data(omicsData=filtered_data(), norm_fn=input$normFunc, normalize=TRUE))
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
  output$norm_plot <- renderPlotly({
    #plot(normalized_data(), class="Phylum")
    #print(norm_plot_obj())
    #plotly::ggplotly( plot(normalized_data(), class=input$norm_class))
    p <- plotly::ggplotly(norm_plot_obj(), tooltip = c("x", "y", "fill"))
    p$elementId <- NULL
    p
  })

  # Calculate abundance on normalized data
  abun_norm <- reactive({
    return(suppressWarnings(pmartRseq::abundance_calc(normalized_data())))
  })

  # Calculate abundance on raw data
  abun_raw <- reactive({
    event_data("plotly_selected")
    return(suppressWarnings(pmartRseq::abundance_calc(filtered_data())))
  })

  # Calculate richness on normalized data
  rich_norm <- reactive({
    return(suppressWarnings(pmartRseq::richness_calc(normalized_data(), index="observed")))
  })

  # Calculate richness on raw data
  rich_raw <- reactive({
    event_data("plotly_selected")
    input$remove_outliers
    return(suppressWarnings(pmartRseq::richness_calc(filtered_data(), index="observed")))
  })

  ra_raw_plot <- reactive({
    plot(abun_raw(), rich_raw(), plot_title="Raw Data", samplabel = attr(normalized_data(),"cnames")$fdata_cname)
  })
  # Create a plot of raw abundance vs raw richness
  output$ra_raw <- renderPlotly({
    #plot(abun_raw(), rich_raw(), plot_title="Raw Data")
    p <- plotly::ggplotly(ra_raw_plot(), tooltip = c("x", "y","colour", "label"))
    p$elementId <- NULL
    p
  })

  ra_norm_plot <- reactive({
    plot(abun_norm(), rich_norm(), plot_title="Normalized Data", samplabel = attr(normalized_data(),"cnames")$fdata_cname)
  })
  # Create a plot of normalized abundance vs normalized richness to see if there is a reduction in correlation
  output$ra_norm <- renderPlotly({
    #plot(abun_norm(), rich_norm(), plot_title="Normalized Data")
    p <- plotly::ggplotly(ra_norm_plot(), tooltip = c("x", "y", "colour", "label"))
    p$elementId <- NULL
    p
  })

  ################ Community Metrics Tab #################

  # Select what variable to put on x-axis in community metrics plots
  output$xaxis <- renderUI({
    selectInput("xaxis",
                label = "x-axis",
                choices = c(colnames(attr(normalized_data(),"group_DF"))),
                selected = "Group")
  })

  # Select what variable to color by in community metrics plots
  output$color <- renderUI({
    selectInput("color",
                label = "color",
                choices = c(colnames(attr(normalized_data(),"group_DF"))),
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

    return(pmartRseq::alphaDiv_calc(normalized_data(), index=input$adiv_index))
  })

  adiv_plot_obj <- reactive({
    plot(a_div(), x_axis=input$xaxis, color=input$color, samplabel=input$new_f_data_cname)
  })
  # Show alpha diversity plot
  output$adiv_plot <- renderPlotly({
    #plot(a_div(), x_axis=input$xaxis, color=input$color)
  p <-  plotly::ggplotly(plot(a_div(), x_axis=input$xaxis, color=input$color, scales = 'free', samplabel = attr(normalized_data(),"cnames")$fdata_cname),
                         tooltip = c("x", "y", "label")) #ggplotly bugs without free scale
  p$elementId <- NULL
  p
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

    return(pmartRseq::richness_calc(filtered_data(), index=input$rich_index))
  })

  rich_plot_obj <- reactive({
    plot(rich(), x_axis=input$xaxis, color=input$color, scales = 'free', samplabel = attr(normalized_data(),"cnames")$fdata_cname)#ggplotly bugs without free scale
  })
  # Show richness plot
  output$rich_plot <- renderPlotly({
    #plot(rich(), x_axis=input$xaxis, color=input$color)
    p <- plotly::ggplotly(rich_plot_obj(), tooltip = c("x", "y", "label"))
    p$elementId <- NULL
    p
  })

  output$rich_summary <- renderPrint({
    summary(rich())
  })

  #----------- abundance example ----------#

  # Calculate abundance
  abun <- reactive({
    validate(
      need(length(input$xaxis) > 0, "Must include a valid x-axis for the plot.")
    )

    return(pmartRseq::abundance_calc(normalized_data()))
  })

  abun_plot_obj <- reactive({
    plot(abun(), x_axis=input$xaxis, color=input$color,  samplabel = attr(normalized_data(),"cnames")$fdata_cname)
  })
  # Show evenness plot
  output$abun_plot <- renderPlotly({
    #plot(even(), x_axis=input$xaxis, color=input$color)
    p <- plotly::ggplotly(abun_plot_obj(), tooltip = c("x", "y", "label"))
    p$elementId <- NULL
    p
  })

  output$abun_summary <- renderPrint({
    summary(abun())
  })

  #----------- evenness example ----------#

  # Which evenness indices to calculate
  output$even_index <- renderUI({
    checkboxGroupInput("even_index",
                       label = "Evenness Index",
                       choices = list("Shannon"="shannon","Simpson"="simpson"),
                       selected = c("shannon","simpson"))
  })

  # Calculate evenness
  even <- reactive({
    validate(
      need(length(input$even_index) > 0, "There needs to be at least one evenness index")
    )

    return(pmartRseq::evenness_calc(normalized_data(), index=input$even_index))
  })

  even_plot_obj <- reactive({
    plot(even(), x_axis=input$xaxis, color=input$color, scales = 'free',  samplabel = attr(normalized_data(),"cnames")$fdata_cname)#ggplotly bugs without free scale
  })
  # Show evenness plot
  output$even_plot <- renderPlotly({
    #plot(even(), x_axis=input$xaxis, color=input$color)
    p <- plotly::ggplotly(even_plot_obj(), tooltip = c("x", "y", "label"))
    p$elementId <- NULL
    p
  })

  output$even_summary <- renderPrint({
    summary(even())
  })

  #----------- effective species example ----------#

  # Calculate effsp
  effsp <- reactive({
    validate(
      need(length(input$xaxis) > 0, "Must include a valid x-axis for the plot.")
    )

    return(pmartRseq::effsp_calc(normalized_data()))
  })

  effsp_plot_obj <- reactive({
    plot(effsp(), x_axis=input$xaxis, color=input$color, samplabel = attr(normalized_data(),"cnames")$fdata_cname)
  })
  # Show evenness plot
  output$effsp_plot <- renderPlotly({
    #plot(effsp(), x_axis=input$xaxis, color=input$color)
    p <- plotly::ggplotly(effsp_plot_obj(), tooltip = c("x", "y", "label"))
    p$elementId <- NULL
    p
  })

  output$effsp_summary <- renderPrint({
    summary(effsp())
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

  # observeEvent(input$submit_goe, {
  #   dimcheck_obj <<- reactive({
  #     goeveg::dimcheckMDS(vegdata(), distance = input$beta_index, autotransform = FALSE)
  #   })
    output$dimcheck <- renderPlot({
      #goeveg::dimcheckMDS(vegdata(), distance = input$beta_index, autotransform = FALSE)
      req(input$submit_goe)
      Sys.sleep(5)
      # print(dimcheck_obj())
      print(goeveg::dimcheckMDS(vegdata(), distance = input$beta_index, autotransform = FALSE))
    })
  # })

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
                 value = 2)
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
  })
    # Plot showing beta diversity
    output$ord_plot <- renderPlotly({
      #if(input$ord_method == "NMDS"){
      # pmartRseq::pmartRseq_NMDS(res = vegmds(),
      #           grp = as.factor(attr(normalized_data(),"group_DF")[match(rownames(vegdata()), attr(normalized_data(),"group_DF")[,attr(normalized_data(),"cnames")$fdata_cname]),input$ord_colors]),ellipses=input$ellipses)
      #pmartRseq::pmartRseq_NMDS(res = vegmds(), omicsData = normalized_data(), grp = input$ord_colors, k = input$k,
      #                         x_axis = input$ord_x, y_axis = input$ord_y, ellipses=input$ellipses)
      req(input$submit_ord)
      #Sys.sleep(5)
      #print(ord_plot_obj())
      m <- ord_plot_obj() + theme(aspect.ratio=NULL)
      p <- plotly::ggplotly(m, height = "100%", width = "100%")
      p$elementId <- NULL
      p
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

  output$pval_thresh <- renderUI({
    sliderInput("pval_thresh",
                label = "P-value significance threshold",
                min = 0,
                max = 0.5,
                step = 0.001,
                value = 0.05)
  })

  output$comparisons <- renderUI({
    checkboxGroupInput("comparisons",
                       label = "Differential abundance pairwise comparisons",
                       choices = lapply(c(1:(factorial(length(unique(attr(filtered_data(),"group_DF")$Group)))/(factorial(2)*factorial(length(unique(attr(filtered_data(),"group_DF")$Group))-2)))), function(x) paste(combn(unique(attr(filtered_data(),"group_DF")$Group),2)[,x], collapse="  VS  ")))
  })

  observe({
    if(input$selectall == 0){
      return(NULL)
    }else if(input$selectall%%2 == 0){
      updateCheckboxGroupInput(session, "comparisons","Differential abundance pairwise comparisons", choices=lapply(c(1:(factorial(length(unique(attr(filtered_data(),"group_DF")$Group)))/(factorial(2)*factorial(length(unique(attr(filtered_data(),"group_DF")$Group))-2)))), function(x) paste(combn(unique(attr(filtered_data(),"group_DF")$Group),2)[,x], collapse="  VS  ")))
    }else{
      updateCheckboxGroupInput(session, "comparisons","Differential abundance pairwise comparisons", choices=lapply(c(1:(factorial(length(unique(attr(filtered_data(),"group_DF")$Group)))/(factorial(2)*factorial(length(unique(attr(filtered_data(),"group_DF")$Group))-2)))), function(x) paste(combn(unique(attr(filtered_data(),"group_DF")$Group),2)[,x], collapse="  VS  ")), selected=sapply(c(1:(factorial(length(unique(attr(filtered_data(),"group_DF")$Group)))/(factorial(2)*factorial(length(unique(attr(filtered_data(),"group_DF")$Group))-2)))), function(x) paste(combn(unique(attr(filtered_data(),"group_DF")$Group),2)[,x], collapse="  VS  ")))
    }
  })

  #observeEvent(input$submit_da, {
    # Calculate normalization factors to use in differential abundance test - will use the same that was used on normalization tab
    norm_factors <- reactive({
      validate(need(length(input$normFunc) == 1, "Need to specify a normalization function."))
      validate(need(input$normFunc %in% c("percentile","tss","rarefy","poisson","deseq","tmm","css","log","clr","none"), "Normalization function must be one of the options specified."))

      return(pmartRseq::normalize_data(omicsData=filtered_data(), norm_fn=input$normFunc, normalize=FALSE))
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
      req(input$submit_da)
      validate(need(length(input$da_index) == 1, "Need to specify a differential abundance test"))
      validate(need(length(input$pval_adjust) == 1, "Need to specify a p-value adjustment method"))

      if(input$normFunc %in% c("rarefy","log","clr")){
        return(pmartRseq::countSTAT(omicsData = normalized_data(), norm_factors = norm_factors()$scale_param, comparisons = comps(), control = NULL, test = input$da_index, pval_adjust = input$pval_adjust, pval_thresh = 0.05))
      }else{
        return(pmartRseq::countSTAT(omicsData = filtered_data(), norm_factors = norm_factors()$scale_param, comparisons = comps(), control = NULL, test = input$da_index, pval_adjust = input$pval_adjust, pval_thresh = 0.05))
      }

    })

    # Look at the results - this is hard to look at, maybe remove?
    output$da_res <- renderDataTable({
      req(input$submit_da)
      diffabun_res()$allResults
      })

    output$da_summary <- renderPrint({
      req(input$submit_da)
      summary(diffabun_res())
    })

    da_flag_plot_obj <<- reactive({
      plot(diffabun_res(), type = "flag")
    })
    # Plot showing number differentially abundant in each comparison and direction of change
    output$flag_plot <- renderPlot({
      req(input$submit_da)

      #plot(diffabun_res(), type = "flag")
      print(da_flag_plot_obj())
    })

    da_logfc_plot_obj <<- reactive({
      plot(diffabun_res(), type = "logfc")
    })
    # Heatmap showing the log2foldchanges of differentially abundant features
    output$logfc_plot <- renderPlot({
      req(input$submit_da)
      #plot(diffabun_res(), type = "logfc")
      print(da_logfc_plot_obj())
    })

    plot_all_da_obj <<- reactive({
      pmartRseq::plot_all_diffabun(countSTAT_results = diffabun_res(), omicsData = normalized_data(), x_axis = "taxonomy2", x_lab = "Phylum")
    })
    # Plot showing log fold changes and p-values of all features, grouped by taxonomy
    output$plot_all_da <- renderPlot({
      req(input$submit_da)
      #pmartRseq::plot_all_diffabun(countSTAT_results = diffabun_res(), omicsData = normalized_data(), x_axis = "taxonomy2", x_lab = "Phylum")
      print(plot_all_da_obj())
    })
  #}, autoDestroy = FALSE)

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

  #observeEvent(input$submit_is, {
    indsp_res <<- reactive({
      req(input$submit_is)
      if(input$within == "NA"){
        return(pmartRseq::indsp_calc(omicsData = normalized_data(), within = NULL, pval_thresh = input$is_pval_thresh))
      }else{
        return(pmartRseq::indsp_calc(omicsData = normalized_data(), within = input$within, pval_thresh = input$is_pval_thresh))
      }
    })

    output$indsp_results <- renderDataTable({
      req(input$submit_is)
      indsp_res()
      })

    output$indsp_summary <- renderPrint({
      req(input$submit_is)
      summary(indsp_res())
    })


    output$indsp_xaxis <- renderUI({
      req(input$submit_is)
      selectInput("indsp_xaxis",
                  label = "X-Axis for Indicator Species Plot",
                  choices = colnames(attr(normalized_data(),"group_DF"))[-which(colnames(attr(normalized_data(),"group_DF"))==attr(normalized_data(),"cnames")$fdata_cname)],
                  selected = "Group")
    })

    output$indsp_group <- renderUI({
      req(input$submit_is)
      selectInput("indsp_group",
                  label = "Fill Variable for Indicator Species Plot",
                  choices = colnames(normalized_data()$e_meta),
                  selected = colnames(normalized_data()$e_meta)[3])
    })

    indsp_plot_obj <<- reactive({
      validate(need(!is.null(indsp_res()), message = "Sumbit analysis"),
               need(!is.null(input$indsp_xaxis), message = "Please wait...getting x-axis together"))
      pmartRseq::plot_indsp(indsp = indsp_res(), omicsData = normalized_data(), x_axis = , group = input$indsp_group)
    })
    output$indsp_plot <- renderPlotly({
      validate(need(!is.null(indsp_plot_obj()), message = "Sumbit analysis"))
      req(input$submit_is)
      #pmartRseq::plot_indsp(indsp = indsp_res(), omicsData = normalized_data(), x_axis = input$indsp_xaxis, group = input$indsp_group)
      m <- indsp_plot_obj()+theme(axis.text.y = element_blank(), aspect.ratio=NULL)
      p <- plotly::ggplotly(m, tooltip = c("x", "y", "fill"), width = "100%", height = "100%")
      p$elementId <- NULL
      p
    })
  #}, autoDestroy = FALSE)


  # ################ Stats Results Tab #################
  # #if(exists(indsp_res()) & exists(diffabun_res())){
  # #----------- differential abundance ----------#
  # diffres <- reactive({
  #   t1 <- diffabun_res()$allResults[,grep("Flag",colnames(diffabun_res()$allResults))]
  #   if(is.null(ncol(t1))){
  #     t1 <- data.frame(OTU=rownames(diffabun_res()$allResults), t1)
  #     rownames(t1) <- t1$OTU
  #     colnames(t1)[1] <- attr(normalized_data(), "cnames")$edata_cname
  #     colnames(t1)[2] <- colnames(diffabun_res()$allResults)[grep("Flag",colnames(diffabun_res()$allResults))]
  #   }else{
  #     t1 <- t1[-which(rowSums(abs(t1), na.rm=TRUE) == 0),]
  #   }
  #   return(t1)
  # })
  #
  # #----------- indicator species ----------#
  # isres <- reactive({
  #   t1 <- indsp_res()[,grep("Flag", colnames(indsp_res()))]
  #   if(is.null(ncol(t1))){
  #     t1 <- data.frame(OTU=rownames(indsp_res()), t1)
  #     rownames(t1) <- t1$OTU
  #     colnames(t1)[1] <- attr(normalized_data(), "cnames")$edata_cname
  #     colnames(t1)[2] <- colnames(indsp_res())[grep("Flag",colnames(indsp_res()))]
  #   }else{
  #     t1 <- t1[-which(rowSums(abs(t1), na.rm=TRUE) == 0),]
  #   }
  #   return(t1)
  # })
  #
  # #----------- combine results ----------#
  # statsres <- reactive({
  #   if(is.null(attr(indsp_res(), "within"))){
  #     res <- intersect(rownames(diffres()), rownames(isres()))
  #     res <- data.frame(res)
  #     colnames(res)[1] <- attr(normalized_data(), "cnames")$edata_cname
  #     return(res)
  #   }else{
  #     res <- lapply(unique(normalized_data()$f_data[,attr(indsp_res(), "within")]), function(x){
  #       da <- diffres()[,grep(x, colnames(diffres()))]
  #       da <- da[-which(rowSums(abs(da)) == 0),]
  #       is <- isres()[,grep(x, colnames(isres()))]
  #       is <- data.frame(OTU=rownames(isres()),is)
  #       rownames(is) <- rownames(isres())
  #       is <- is[which(abs(is[,-1]) > 0),]
  #       res <- intersect(rownames(da), rownames(is))
  #       res <- data.frame(Within=x, OTU=res)
  #       return(res)
  #     })
  #     res <- do.call(rbind, res)
  #     colnames(res)[which(colnames(res) == "Within")] <- attr(indsp_res(), "within")
  #     colnames(res)[which(colnames(res) == "OTU")] <- attr(normalized_data(), "cnames")$edata_cname
  #     return(res)
  #   }
  # })
  #
  # taxares <- reactive({
  #   res <- merge(statsres(), normalized_data()$e_meta, by=attr(normalized_data(), "cnames")$edata_cname)
  #   return(res)
  # })
  #
  # output$stats_res <- DT::renderDataTable(taxares())
  # output$newisplot <- renderPlot({
  #   pmartRseq::plot_indsp(indsp = indsp_res(), omicsData = normalized_data(), x_axis = input$indsp_xaxis, group = input$indsp_group)
  # })
  # output$newdaplot <- renderPlot({
  #   pmartRseq::plot_all_diffabun(countSTAT_results = diffabun_res(), omicsData = normalized_data(), x_axis = "taxonomy2", x_lab = "Phylum")
  # })
  # #}

  ################ ALDEx2 Tab #################
  #----------- aldex2 ----------#

  # Select which variables to use as main effects
  output$pa_mainEffects <- renderUI({
      selectInput("pa_mainEffects",
                  label = "Main Effect(s) to use in Model",
                  choices = c("NA"="NULL",colnames(filtered_data()$f_data)[group_vars()]),
                  selected = "NULL",
                  multiple = TRUE)
  })

  # Select which variables to use as random effects
  output$pa_randomEffect <- renderUI({
    selectInput("pa_randomEffect",
                label = "Optional, Random Effect to use in Model",
                choices = c("NA"="NULL",colnames(filtered_data()$f_data)[group_vars()]),
                selected = "NULL",
                multiple = TRUE)
  })

  output$pa_Interactions <- renderUI({
    checkboxInput("pa_Interactions",
                  label = "Include Interactions",
                  value=FALSE)
  })

  output$mcsamples <- renderUI({
    numericInput("mcsamples",
                 label="Number of Monte Carlo Samples",
                 value=128)
  })


  #observeEvent(input$submit_pa, {

    # Perform differential abundance analysis
    pa_results <<- reactive({
      req(input$submit_pa)
      validate(need(length(input$mcsamples) == 1, "Need to specify number of Monte Carlo samples"))
      #validate(need(length(input$pval_adjust) == 1, "Need to specify a p-value adjustment method"))

      if("NULL" %in% input$pa_mainEffects){
          pa_mE <- NULL
      }else{
        pa_mE <- input$pa_mainEffects
      }

      if("NULL" %in% input$pa_randomEffect){
        pa_rE <- NULL
      }else{
        pa_rE <- input$pa_randomEffect
      }

      return(pmartRseq::pmartRseq_aldex2(omicsData = filtered_data(), mainEffects = pa_mE, randomEffect = pa_rE, interactions = input$pa_Interactions, mc.samples = input$mcsamples))

    })

    # Look at the results - this is hard to look at, maybe remove?
    output$pa_res <- renderDataTable({
      req(input$submit_pa)
      pa_results()$results
      })

    output$pa_summary <- renderPrint({
      req(input$submit_pa)
      summary(pa_results())
    })

    pa_pval_plot_obj <<- reactive({
      plot(pa_results(), type = "pvals")
    })
    # Heatmap showing the log2foldchanges of differentially abundant features
    output$pa_pval_plot <- renderPlot({
      req(input$submit_pa)
      #plot(diffabun_res(), type = "logfc")
      print(pa_pval_plot_obj())
    })

    pa_flag_plot_obj <<- reactive({
      plot(pa_results(), type = "flag")
    })
    # Plot showing number differentially abundant in each comparison and direction of change
    output$pa_flag_plot <- renderPlot({
      req(input$submit_pa)
      #plot(diffabun_res(), type = "flag")
      print(pa_flag_plot_obj())
    })

  #}, autoDestroy = FALSE)

  ################ Network Analysis Tab #################
  #----------- network ----------#

  # # Select which variables to use as main effects
  # output$na_corrtype <- renderUI({
  #   selectInput("na_corrtype",
  #               label = "Correlation function to use for network",
  #               choices = c("spearman","pearson"),
  #               selected = "spearman",
  #               multiple = FALSE)
  # })

  # Select if should run network analysis in groups
  output$na_group <- renderUI({
    checkboxInput("na_group",
                  label = "Run separate networks for each group?",
                  value = FALSE)
  })

  # Select which variable to use for grouping
  output$na_group_var <- renderUI({
    selectInput("na_group_var",
                label = "Which variable to use for grouping",
                choices = c("NA",colnames(normalized_data()$f_data)[group_vars()]),
                selected = "NA",
                multiple = FALSE)
  })

  # Select if should use 0 or NA for missing values
  output$na_missingval <- renderUI({
    selectInput("na_missingval",
                  label = "Which to use for 'missing' values",
                  choices = c("NA",0),
                  selected = 0,
                  multiple = FALSE)
  })

  # igraph corr coeff cutoff
  output$na_coeff <- renderUI({
    sliderInput("na_coeff",
                label = "Absolute value of correlation coefficient must be higher than this value to be included in plot",
                min = 0,
                max = 1,
                step = 0.01,
                value = 0.3)
  })

  # # igraph corr coeff cutoff
  # output$na_pval <- renderUI({
  #   sliderInput("na_pval",
  #               label = "Correlation pvalue must be lower than this value to be included in plot",
  #               min = 0,
  #               max = 1,
  #               step = 0.01,
  #               value = 0.05)
  # })

  # igraph corr coeff cutoff
  output$na_qval <- renderUI({
    sliderInput("na_qval",
                label = "Correlation qvalue must be lower than this value to be included in plot",
                min = 0,
                max = 1,
                step = 0.01,
                value = 0.05)
  })

  # colour of vertices in graph
  output$na_colour <- renderUI({
    selectInput("na_colour",
                label = "Which variable to use for vertex colour",
                choices = c("NA",colnames(normalized_data()$e_meta)),
                selected = "NA")
  })

  # Select if should scale vertex size by abundance
  output$na_size <- renderUI({
    checkboxInput("na_size",
                  label = "Scale vertex size by abundance?",
                  value = FALSE)
  })

  # Select module algorithm to use
  output$na_cluster <- renderUI({
    selectInput("na_cluster",
                label = "Which clustering method to use",
                choices = c("edge_betweenness", "fast_greedy", "infomap", "label_prop", "leading_eigen", "louvain", "optimal", "spinglass", "walktrap"),
                selected = "louvain")
  })

  # Select min module size
  output$na_mod_size <- renderUI({
    numericInput("na_mod_size",
                 label = "Modules must have a minimum of this many features, otherwise will be grouped into 'non-modular'",
                 value = 5)
  })

  # Select environmental variables to correlate to modules
  output$na_envvars <- renderUI({
    selectInput("na_envvars",
                label = "Which environmental variables to correlate to modules?",
                choices = c("NA",colnames(normalized_data()$f_data)),
                selected = "NA",
                multiple = TRUE)
  })

  # Select module envvars correlation p-value
  output$env_pval <- renderUI({
    sliderInput("env_pval",
                label = "Highlight p-values lower than this value in plot",
                min = 0,
                max = 1,
                step = 0.01,
                value = 0.05)
  })


  #observeEvent(input$submit_na, {

    # Perform differential abundance analysis
    na_network <<- reactive({
      req(input$submit_na)
      validate(need(length(input$na_missingval) == 1, "Need to specify what value to use for missing values"))
      #validate(need(length(input$pval_adjust) == 1, "Need to specify a p-value adjustment method"))

      if(!input$na_group){
        na_group_var <- NULL
        na_group <- FALSE
      }else{
        na_group <- TRUE
        if(input$na_group_var == "NA"){
          na_group_var <- NULL
        }else{
          na_group_var <- input$na_group_var
        }
      }

      return(pmartRseq::network_calc(omicsData = normalized_data(), type="spearman", group=na_group, group_var=na_group_var, fdr_method="fndr", missing_val=input$na_missingval))
    })

    na_igraph <<- reactive({
      req(input$submit_na)
      validate(need(length(input$na_missingval) == 1, "Need to specify what value to use for missing values"))

      return(pmartRseq::pmartRseq_igraph(netData = na_network(), coeff = input$na_coeff, qval = input$na_qval, pval = NULL))
    })

    na_network_plot <<- reactive({
      req(input$submit_na)
      if(input$na_colour != "NA"){
        na_colour <- input$na_colour
      }else{
        na_colour <- NULL
      }
      return(pmartRseq::network_plot(netGraph=na_igraph(), omicsData=normalized_data(), modData=NULL, colour=na_colour, vsize=input$na_size, legend.show=TRUE, legend.pos="bottomleft"))
    })

    # Network Plot
    output$na_network_plot <- renderPlot({
      req(input$submit_na)
      #plot(diffabun_res(), type = "logfc")
      print(na_network_plot())
    })
 # }, autoDestroy = FALSE)

  #observeEvent(input$submit_modules, {
    # Network indices
    na_net_indc <<- reactive({
      req(input$submit_modules)
      validate(need(exists(na_igraph()), "network graph object must be created first"))

      return(pmartRseq::network_indices(netGraph = na_igraph()))
    })

    # Modules
    na_mods <<- reactive({
      req(input$submit_modules)
      validate(need(length(input$na_cluster) == 1, "Need to specify a clustering algorithm to use"))
      validate(need(length(input$na_mod_size) == 1, "Need to specify a minimum module size"))

      return(pmartRseq::detect_modules(netGraph = na_igraph(), cluster = input$na_cluster, cutoff = input$na_mod_size))
    })

    na_env <<- reactive({
      req(input$submit_modules)
      validate(need(length(input$na_envvars) >= 1, "Need to specify environmental variables"))

      return(pmartRseq::mod_env(omicsData = normalized_data(), modData = na_mods(), envVars = input$na_envvars, pca.method="svd", cor.method="spearman", use="pairwise", padjust="BH"))
    })

    mod_plot <<- reactive({
      req(input$submit_modules)
      validate(need(length(input$na_size) == 1, "Need to specify vertex size"))

      return(pmartRseq::network_plot(netGraph = na_igraph(), omicsData = normalized_data(), modData = na_mods(), colour = "Module", vsize = input$na_size, legend.show=TRUE, legend.pos = "bottomleft"))
    })

    output$na_mod_plot <- renderPlot({
      req(input$submit_modules)
      print(mod_plot())
    })
  #},autoDestroy = FALSE)

  #observeEvent(input$submit_envvars, {
    mod_env_plot <<- reactive({
      req(input$submit_envvars)
      validate(need(length(input$env_pval) == 1, "Need to specify a p-value cutoff for module and environmental variable correlations"))

      return(plot(na_env(), pval.thresh=input$env_pval))
    })

    # EnvVars Plot
    output$na_envvars_plot <- renderPlot({
      req(input$submit_envvars)
      #plot(diffabun_res(), type = "logfc")
      print(mod_env_plot())
    })

 # },autoDestroy = FALSE)

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
        rep$data <- rRNA_agg()
        write.csv(rRNA_agg()$e_data, file="raw.csv")
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
