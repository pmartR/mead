
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
  #---------- Outputs -------------------------#
 
#   output$lib_size_hist = reactive({
#     xlab = "Number of Reads (Counts)"
#     ylab = "Number of Libraries"
#     if (!is.null(full_data())){
#       plot_data <- taxa_sums(full_data())
#       p <- sums_hist(plot_data, xlab, ylab)
#     }
#     if (!is.null(full_data())){
#       p <- NULL
#     }
#     
#     return(p)
#   })
#   
#   output$content <- renderText({
#     if (!is.null(full_data())){
#       print(nsamples(full_data())) 
#     }
#     else (NULL)
#   })
# 
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
  
  output$downloadRichnessEstimates <- downloadHandler(
    filename = "alpha_diversity_estimates.csv",
    content = function(file) {
    write.csv(estimate_richness(pruned_data(), split = as.logical(input$split), measures =input$dist_measures), file)
  })
#   
  output$richness_plot <- renderPlot({
    if (!is.null(pruned_data())) {
      plot_data <- pruned_data()
    return(plot_richness(plot_data, measures = input$dist_measures))
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