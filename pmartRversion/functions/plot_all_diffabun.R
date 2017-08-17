#' Plot all differential abundance results
#'
#' This function takes the differential expression results of a mintR data object and creates a plot showing all the log2 fold changes of all features and the corresponding p-values. This looks at the entire dataset to help determine shifts in the whole community.
#'
#' @param countSTAT_results Required, an object of the class \code{countSTAT} created by \code{\link{countSTAT}}.
#' @param omicsData Required, omicsData ,an object of the class \code{cDNAdata}, \code{gDNAdata}, or \code{rRNAdata} usually created by \code{\link{as.cDNAdata}}, \code{\link{as.gDNAdata}}, or \code{\link{as.rRNAdata}}, respectively.
#' @param x_axis Required, a character vector stating which variable to group data by and put on the x-axis, must be one of "Comparison" or one of the column names in omicsData$e_meta. Default is "Phylum".
#'@param facet Optional, a character vector depicting the formula to use when faceting the plots. Default is "~Comparison".
#'@param scales Optional, a character vector to describe if any of the axes should be free when faceting. Default is "fixed".
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
#'
#' @return A plot showing log2 fold change and p-values for all comparisons and all features
#'
#' @examples
#' library(mintR)
#' data("rRNAdata")
#'
#' rRNAdata <- group_designation(omicsData=rRNA_data, main_effects="Treatment")
#' norm_factors <- normalize_data(omicsData=rRNAdata, subset_fn="none", norm_fn="percentile", normalize=FALSE)
#' norm_data <- normalize_data(omicsData=rRNAdata, subset_fn="none", norm_fn="percentile", normalize=TRUE)
#' diffexp <- countSTAT(omicsData=rRNAdata, norm_factors=norm_factors, comparisons="all", test="dw", pval_adjust="none", pval_thresh=0.05)
#'
#' plot_all_diffabun(countSTAT_results=diffexp, omicsData=norm_data, x_axis="Phylum", facet="~Comparison", scales="fixed")
#'
#' @author Allison Thompson
#'
#' @export
#'
plot_all_diffabun <- function(countSTAT_results, omicsData, x_axis="Phylum", facet="~Comparison", scales="fixed",
                                x_lab=NULL, y_lab=NULL, leglab=NULL, plot_title=NULL){

  library(ggplot2)
  library(reshape2)

  if(is.null(x_axis)){
    stop("x_axis cannot be NULL")
  }

  if(!(x_axis %in% c("Comparison", colnames(omicsData$e_meta)))){
    stop("x_axis must be either 'Comparison' or one of the column names in omicsData$e_meta")
  }

  if(!is.null(scales)){
    if(!(scales %in% c("free_y","free_x","free","fixed"))){
      stop("scales must be one of 'free_y', 'free_x', 'free', or 'fixed'.")
    }
  }

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }

  # Extract differential expression results
  data <- countSTAT_results$allResults

  # Create a new plot for every test that was run
 lapply(attr(countSTAT_results, "Tests")$Test, function(t){
    # Extract only data for specific test
    data <- data[,grep(t, colnames(data))]

    # Extract and format adjusted p-values
    data.padj <- data[,grep("padj", colnames(data))]
    data.padj <- data.frame(Temp=rownames(data), data.padj)
    colnames(data.padj)[1] <- attr(omicsData, "cnames")$edata_cname
    colnames(data.padj)[-1] <- colnames(data)[grep("padj",colnames(data))]
    data.padj <- melt(data.padj)
    colnames(data.padj)[3] <- "padj"
    data.padj$Comparison <- unlist(lapply(data.padj$variable, function(x) strsplit(as.character(x),"padj_")[[1]][2]))

    # Extract and format log2 fold change values
    data.logfc <- data[,grep("logFC", colnames(data))]
    data.logfc <- data.frame(Temp=rownames(data), data.logfc)
    colnames(data.logfc)[1] <- attr(omicsData, "cnames")$edata_cname
    colnames(data.logfc)[-1] <- colnames(data)[grep("logFC",colnames(data))]
    data.logfc <- melt(data.logfc)
    colnames(data.logfc)[3] <- "logFC"
    data.logfc$Comparison <- unlist(lapply(data.logfc$variable, function(x) strsplit(as.character(x),"logFC_")[[1]][2]))

    # Merge adjusted p-values and log2fc data
    plot.data <- merge(data.padj, data.logfc, by=c(attr(omicsData, "cnames")$edata_cname, "Comparison"))
    # Merge data with e_meta
    plot.data <- merge(plot.data, omicsData$e_meta, by=attr(omicsData, "cnames")$edata_cname)

    map <- aes_string(x=x_axis, y="logFC", colour="padj")
    p <- ggplot(plot.data, map) +
      geom_point()+
      geom_hline(yintercept=log2(2),lty="dashed")+
      geom_hline(yintercept=-log2(2),lty="dashed")+
      geom_hline(yintercept=log2(10),lty="dotted")+
      geom_hline(yintercept=-log2(10),lty="dotted")+
      theme_bw()+
      theme(axis.text.x=element_text(angle=90),
            axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour="black"),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank())

    if(!is.null(x_lab)){
      p <- p + labs(x=x_lab)
    }else{
      p <- p + labs(x=x_axis)
    }

    if(!is.null(y_lab)){
      p <- p + labs(y=y_lab)
    }else{
      p <- p + labs(y="Log2 Fold Change")
    }

    if(!is.null(plot_title)){
      p <- p + ggtitle(plot_title)
    }else{
      p <- p + ggtitle(paste(t," Results"))
    }

    if(!is.null(leglab)){
      p <- p + guides(colour=guide_legend(title=leglab))
    }

    if(!is.null(facet)){
      p <- p + facet_wrap(as.formula(facet), scales=scales)
    }

    return(p)

  })

}
