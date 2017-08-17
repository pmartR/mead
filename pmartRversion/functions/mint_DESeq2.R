#' DESeq2 analysis of mintR objects
#'
#' Differential abundance analysis of count data using DESeq2
#'
#' @param omicsData an object of the class 'gDNAdata', 'cDNAdata', or 'rRNAdata' usually created by \code{\link{as.gDNAdata}}, \code{\link{as.cDNAdata}},or \code{\link{as.rRNAdata}}, respectively.
#' @param norm_factors Named vector of normalization parameters to put into DESeq2. If NULL, will use DESeq2's inbuilt normalization. Default is NULL.
#' @param pairs a matrix dictating which pairwise comparisons to make
#' @param test name of which differential expression test to use, options are "wald" or "paired". Default is "wald". See details for further explanation.
#' @param adj multiple comparison adjustment method to use.
#' @param thresh p-value threshold for significance.
#'
#' @details Performs differential abundance testing on two groups using DESeq2.
#'
#' @return DESeqREsults object, which is a simple subclass of DataFrame. Columns include baseMean, log2FoldChange, lfcSE (standard error of log2FoldChange), stat (Wald statistic), pvalue, and padj (BH adjusted p-values). Use mcols(res)$description.
#'
#' @examples
#' library(mintJansson)
#' data(cDNA_hiseq_data)
#' mycdnadata <- group_designation(omicsData = cDNA_hiseq_data, main_effects = c("treatment"), time_course=NULL)
#' mycdnadata_norm <- normalize_data(omicsData = mycdnadata, norm_fn = "percentile")
#' mycdnadata_DESeq2 <- mint_DESeq2(omicsData = mycdnadata_norm, test="wald", pairs = cbind(list("Neg","Plus")), adj = "BH", thresh = 0.05)
#'
#' @author Allison Thompson
#'
#' @references Love MI, Huber W and Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, pp. 550.
#'

mint_DESeq2 <- function(omicsData, norm_factors=NULL, test="wald", pairs, adj, thresh){

  library(DESeq2)

  if(!(class(omicsData) %in% c("gDNAdata","cDNAdata","rRNAdata"))){
    stop("Data must be of class gDNAdata, cDNAdata, or rRNAdata")
  }

  # Make features be rownames and remove that column from the data, so the data is all numeric
  rownames(omicsData$e_data) <- omicsData$e_data[,attr(omicsData,"cnames")$edata_cname]
  omicsData$e_data <- omicsData$e_data[,-which(colnames(omicsData$e_data) == attr(omicsData,"cnames")$edata_cname)]

  # Change NA's to 0 for DESeq2
  omicsData$e_data[is.na(omicsData$e_data)] <- 0

  # Check if group1 and/or group2 is null
  if(is.null(attributes(omicsData)$group_DF)){
    stop("run group designation function first")
  }


  if(!("paired" %in% test)){
    # Set the metadata (group indication) object for DESeq2
    colData <- attributes(omicsData)$group_DF
    colData <- colData[match(colnames(omicsData$e_data), colData[,attr(omicsData,"cnames")$fdata_cname]),]

    # Create the appropriate data object
    DESeq_data_set_object <- DESeqDataSetFromMatrix(countData = omicsData$e_data, colData = colData, design = ~Group)
  }else{
    colData <- attributes(omicsData)$group_DF
    colData <- colData[match(colnames(omicsData$e_data), colData[,attr(omicsData,"cnames")$fdata_cname]),]
    DESeq_data_set_object <- DESeqDataSetFromMatrix(countData=omicsData$e_data,colData=colData,design=~1)
  }

  # Set own sizeFactors(normalization factors), after running normalization
  if(is.null(norm_factors)){
    warning("Using DESeq2's inbuilt normalization - if want to use own normalization, input norm_factors")
    DESeq_data_set_object <- estimateSizeFactors(DESeq_data_set_object)
  }else{
    # check that normalization (size) factors are present for every sample
    if(all(colnames(omicsData$e_data) %in% names(norm_factors))){
      norm_factors <- norm_factors[match(colnames(omicsData$e_data), names(norm_factors))]
      sizeFactors(DESeq_data_set_object) <- norm_factors
    }else{
      stop("Must have a normalization (size) factor for every sample in the data")
    }
  }

  # DESeq2 analysis
  if(tolower(test)=="wald"){
    DESeq_data_set_object <- DESeq(DESeq_data_set_object, test="Wald", quiet=TRUE)
    # extract results for each pairwise comparison
    res <- apply(pairs, 2, function(x){
      x <- as.character(x)
      myContrast <- c("Group", x[1], x[2])
      res <- results(DESeq_data_set_object, contrast=myContrast, pAdjustMethod=adj)
      res$Flag <- ifelse(res$padj <= thresh, 1, 0)
      res$Flag[which(res$log2FoldChange < 0)] <- res$Flag[which(res$log2FoldChange < 0)] * -1
      names(res) <- unlist(lapply(names(res), function(y) paste(y,"_",x[1],"_vs_",x[2],sep="")))
      res <- res[order(rownames(res)), ]
      return(res)
    })
    # combine results from all pairwise comparisons
    res2 <- do.call(cbind, res)
  }else if(tolower(test)=="lrt"){
    DESeq_data_set_object <- DESeq(DESeq_data_set_object, test="LRT", quiet=TRUE, reduced = ~ 1)
    # extract results for each pairwise comparison
    res <- apply(pairs, 2, function(x){
      x <- as.character(x)
      myContrast <- c("Group", x[1], x[2])
      res <- results(DESeq_data_set_object, contrast=myContrast, pAdjustMethod=adj)
      res$Flag <- ifelse(res$padj <= thresh, 1, 0)
      res$Flag[which(res$log2FoldChange < 0)] <- res$Flag[which(res$log2FoldChange < 0)] * -1
      #names(res) <- unlist(lapply(names(res), function(y) paste(y,"_",x[1],"_vs_",x[2],sep="")))
      res <- res[order(rownames(res)), ]
      return(res)
    })
    # combine results from all pairwise comparisons
    res2 <- do.call(cbind, res)
  }else if(test=="paired"){
    paired <- paste("~",pairs[2],"+",pairs[1],sep="")
    design(DESeq_data_set_object) <- formula(paired)
    DESeq_data_set_object <- DESeq(DESeq_data_set_object, quiet=TRUE)
    # extract results for paired condition
    res <- results(DESeq_data_set_object, pAdjustMethod=adj)
    res$Flag <- ifelse(res$padj <= thresh, 1, 0)
    res$Flag[which(res$log2FoldChange < 0)] <- res$Flag[which(res$log2FoldChange < 0)] * -1
    names(res) <- unlist(lapply(names(res), function(y) paste(y,"_",levels(attributes(omicsData)$group_DF[,pairs[[1]][1]])[2],"_vs_",levels(attributes(omicsData)$group_DF[,pairs[[1]][1]])[1],sep="")))
    res2 <- res[order(rownames(res)), ]
#     res <- apply(pairs, 2, function(x){
#       myContrast <- c("Group", x[1], x[2])
#       res <- results(analysis_data_set_object, contrast=myContrast)
#       names(res) <- unlist(lapply(names(res), function(y) paste(y,"_",x[1],"_",x[2],sep="")))
#       res <- res[order(rownames(res)), ]
#       return(res)
#     })
  }else{
    stop("test must be one of 'wald', 'lrt', or 'paired")
  }
  # Wald Test - uses the estimated standard error of a log2 fold change to test if it is equal to zero, this is the default in DESeq2
  # LRT - examines both a full model and a reduced model (where some terms in full model are removed and determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero - more useful for testing multiple terms at once, e.g. 3 or more levels of a factor at once or all interactiosn between two variables. The LRT for count data is conceptually similar to an analysis of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a reduced model.
  # Example - DESeq(dds, test="LRT", reduced=~1)
  # Example - DESeq(dds, test="LRT", reduced=~batch)

#   DESeq_data_set_object <- estimateDispersions(DESeq_data_set_object, quiet=TRUE)
#   DESeq_data_set_object <- nbinomWaldTest(DESeq_data_set_object, quiet=TRUE)

  return(res2)
}
