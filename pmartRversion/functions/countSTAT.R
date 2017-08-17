# Wrapper for count data statistical tests

# Inputs - mintRdata, comparisons, test (DESeq2/Wald, DESeq2/LRT, DESeq2/Paired, edgeR/QCML, edgeR/LRT, edgeR/ANOVA, edgeR/Paired), control (if comparing against control), p-value threshold, pvalue adjustment method,

#' Statistical tests for count data
#'
#' Differential expression tests for count data
#'
#' @param omicsData an object of the class 'gDNAdata', 'cDNAdata', or 'rRNAdata' usually created by \code{\link{as.gDNAdata}}, \code{\link{as.cDNAdata}},or \code{\link{as.rRNAdata}}, respectively.
#' @param norm_factors Named vector of normalization parameters to put into DESeq2 or edgeR. If NULL, will use DESeq2's or edgeR's inbuilt normalization. Default is NULL.
#' @param comparisons dictates which pairwise comparisons to make. 'all' will give the results for all pairwise comparisons, 'control' will give the results for all comparisons against a specified control group, or a list of specific comparisons to be made or terms for a paired test can be given.
#' @param control only necessary when performing comparisons against control, name of the control group. default is NULL.
#' @param test names of which differential expression test(s) to use, options are "dw" (DESeq2 analysis with Wald test), "dl" (DESeq2 analysis with LRT test), "dp" (DESeq2 analysis on paired data), "eq" (edgeR analysis with QCML test), "el" (edgeR analysis with LRT test), "ef" (edgeR analysis with QL F-test), and/or "ep" (edgeR analysis on paired data). At this point, cannot perform a paired analysis at the same time as the other tests.
#' @param pval_adjust Name of which multiple comparisons adjustment to use. Options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "FDR", and "none". See ?p.adjust for more information. The default is "BH".
#' @param pval_thresh P-value threshold for significance. The default is 0.05.
#'
#' @details Perform pairwise differential abundance analysis of mintR count data using DESeq2 and/or edgeR.
#'
#' @return A results object containing the log2 fold change (logFC), log-average concentration/abundance (logCPM), likelihood ratio (LR), exact p-value for differential expression using the negative binomial model (PValue), and the p-value adjusted for multiple testing (FDR) for every pairwise comparison and every feature.
#'
#' @examples
#' library(mintJansson)
#' data(cDNA_hiseq_data)
#' mycdnadata <- group_designation(omicsData = cDNA_hiseq_data, main_effects = c("treatment"), time_course=NULL)
#' mycdnadata_norm <- normalize_data(omicsData = mycdnadata, norm_fn = "percentile")
#' mycdnadata_dw_results <- countSTAT(omicsData = mycdnadata_norm, comparisons = list(c("Neg","Plus")), control = NULL, test = "dw", pval_adjust = "none", pval_thresh = 0.05)
#' mycdnadata_all_results <- countSTAT(omicsData = mycdnadata_norm, comparisons = "all", control = NULL, test = c("dw","eq","el","ef"), pval_adjust = "BH", pval_thresh = 0.05)
#'
#' @author Allison Thompson
#'
#' @export
countSTAT <- function(omicsData, norm_factors=NULL, comparisons, control = NULL, test, pval_adjust = "BH", pval_thresh = 0.05 ){

  library(reshape2)

  if(!(class(omicsData) %in% c("gDNAdata","cDNAdata","rRNAdata"))){
    stop("Data must be of class gDNAdata, cDNAdata, or rRNAdata")
  }

  if(!(any((tolower(test)) %in% c("dw","dl","dp","eq","el","ef","ep")))){
    stop("test must only contain dw, dl, dp, eq, el, ef, ep")
  }

  # determine which pairwise comparisons to make
  if(tolower(comparisons) == "all"){
    pairs <- combn(levels(attributes(omicsData)$group_DF$Group),2)
  }else if(tolower(comparisons) == "control"){
    pairs <- combn(levels(attributes(omicsData)$group_DF$Group),2)
    if(ncol(pairs) > 1){
      if(length(unique(c(grep(control,pairs[1,]),grep(control,pairs[2,])))) > 1){
        pairs <- pairs[,unique(c(grep(control, pairs[1,]),grep(control, pairs[2,])))]
      }
      if(any(pairs[1,]==control)){
        pairs[,which(pairs[1,] == control)] <- apply(pairs[,which(pairs[1,] == control)], 2, rev)
      }
    }
  }else if(is.list(comparisons)){
    pairs <- do.call(cbind, comparisons)
  }else{
    stop("check that comparisons argument is either 'all', 'control', or a list of specific comparisons or terms for a paired test")
  }

  res <- list()
  test_text <- "Differential expression analysis of the data was performed using the following method(s)"

  if("dw" %in% tolower(test)){
    # Run DESeq2 function, using the Wald test
    res[["dw"]] <- mint_DESeq2(omicsData=omicsData, norm_factors=norm_factors, test="wald", pairs=pairs, adj=pval_adjust, thresh=pval_thresh)
    names(res[["dw"]]) <- gsub("log2FoldChange","logFC",names(res[["dw"]]))
    names(res[["dw"]]) <- gsub("pvalue","PValue",names(res[["dw"]]))
    #names(res[["dw"]]) <- unlist(lapply(names(res[["dw"]]), function(x) paste("dw_", x, sep="")))
    res[["dw"]] <- as.data.frame(res[["dw"]])
    test_text <- c(test_text,", DESeq2 using the Wald test")
  }
#   if("dl" %in% tolower(test)){
#     # Run DESeq2 function, using the LRT test
#     # DOESN'T WORK, PROBLEM WITH REDUCED MODEL
#     res[["dl"]] <- mint_DESeq2(omicsData=omicsData, test="lrt", pairs=pairs, adj=pval_adjust, thresh=pval_thresh)
#     names(res[["dl"]]) <- gsub("log2FoldChange","logFC",names(res[["dl"]]))
#     #names(res[["dl"]]) <- unlist(lapply(names(res[["dl"]]), function(x) paste("dl_", x, sep="")))
#     res[["dl"]] <- as.data.frame(res[["dl"]])
#     test_text <- c(test_text,", DESeq2 using the likelihood ratio test")
#   }
  if("dp" %in% tolower(test)){
    # Run DESeq2 function, performing paired test
    res[["dp"]] <- mint_DESeq2(omicsData=omicsData, norm_factors=norm_factors, test="paired", pairs=pairs, adj=pval_adjust, thresh=pval_thresh)
    names(res[["dp"]]) <- gsub("log2FoldChange","logFC",names(res[["dp"]]))
    names(res[["dp"]]) <- gsub("pvalue","PValue",names(res[["dp"]]))
    #names(res[["dp"]]) <- unlist(lapply(names(res[["dp"]]), function(x) paste("dp_", x, sep="")))
    res[["dp"]] <- as.data.frame(res[["dp"]])
    test_text <- c(test_text,", DESeq2 using the paired test")
  }
  if("eq" %in% tolower(test)){
    # Run edgeR function, using QCML test
    res[["eq"]] <- mint_edgeR(omicsData=omicsData, norm_factors=norm_factors, test="qcml", pairs=pairs, adj=pval_adjust, thresh=pval_thresh)
    names(res[["eq"]]) <- gsub("FDR","padj",names(res[["eq"]]))
    #names(res[["eq"]]) <- unlist(lapply(names(res[["eq"]]), function(x) paste("eq_", x, sep="")))
    test_text <- c(test_text,", edgeR using the quantile-adjusted conditional maximum likelihood test")
  }
  if("el" %in% tolower(test)){
    # Run edgeR function, using LRT test
    res[["el"]] <- mint_edgeR(omicsData=omicsData, norm_factors=norm_factors, test="lrt", pairs=pairs, adj=pval_adjust, thresh=pval_thresh)
    names(res[["el"]]) <- gsub("FDR","padj",names(res[["el"]]))
    #names(res[["el"]]) <- unlist(lapply(names(res[["el"]]), function(x) paste("el_", x, sep="")))
    test_text <- c(test_text,", edgeR using the likelihood ratio test")
  }
  if("ef" %in% tolower(test)){
    # Run edgeR function, using ANOVA test
    res[["ef"]] <- mint_edgeR(omicsData=omicsData, norm_factors=norm_factors, test="qlftest", pairs=pairs, adj=pval_adjust, thresh=pval_thresh)
    names(res[["ef"]]) <- gsub("FDR","padj",names(res[["ef"]]))
    #names(res[["ef"]]) <- unlist(lapply(names(res[["ef"]]), function(x) paste("ef_", x, sep="")))
    test_text <- c(test_text,", edgeR using the quasi-likelihood F-test")
  }
  if("ep" %in% tolower(test)){
    # Run edgeR function, performing paired test
    res[["ep"]] <- mint_edgeR(omicsData=omicsData, norm_factors=norm_factors, test="paired", pairs=pairs, adj=pval_adjust, thresh=pval_thresh)
    names(res[["ep"]]) <- gsub("FDR","padj",names(res[["ep"]]))
    #names(res[["ep"]]) <- unlist(lapply(names(res[["ep"]]), function(x) paste("ep_", x, sep="")))
    test_text <- c(test_text,", edgeR using the paired test")
  }

  # Check that differential expression results are returned for all biomolecules
  temp.rn <- lapply(res, function(x) nrow(x))
  if(!(all(temp.rn == nrow(omicsData$e_data)))){
    stop("Something weird happened, why aren't we returning all results?")
  }

  # Put all results in the same order in order to bind them together
  if(length(res) > 1){
    allRes <- lapply(res, function(x) x[match(rownames(res[[1]]),rownames(x)),])
    allRes <- do.call(cbind, allRes)
  }else{
    allRes <- as.data.frame(res)
  }

  allRes <- allRes[,-grep("PValue", colnames(allRes))]

  if(pval_adjust == "none"){
    test_text <- c(test_text," and without any adjustment.")
  }else{
    test_text <- c(test_text," and with ", pval_adjust, " adjustment.")
  }

  test_text <- paste(test_text, collapse="")

  results <- list()
  results$allResults <- allRes

  if("dp" %in% tolower(test) | "ep" %in% tolower(test)){
    compars <- paste(levels(attributes(omicsData)$group_DF[,pairs[[1]][1]])[2],"_vs_",levels(attributes(omicsData)$group_DF[,pairs[[1]][1]])[1],sep="")
    attr(results, "comparisons") = list(comparison = compars, terms = paste("~",pairs[2],"+",pairs[1],sep=""))
    #attr(results, "comparisons")$terms =  paste("~",pairs[2],"+",pairs[1],sep="")
  }else{
    compars <- apply(pairs, 2, function(x) paste(x[1],"_vs_",x[2],sep=""))
    attr(results, "comparisons") = list(comparison = compars, term = "Group")
    #attr(results, "comparisons")$terms = "Group"
  }

  attr(results,"Adjustment")=pval_adjust
  attr(results,"Threshold")=pval_thresh
  attr(results,"Tests")$Test=test
  attr(results,"Tests")$Text=test_text

  attr(results,"group_DF")=attr(omicsData,"group_DF")

  class(results) <- c("countSTAT_results", class(omicsData))

  #results <- list(allResults=allRes, Comparisons=pairs, Adjustment.Method=pval_adjust, Threshold=pval_thresh, Tests=test, Text=test_text)

  return(results)

}
