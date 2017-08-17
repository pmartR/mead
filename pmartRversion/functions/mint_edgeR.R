#' edgeR analysis of mintR objects
#'
#' Differential abundance analysis of count data using edgeR
#'
#' @param omicsData an object of the class 'gDNAdata', 'cDNAdata', or 'rRNAdata' usually created by \code{\link{as.gDNAdata}}, \code{\link{as.cDNAdata}},or \code{\link{as.rRNAdata}}, respectively.
#' @param norm_factors Named vector of normalization parameters to put into DESeq2. If NULL, will use DESeq2's inbuilt normalization. Default is NULL.
#' @param pairs a matrix dictating which pairwise comparisons to make
#' @param test name of which differential expression test to use, options are "qcml", "lrt", "qlftest", or "paired". Default is "qcml". See details for further explanation.
#' @param adj multiple comparison adjustment method to use.
#' @param thresh p-value threshold for significance.
#'
#' @details Perform pairwise differential abundance analysis of count data using edgeR. For test parameter: "qcml" refers to quantile-adjusted conditional maximum likelihood method, which is the most reliable in terms of bias on a wide range of conditions and specifically performs best in the situation of many small samples with a common dispersion. "lrt" refers to generalized linear model likelihood ratio test, best used for cases where there a multiple treatment groups and provides inferences with GLMs. "qlftest" refers to the QL F-test, preferred as it reflects the uncertainty in estimating the dispersion for each gene and provides a more robust and reliable error rate control when the number of replicates is small. "paired" refers to a paired test and performs the generalized liner model likelihood ratio test for paired data (this should also be used for experiments with a block or batch effect). In the "paired" test, 'group1' should be the main effect of interest and 'group2' should be the covariate/block/batch. Currently, this can be used only if there are 2 main effects of interest.
#'
#' @return A matrix containing the log2 fold change (logFC), log-average concentration/abundance (logCPM), likelihood ratio (LR), exact p-value for differential expression using the negative binomial model (PValue), and the p-value adjusted for multiple testing (FDR) for every pairwise comparison and every feature.
#'
#' @examples
#' library(mintJansson)
#' data(cDNA_hiseq_data)
#' mycdnadata <- group_designation(omicsData = cDNA_hiseq_data, main_effects = c("treatment"), time_course=NULL)
#' mycdnadata_norm <- normalize_data(omicsData = mycdnadata, norm_fn = "percentile")
#' mycdnadata_edgeR <- mint_edgeR(omicsData = mycdnadata_norm, test="qcml", pairs = cbind(list("Neg","Plus")), adj = "BH", thresh = 0.05)
#'
#' @author Allison Thompson
#'
#' @references Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26, pp. -1.
#' McCarthy, J. D, Chen, Yunshun, Smyth and K. G (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research, 40(10), pp. -9.
#' Robinson MD and Smyth GK (2007). Moderated statistical tests for assessing differences in tag abundance. Bioinformatics, 23, pp. -6.
#' Robinson MD and Smyth GK (2008). Small-sample estimation of negative binomial dispersion, with applications to SAGE data. Biostatistics, 9, pp. -11.
#' Zhou X, Lindsay H and Robinson MD (2014). Robustly detecting differential expression in RNA sequencing data using observation weights. Nucleic Acids Research, 42, pp. e91.
#'

mint_edgeR <- function(omicsData, norm_factors=NULL, test="qcml", pairs, adj, thresh){
  library(edgeR)
  library(dplyr)

  if(!(class(omicsData) %in% c("gDNAdata","cDNAdata","rRNAdata"))){
    stop("Data must be of class gDNAdata, cDNAdata, or rRNAdata")
  }

  if(!(tolower(test) %in% c("qcml","lrt","qlftest","paired"))){
    stop("test must be qcml, lrt, qlftest, or paired")
  }

  # Make features be rownames and remove that column from the data, so the data is all numeric
  rownames(omicsData$e_data) <- omicsData$e_data[,attr(omicsData,"cnames")$edata_cname]
  omicsData$e_data <- omicsData$e_data[,-which(colnames(omicsData$e_data) == attr(omicsData,"cnames")$edata_cname)]

  # Change NA's to 0 for edgeR
  omicsData$e_data[is.na(omicsData$e_data)] <- 0

  # Check if group1 and/or group2 is null
  if(is.null(attributes(omicsData)$group_DF)){
    stop("run group designation function first")
  }

  # format metadata
  metadata <- attributes(omicsData)$group_DF
  rownames(metadata) <- metadata[,attr(omicsData,"cnames")$fdata_cname]
  metadata <- metadata[which(rownames(metadata)%in%colnames(omicsData$e_data)),]
  metadata <- metadata[match(colnames(omicsData$e_data), metadata[,attr(omicsData,"cnames")$fdata_cname]),]

  # Creates object to be used in differential expression tests
  y <- DGEList(counts = omicsData$e_data, group = metadata$Group)

  # Set own sizeFactors(normalization factors), after running normalization
  if(is.null(norm_factors)){
    warning("using edgeR's inbuilt TMM normalization - if want to use own normalization, run normalization function first")
    y <- calcNormFactors(y)
  }else{
    # check that normalization (size) factors are present for every sample
    if(all(colnames(omicsData$e_data) %in% names(norm_factors))){
      norm_factors <- norm_factors[match(colnames(omicsData$e_data), names(norm_factors))]
      y$samples$norm.factors <- norm_factors
    }else{
      stop("Must have a normalization (size) factor for every sample in the data")
    }
  }

  if(test == "qcml"){
    # quantile adjusted conditional maximum likelihood test
    group <- metadata$Group
    design <- model.matrix(~0+group, data=y$samples)

    y <- estimateDisp(y, design)

    # get results from all pairwise comparisons
    results <- apply(pairs, 2, function(x){
      x <- as.character(x)
      eres <- topTags(exactTest(y, pair=c(x[2],x[1])), n=nrow(omicsData$e_data), adjust.method=adj)
    })

    if(adj == "none"){
      results <- lapply(results, function(x){
        x@.Data[[1]]$FDR=x@.Data[[1]]$PValue
        return(x)
      })
    }

    # format results for all comparisons
    res <- lapply(results, function(x){
      res <- as.data.frame(x@.Data[[1]])
      res$Flag <- ifelse(res$FDR <= thresh, 1, 0)
      res$Flag[which(res$logFC < 0)] <- res$Flag[which(res$logFC < 0)] * -1
      names(res) <- unlist(lapply(names(res), function(y) paste(y,"_",x@.Data[[3]][2],"_vs_",x@.Data[[3]][1],sep="")))
      res <- res[order(rownames(res)),]
      return(res)
    })

    res2 <- do.call(cbind, res)

  }else if(test == "lrt"){
    # generalized linear model likelihood ratio test
    group <- metadata$Group
    design <- model.matrix(~0+group, data=y$samples)

    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)

    # get results from all comparisons
    results <- apply(pairs, 2, function(x){
      x <- as.character(x)
      cntrst <- rep(0, length(levels(metadata$Group)))
      cntrst[which(levels(metadata$Group) == x[1])] <- 1
      cntrst[which(levels(metadata$Group) == x[2])] <- -1
      lrt <- glmLRT(fit, contrast=cntrst)
      results <- topTags(lrt, n=nrow(omicsData$e_data), adjust.method=adj)
      return(results)
    })

    if(adj == "none"){
      results <- lapply(results, function(x){
        x@.Data[[1]]$FDR=x@.Data[[1]]$PValue
        return(x)
      })
    }

    # format results for all comparisons
    res <- lapply(results, function(x){
      res <- as.data.frame(x@.Data[[1]])
      res$Flag <- ifelse(res$FDR <= thresh, 1, 0)
      res$Flag[which(res$logFC < 0)] <- res$Flag[which(res$logFC < 0)] * -1
      grps <- unlist(strsplit(x@.Data[[3]], " "))
      grp2 <- strsplit(grps[grep("-1",grps)],"group")[[1]][2]
      grp1 <- strsplit(grps[-grep("-1",grps)],"group")[[1]][2]
      names(res) <- unlist(lapply(names(res), function(y) paste(y,"_",grp1,"_vs_",grp2,sep="")))
      res <- res[order(rownames(res)),]
      return(res)
    })

    res2 <- do.call(cbind, res)

  }else if(test == "qlftest"){
    # QL F-test
    group <- metadata$Group
    design <- model.matrix(~0+group, data=y$samples)

    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)

    # get results from all comparisons
    results <- apply(pairs, 2, function(x){
      x <- as.character(x)
      cntrst <- rep(0, length(levels(metadata$Group)))
      cntrst[which(levels(metadata$Group) == x[1])] <- 1
      cntrst[which(levels(metadata$Group) == x[2])] <- -1
      qlf <- glmQLFTest(fit, contrast=cntrst)
      results <- topTags(qlf, n=nrow(omicsData$e_data), adjust.method=adj)
      return(results)
    })

    if(adj == "none"){
      results <- lapply(results, function(x){
        x@.Data[[1]]$FDR=x@.Data[[1]]$PValue
        return(x)
      })
    }

    # format results for all comparisons
    res <- lapply(results, function(x){
      res <- as.data.frame(x@.Data[[1]])
      res$Flag <- ifelse(res$FDR <= thresh, 1, 0)
      res$Flag[which(res$logFC < 0)] <- res$Flag[which(res$logFC < 0)] * -1
      grps <- unlist(strsplit(x@.Data[[3]], " "))
      grp2 <- strsplit(grps[grep("-1",grps)],"group")[[1]][2]
      grp1 <- strsplit(grps[-grep("-1",grps)],"group")[[1]][2]
      names(res) <- unlist(lapply(names(res), function(y) paste(y,"_",grp1,"_vs_",grp2,sep="")))
      res <- res[order(rownames(res)),]
      return(res)
    })

    res2 <- do.call(cbind, res)

  }else if(test == "paired"){
    # paired generalized linear model likelihood ratio test
    grp1 <- factor(attributes(omicsData)$group_DF[,pairs[1]])
    grp2 <- factor(attributes(omicsData)$group_DF[,pairs[2]])
    design <- model.matrix(~grp2+grp1)

    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit)

    results <- topTags(lrt, n=nrow(omicsData$e_data), adjust.method=adj)

    if(adj == "none"){
      results@.Data[[1]]$FDR = results@.Data[[1]]$PValue
#       results <- lapply(results, function(x){
#         x@.Data[[1]]$FDR=x@.Data[[1]]$PValue
#         return(x)
#       })
    }

    res <- as.data.frame(results@.Data[[1]])
    res$Flag <- ifelse(res$FDR <= thresh, 1, 0)
    res$Flag[which(res$logFC < 0)] <- res$Flag[which(res$logFC < 0)] * -1
    names(res) <- unlist(lapply(names(res), function(x) paste(x,"_",levels(grp1)[2],"_vs_",levels(grp1)[1],sep="")))
    res2 <- res[order(rownames(res)),]

  }else{
    stop("test must be one of 'qcml', 'lrt', 'qlftest', or 'paired")
  }

  return(res2)
}
