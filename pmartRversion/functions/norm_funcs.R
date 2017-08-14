#' Cumulative sum scaling normalization of count data
#'
#' The method normalizes count data by cumulative sum scaling, using a specified quantile (e.g., 75th quantile)
#'
#' @param e_data a \eqn{p \times n} data.frame of count data, where \eqn{p} is the number of features and \eqn{n} is the number of samples. Each row corresponds to data for a feature, with the first column giving the feature name.
#' @param edata_id character string indicating the name of the feature identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param q a number to indicate which quantile to normalize with. Default is 0.75.
#' @param qg a number to use for scaling all samples. Default is 1000 (as in reference). Can also specify "median" to use the median value scaling value across all samples.
#'
#' @details Count data is normalized by a given quantile, dividing by the sum of all values up to and including the given quantile of the sample and multiplying by the given scaling value (either 1000 or the median scaling value across all samples).
#'
#' @return List containing 3 elements: norm_data is a data.frame with same structure as e_data that contains the quantile-normalized data, location_param is NULL, and scale_param is a numeric vector containing, for every sample, the value of the sample at the designated quantile (q) divided by the value of the global quantile (q).
#'
#' @references Paulson, Joseph N, O Colin Stine, Hector Corrada Bravo, and Mihai Pop. "Differential abundance analysis for microbial marker-gene surveys." Nature Methods. 10.12 (2013)
#'
#' @examples
#' library(mintJansson)
#' data(rRNA_data)
#' rRNA_CSS <- CSS_Norm(e_data = rRNA_data$e_data, edata_id = attr(rRNA_data, "cnames")$edata_cname)
#' norm_factors <- attr(rRNA_CSS,"data_info")$scale_param
#'
#' @author Allison Thompson, Lisa Bramer
#'

CSS_Norm <- function(e_data, edata_id, q=0.75, qg="median"){
  e_data_norm <- e_data[,-which(colnames(e_data)==edata_id)]

  # calculate the q quantile of data, per sample and globally #
  col.q = apply(e_data_norm, 2, function(x) sum(x[x<=quantile(x[x!=0], probs = q, na.rm=TRUE)], na.rm=TRUE))
  #g.q = sum(e_data_norm[e_data_norm <= quantile(e_data_norm, probs = q, na.rm=TRUE)], na.rm=TRUE)

  if(qg == 1000){
    g.q = 1000
  }else if(qg == "median"){
    g.q = median(col.q, na.rm=TRUE)
  }else{
    stop("Invalid value for qg")
  }
  #g.q = 1000
  #g.q = median(col.q, na.rm=TRUE)

  # normalize omics_data data by q quantile and transform back to count data #
  for(i in 1:ncol(e_data_norm)){
    e_data_norm[,i] = (e_data_norm[,i] / col.q[i]) * g.q
  }

  e_data <- data.frame(e_data[,which(colnames(e_data)==edata_id)],e_data_norm)
  colnames(e_data)[1] <- edata_id
  return(list(normed_data=e_data, location_param=NULL, scale_param=col.q/g.q))
}


#' Normalization of count data via DESeq scale factors
#'
#' The method normalizes count data by scaling the raw counts by the median scaled counts.  The normalized factors used by
#' this procedure are also called DEseq size scores (Anders et al. 2010).  The same normalization technique is used in DESeq2 as well (Love et al. 2014).
#'
#' @param e_data a \eqn{p \times n} data.frame of count data, where \eqn{p} is the number of features and \eqn{n} is the number of samples. Each row corresponds to data for a feature, with the first column giving the feature name.
#' @param edata_id character string indicating the name of the peptide, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#'
#' @details Count data is normalized by the median scaled score
#'
#' @return List containing 3 elements: norm_data is a data.frame with same structure as e_data that contains the normalized data, location_param is NULL, scale_param is a numeric vector of DESeq scores.
#'
#' @author Bryan Stanfill
#'
#' @references
#' Anders, Simon, and Wolfgang Huber. "Differential expression analysis for sequence count data." Genome biol 11.10 (2010): R106.
#' Love, Michael I., Wolfgang Huber, and Simon Anders. "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome biology 15.12 (2014): 1-21.
#' @examples
#' library(mintJansson)
#' e_data_id <- attr(rRNA_data, "cnames")$edata_cname
#' deseq_n <- med_scounts_norm(e_data = rRNA_data$e_data, e_data_id)
#'

med_scounts_norm <- function(e_data, edata_id){

  #Which column has the ID variable we want to remove?
  col_to_omit <- which(colnames(e_data)==edata_id)
  if(length(col_to_omit)==0){
    stop("The supplied edata_id name couldn't be found in e_data.")
  }

  log_scale_mean <- rowMeans(log(e_data[,-col_to_omit]),na.rm=TRUE)

  #Compute size factors: sjs
  sjs <- as.numeric(apply(e_data[,-col_to_omit], 2, function(cnts) median(exp(log(cnts) - log_scale_mean),na.rm=T)))

  e_data[,-col_to_omit] <- data.matrix(e_data[,-col_to_omit])%*%diag(1/sjs)

  #if(nas){
  #  edata_nona[which_na] <- NA
  #}
  return(list(normed_data=e_data, location_param=NULL, scale_param=sjs))
}


#' Normalization of count data via Poisson Sampling
#'
#' The method normalizes count data by resampling the data from a Poisson distribution with parameter estimated from the raw counts.  More specifically,
#' the quantity that is returned is given by equation 2.4 of Li and Tibshirani (2013).
#'
#' @param e_data a \eqn{p \times n} data.frame of count data, where \eqn{p} is the number of features and \eqn{n} is the number of samples. Each row corresponds to data for a feature, with the first column giving the feature name.
#' @param edata_id character string indicating the name of the peptide, lipid, or metabolite identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#'
#' @details Count data resampled from the Poisson distribution
#'
#' @return List containing 3 elements: norm_data is a data.frame with same structure as e_data that contains the normalized data, location_param is NULL, scale_param is a numeric vector of DESeq scores.
#'
#' @author Bryan Stanfill
#'
#' @references
#'  Li, Jun, and Robert Tibshirani. "Finding consistent patterns: a nonparametric approach for identifying differential expression in RNA-Seq data." Statistical methods in medical research 22.5 (2013): 519-536.
#'
#' @examples
#' library(mintJansson)
#' e_data <- rRNA_data$e_data
#' e_data_id <- attr(rRNA_data, "cnames")$edata_cname
#' pois_ndata <- poisson_norm(rRNA_data$e_data, e_data_id)
#'

poisson_norm <- function(e_data, edata_id){

  #Compute sequence depts for each experiment
  ds <- colSums(e_data[,-1],na.rm=T)

  #Compute geometric mean of counts by computing the mean on the log scale and exponentiating
  dbar <- exp(mean(log(ds)))

  #Which column has the ID variable we want to remove?
  col_to_omit <- which(colnames(e_data)==edata_id)
  if(length(col_to_omit)==0){
    stop("The supplied edata_id name couldn't be found in e_data.")
  }

  #For each non-nNA cell, sample from a poisson dist with parameter (dbar/ds[i])N_{ij}
  e_data[,-col_to_omit] <- apply(e_data[,-col_to_omit],2,pois_samp,dbar=dbar)

  return(list(normed_data=e_data, location_param=NULL, scale_param=dbar/ds))

}


pois_samp <- function(Nij,dbar){
  rmna <- which(is.na(Nij))
  di <- sum(Nij,na.rm=T)
  Nij[-rmna] <- rpois(length(Nij)-length(rmna),Nij[-rmna]*dbar/di)
  return(Nij)
}


#' Quantile normalization of count data
#'
#' The method normalizes count data by a specified quantile (e.g., 75th quantile)
#'
#' @param e_data a \eqn{p \times n} data.frame of count data, where \eqn{p} is the number of features and \eqn{n} is the number of samples. Each row corresponds to data for a feature, with the first column giving the feature name.
#' @param edata_id character string indicating the name of the feature identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param q a number to indicate which quantile to normalize with. Default is 0.75.
#'
#' @details Count data is normalized by a given quantile, dividing by the given quantile of the sample and multiplying by the given quantile of the entire dataset.
#'
#' @return List containing 3 elements: norm_data is a data.frame with same structure as e_data that contains the quantile-normalized data, location_param is NULL, and scale_param is a numeric vector containing, for every sample, the value of the sample at the designated quantile (q) divided by the value of the global quantile (q).
#'
#' @examples
#' library(mintJansson)
#' data(cDNA_hiseq_data)
#' cDNA_quant <- Quant_Norm(e_data = cDNA_hiseq_data$e_data, edata_id = attr(cDNA_hiseq_data, "cnames")$edata_cname)
#' norm_factors <- attr(cDNA_quant,"data_info")$scale_param
#'
#' @author Allison Thompson, Lisa Bramer
#'

Quant_Norm <- function(e_data, edata_id, q=0.75){
  e_data_norm <- e_data[,-which(colnames(e_data)==edata_id)]

  # calculate the q quantile of data, per sample and globally #
  col.q = apply(e_data_norm, 2, function(x) quantile(x, probs = q, na.rm=TRUE))
  g.q = quantile(e_data_norm, probs = q, na.rm=TRUE)

  # normalize omics_data data by q quantile and transform back to count data #
  for(i in 1:ncol(e_data_norm)){
    e_data_norm[,i] = (e_data_norm[,i] / col.q[i]) * g.q
  }

  e_data <- data.frame(e_data[,which(colnames(e_data)==edata_id)],e_data_norm)
  colnames(e_data)[1] <- edata_id
  return(list(normed_data=e_data, location_param=NULL, scale_param=col.q/g.q))
}


#' Rarefying of count data
#'
#' The method normalizes count data by rarefying
#'
#' @param e_data a \eqn{p \times n} data.frame of count data, where \eqn{p} is the number of features and \eqn{n} is the number of samples. Each row corresponds to data for a feature, with the first column giving the feature name.
#' @param edata_id character string indicating the name of the feature identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param size the library size to rarefy down to. Default uses the minimum sample size.
#'
#' @details Count data is normalized by rarefying, subsampling samples down to a specified library size. If the specified library size is larger than a sample's library size, the sample will be discarded. A warning message will display which samples are discarded. This normalization method is likely not the best course of action.
#'
#' @return List containing 4 elements: norm_data is a data.frame with same structure as e_data that contains the rarefied data, location_param is NULL, and scale_param is the library size counts were rarefied to.
#'
#' @examples
#' library(mintJansson)
#' data(cDNA_hiseq_data)
#' cDNA_rarefy <- Rarefy(e_data = cDNA_hiseq_data$e_data, edata_id = attr(cDNA_hiseq_data, "cnames")$edata_cname)
#' library_size <- attr(cDNA_rarefy,"data_info")$scale_param
#'
#' @author Allison Thompson
#'
#' @references McMurdie, Paul J., and Susan Holmes. "Waste Not, Want Not: Why Rarefying Microbiome Data Is Inadmissible." PLOS Computational Biology. 10.4 (2014)
#'

Rarefy <- function(e_data, edata_id, size=NULL){
  e_data_norm <- e_data[,-which(colnames(e_data)==edata_id)]
  # check library size
  if(is.null(size)){
    size <- min(colSums(e_data_norm,na.rm=TRUE))
  }else{
    size <- size
  }

  # if input library size is larger than any library size for any sample, discard that sample, check group designations first  #
  samples <- NA
  if(any(colSums(e_data_norm,na.rm=TRUE) < size)){
    samples <- names(which(colSums(e_data_norm,na.rm=TRUE) < size))
    # write a warning stating which samples are being removed #
    warning(paste("Input size is larger than the library size of samples ", paste(samples,collapse=", "), ". These samples will be removed from the data.",sep=""))
  }

  # Change NA to 0 in order to extract features
  e_data_norm[is.na(e_data_norm)]<- 0

  # Create a temporary empty matrix to fill with rarefied data
  temp <- matrix(nrow=nrow(e_data_norm),ncol=ncol(e_data_norm))
  rownames(temp) <- e_data[,edata_id]
  colnames(temp) <- colnames(e_data_norm)

  for(i in 1:ncol(e_data_norm)){
    # Create a vector of features
    features <- rep(e_data[,edata_id],e_data_norm[,i])

    # Subsample features to a predetermined library size
    subsamp <- sample(features, ifelse(length(features)>=size,size,0),replace=FALSE)
    data <- table(subsamp)

    # Put rarefied counts back into temporary matrix
    for(j in 1:length(data)){
      temp[which(rownames(temp)==names(data)[j]),i] <- data[j]
    }
    colnames(temp)[i] <- colnames(e_data_norm)[i]
  }

  if(any(rowSums(temp, na.rm=TRUE) == 0)){
    # Remove any extra features
    temp <- temp[-which(rowSums(temp,na.rm=TRUE)==0),]
  }

  # Set 0 back to NA
  temp[temp==0] <- NA

  # Reorganize so this resembles e_data
  temp <- as.data.frame(temp)
  temp <- data.frame(rownames(temp),temp)
  names(temp)[1] <- edata_id
  rownames(temp) <- NULL

  # Replace e_data with rarefied counts
  e_data <- temp

  return(list(normed_data=e_data, location_param=NULL, scale_param=size))
}


#' Trimmed Mean of M Values normalization of count data
#'
#' The method normalizes count data by the trimmed mean of m values in each sample
#'
#' @param e_data a \eqn{p \times n} data.frame of count data, where \eqn{p} is the number of features and \eqn{n} is the number of samples. Each row corresponds to data for a feature, with the first column giving the feature name.
#' @param edata_id character string indicating the name of the feature identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#' @param reference which column in e_data should be used as the reference, default is to use the sample with the least amount of missing data.
#' @param qm percentage by which to trim M values (gene-wise log-fold-changes), default is 0.30 (30\%)
#' @param qa percentage by which to trim A values (absolute expression levels), default is 0.05 (5\%)
#'
#' @details Count data is normalized by the trimmed mean of m values.
#'
#' @return List containing 3 elements: norm_data is a data.frame with same structure as e_data that contains the TMM-normalized data, location_param is a numeric vector of the TMM values for each sample, and scale_param is NULL.
#'
#' @examples
#' library(mintJansson)
#' data(cDNA_hiseq_data)
#' cDNA_TMM <- Quant_Norm(e_data = cDNA_hiseq_data$e_data, edata_id = attr(cDNA_hiseq_data, "cnames")$edata_cname)
#' norm_factors <- attr(cDNA_TMM,"data_info")$scale_param
#'
#' @author Allison Thompson, Lisa Bramer
#'

TMM_Norm <- function(e_data,edata_id, reference=NULL,qm=0.30,qa=0.05){
  e_data[which(e_data==0)] <- NA
  e_data_norm <- e_data[,-which(colnames(e_data)==edata_id)]
  # determine which sample should be used as the reference
  if(is.null(reference)){
    # determine which sample has the least amount of missing data
    ref = which(apply(e_data_norm, 2, function(x) length(which(is.na(x)))) == min(apply(e_data_norm, 2,function(x) length(which(is.na(x))))))
    if(length(ref) > 1){
      ref = which(apply(e_data_norm[,ref],2,function(x) mean(x,na.rm=TRUE)) == max(apply(e_data_norm[,ref],2,function(x) mean(x,na.rm=TRUE))))
    }
    if(length(ref) > 1){
      ref = which(colSums(e_data_norm[,ref], na.rm=TRUE) == max(colSums(e_data_norm[,ref], na.rm=TRUE)))
    }
  }else{
    ref = reference
  }

  # remove genes where the reference has a value of NA
  data <- e_data_norm[-which(is.na(e_data_norm[,ref])),]
  names <- e_data[,edata_id]
  names <- names[-which(is.na(e_data_norm[,ref]))]
  l2TMM <- vector()
  TMM <- vector()

  for(k in 1:ncol(data)){
    # remove genes where the sample of interest has a value of NA
    if(length(which(is.na(data[,k]))) > 0){
      trimmed.data = data[-which(is.na(data[,k])),]
    }else{
      trimmed.data = data
    }

    # determine the total number of reads in the reference and interest samples
    Nr = colSums(trimmed.data,na.rm=TRUE)[ref]
    Nk = colSums(trimmed.data,na.rm=TRUE)[k]

    names2 <- names[-which(is.na(data[,k]))]

    # create empty matrices for M and A values
    Mgkr <- matrix(nrow=nrow(trimmed.data),ncol=1)
    rownames(Mgkr) <- names2
    Agkr <- matrix(nrow=nrow(trimmed.data),ncol=1)
    rownames(Agkr) <- names2

    # create empty vectors for counts
    Ygr <- vector()
    Ygk <- vector()

    for(g in 1:nrow(trimmed.data)){
      Ygr[g] = trimmed.data[g,ref]
      Ygk[g] = trimmed.data[g,k]

      # calculate M and A for each gene and each sample
      #Mgkr[g,1] = log2(Ygk[g]/Nk) / log2(Ygr[g]/Nr)
      Mgkr[g,1] = log2((Ygk[g]/Nk) / (Ygr[g]/Nr))
      Agkr[g,1] = 1/2 * log2(Ygk[g]/Nk * Ygr[g]/Nr)
    }

    # trim M values by qm
    Qm = quantile(Mgkr[,1],c(qm,1-qm),na.rm=TRUE)
    Mgkr = Mgkr[which(Qm[1] < Mgkr[,1] & Mgkr[,1] < Qm[2]),]
    Agkr = Agkr[which(rownames(Agkr) %in% names(Mgkr)),]

    # trim A values by qa
    Qa = quantile(Agkr,c(qa,1-qa),na.rm=TRUE)
    Agkr = Agkr[which(Qa[1] < Agkr & Agkr < Qa[2])]
    Mgkr = Mgkr[which(names(Mgkr) %in% names(Agkr))]

    # only keep counts that remain after trimming
    Ygr = Ygr[which(names2 %in% names(Mgkr))]
    Ygk = Ygk[which(names2 %in% names(Mgkr))]

    wgkr <- vector()
    for(g in 1:length(Mgkr)){
      # calculate weights
      wgkr[g] = (Nk-Ygk[g])/(Nk*Ygk[g]) + (Nr-Ygr[g])/(Nr*Ygr[g])
    }

    # calculate log2(TMM) values for sample
    l2TMM[k] = sum(wgkr*Mgkr) / sum(wgkr)

    # normalize sample by log2(TMM)
    #e_data_norm[,k] <- log2(e_data_norm[,k]) - l2TMM[k]
    TMM[k] <- 2^l2TMM[k]
    e_data_norm[,k] <- e_data_norm[,k] / TMM[k]
  }

  e_data <- data.frame(e_data[,edata_id],e_data_norm)
  names(e_data)[1] <- edata_id
  names(TMM) <- names(data)
  TMM[is.na(TMM)] <- 1

  return(list(normed_data=e_data, location_param=NULL, scale_param=TMM))
}


#' Total sum scaling normalization of count data
#'
#' The method normalizes count data by the total number of counts in each sample
#'
#' @param e_data a \eqn{p \times n} data.frame of count data, where \eqn{p} is the number of features and \eqn{n} is the number of samples. Each row corresponds to data for a feature, with the first column giving the feature name.
#'@param edata_id character string indicating the name of the feature identifier. Usually obtained by calling \code{attr(omicsData, "cnames")$edata_cname}.
#'
#' @details Count data is normalized by the total count, dividing by the total count of each sample
#'
#' @return List containing 3 elements: norm_data is a data.frame with same structure as e_data that contains the TSS-normalized data, location_param is NULL, scale_param is a numeric vector of the total count in each sample.
#'
#' @examples
#' library(mintJansson)
#' data(cDNA_hiseq_data)
#' cDNA_TSS <- Quant_Norm(e_data = cDNA_hiseq_data$e_data, edata_id = attr(cDNA_hiseq_data, "cnames")$edata_cname)
#' norm_factors <- attr(cDNA_TSS,"data_info")$scale_param
#'
#' @author Allison Thompson, Lisa Bramer
#'

TSS_Norm <- function(e_data, edata_id){
  e_data_norm <- e_data[,-which(colnames(e_data)==edata_id)]

  # normalize by the total sum in each sample
  scale_param <- colSums(e_data_norm, na.rm=T)
  e_data_norm <- apply(e_data_norm, 2, function(x) x / sum(x, na.rm = TRUE))

  e_data <- data.frame(e_data[,which(colnames(e_data)==edata_id)],e_data_norm)
  colnames(e_data)[1] <- edata_id

  return(list(normed_data=e_data, location_param=NULL, scale_param=scale_param))
}
