#' Normalize data
#'
#' Normalizes the data using the specified normalization function
#'
#' For count data (16S data), the default normalization is Cumulative Sum Scaling \code{norm_fn="css"}.
#' The choices for normalization currently available, \code{norm_fn}, are:
#' \tabular{ll}{
#'  \code{"percentile"}  \tab Standardize the data by dividing each feature in \code{e_data} by the sample-wide qth percentile and multiply by the gloabal qth percentile\cr
#'  \tab \cr
#'  \code{"tss"}  \tab Standardize the data by dividing each feature in \code{e_data} by sum of the sample-wide counts\cr
#'  \tab \cr
#'  \code{"rarefy"} \tab Normalize the data by subsampling down to specified library size\cr
#'  \tab \cr
#'  \code{"poisson"} \tab Normalize the data by sampling from a Poisson distribution with appropriate mean value, see Section 2.2 of Li et al. (2013)\cr
#'  \tab \cr
#'  \code{"deseq"} \tab Normalization method used in DESeq and DESeq2, which uses size factors to standardize sequencing depths across samples\cr
#'   \tab \cr
#'  \code{"css"}  \tab Normalize the data by dividing each feature in \code{e_data} by the sum of the counts up to a specified quantile and multiplying by a global scaling factor\cr
#'  \tab \cr
#'  \code{"tmm"}  \tab Normalize the data using the trimmed mean of M values\cr
#'  \tab \cr
#'  }
#'
#' @param omicsData an object of the class 'gDNAdata', 'cDNAdata', or 'rRNAdata' usually created by \code{\link{as.gDNAdata}}, \code{\link{as.cDNAdata}}, or \code{\link{as.rRNAdata}}, respectively.
#' @param norm_fn character vector indicating the normalization function to use for normalization. See details for valid options.
#' @param normalize For count data, this function will only return the scale and location parameters - will not return normalized data unless this parameter is set to TRUE. This is due to later statistics requiring raw data for count data analyses. Default is FALSE.
#' @param ... additional arguments passed to the chosen normalization function.
#'
#' @return If normalize=FALSE, a list containing the location and scale parameters to use when normalizing the data. If normalize=TRUE, returns the omicsData object, where e_data has been normalized with the appropriate parameters and the scale and location parameters are returned as an attribute of the data.
#'
#' @author Kelly Stratton, Lisa Bramer, Bryan Stanfill, Allison Thompson
#' @references
#'  Li, Jun, and Robert Tibshirani. Finding consistent patterns: a nonparametric approach for identifying differential expression in RNA-Seq data. Statistical methods in medical research 22.5 (2013): 519-536.
#'  Anders, Simon and Wolfgang Huber. Differential expression analysis for sequence count data. Genome Biology 11:R106 (2010).
#'  Paulson, Joseph N, O Colin Stine, Hector Corrada Bravo, and Mihai Pop. Differential abundance analysis for microbial marker-gene surveys. Nature Methods. 10.12 (2013)
#'  Robinson, Mark D and Alicia Oshlack. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biology 11:R25 (2010).
#'
#' @examples
#' library(mintJansson)
#'
#' #Count data are passed to a quantile normalization by default
#' normalized_rRNAdata <- normalize_data(rRNA_data)
#' normalized_rRNAdata <- normalize_data(rRNA_data, norm_fn = 'percentile')
#'
#' #One can also use the TSS normalization for count data, the normalization function used by DESeq/DESeq2, normalization from SAMSeq (aka Poisson), another from metagenomeSeq (cumulative sum scaling normalization (CSS)), or yet another competitor edgeR (TMM normalization)
#' normalized_rRNAdata <- normalize_data(rRNA_data, norm_fn = "tss")
#' normalized_rRNAdata <- normalize_data(rRNA_data, norm_fn = "deseq")
#' normalized_rRNAdata <- normalize_data(rRNA_data, norm_fn = "poisson")
#' normalized_rRNAdata <- normalize_data(rRNA_data, norm_fn = "css")
#' normalized_rRNAdata <- normalize_data(rRNA_data, norm_fn = "tmm")
#'
#' #One could also rarefy the data - though this is highly NOT recommended
#' normalized_rRNAdata <- normalize_data(rRNA_data, norm_fn = "rarefy")
#'
#' @export
normalize_data <- function(omicsData, norm_fn, normalize=FALSE, ...){
  ## initial checks ##

  #Store data class as it will be referred to often
  dat_class <- class(omicsData)

  # check that omicsData is of the appropriate class
  if(!(dat_class%in% c("cDNAdata", "gDNAdata","rRNAdata"))) stop("omicsData is not an object of appropriate class")

  edata_id <- attr(omicsData, "cnames")$edata_cname
  samp_id <- attr(omicsData, "cnames")$fdata_cname

  #Use default normalization if not specified
  if(missing(norm_fn)){
    warning("norm_fn wasn't specified so the default for this data type was used, see 'help(normalize_data)'.")
    #For count data, cumulative sum scaling is used
    norm_fn <- "css"
    qg <- "median"
  }

  #Peel off the extra arguments passed
  extra_args <- list(...)

  ## Normalize ##
  edata <-  omicsData$e_data
  mintR_groupDF <- attr(omicsData, "group_DF")

  norm_fn <- try(match.arg(tolower(norm_fn),c("percentile","tss","rarefy","poisson","deseq","tmm","css")),silent=TRUE)

  # apply normalization scheme

  if(norm_fn%in%c("percentile","tss","poisson","deseq","tmm","css")){

    fn_to_use <- switch(norm_fn,percentile=Quant_Norm,tss=TSS_Norm,poisson=poisson_norm,deseq=med_scounts_norm,tmm=TMM_Norm,css=CSS_Norm)
    temp <- fn_to_use(e_data = edata, edata_id=edata_id, ...)
    norm_results <- list(norm_data=temp$normed_data, location_param=temp$location_param, scale_param=temp$scale_param)

  }else if(norm_fn=="rarefy"){
    message("The size of the subsample used in the Rarefy function is passed as the 'scale_param'.")
    temp <- Rarefy(e_data = edata,edata_id = edata_id,...)
    norm_results <- list(norm_data=temp$normed_data, location_param=temp$location_param, scale_param=temp$scale_param)

  }else{
    stop("The specified 'norm_fn' option provided is not currently available (or you made a typo) so no normalization was done.")
    norm_results <- list(norm_data=omicsData$e_data , location_param=NULL, scale_param=NULL)
  }

  # change attributes of omics_data to indicate that data has been normalized #
  attributes(omicsData)$data_info$data_norm = TRUE
  attributes(omicsData)$data_info$norm_method = norm_fn

  if(attributes(omicsData)$data_info$data_scale == "count" & !normalize){
    # add attributes with normalization parameters #
    data <- list()
    data$scale_param <- norm_results$scale_param
    data$location_param <- norm_results$location_param
    attr(data, "cnames") <- attr(omicsData, "cnames")
    attr(data, "group_DF") <- attr(omicsData, "group_DF")
    class(data) <- c("NormFactors",class(norm_results$scale_param),class(norm_results$location_param))
    #attributes(omicsData)$data_info$scale_param = norm_results$scale_param
    #attributes(omicsData)$data_info$location_param = norm_results$location_param
    attr(data, "raw_data") <- omicsData$e_data
    return(data)
  }else{
    # Replace the "e_data" in the original data with the normalized data
    omicsData$e_data <- norm_results$norm_data
    attributes(omicsData)$data_info$scale_param = norm_results$scale_param
    attributes(omicsData)$data_info$location_param = norm_results$location_param
    return(omicsData)
  }
}





