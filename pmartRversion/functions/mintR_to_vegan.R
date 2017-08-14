#' Transforms a mintR data object into a vegan data object
#'
#' This function takes a mintR data object and transforms it into a data object suitable for use in the vegan package
#'
#' @param omicsData an object of one of the classes "cDNAdata", "rRNAdata", or "gDNAdata"
#'
#' @return A data object for use with the package 'vegan'
#'
#' @author Allison Thompson
#'
#' @references
#'   Jari Oksanen, F. Guillaume Blanchet, Michael Friendly, Roeland Kindt, Pierre Legendre, Dan McGlinn, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Eduard Szoecs and Helene Wagner (2016). vegan: Community Ecology Package. R package version 2.4-1. https://CRAN.R-project.org/package=vegan
#'
#' @examples
#'
#' library(vegan)
#' library(mintJansson)
#' data(rRNA_data)
#'
#' data <- mintR_to_vegan(rRNA_data)
#' ord <- metaMDS(data)
#' plot(ord, display="sites", type="points")
#'
#' @export
mintR_to_vegan <- function(omicsData){
  # Extract e_data
  e_data <- omicsData$e_data
  
  # Format so this can be a matrix
  rownames(e_data) <- e_data[,which(colnames(e_data) == attr(omicsData, "cnames")$edata_cname)]
  e_data <- e_data[,-which(colnames(e_data) == attr(omicsData, "cnames")$edata_cname)]
  
  # Change NA to 0
  e_data[is.na(e_data)] <- 0
  
  # Transpose matrix, so samples are in rows (not columns)
  e_data <- t(e_data)
  
  # Return data
  return(e_data)
  
}