#' Create an object of class markerMap.
#'
#' Create a markerMap object, required input for \code{\link{sim_RVstudy}} and \code{\link{sim_RVseq}} functions.
#'
#' The named columns in the data frame, \code{markerDF}, are described as follows:
#' \describe{
#'  \item{marker}{(OPTIONAL) An identifying name/number for the variant. (If missing this field will be automatically generated from the columns named "chrom" and "position" as "chrom_position"). NOTE: These names must match \emph{exactly} the column names (i.e. the marker names) in the haplotype_dist dataframe.}
#'  \item{chrom}{The number of the chromosome on which the marker resides.}
#'  \item{position}{The location of the variant, in basepairs.}
#'  \item{pathwayID}{A pathwayID for the variant. (**Should be optional, I don't think we use this yet**)}
#'  \item{possibleRV}{(OPTIONAL) Logical variable, indicates whether the variant may be chosen as a causal familial variant.  If this field is missing, it is assumed that all variants listed are candidates for selection.}
#'  \item{probCausal}{(OPTIONAL) The probability that the variant is selected as a familial RV.  If this field is missing, it will be assumed that all variants for which possibleRV is TRUE are equally likely.}
#' }
#'
#'
#' @param markerDF A data frame containing pertinent information on markers.   See details.
#'
#' @return An object of class markerMap.
#' @export
#'
#' @examples
#' data(mark_map)
#'
#' head(mark_map)
#' class(mark_map)
#'
#' markObj <- markerMap(mark_map)
#' head(markObj)
#' class(markObj)
#'
markerMap <- function(markerDF) {
  # Set up a new object of class markerMap

  if (!"chrom" %in% colnames(markerDF) |
      !"position" %in% colnames(markerDF) |
      !"pathwayID" %in% colnames(markerDF)) {
    stop('please provide a data.frame with the following variables: chrom, position, and pathwayID')
  }

  #add marker name, if not provided
  if (is.null(markerDF$marker)) {
    markerDF$marker <- paste0(markerDF$chrom, sep = "_", markerDF$position)
  }

  #add possibleRV, if not provided
  if (is.null(markerDF$possibleRV)) markerDF$possibleRV <- TRUE

  #add probCausal, if not provided
  if (is.null(markerDF$probCausal)) {
    markerDF$probCausal <- markerDF$possibleRV/sum(markerDF$possibleRV)
  }

  if (any(is.na(markerDF))) {
    stop('markerDF contains missing values')
  }

  if (any(!is_int(markerDF$chrom))) {
    stop('chromosome numbers must be expressed as whole numbers')
  }


  obj <- markerDF
  class(obj) <- c("markerMap", class(markerDF))
  return(obj)
}


#' Check to see if object is of class markerMap
#'
#' @param x An R object.
#'
#' @return Logical. Indicates if \code{x} is of class \code{markerMap}.
#' @keywords internal
#'
is.markerMap <- function(x) {
  return(inherits(x, "markerMap"))
}

