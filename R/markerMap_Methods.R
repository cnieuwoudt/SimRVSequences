#' Create an object of class markerMap.
#'
#' Create a markerMap object, required input for \code{\link{sim_RVstudy}} and \code{\link{sim_RVseq}} functions.
#'
#' The argument \code{markerDF} is a data frame consisting of SNP data.  If \code{\link{read_slim}} was used to read a slim output file, users may supply the first item returned by \code{read_slim} as \code{markerDF}.
#'
#' The named columns in the data frame, \code{markerDF}, are described as follows:
#' \describe{
#'  \item{marker}{(OPTIONAL) An identifying name/number for the variant. (If missing this field will be automatically generated from the columns named "chrom" and "position" as "chrom_position"). NOTE: These names must match \emph{exactly} the column names (i.e. the marker names) in the haplotype_dist dataframe.}
#'  \item{chrom}{The number of the chromosome on which the marker resides.}
#'  \item{position}{The location of the variant, in basepairs.}
#'  \item{possibleRV}{(OPTIONAL) Logical variable, indicates whether the variant may be chosen as a causal familial variant.  If this field is missing and \code{pathwayDF} is not provided it is assumed that all variants listed are candidates for selection.}
#'  \item{probCausal}{(OPTIONAL) The probability that the variant is selected as a familial RV.  If this field is missing, it will be assumed that all variants for which possibleRV is TRUE are equally likely.}
#' }
#'
#' TO FIX: CURRENTLY NOT ACCOUNTING FOR CARRIER PROB WHEN CHOOSING VARIANTS
#'
#' Describe rare variant selection and pathwayDF.
#'
#' @param markerDF A data frame containing pertinent information on markers.   See details.
#' @param pathwayDF A data frame containing pathway data.  See details.
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
markerMap <- function(markerDF, pathwayDF = NULL) {
  # Set up a new object of class markerMap

  if (!"chrom" %in% colnames(markerDF) |
      !"position" %in% colnames(markerDF)) {
    stop('please provide a data.frame with named variables: chrom and position')
  }

  if (any(is.na(markerDF))) {
    stop('markerDF contains missing values')
  }

  if (any(!is_int(markerDF$chrom))) {
    stop('chromosome numbers must be expressed as integers')
  }

  #add marker name, if not provided
  if (is.null(markerDF$marker)) {
    markerDF$marker <- paste0(markerDF$chrom, sep = "_", markerDF$position)
  }

  #create possibleRV
  if (is.null(pathwayDF) & is.null(markerDF$possibleRV)) {
    markerDF$possibleRV <- TRUE
  } else if (!is.null(pathwayDF)) {
    if (!is.null(pathwayDF) & !is.null(markerDF$possibleRV)){
      warning('Redefining possibleRV variable based on the data provided in pathwayDF')
    }

    markerDF <- do.call(rbind, lapply(unique(markerDF$chrom), function(x){
      identify_possibleRVs(pathwayDF[pathwayDF$chrom == x, ],
                           markerDF[markerDF$chrom == x, ])
    }))
  }

  #add probCausal, if not provided
  if (is.null(markerDF$probCausal)) {
    markerDF$probCausal <- markerDF$possibleRV/sum(markerDF$possibleRV)
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

#' Identify possible variants based on defined pathway
#'
#' @param path_by_chrom The pathway data for the chromosome under consideration
#' @param marker_map_by_chrom The marker_map for the chromosome under consideration
#'
#' @return marker_map_by_chrom with possibleRV identified
#' @keywords internal
identify_possibleRVs <- function(path_by_chrom, marker_map_by_chrom){
  if(nrow(path_by_chrom) == 0){
    marker_map_by_chrom$possibleRV <- FALSE
  } else {
    #since we will want to include variants that occur at the first base pair
    #location, subtracting 1 from start positions
    #Similarly, since we want the first mutation in marker_map_by_chrom
    #to have a bin, subtracting 1 from this position
    if(min(marker_map_by_chrom$position) != min(path_by_chrom$exonStart)){
      cbreaks <- sort(c((min(marker_map_by_chrom$position) - 1),
                        (path_by_chrom$exonStart - 1),
                        unique(c(path_by_chrom$exonEnd,
                                 max(marker_map_by_chrom$position)))))

      keep_bins <- seq(2, length(cbreaks) - 1, by = 2)
    } else {
      cbreaks = sort(c((path_by_chrom$exonStart - 1),
                       unique(c(path_by_chrom$exonEnd,
                                max(marker_map_by_chrom$position)))))
      keep_bins <- seq(1, length(cbreaks) - 1, by = 2)
    }

    if(any(duplicated(cbreaks))){
      stop("Expecting non-overlapping exonStart and ExonEnd positions. \n Please combine overlapping segments into a single entry.")}

    marker_map_by_chrom$possibleRV <- cut(marker_map_by_chrom$position,
                                          breaks = cbreaks,
                                          labels = FALSE)
    marker_map_by_chrom$possibleRV <- marker_map_by_chrom$possibleRV %in% keep_bins
  }

  return(marker_map_by_chrom)
}
