#' Combine overlapping exons into a single observation
#'
#' @param exon_data The exon data to combine.  This data must unclude named variables: 'chrom', a chromosome identifer; 'exonStart', the first position of the exon in base pairs; and 'exonEnd', the last position of the exon in base pairs.
#' @param zero_start. Logical.  \code{zero_start = TRUE} if the first position is 0-based, and \code{zero_start = FALSE} if the first position is 1-based.
#'
#' @return A matrix of exon unions, i.e. non-overlapping observations, by chromosome.
#' @export
#'
#' @examples
#' # create an example data frame that contains the
#' # the variables: chrom, exonStart, and exonEnd
#' exDat <- data.frame(chrom     = c(1, 1, 1, 2, 2, 2),
#'                     exonStart = c(1, 2, 5, 1, 3, 3),
#'                     exonEnd   = c(3, 4, 7, 4, 5, 6))
#'
#' exDat
#'
#' # supply exDat to combine_exons to combine overlapping segments
#' # into a single observation.
#' combine_exons(t1)
#'
combine_exons <- function(exon_data, zero_start = FALSE){

  #check to see if exon_data contains the three named
  #columns we require to combine exons
  if(any(!(c("chrom", "exonStart", "exonEnd") %in% colnames(exon_data)))){
    stop("exon_data does not include named columns: 'chrom', 'exonStart', and 'exonEnd'.")
  }

  exstart = which(colnames(exon_data) == "exonStart")
  exend = which(colnames(exon_data) == "exonEnd")

  cexons <- do.call(rbind, lapply(unique(exon_data$chrom), function(x){
    combine_exons_by_chrom(chrom = x,
                           start_stop_dat = exon_data[exon_data$chrom == x, c(exstart, exend)])
  }))

 return(cexons)

}


#' Combine exons within a chromosome
#'
#' @param chrom the chromosome number
#' @param start_stop_dat the exon start and stop data, i.e. two columns of a dataframe or matrix.  Start positons should be contained in column 1 and stop positions in column 2.
#'
#' @importFrom intervals interval_union
#' @importFrom intervals Intervals
#'
#' @return a matrix with combined exons for a single chromosome
#' @keywords internal
combine_exons_by_chrom <- function(chrom, start_stop_dat){
  #combine exons in this chromosome using the interval_union
  #and Intervals functions provided by the intervals package
  ex_mat <- interval_union(Intervals(start_stop_dat))@.Data

  #add chromosome variable
  ex_mat <- cbind(rep(chrom, nrow(ex_mat)), ex_mat)

  #add column names
  colnames(ex_mat) <- c('chrom', 'exonStart', 'exonEnd')

  return(ex_mat)
}
