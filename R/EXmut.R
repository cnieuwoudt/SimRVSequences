#' Example Mutations Dataset
#'
#' This data set contains mutations along segments of length 20,000 base pairs from human chromosomes 1-5. This dataset is intended to accompany the EXgen dataset and is used for demonstration purposes.
#'
#' @docType data
#'
#' @seealso EXgen
#'
#' @format A data set with 250 rows and 4 variables:
#' \describe{
#'   \item{colID}{Numeric. The corresponding column number of the SVN in the EXgen dataset.}
#'   \item{chrom}{Numeric. The chromosome number.}
#'   \item{position}{Numeric. The location of the SNV, in base pairs.}
#'   \item{afreq}{Numeric. The derived allele frequency of the SNV.}
#'   \item{marker}{Character. The names of the genes contained in the combined exon.}
#'   \item{pathwaySNV}{Logical. Indicates if the SNV is located within the pathway.}
#' }
"EXmut"
