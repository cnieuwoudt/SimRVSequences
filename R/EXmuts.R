#' Example Mutations Dataset
#'
#' This data set catalogs the mutations in the \code{\link{EXhaps}} data set.  Each row of \code{EXmuts} describes a column (i.e. SNV) in \code{EXhaps}. This toy data set is used primarily for demonstration, and does not represent a complete sample of all mutations that would accompany a full exon-only simulation.
#'
#' @docType data
#'
#' @seealso EXgen
#'
#' @format A data set with 450 rows and 6 variables:
#' \describe{
#'   \item{colID}{Numeric. The corresponding column number of the SVN in the EXgen dataset.}
#'   \item{chrom}{Numeric. The chromosome number.}
#'   \item{position}{Numeric. The location of the SNV, in base pairs.}
#'   \item{afreq}{Numeric. The derived allele frequency of the SNV.}
#'   \item{marker}{Character. The names of the genes contained in the combined exon.}
#'   \item{pathwaySNV}{Logical. Indicates if the SNV is located within the pathway.}
#' }
"EXmuts"
