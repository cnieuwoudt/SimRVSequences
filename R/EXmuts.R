#' Example Mutations dataset
#'
#' This data set catalogs the 500 single-nucleotide variants contained in the \code{EXhaps} dataset.  This dataset is intended to accompany the \code{EXhaps} dataset; each row of \code{EXmuts} describes a column (i.e. SNV) in \code{EXhaps}.
#'
#' Together, the \code{EXmuts} and \code{EXhaps} datasets represent example output of the \code{read_slim} function.  The \code{EXhaps} data set represents the sparse matrix \code{Haplotypes} returned by \code{read_slim}, and the \code{EXmuts} data set represents the \code{Mutations} data frame returned by \code{read_slim}.  This toy data set is used primarily for demonstration, and does not represent a complete sample of all mutations that would accompany a genome-wide, exon-only sequence data.
#'
#'
#' @docType data
#'
#' @seealso \code{\link{EXhaps}}, \code{\link{read_slim}}
#'
#' @format A data set with 500 rows and 6 variables:
#' \describe{
#'   \item{colID}{Numeric. The corresponding column number of the SVN in the EXgen dataset.}
#'   \item{chrom}{Numeric. The chromosome number.}
#'   \item{position}{Numeric. The location of the SNV, in base pairs.}
#'   \item{afreq}{Numeric. The derived allele frequency of the SNV.}
#'   \item{marker}{Character. The names of the genes contained in the combined exon.}
#'   \item{pathwaySNV}{Logical. Indicates if the SNV is located within the pathway.}
#' }
"EXmuts"
