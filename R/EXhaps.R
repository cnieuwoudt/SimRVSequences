#' Example Haplotypes Dataset
#'
#' This data set contains 20,000 haplotypes (i.e. rows) spanning 500 single-nucleotide variants.
#'
#'  This dataset is intended to accompany the \code{EXmuts} dataset and is used for demonstration purposes. For additional information regarding the SNVs contained in \code{EXhaps} please refer to the documentation for the \code{\link{EXmuts}} data set.
#'
#' @docType data
#'
#' @seealso \code{\link{EXmuts}}, \code{\link{read_slim}}
#'
#' @format A sparseMatrix of class dgCMatrix with 20000 rows and 500 variables.  Each row represents a haplotype and each column represents an SNV locus.
"EXhaps"
