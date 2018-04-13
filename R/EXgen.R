#' Example Genomes Dataset
#'
#' This data set contains the haplotypes for 10000 diploid individuals along segments of length 20,000 base pairs from human chromosomes 1-5. This dataset is intended to accompany the EXmut dataset and is used for demonstration purposes.
#'
#' @docType data
#'
#' @seealso \code{\link{EXmut}}, \code{\link{read_slim}}
#'
#' @format A sparseMatrix (formal class dgCMatrix) with 20000 rows and 173 variables.  Each row is the haplotype for 1 diploid individual, each column represents a single SNV locus.  Each column in \code{EXgen} is described by a row in \code{\link{EXmut}}.
"EXgen"
