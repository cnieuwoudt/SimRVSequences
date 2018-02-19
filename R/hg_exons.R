#' Human Exon Data
#'
#' This data set catalouges the positions, in base pairs, of the exons residing on chromosomes 1-22.  The positions of the exons were obtained from the UCSC Genome Broswer, using the NCBI RefSeq curated exons.  In this dataset we have combined overlapping exons into a single entry.
#'
#' NOTE: (REMOVE BEFORE RELEASE) CombineExons_NCBI.R contains details
#'
#' @docType data
#'
#' @format A data set with 223289 rows and 4 variables:
#' \describe{
#'   \item{chrom}{Numeric. The chromosome number.}
#'   \item{exonStart}{Numeric. The exon's starting position, in base pairs.}
#'   \item{exonStop}{Numeric. The exon's ending position, in base pairs.}
#'   \item{geneName}{Character. The names of the genes contained in the combined exon.}
#' }
"hg_exons"
