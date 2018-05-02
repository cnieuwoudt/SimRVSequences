#' Human Exon Data
#'
#' This data set catalouges the positions, in base pairs, of the exons residing on human chromosomes 1-22.  The positions of the exons were obtained from the UCSC Genome Broswer, using the hg38 database and ncbiRefSeqCurated table.  In this dataset we have combined overlapping exons into a single entry.
#'
#'
#' In \code{hg_exons} overlapping exons are combined into a single observation. When exons from genes with different NCBI accession numbers have been combined the variable \code{NCBIref} will contain multiple accession numbers, each separated by a comma.  We note that different accession numbers may exist for isoforms or transcript variants of the same gene.
#'
#'
#' @docType data
#'
#' @references Karolchik, D., Hinrichs, A. S., Furey, T. S., Roskin, K. M., Sugnet, C. W., Haussler, D., and Ken, W. J. (2004). The UCSC Table Browser data retrieval tool. \emph{Nucleic Acids Res}.
#' @references Kent, W. J., Sugnet, C. W., Furey, T. S., Roskin, K. M., Pringle, T. H., Zahler, A. M., and Haussler, D. (2002). The human genome browser at UCSC. \emph{Genome Res}, 12(6):996-1006.
#'
#' @format A data set with 223289 rows and 4 variables:
#' \describe{
#'   \item{chrom}{Numeric. The chromosome number.}
#'   \item{exonStart}{Numeric. The exon's starting position, in base pairs.}
#'   \item{exonStop}{Numeric. The exon's ending position, in base pairs.}
#'   \item{NCBIref}{Character. The NCBI reference seqence accession number of the gene(s) in which the exon(s) reside.}
#' }
"hg_exons"
