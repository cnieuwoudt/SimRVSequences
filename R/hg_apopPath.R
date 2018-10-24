#' Pseudo-Apoptosis Pathway Data
#'
#' This data set cataloges the positions of the exons residing in the 25 genes that have the highest interaction with the TNFSF10 gene USCS Genome Browser's Gene Interaction Tool.  We use this data set to model a pseudo-apoptosis sub-pathway in our examples and vignette.
#'
#'  In this data set, any overlapping exonshave been combined into a single observation.  When exons from genes with different NCBI accession numbers have been combined the variable \code{NCBIref} will contain multiple NCBI reference sequence accession numbers, each separated by a comma.  We note that different NCBI reference sequence accession numbers may exist for isoforms or transcript variants of the same gene.
#'
#' @docType data
#'
#' @references Karolchik, D., Hinrichs, A. S., Furey, T. S., Roskin, K. M., Sugnet, C. W., Haussler, D., and Ken, W. J. (2004). The UCSC Table Browser data retrieval tool. \emph{Nucleic Acids Res}. Accessed on 20 February 2018.
#' @references Kent, W. J., Sugnet, C. W., Furey, T. S., Roskin, K. M., Pringle, T. H., Zahler, A. M., and Haussler, D. (2002). The human genome browser at UCSC. \emph{Genome Res}, 12(6):996-1006.
#' @references Poon, H., Quirk, C., DeZiel, C., and Heckerman, D. (2014). Literome: Pubmed-scale genomic knowledge base in the cloud. \emph{Bioinformatics}, 30:2840-2842.
#'
#' @format A data set with 253 rows and 5 variables:
#' \describe{
#'   \item{chrom}{Numeric. The chromosome number.}
#'   \item{exonStart}{Numeric. The exon's starting position, in base pairs.}
#'   \item{exonStop}{Numeric. The exon's ending position, in base pairs.}
#'   \item{NCBIref}{Character. The NCBI reference seqence accession number of the gene(s) in which the exon(s) reside.}
#'   \item{gene}{Character. The name(s) of the gene.}
#' }
"hg_apopPath"
