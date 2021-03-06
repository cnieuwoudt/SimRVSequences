% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_exons.R
\name{combine_exons_by_chrom}
\alias{combine_exons_by_chrom}
\title{Combine exons within a chromosome}
\usage{
combine_exons_by_chrom(chrom, start_stop_dat)
}
\arguments{
\item{chrom}{the chromosome number}

\item{start_stop_dat}{the exon start and stop data, i.e. two columns of a dataframe or matrix.  Start positons should be contained in column 1 and stop positions in column 2.}
}
\value{
a matrix with combined exons for a single chromosome
}
\description{
Combine exons within a chromosome
}
\keyword{internal}
