% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SlimFunctions.R
\name{combine_identicalmutations}
\alias{combine_identicalmutations}
\title{Combine identical mutations}
\usage{
combine_identicalmutations(mutmap, hapmat, mut_type)
}
\arguments{
\item{mutmap}{data.frame. The SNV_map with identical mutations}

\item{hapmat}{sparseMatrix. The sparseMatrix of haplotypes}

\item{mut_type}{character. the name of the mutation type for which to re-code identical mutations}
}
\value{
A list containing:

\item{\code{hapmat}}{A sparse matrix of haplotypes. See details.}

\item{\code{mutmap}}{A data frame cataloging SNVs in \code{hapmap}.}
}
\description{
Assuming that all mutations are of the same type (in SLIM simulation), different mutations at the same site are actually identical mutations from different lineages. This function re-codes these mutations so that they are only cataloged once.
}
\keyword{internal}
