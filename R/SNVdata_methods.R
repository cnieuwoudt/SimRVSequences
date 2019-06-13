#' Constructor function for an object of class famStudy
#'
#' @param SNV_data The list of items returned by the \code{read_slim} and \code{import_exons} function.
#'
#' @return an object of class \code{SNVdata}.
#' @export
SNVdata <- function(SNV_data) {
  class(SNV_data) <- c("SNVdata", class(SNV_data))
  return(SNV_data)
}

#' Check to see if object is of class ped
#'
#' @param x An R object.
#'
#' @return Logical. Indicates if \code{x} is of class \code{ped}.
#'
#' @export
is.SNVdata <- function(x) {
  return(inherits(x, "SNVdata"))
}
