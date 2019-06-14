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

#' Import 1000 genomes exon data formatted for sim_RVstudy
#'
#' @param chrom The chromosome number.  Can be a numeric list of chromosome numbers or a character argument set to \code{"ALL"} if data from all chromosomes should be imported.
#'
#' @return An object of class \code{SNVdata} or a list of objects of class \code{SNVdata}, i.e. one for each chomosome that was imported.
#' @export
#'
#' @examples
#' exdata = import_SNVdata(0)
#'
#' head(exdata$Mutations)
#' exdata$Haplotypes[1:20, 1:10]
import_SNVdata <- function(chrom){

  if (any(chrom == 0)){
    chrom = seq(1:22)
  } else if (!(class(chrom) %in% c("numeric", "integer"))){
    stop("\n chrom must be a numeric list of chromosome numbers.\n
           Note if chrom = 0, all autosomes are imported (i.e. chroms 1 - 22 )")
  }

  #store the object names for the imported data by chromosome number
  object_names <- paste0("SNVdata_chrom", chrom)


  #import the formatted date from github
  for (i in 1:length(chrom)){
    load(url(paste0("https://github.com/cnieuwoudt/1000-Genomes-Exon-Data/raw/master/Formatted-SNVdata/SNVdata_chrom",
                    chrom[i], ".rda", sep = "")))
  }

  #store each to SNVdata object to a list
  chrom_dat = lapply(object_names, function(x){get(x)})

  if (length(chrom) > 1){
    #combine the data from the different chromosomes, for ease of use
    Haplotypes <- do.call(cbind, lapply(chrom_dat, `[[`, 1))
    Mutations <- do.call(rbind, lapply(chrom_dat, `[[`, 2))
  } else {
    Haplotypes <- chrom_dat[[1]]$Haplotypes
    Mutations <- chrom_dat[[1]]$Mutations
  }

  return(SNVdata(list(Haplotypes = Haplotypes, Mutations = Mutations)))
}
