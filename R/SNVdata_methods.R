#' Constructor function for an object of class famStudy
#'
#' @param Haplotypes sparseMatrix. A sparse matrix of haplotype data, which contains the haplotypes for unrelated individuals representing the founder population.  Rows are assumed to be haplotypes, while columns represent SNVs.  If the \code{\link{read_slim}} function was used to import SLiM data to \code{R}, users may supply the sparse matrix \code{Haplotypes} returned by \code{read_slim}.
#' @param Mutations Data frame. A data frame that catalogs the SNVs in \code{haplos}.  If the \code{\link{read_slim}} function was used to import SLiM data to \code{R}, the data frame \code{Mutations} is of the proper format for \code{SNV_map}.  However, users must add the variable \code{is_CRV} to this data frame, see details.
#' @param Samples Format-free.  Describe me...
#' @param MetaData Format-free.  Describe me...
#'
#' @return an object of class \code{SNVdata}.
#' @export
SNVdata <- function(Haplotypes, Mutations, Samples = NULL, MetaData = NULL) {

  #check SNV_map for possible issues
  check_SNV_map(Mutations)

  if (!"marker" %in% colnames(Mutations)) {
    Mutations$marker <- make.unique(paste0(Mutations$chrom, sep = "_", Mutations$position))
  }

  if (nrow(Mutations) != ncol(Haplotypes)) {
    stop("\n nrow(Mutations) != ncol(Haplotypes). \n Mutations must catalog every SNV in Haplotypes.")
  }



  #create list containing all relavant of SNVdata information
  SNV_data = list(Haplotypes = Haplotypes, Mutations = Mutations,
                  Samples = Samples, MetaData = MetaData)

  class(SNV_data) <- c("SNVdata", class(SNV_data))
  return(SNV_data)
}

#' Check to see if object is of class ped
#'
#' @param x An R object.
#'
#' @return Logical. Indicates if \code{x} is of class \code{ped}.
#'
#' @keywords internal
is.SNVdata <- function(x) {
  return(inherits(x, "SNVdata"))
}

#' Import 1000 genomes exon data formatted for sim_RVstudy
#'
#' Import 1000 genomes exon data formatted for sim_RVstudy
#'
#' We expect that `pathwayDF` does not contain any overlapping segments.  Users may combine overlapping exons into a single observation with our `combine_exons` function.  For additional information regarding the `combine_exons` function please  execute `help(combine_exons)` in the console.
#'
#' @param chrom Numeric.  The chromosome number(s).  A numeric list of chromosome numbers representing the 1000 genomes exon-data to be imported.
#' @param pathway_df Data frame. (Optional) A data frame that contains the positions for each exon in a pathway of interest.  This data frame must contain the variables `chrom`, `exonStart`, and `exonEnd`. See Details.
#'
#' @return An object of class \code{SNVdata} or a list of objects of class \code{SNVdata}, i.e. one for each chomosome that was imported.
#' @export
#'
#' @examples
#' exdata = import_SNVdata(21)
#'
#' head(exdata$Mutations)
#' exdata$Haplotypes[1:20, 1:10]
import_SNVdata <- function(chrom, pathway_df = NULL){

  if (any(!chrom %in% seq(1:22))){
    stop("\n We expect 'chrom' to be a numeric list of automosome numbers.\n
           Note: We do not provide sex chromosomes (i.e. chromosomes X and Y)")
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

  #----------------------#
  # Identify Pathway RVs #
  #----------------------#
  #if pathway data has been supplied, identify pathway SNVs
  if (!is.null(pathway_df)) {
    Mutations <- identify_pathwaySNVs(markerDF = Mutations,
                                      pathwayDF = pathway_df)
  }

  return(SNVdata(Haplotypes = Haplotypes, Mutations = Mutations))
}
