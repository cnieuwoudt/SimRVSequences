#' Summary function for object of class fam_study
#'
#' @param fam_study The object returned by sim_RVstudy
#'
#' @return Not sure yet..
#' @export
#'
#' @examples
#' #create example
summarize_study <- function(fam_study){
  Fids <- sort(unique(fam_study$ped_files$FamID))

  aff_allele_counts <- lapply(Fids, function(x){
    affected_allele_count(ped_haps = fam_study$ped_haplos[fam_study$haplo_map$FamID == x, ],
                          hap_map  = fam_study$haplo_map[fam_study$haplo_map$FamID == x, ],
                          ped_file = fam_study$ped_files[fam_study$ped_files$FamID == x, ])
  })

  allele_count_dat <- do.call(rbind, aff_allele_counts)
  allele_count_dat <- cbind(Fids, allele_count_dat)
  colnames(allele_count_dat) <- c("FamID", paste0("SNV", sep = "_", seq(1:ncol(fam_study$ped_haplos))))
  return(allele_count_dat)
}



#' Determine total number of like-alleles shared by affecteds in a family
#'
#' @param ped_haps sparse matrix.  The familial haplotype data.
#' @param hap_map data frame. Mapping data: maps individuals in \code{ped_file} to the haplotypes in \code{ped_haps}.
#' @param ped_file data frame.  Pedigree file, a for single family; i.e. not an entire study of ped files.
#'
#' @importFrom Matrix colSums
#'
#' @return A list of sharing counts, in the same order as the SNVs in \code{ped_haps}.
#' @export
#'
#' @examples
#' #no examples yet
affected_allele_count <- function(ped_haps, hap_map, ped_file){
  #determine the locations (rows) of the affecteds in ped_haps
  aff_IDs <- ped_file$ID[ped_file$affected]
  aff_rows <- which(hap_map$ID %in% aff_IDs)

  total_count <- colSums(ped_haps[aff_rows, ])
  return(total_count)
}
