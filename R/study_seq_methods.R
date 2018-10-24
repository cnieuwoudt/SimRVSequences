#' Summary function for object of class fam_study
#'
#' @param fam_study The object returned by sim_RVstudy
#'
#' @return Not sure yet..
#' @export
#'
#' @examples
#' library(SimRVSequences)
#'
#' #load pedigree, haplotype, and mutation data
#' data(study_peds)
#' data(EXmuts)
#' data(EXhaps)
#'
#' #create variable is_CRV in EXmuts to identify causal
#' #rare variants, from which to sample familial variants.
#' EXmuts$is_CRV = FALSE
#' EXmuts$is_CRV[c(2, 3, 12, 24)] = TRUE
#'
#' #supply required inputs to the sim_RVstudy function
#' seqDat = sim_RVstudy(ped_files = study_peds,
#'                      SNV_map = EXmuts,
#'                      haplos = EXhaps)
#'
#' #to count the number of SNVs shared by the disease-affected
#' #relatives in each pedigree, supply the output returned by
#' #sim_RVstudy to count_affectedRV
#' count_affectedRV(seqDat)
#'
#'
count_affectedRV <- function(fam_study){
  Fids <- sort(unique(fam_study$ped_files$FamID))

  aff_allele_counts <- lapply(Fids, function(x){
    affected_allele_count(ped_haps = fam_study$ped_haplos[fam_study$haplo_map$FamID == x, ],
                          hap_map  = fam_study$haplo_map[fam_study$haplo_map$FamID == x, ],
                          ped_file = fam_study$ped_files[fam_study$ped_files$FamID == x, ])
  })

  allele_count_dat <- do.call(rbind, aff_allele_counts)
  allele_count_dat <- cbind(Fids, allele_count_dat)
  colnames(allele_count_dat) <- c("FamID", fam_study$SNV_map$marker)

  #remove SNVs not carried by affecteds
  allele_count_dat <- allele_count_dat[, which(colSums(allele_count_dat) != 0)]
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
