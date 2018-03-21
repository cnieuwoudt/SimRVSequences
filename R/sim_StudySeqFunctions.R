#' Draw Founder Genotypes from Haplotype Distribution Given Familial Risk Variant
#'
#'
#' @return list of familial founder genotypes
#' @export
#'
sim_FGenos <- function(founder_ids, RV_founder, FamID,
                       haplotype_dist, FamRV, marker_map) {

  #store chromosomal postion (i.e. list position in haplotype_dist) of risk variant
  RV_chrom_pos <- which(unique(marker_map$chrom) == marker_map$chrom[marker_map$marker == FamRV])

  fam_genos <- list()
  #NOTE: If we want to change from SNVs to another variant
  #type we will need to change the functionality of RV_status below

  for (k in 1:length(unique(marker_map$chrom))) {
    if (k == RV_chrom_pos) {
      #Determine which haplotypes carry the familial rare variant and which so not
      RV_haps <- which(haplotype_dist[[k]][, marker_map$colID[marker_map$marker == FamRV]] == 1)
      noRV_haps <- which(haplotype_dist[[k]][, marker_map$colID[marker_map$marker == FamRV]] == 0)

      #for the seed founder sample one haplotype from those that carry the RV
      # and one haplotype from those that DO NOT carry the RV
      RV_founder_dat = haplotype_dist[[k]][c(sample(x = RV_haps, size = 1),
                                             sample(x = noRV_haps, size = 1)), ]

      #sample 1 haplotype that does not carry the RV for the seed founder and
      #2 haplotypes that do not carry the RV for all other founders
      NRV_founder_dat <- haplotype_dist[[k]][sample(x = noRV_haps,
                                                    size = 2*length(founder_ids),
                                                    replace = TRUE), ]

      fam_genos[[k]] <- rbind(RV_founder_dat, NRV_founder_dat)
    } else {

      #since the familial RV does not reside on the kth chromosome we
      #sample founder haplotypes according to from the haplotype distribution

      fam_genos[[k]] = haplotype_dist[[k]][sample(x = c(1:nrow(haplotype_dist[[k]])),
                                                  size = 2*(length(founder_ids) + 1),
                                                  replace = TRUE), ]
    }
  }

  founder_genos <- do.call("cbind", fam_genos)
  #ramdomly permute RV founders rows so that the RV is not always paternally inherited.
  founder_genos[c(1,2), ] <- founder_genos[sample(x = c(1, 2), size = 2, replace = F), ]

  founder_genos_ID <- rep(c(RV_founder, founder_ids), each = 2)

  return(list(founder_genos, founder_genos_ID))
}

#' Simulate sequence data for a study
#'
#' @inheritParams sim_RVseq
#' @param ped_files Data frame. Must match format of pedigree simulated by sim_RVped
#' @param marker_map Data.frame. Must contain three columns with: column 1: marker names, must be listed in the same order as in the founder genotype file, column 2: the chromosomal position of the marker, column 3: the position of the marker in cM.
#' @param haplotype_dist. List of dataframes.  Each list entry should contain a dataframe containing the haplotype distribution for that chromosome (in same order as marker_map).
#' @param founder_genotypes list of dataframes. Should contain, in the order listed in chrom_map, 1 dataframe for every chromosome contained in chrom_map.  Each data frame should contain two rows for each founder: the first to correspond to the paternally inherited gamete and the second to correspond to the maternally inherited gamete.
#' @param RV_markers character. A list of possible RV markers. If missing all markers in \code{marker_map} are assumed to be possible RV markers.
#'
#' @return study_sequences
#' @export
#'
#' @examples
#' library(SimRVPedigree)
#' data(EgPeds)
#'
#' ex_study_peds <- EgPeds[EgPeds$FamID == 1, ]
#'
#' data(hg_chrom)
#' my_chrom_map = hg_chrom
#' my_chrom_map
#'
#' head(mark_map)
#' mm_obj <- markerMap(mark_map)
#' head(mm_obj)
#'
#' h_dist <- list()
#' h_dist[[1]]<- SNP_dat
#'
#' set.seed(6)
#' ped_seq <- sim_RVstudy(ped_files = ex_study_peds,
#'                        haplotype_dist = h_dist,
#'                        marker_map = mm_obj,
#'                        chrom_map = my_chrom_map)
#' ped_seq
#'
sim_RVstudy <- function(ped_files, marker_map, chrom_map,
                        haplotype_dist,
                        affected_only = TRUE,
                        convert_to_cM = TRUE,
                        burn_in = 1000,
                        gamma_params = c(2.63, 2.63/0.5)){

  if(!is.markerMap(marker_map)) stop("Expecting class(marker_map) to include markerMap")

  #convert from base pairs to centiMorgan
  if (convert_to_cM) {
    options(digits = 9)
    chrom_map$start_pos <- convert_BP_to_cM(chrom_map$start_pos)
    chrom_map$end_pos <- convert_BP_to_cM(chrom_map$end_pos)
    chrom_map$center <- convert_BP_to_cM(chrom_map$center)

    marker_map$position <- convert_BP_to_cM(marker_map$position)
  }

  FamIDs <- unique(ped_files$FamID)
  #reduce to affected only pedigrees unless otherwise specified
  if (affected_only) {
    Afams <- lapply(c(1:length(FamIDs)), function(x){
      affected_onlyPed(ped_file = ped_files[which(ped_files$FamID == FamIDs[x]),])
      })

    ped_files <- do.call("rbind", Afams)
  }

  # WILL PROBABLY WANT TO CHANGE PROCESS HERE
  # TO CHOOSE PATHWAY THEN RARE VARIANT


  #sampling from RV markers (with probability
  #probCausal)to determine familial RV locus
  Fam_RVs <- sample(x = marker_map$marker,
                    prob = marker_map$probCausal,
                    size = length(FamIDs),
                    replace = TRUE)

  #Given the location of familial risk variants,
  #sample familial founder haplotypes from
  #conditional haplotype distribution
  f_genos <- lapply(c(1:length(FamIDs)), function(x){
    sim_FGenos(founder_ids = ped_files$ID[which(ped_files$FamID == FamIDs[x]
                                                & is.na(ped_files$dadID)
                                                & (ped_files$DA1 + ped_files$DA2) == 0)],
               RV_founder = ped_files$ID[which(ped_files$FamID == FamIDs[x]
                                               & is.na(ped_files$dadID)
                                               & (ped_files$DA1 + ped_files$DA2) == 1)],
               FamID = FamIDs[x], haplotype_dist, FamRV = Fam_RVs[x], marker_map)
  })

  #simulate non-founder haploypes via conditional gene drop
  ped_seqs <- lapply(c(1:length(FamIDs)), function(x){
    sim_RVseq(ped_file = ped_files[which(ped_files$FamID == FamIDs[x]), ],
              founder_genos = f_genos[[x]],
              marker_map, chrom_map,
              RV_marker = Fam_RVs[x],
              burn_in, gamma_params)
    })

  ped_genos <- do.call("rbind", lapply(ped_seqs, function(x){x$ped_genos}))
  geno_map <- do.call("rbind", lapply(ped_seqs, function(x){x$geno_map}))


  # #return data for affecteds only
  # if (affected_only) {
  #   study_sequenceDat <- study_sequenceDat[which(study_sequenceDat$affected == 1), ]
  #   #study_sequenceDat <- study_sequenceDat[, -which(colnames(study_sequenceDat) == "affected")]
  # }

  return(list(ped_genos = ped_genos, geno_map = geno_map))
}
