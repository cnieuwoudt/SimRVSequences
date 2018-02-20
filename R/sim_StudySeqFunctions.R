#' Condition Haplotype distribution
#'
#' @param Chaplo_dist haplotype distribution for a single chromosome
#' @param RV_marker familial RV marker
#' @param RV_status 0 or 1, 1 if RV is inherited
#'
#' @return The conditioned haplotype distribution.
#' @export
#'
condition_haploDist <- function(Chaplo_dist, RV_marker, RV_status){
  hap_prob <- rep(1/nrow(Chaplo_dist), nrow(Chaplo_dist))

  #determine which rows of Chaplo_dist are appropriate given rv status
  keep_rows <- which(Chaplo_dist[, which(colnames(Chaplo_dist) == RV_marker)] == RV_status)
  if(length(keep_rows) == 1){
    cond_haploDist <- as.data.frame(t(as.matrix(Chaplo_dist[keep_rows, ])))
    cond_haploDist$prob <- hap_prob[keep_rows]/sum(hap_prob[keep_rows])
  } else {
    cond_haploDist <- as.data.frame(Chaplo_dist[keep_rows, ])
    cond_haploDist$prob <- hap_prob[keep_rows]/sum(hap_prob[keep_rows])
  }

  return(cond_haploDist)
}


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
      #Since the risk variant is located on this chromosome, we must condition the
      #haplotype distribution on RV status before drawing haplotype
      RVdist <- condition_haploDist(haplotype_dist[[k]],
                                    RV_marker = FamRV,
                                    RV_status = 1)

      RV_founder_dat = RVdist[sample(x = c(1:nrow(RVdist)),
                                     size = 1,
                                     prob = RVdist$prob),
                              -ncol(RVdist)]

      NRVdist <- condition_haploDist(haplotype_dist[[k]],
                                     RV_marker = FamRV,
                                     RV_status = 0)

      NRV_founder_dat <- NRVdist[sample(x = c(1:nrow(NRVdist)),
                                        size = (2*length(founder_ids) + 1),
                                        replace = TRUE,
                                        prob = NRVdist$prob),
                                 -ncol(NRVdist)]

      fam_genos[[k]] <- rbind(RV_founder_dat, NRV_founder_dat)
    } else {

      #since the familial RV does not reside on the kth chromosome we
      #sample founder haplotypes according to from the haplotype distribution

      fam_genos[[k]] = haplotype_dist[[k]][sample(x = c(1:nrow(haplotype_dist[[k]])),
                                                  size = 2*(length(founder_ids) + 1),
                                                  replace = TRUE,
                                                  prob = haplotype_dist[[k]]$prob),
                                           -ncol(haplotype_dist[[k]])]
    }
  }

  founder_genos <- do.call("cbind", fam_genos)
  founder_genos$ID <- rep(c(RV_founder, founder_ids), each = 2)
  rownames(founder_genos) <- NULL

  #ramdomly permute RV founders rows so that the RV is not always paternally inherited.
  founder_genos[c(1,2), ] <- founder_genos[sample(x = c(1, 2), size = 2, replace = F), ]

  return(founder_genos)
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
#' data(hg_chrom)
#' data(mark_map)
#' data(SNP_dat)
#'
#' library(SimRVPedigree)
#' data(EgPeds)
#'
#' ex_study_peds <- EgPeds
#'
#' my_chrom_map = hg_chrom[17, ]
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

  #simulate non-founder haploypes from
  #founder haplotypes
  ped_seqs <- lapply(c(1:length(FamIDs)), function(x){
    sim_RVseq(ped_file = ped_files[which(ped_files$FamID == FamIDs[x]), ],
              founder_genos = f_genos[[x]],
              marker_map, chrom_map,
              RV_marker = Fam_RVs[x],
              burn_in, gamma_params)
    })

  study_sequenceDat <- do.call("rbind", ped_seqs)

  if (affected_only) {
    study_sequenceDat <- study_sequenceDat[which(study_sequenceDat$affected == 1), ]
    #study_sequenceDat <- study_sequenceDat[, -which(colnames(study_sequenceDat) == "affected")]
  }
  rownames(study_sequenceDat) = NULL

  return(study_sequenceDat)
}
