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
#' library(kinship2)
#' #Read in age-specific hazards
#' data(EgPeds)
#'
#' EgPeds5 <- EgPeds4 <- EgPeds3 <- EgPeds2 <- EgPeds
#' EgPeds2$FamID <- EgPeds2$FamID + 5
#' EgPeds3$FamID <- EgPeds3$FamID + 10
#' EgPeds4$FamID <- EgPeds4$FamID + 15
#' EgPeds5$FamID <- EgPeds5$FamID + 20
#' ex_study_peds <- rbind(EgPeds, EgPeds2, EgPeds3, EgPeds4, EgPeds5)
#' ex_study_peds <- EgPeds
#'
#' nrow(ex_study_peds[which(is.na(ex_study_peds$dadID)), ])
#' length(unique(ex_study_peds$FamID))
#'
#' my_chrom_map = data.frame(chrom     = c(1, 2),
#'                           start_pos = c(0, 0),
#'                           end_pos   = c(247199719, 242751149),
#'                           center = c(98879888, 97100460))
#' my_chrom_map
#'
#' mark_map <- data.frame(chrom = c(1, 1, 1, 1, 1, 1, 1, 1,
#'                                       2, 2, 2, 2, 2, 2),
#'                        position = c(5012368, 5012369, 5012370,
#'                                     78541008, 78541009, 78541010,
#'                                     247199219, 247199220,
#'                                     11330, 11332,
#'                                     234577, 234578, 234579,
#'                                     18799180),
#'                        pathwayID = c(1, 1, 1, 2, 2, 2, 3, 3,
#'                                      2, 2, 3, 3, 3, 1),
#'                        possibleRV = c(0, 1, 1, 0, 0, 0, 1, 1,
#'                                       0, 0, 0, 1, 1, 1))
#' mark_map$marker <- paste0(mark_map$chrom, sep = "_", mark_map$position)
#' mark_map <- mark_map[, c(5, 1:4)]
#' mark_map
#'
#' set.seed(1)
#' founder_seq <- as.data.frame(matrix(sample(2*nrow(mark_map)*1000,
#'                                     x = c(0, 1),
#'                                     replace = T,
#'                                     prob = c(0.95, 0.05)),
#'                              nrow = 2*1000))
#' colnames(founder_seq) = as.character(mark_map$marker)
#' hdist = estimate_haploDist(founder_seq, mark_map)
#'
#' hdist2 = estimate_haploDist(founder_seq[which(founder_seq[, 2] == 1), ], mark_map)
#' hdist2[[1]]
#'
#' condition_haploDist(hdist[[1]], "1_5012369", 1)
#'
#' set.seed(6)
#' ped_seq <- sim_RVstudy(ped_files = ex_study_peds,
#'                        haplotype_dist = estimate_haploDist(founder_seq, mark_map),
#'                        marker_map = mark_map,
#'                        chrom_map = my_chrom_map)
#' ped_seq
#'
#' set.seed(6)
#' system.time(sim_RVstudy(ped_files = ex_study_peds,
#'                         haplotype_dist = estimate_haploDist(founder_seq, mark_map),
#'                         marker_map = mark_map,
#'                         chrom_map = my_chrom_map))
#'
sim_RVstudy <- function(ped_files, marker_map, chrom_map,
                        haplotype_dist,
                        affected_only = TRUE,
                        convert_to_cM = TRUE,
                        burn_in = 1000,
                        gamma_params = c(2.63, 2.63/0.5)){

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

  # WILL PROBABLY WANT TO CHANGE CODE HERE
  # TO CHOOSE PATHWAY THEN RARE VARIANT

  if (is.null(marker_map$probCausal)) {
    marker_map$probCausal <- marker_map$possibleRV/sum(marker_map$possibleRV)
  }

  #determine familial RV locus
  Fam_RVs <- sample(x = marker_map$marker, prob = marker_map$probCausal,
                    size = length(FamIDs), replace = TRUE)

  #simulate founder haplotypes from haplotype distribution given the location
  #of the familial risk variant
  f_genos <- lapply(c(1:length(FamIDs)), function(x){
    sim_FGenos(founder_ids = ped_files$ID[which(ped_files$FamID == FamIDs[x]
                                                & is.na(ped_files$dadID)
                                                & (ped_files$DA1 + ped_files$DA2) == 0)],
               RV_founder = ped_files$ID[which(ped_files$FamID == FamIDs[x]
                                               & is.na(ped_files$dadID)
                                               & (ped_files$DA1 + ped_files$DA2) == 1)],
               FamID = FamIDs[x], haplotype_dist, FamRV = Fam_RVs[x], marker_map)
  })

  #simulate non-founder haploypes from haplotype distribution given location
  #of the familial risk variant
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
