#' Simulate sequence data for a study
#'
#' @inheritParams sim_RVseq
#' @param ped_files Data frame. Must match format of pedigree simulated by sim_RVped
#' @param founder_genotypes list of dataframes. Should contain, in the order listed in chrom_map, 1 dataframe for every chromosome contained in chrom_map.  Each data frame should contain two rows for each founder: the first to correspond to the paternally inherited gamete and the second to correspond to the maternally inherited gamete.
#' @param linkage_map Data.frame. Must contain three columns with: column 1: marker names, must be listed in the same order as in the founder genotype file, column 2: the chromosomal position of the marker, column 3: the position of the marker in cM.
#' @param RV_markers character. A list of possible RV markers. If missing all markers in \code{linkage_map} are assumed to be possible RV markers.
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
#' my_study_peds = EgPeds
#' my_chrom_map = data.frame(chrom     = c(1, 2),
#'                           start_pos = c(0, 0),
#'                           end_pos   = c(270, 200),
#'                           center = c(55, 40))
#' my_chrom_map
#'
#' link_map <- data.frame(chrom = c(1, 1, 1, 1, 1, 1, 1,
#'                                       2, 2, 2, 2, 2, 2),
#'                        position = c(20, 50, 100, 125, 175, 200, 250,
#'                                      10, 25,  75, 125, 150, 200))
#' link_map$marker <- paste0(link_map$chrom, sep = "_", link_map$position)
#' link_map <- link_map[, c(3, 1, 2)]
#' link_map
#'
#' my_RV_markers <- c('1_50', '2_75', '2_150')
#' my_RV_markers <- c('1_50')
#'
#'
#' set.seed(1)
#' founder_seq <- as.data.frame(matrix(sample(2*nrow(link_map)*1000,
#'                                     x = c(0, 1),
#'                                     replace = T,
#'                                     prob = c(0.95, 0.05)),
#'                              nrow = 2*1000))
#' colnames(founder_seq) = as.character(link_map$marker)
#' length(which(founder_seq$`1_50` == 1))
#'
#'
#' set.seed(6)
#' ped_seq <- sim_RVstudy(ped_files = my_study_peds,
#'                        founder_genotypes = founder_seq,
#'                        linkage_map = link_map,
#'                        chrom_map = my_chrom_map,
#'                        RV_markers = my_RV_markers)
#' ped_seq
#'
#' set.seed(6)
#' system.time(sim_RVstudy(ped_files = my_study_peds,
#'                         founder_genotypes = founder_seq,
#'                         linkage_map = link_map,
#'                         chrom_map = my_chrom_map,
#'                         RV_marker = my_RV_marker))
#'
sim_RVstudy <- function(ped_files, founder_genotypes,
                        linkage_map, chrom_map, RV_markers,
                        RV_marker_probs,
                        affected_only = TRUE,
                        burn_in = 1000, gamma_params = c(2.63, 2.63/0.5)){

  #initialize founder_genotypes ID and FamID Variables
  founder_genotypes$ID <- NA

  FamIDs <- unique(ped_files$FamID)
  #reduce to affected only pedigrees unless otherwise specified
  if (affected_only) {
    Afams <- lapply(c(1:length(FamIDs)), function(x){
      affected_onlyPed(ped_file = ped_files[which(ped_files$FamID == FamIDs[x]),])
      })

    ped_files <- do.call("rbind", Afams)
  }

  if (missing(RV_marker_probs) & !missing(RV_markers)) {
    RV_marker_probs <- rep(1/length(RV_markers), length(RV_markers))
  } else if (missing(RV_marker_probs) & missing(RV_markers)) {
    RV_markers <- linkage_map$marker
    RV_marker_probs <- rep(1/length(RV_markers), length(RV_markers))
  }

  #determine familial RV locus
  Fam_RVs <- sample(x = RV_markers, prob = RV_marker_probs,
                    size = length(FamIDs), replace = TRUE)

  f_genos <- list()
  open_genos <- founder_genotypes
  for(i in 1:length(Fam_RVs)){
    loop_genos <- get_FGenos(founder_ids = ped_files$ID[which(ped_files$FamID == FamIDs[i]
                                                & is.na(ped_files$dad_id)
                                                & (ped_files$DA1 + ped_files$DA2) == 0)],
                             RV_founder = ped_files$ID[which(ped_files$FamID == FamIDs[i]
                                                             & is.na(ped_files$dad_id)
                                                             & (ped_files$DA1 + ped_files$DA2) == 1)],
                             FamID = FamIDs[i], founder_genotypes = open_genos, FamRV = Fam_RVs[i])
    f_genos[[i]] <- loop_genos[[1]]
    open_genos <- loop_genos[[2]]
  }

  # length(unique(rownames(do.call("rbind", f_genos))))
  # length(rownames(do.call("rbind", f_genos)))

  ped_seqs <- lapply(c(1:length(FamIDs)), function(x){
    sim_RVseq(ped_file = ped_files[which(ped_files$FamID == FamIDs[x]), ],
              founder_genotypes = f_genos[[x]],
              linkage_map, chrom_map,
              RV_marker = Fam_RVs[x],
              burn_in, gamma_params)
  })

  study_sequenceDat <- do.call("rbind", ped_seqs)

  return(study_sequenceDat)
}
