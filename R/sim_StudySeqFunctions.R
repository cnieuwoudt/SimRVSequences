#' Draw Founder Genotypes from Haplotype Distribution Given Familial Risk Variant
#'
#' \strong{For internal use.}
#'
#' @param founder_ids Numeric list. The ID numbers of all non-seed founders.
#' @param RV_founder Numeric. The ID number of the seed founder.
#' @param RV_founder_pat Numeric. RV_founder_pat == 1 if RV founder inherited the RV from dad, and 0 if inherited RV from mom.
#' @param haplotype_dist sparseMatrix.  The sparseMatrix of genomes returned by \code{read_slimOut}.
#' @param RV_col_loc Nueric. The column location of the familial RV in haplotype_dist.
#'
#' @return list of familial founder genotypes
#' @export
#'
sim_FGenos <- function(founder_ids, RV_founder, RV_founder_pat,
                       haplotype_dist, RV_col_loc) {

  #Determine which haplotypes carry the familial rare variant and which so not
  RV_haps <- which(haplotype_dist[, RV_col_loc] == 1)
  noRV_haps <- which(haplotype_dist[, RV_col_loc] == 0)

  #for the seed founder: sample one haplotype from those that carry the RV
  # and one haplotype from those that DO NOT carry the RV
  #for all other founders: sample 2 haplotypes that do not carry the RV
  founder_genos = haplotype_dist[c(sample(x = RV_haps, size = 1),
                                   sample(x = noRV_haps,
                                          size = (2*length(founder_ids) + 1),
                                          replace = TRUE)), ]


  #Asscociate RV to row 1, if paternally inherited OR row 2 if maternally inherited.
  if(RV_founder_pat == 0){
    founder_genos[c(1,2), ] <- founder_genos[c(2, 1), ]
  }


  #create IDs to associate founders to rows in founder_genos
  founder_genos_ID <- rep(c(RV_founder, founder_ids), each = 2)

  return(list(founder_genos, founder_genos_ID))
}

#' Simulate sequence data for a study
#'
#' @inheritParams sim_RVseq
#' @param ped_files Data frame. Must match format of pedigree simulated by sim_RVped
#' @param marker_map Data.frame. Must contain three columns with: column 1: marker names, must be listed in the same order as in the founder genotype file, column 2: the chromosomal position of the marker, column 3: the position of the marker in cM.
#' @param haplotype_dist sparseMatrix. The genomes matrix returned by \code{read_slim}.  Mutations in haplotype_dist are described in \code{marker_map}.
#' @param affected_only Logical. When \code{affected_only = TRUE} pedigrees are reduced contain only diesease-affected relatives, their parents, and any obligate carriers or founders; sequence data is simulated only for retained family memebres. When \code{affected_only = FALSE} sequence data is simulated for all.  By default, \code{affected_only = TRUE}.
#' @param convert_to_cM Logical. The setting \code{convert_to_cM = TRUE} must be used if genomic positions are given in base pairs.  If the genomic postions in \code{marker_map} are listed in centiMorgan, please set \code{convert_to_cM = FALSE}.  By default, \code{convert_to_cM = TRUE}.
#'
#' @return study_sequences
#' @export
#'
#' @examples
#'
#' #FIND SHORT WORKING EXAMPLE
#'
sim_RVstudy <- function(ped_files, marker_map, chrom_map,
                        haplotype_dist,
                        affected_only = TRUE,
                        convert_to_cM = TRUE,
                        burn_in = 1000,
                        gamma_params = c(2.63, 2.63/0.5)){

  #if(!is.markerMap(marker_map)) stop("Expecting class(marker_map) to include markerMap")

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
    Afams <- lapply(FamIDs, function(x){
      affected_onlyPed(ped_file = ped_files[which(ped_files$FamID == x),])
      })

    ped_files <- do.call("rbind", Afams)
  }

  #sampling from RV markers
  #to determine familial RV locus
  Fam_RVs <- sample(x = marker_map$marker[marker_map$possibleRV],
                    #prob = marker_map$probCausal,
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
               RV_founder_pat = ped_files$DA1[which(ped_files$FamID == FamIDs[x]
                                                   & is.na(ped_files$dadID)
                                                   & (ped_files$DA1 + ped_files$DA2) == 1)],
               haplotype_dist, RV_col_loc = which(marker_map$marker == Fam_RVs[x]))
  })

  #simulate non-founder haploypes via conditional gene drop
  ped_seqs <- lapply(c(1:length(FamIDs)), function(x){
    sim_RVseq(ped_file = ped_files[ped_files$FamID == FamIDs[x], ],
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
