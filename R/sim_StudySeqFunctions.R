#' Reduce haplos to contain non-cSNV data
#'
#' @inheritParams sim_FGenos
#'
#' @return The reduced haplotype matrix.
#'
#' @export
condition_haplos_no_cSNV <- function(haplos, RV_pool_loc){
  #determine which of the haplotypes carry a causal SNV
  RV_haps <- lapply(RV_pool_loc, function(x){
    which(haplos[, x] == 1)
  })

  RV_rows <- Reduce(union, RV_haps)

  return(haplos[-RV_rows, ])
}

#' Draw Founder Genotypes from Haplotype Distribution Given Familial Risk Variant
#'
#' \strong{For internal use.}
#'
#' @param founder_ids Numeric list. The ID numbers of all non-seed founders.
#' @param RV_founder Numeric. The ID number of the seed founder.
#' @param RV_founder_pat Numeric. RV_founder_pat == 1 if RV founder inherited the RV from dad, and 0 if inherited RV from mom.
#' @param haplos sparseMatrix.  The sparseMatrix of genomes returned by \code{read_slimOut}.
#' @param RV_col_loc Nueric. The column location of the familial RV in haplos.
#' @param RV_pool_loc The column locations of each SNV in the pool of candidate SNVs.
#'
#' @return list of familial founder genotypes
#' @export
#'
sim_FGenos <- function(founder_ids, RV_founder, RV_founder_pat,
                       haplos, RV_col_loc, RV_pool_loc) {

  #here we handle the fully sporatic families
  #i.e. families that do not segregate any cSNVs
  if(is.null(RV_founder)){

    #reduce haplos to contain only haplotypes that do
    #not carry any of the cRVs in our pool of possible SNVs
    noRV_haps <- condition_haplos_no_cSNV(haplos, RV_pool_loc)

    #sample all founder data from this pool
    founder_genos <- noRV_haps[c(sample(x = 1:nrow(noRV_haps),
                                       size = (2*length(founder_ids) + 2),
                                       replace = TRUE)), ]

    #create IDs to associate founders to rows in founder_genos
    founder_genos_ID <- rep(c(RV_founder, founder_ids), each = 2)

  } else {

    #Determine which haplotypes carry the familial rare variant and which so not
    RV_hap_loc <- which(haplos[, RV_col_loc] == 1)
    noRV_hap_loc <- which(haplos[, RV_col_loc] == 0)

    #NOTE: Under this scheme, marry-ins may NOT introduce any SNV
    #from our pool of causal rare variants

    #NOTE: Under this scheme, marry-ins may introduce SNVs from our pool of causal
    #rare variants, just not the one that is the familial cRV


    #for the seed founder: sample one haplotype from those that carry the RV
    #and one haplotype from those that DO NOT carry the RV
    #for all other founders: sample 2 haplotypes that do not carry the RV
    if(length(RV_hap_loc) == 1){
      founder_genos <- haplos[c(RV_hap_loc,
                               sample(x = noRV_hap_loc,
                                      size = (2*length(founder_ids) + 1),
                                      replace = TRUE)), ]
    } else {
      founder_genos <- haplos[c(sample(x = RV_hap_loc, size = 1),
                               sample(x = noRV_hap_loc,
                                      size = (2*length(founder_ids) + 1),
                                      replace = TRUE)), ]
    }



    #Asscociate RV to row 1, if paternally inherited OR row 2 if maternally inherited.
    if(RV_founder_pat == 0){
      founder_genos[c(1,2), ] <- founder_genos[c(2, 1), ]
    }


    #create IDs to associate founders to rows in founder_genos
    founder_genos_ID <- rep(c(RV_founder, founder_ids), each = 2)
  }

  return(list(founder_genos, founder_genos_ID))
}

#' Remove unmutated markers from data
#'
#' Remove any markers for which all founders, in the study, are homozygous for the wild-type allele.  Since we do not model de novo mutations, it is not possible for non-founders to develop mutations at these loci.
#'
#' @param f_haps The founder haplotypes data. This is a list of lists (by family). By family, this contains the haplotypes for each founder (first item), and a list of ID numbers (second item) which is used to map the haplotype to the person to whom it belongs.
#' @param SNV_map data.frame. Catalogs the SNV data contained in the familial haplotypes.
#'
#' @return A list (by family) of haplotype matrices and ID vectors and the reduce marker data set.
#' @importFrom Matrix colSums
#' @export
#'
#' @examples
#' #probably no examples
#'
remove_allWild <- function(f_haps, SNV_map){

  #determine which columns are all zero in founder haplotype data.
  #These are markers at which no one in ped_files will carry a SNV
  #(i.e. all family members will carry the wild type allele) since no
  #founders have mutated alleles to pass on.
  #We will reduce the size of the data by removing these superfluous

  #Determine all wild columns by family
  wild_col <- lapply(f_haps, function(x){
    which(colSums(x[[1]]) == 0)
  })

  #determine all wild columns overall - by study
  remove_cols <- Reduce(intersect, wild_col)

  #remove all wild columns from founder haplotypes
  red_haps <- lapply(f_haps, function(x){
    list(x[[1]][, -remove_cols],
         x[[2]])
  })

  #remove all wild columns from SNV_map
  #RECALL: rows of SNV_map are columns in haplotype data
  SNV_map <- SNV_map[-remove_cols, ]
  SNV_map$colID <- seq(1:nrow(SNV_map))
  row.names(SNV_map) = NULL

  return(list(red_haps, SNV_map))
}

#' Simulate sequence data for a study
#'
#' @inheritParams sim_seq
#' @param ped_files Data frame. Must match format of pedigree simulated by sim_RVped
#' @param SNV_map Data.frame. Must contain three columns with: column 1: marker names, must be listed in the same order as in the founder genotype file, column 2: the chromosomal position of the marker, column 3: the position of the marker in cM.
#' @param haplos sparseMatrix. The genomes matrix returned by \code{read_slim}.  Mutations in haplos are described in \code{SNV_map}.
#' @param affected_only Logical. When \code{affected_only = TRUE} pedigrees are reduced contain only diesease-affected relatives, their parents, and any obligate carriers or founders; sequence data is simulated only for retained family memebres. When \code{affected_only = FALSE} sequence data is simulated for all.  By default, \code{affected_only = TRUE}.
#' @param pos_in_bp Logical. The setting \code{pos_in_bp = TRUE} must be used if genomic positions are given in base pairs.  If the genomic postions in \code{SNV_map} are listed in centiMorgan, please set \code{pos_in_bp = FALSE}.  By default, \code{pos_in_bp = TRUE}.
#' @param remove_wild Logical. Should markers at which no member of study carry a mutated allele be removed from the data. By default, \code{remove_wild = TRUE}.
#'
#' @return study_sequences
#' @export
#'
#' @examples
#' library(SimRVPedigree)
#' data(EgPeds)
#'
#' library(SimRVSequences)
#' data(EXmuts)
#' data(EXhaps)
#'
#' markers = EXmuts
#' markers$is_CRV = FALSE
#' markers$is_CRV[1] = TRUE
#'
#' seqDat = sim_RVstudy(ped_files = EgPeds,
#'                      SNV_map = markers,
#'                      haplos = EXhaps)
#'
#' summary(seqDat)
#'
sim_RVstudy <- function(ped_files, SNV_map, haplos,
                        affected_only = TRUE,
                        remove_wild = TRUE,
                        pos_in_bp = TRUE,
                        burn_in = 1000,
                        gamma_params = c(2.63, 2.63/0.5)){

  #check SNV_map for possible issues
  check_SNV_map(SNV_map)

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
  #NOTE: IF POSSIBLE SNV NOT SPECIFIED A SINGLE SNV IS CHOSEN AS
  #      CAUSAL FOR ALL FAMILIES IN STUDY.  DOES NOT CONSIDER ALLELE FREQUENCY.
  if (is.null(SNV_map$is_CRV)) {
    SNV_map$is_CRV = FALSE
    SNV_map$is_CRV[sample(1:nrow(SNV_map), size = 1)] = TRUE
  }

  #--------------------#
  #      FIX THIS      #
  #--------------------#
  #NOTE WILL BREAK IF ONLY ONE CRV IN POOL#
  #sample the familial cSV from the pool of potential cRVs with replacement.
  Fam_RVs <- sample(x = SNV_map$marker[SNV_map$is_CRV],
                    size = length(FamIDs),
                    replace = TRUE)

  #Given the location of familial risk variants, sample familial founder
  #haplotypes from conditional haplotype distribution
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
               haplos, RV_col_loc = which(SNV_map$marker == Fam_RVs[x]),
               RV_pool_loc = SNV_map$colID[SNV_map$is_CRV])
  })

  #If desired by user, we now reduce the size of the data by removing
  #markers at which no member of our study carries a mutated allele.
  if (remove_wild) {
    reduced_dat <- remove_allWild(f_haps = f_genos, SNV_map)
    f_genos <- reduced_dat[[1]]
    SNV_map <- reduced_dat[[2]]
  }

  #create chrom_map, this is used to determine the segments over
  #which we will simulate genetic recombination
  chrom_map <- create_chrom_map(SNV_map)

  #convert from base pairs to centiMorgan
  if (pos_in_bp) {
    options(digits = 9)
    chrom_map$start_pos <- convert_BP_to_cM(chrom_map$start_pos)
    chrom_map$end_pos <- convert_BP_to_cM(chrom_map$end_pos)

    SNV_map$position <- convert_BP_to_cM(SNV_map$position)
  }

  #simulate non-founder haploypes via conditional gene drop
  ped_seqs <- lapply(c(1:length(FamIDs)), function(x){
    sim_seq(ped_file = ped_files[ped_files$FamID == FamIDs[x], ],
              founder_genos = f_genos[[x]],
              SNV_map, chrom_map,
              RV_marker = Fam_RVs[x],
              burn_in, gamma_params)
    })

  ped_haplos <- do.call("rbind", lapply(ped_seqs, function(x){x$ped_genos}))
  haplo_map <- do.call("rbind", lapply(ped_seqs, function(x){x$geno_map}))

  #convert back to base pairs if we converted to CM
  if (pos_in_bp) {
    options(digits = 9)
    SNV_map$position <- convert_CM_to_BP(SNV_map$position)
  }

  return(list(ped_haplos = ped_haplos, haplo_map = haplo_map, SNV_map = SNV_map, ped_files = ped_files))
}
