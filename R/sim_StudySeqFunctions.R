#' Draw Founder Genotypes from Haplotype Distribution Given Familial Risk Variant
#'
#' \strong{For internal use.}
#'
#' @param founder_ids Numeric list. The ID numbers of all non-seed founders.
#' @param RV_founder Numeric. The ID number of the seed founder.
#' @param RV_founder_pat Numeric. RV_founder_pat == 1 if RV founder inherited the RV from dad, and 0 if inherited RV from mom.
#' @param haplos sparseMatrix.  The sparseMatrix of genomes returned by \code{read_slimOut}.
#' @param RV_col_loc Nueric. The column location of the familial RV in haplos.
#'
#' @return list of familial founder genotypes
#' @export
#'
sim_FGenos <- function(founder_ids, RV_founder, RV_founder_pat,
                       haplos, RV_col_loc) {

  #Determine which haplotypes carry the familial rare variant and which so not
  RV_haps <- which(haplos[, RV_col_loc] == 1)
  noRV_haps <- which(haplos[, RV_col_loc] == 0)

  #for the seed founder: sample one haplotype from those that carry the RV
  # and one haplotype from those that DO NOT carry the RV
  #for all other founders: sample 2 haplotypes that do not carry the RV
  if(length(RV_haps) == 1){
    founder_genos = haplos[c(RV_haps,
                                     sample(x = noRV_haps,
                                            size = (2*length(founder_ids) + 1),
                                            replace = TRUE)), ]
  } else {
    founder_genos = haplos[c(sample(x = RV_haps, size = 1),
                                     sample(x = noRV_haps,
                                            size = (2*length(founder_ids) + 1),
                                            replace = TRUE)), ]
  }



  #Asscociate RV to row 1, if paternally inherited OR row 2 if maternally inherited.
  if(RV_founder_pat == 0){
    founder_genos[c(1,2), ] <- founder_genos[c(2, 1), ]
  }


  #create IDs to associate founders to rows in founder_genos
  founder_genos_ID <- rep(c(RV_founder, founder_ids), each = 2)

  return(list(founder_genos, founder_genos_ID))
}

#' Remove markers at which all individuals carry wild type allele
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
#' @inheritParams sim_RVseq
#' @param ped_files Data frame. Must match format of pedigree simulated by sim_RVped
#' @param SNV_map Data.frame. Must contain three columns with: column 1: marker names, must be listed in the same order as in the founder genotype file, column 2: the chromosomal position of the marker, column 3: the position of the marker in cM.
#' @param haplos sparseMatrix. The genomes matrix returned by \code{read_slim}.  Mutations in haplos are described in \code{SNV_map}.
#' @param affected_only Logical. When \code{affected_only = TRUE} pedigrees are reduced contain only diesease-affected relatives, their parents, and any obligate carriers or founders; sequence data is simulated only for retained family memebres. When \code{affected_only = FALSE} sequence data is simulated for all.  By default, \code{affected_only = TRUE}.
#' @param convert_to_cM Logical. The setting \code{convert_to_cM = TRUE} must be used if genomic positions are given in base pairs.  If the genomic postions in \code{SNV_map} are listed in centiMorgan, please set \code{convert_to_cM = FALSE}.  By default, \code{convert_to_cM = TRUE}.
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
#' data(EXmut)
#' data(EXgen)
#'
#' markers = EXmut
#' markers$possibleSNV = FALSE
#' markers$possibleSNV[1] = TRUE
#'
#' seqDat = sim_RVstudy(ped_files = EgPeds,
#'                      SNV_map = markers,
#'                      haplos = EXgen)
#'
#' summary(seqDat)
#'
sim_RVstudy <- function(ped_files, SNV_map,
                        haplos,
                        affected_only = TRUE,
                        convert_to_cM = TRUE,
                        remove_wild = TRUE,
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
  Fam_RVs <- sample(x = SNV_map$marker[SNV_map$possibleSNV],
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
               haplos, RV_col_loc = which(SNV_map$marker == Fam_RVs[x]))
  })

  #If desired by user, we now reduce the size of the data by removing
  #markers at which no member of our study carries a mutated allele.
  if (remove_wild) {
    reduced_dat <- remove_allWild(f_haps = f_genos, SNV_map)
    f_genos <- reduced_dat[[1]]
    SNV_map <- reduced_dat[[2]]

  }

  #create chrom_map, this is used to determine where we need to
  #simulate recombination
  chrom_map <- create_chrom_map(SNV_map)

  #convert from base pairs to centiMorgan
  if (convert_to_cM) {
    options(digits = 9)
    chrom_map$start_pos <- convert_BP_to_cM(chrom_map$start_pos)
    chrom_map$end_pos <- convert_BP_to_cM(chrom_map$end_pos)

    SNV_map$position <- convert_BP_to_cM(SNV_map$position)
  }

  #simulate non-founder haploypes via conditional gene drop
  ped_seqs <- lapply(c(1:length(FamIDs)), function(x){
    sim_RVseq(ped_file = ped_files[ped_files$FamID == FamIDs[x], ],
              founder_genos = f_genos[[x]],
              SNV_map, chrom_map,
              RV_marker = Fam_RVs[x],
              burn_in, gamma_params)
    })

  ped_genos <- do.call("rbind", lapply(ped_seqs, function(x){x$ped_genos}))
  geno_map <- do.call("rbind", lapply(ped_seqs, function(x){x$geno_map}))

  return(list(ped_genos = ped_genos, geno_map = geno_map, SNV_map = SNV_map))
}
