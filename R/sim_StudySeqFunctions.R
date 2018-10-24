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
#' @param RV_col_loc Numeric. The column location of the familial RV in haplos.
#' @param RV_pool_loc The column locations of each SNV in the pool of candidate SNVs.
#'
#' @return list of familial founder genotypes
#' @export
#'
sim_FGenos <- function(founder_ids, RV_founder, RV_founder_pat,
                       haplos, RV_col_loc, RV_pool_loc) {

  #reduce haplos to contain only haplotypes that do
  #not carry any of the cRVs in our pool of possible SNVs
  no_CRVhaps <- condition_haplos_no_cSNV(haplos, RV_pool_loc)

  #here we handle the fully sporatic families
  #i.e. families that do not segregate any cSNVs
  #In this case, the haplotypes for ALL founders
  #is sampled from no_CRVhaps
  if(length(RV_founder) == 0){

    #sample all founder data from this pool
    founder_genos <- no_CRVhaps[c(sample(x = 1:nrow(no_CRVhaps),
                                       size = 2*length(founder_ids),
                                       replace = TRUE)), ]
  } else {

    #Determine which haplotypes carry the familial rare variant and which do not
    RV_hap_loc <- which(haplos[, RV_col_loc] == 1)

    #NOTE: Under this scheme, marry-ins may NOT introduce any SNV
    #from our pool of causal rare variants

    #for the seed founder: sample one haplotype from those that carry the RV
    #and one haplotype from those that DO NOT carry the RV
    #for all other founders: sample 2 haplotypes that do not carry the RV
    if(length(RV_hap_loc) == 1){
      founder_genos <- rbind(haplos[RV_hap_loc, ],
                             no_CRVhaps[c(sample(x = 1:nrow(no_CRVhaps),
                                                 size = (2*length(founder_ids) + 1),
                                                 replace = TRUE)), ])
    } else {
      founder_genos <- rbind(haplos[sample(x = RV_hap_loc, size = 1), ],
                             no_CRVhaps[c(sample(x = 1:nrow(no_CRVhaps),
                                                 size = (2*length(founder_ids) + 1),
                                                 replace = TRUE)), ])
    }

    #Asscociate CRV to row 1, if paternally inherited OR row 2 if maternally inherited.
    if(RV_founder_pat == 0){
      founder_genos[c(1,2), ] <- founder_genos[c(2, 1), ]
    }

  }

  #create IDs to associate founders to rows in founder_genos
  founder_genos_ID <- rep(c(RV_founder, founder_ids), each = 2)

  return(list(founder_genos, founder_genos_ID))
}

#' Remove unmutated markers from data
#'
#' Remove any markers for which all founders, in the study, are homozygous for the wild-type allele.  Since we do not model de novo mutations, it is not possible for non-founders to develop mutations at these loci.
#'
#' @param f_haps The founder haplotypes data. This is a list of family lists. By family, this contains the haplotypes for each founder (first item), and a list of ID numbers (second item) which is used to map the haplotype to the person to whom it belongs.
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

  #Determine all-wild columns by family
  wild_col <- lapply(f_haps, function(x){
    which(colSums(x[[1]]) == 0)
  })

  #determine all wild columns overall, i.e. for the study
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
#' Simulate SNV data for a sample of ascertained pedigrees
#'
#' NOTE: Due to a forwards-in-time model, certain type of inbreeding/loop may cause sim_RVstudy to crash.
#'
#' @inheritParams sim_seq
#' @param ped_files Data frame. Must match format of pedigree simulated by sim_RVped
#' @param SNV_map Data.frame. Must contain three columns with: column 1: marker names, must be listed in the same order as in the founder genotype file, column 2: the chromosomal position of the marker, column 3: the position of the marker in cM.
#' @param haplos sparseMatrix. The genomes matrix returned by \code{read_slim}.  Mutations in haplos are described in \code{SNV_map}.
#' @param affected_only Logical. When \code{affected_only = TRUE} pedigrees are reduced contain only diesease-affected relatives, their parents, and any obligate carriers or founders; sequence data is simulated only for retained family memebres. When \code{affected_only = FALSE} sequence data is simulated for all.  By default, \code{affected_only = TRUE}.
#' @param pos_in_bp Logical. The setting \code{pos_in_bp = TRUE} must be used if genomic positions are given in base pairs.  If the genomic postions in \code{SNV_map} are listed in centiMorgan, please set \code{pos_in_bp = FALSE}.  By default, \code{pos_in_bp = TRUE}.
#' @param remove_wild Logical. Should markers at which no member of study carry a mutated allele be removed from the data. By default, \code{remove_wild = TRUE}.
#'
#' @references Roeland E. Voorrips and Chris A Maliepaard. (2012). \emph{The simulation of meiosis in diploid and tetraploid organisms using various genetic models}. BMC Bioinformatics, 13:248.
#' @references Christina Nieuwoudt and Jinko Graham. (??) Future bioRxiv article.
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

  #check ped_files for possible issues
  check_peds(ped_files)

  if (nrow(SNV_map) != ncol(haplos)) {
    stop("\n nrow(SNV_map) != ncol(haplos). \n The rows of SNV_map should describe the columns (i.e. mutations) in haplos.")
  }

  #check to see that the sample contains affected relatives when the
  #affected_only setting is used
  if (affected_only & sum(ped_files$affected) == 0) {
    stop("\n There are no disease-affected relatives in this study. \n To simulate data for pedigrees without disease-affected relatives use affected_only = FALSE.")
  }

  #collect list of FamIDs
  FamIDs <- unique(ped_files$FamID)

  #check for pedigree formatting issues
  for (i in FamIDs){
    check_ped(ped_files[ped_files$FamID == i, ])
  }

  #Reduce to affected-only pedigrees
  if (affected_only) {
    #reduce pedigrees to contain only disease-affected relative and
    #the individuals who connect them along a line of descent.
    Afams <- lapply(FamIDs, function(x){
      affected_onlyPed(ped_file = ped_files[which(ped_files$FamID == x),])
      })

    #combine the reduced pedigrees
    ped_files <- do.call("rbind", Afams)

    #check to see if any pedigrees were removed due to lack of
    #disease affected relatives and issue warning for removed pedigrees
    removed_peds <- setdiff(FamIDs, unique(ped_files$FamID))

    if (length(removed_peds) > 0){
      FamIDs <- unique(ped_files$FamID)
      warning("\n There are no disease-affected relatives in the pedigrees with FamID: ",
              paste0(removed_peds, collapse = ", "),
              "\n These pedigrees have been removed from ped_files.")
    }

  }

  #sampling from RV markers
  #to determine familial RV locus
  #NOTE: IF POSSIBLE SNV NOT SPECIFIED A SINGLE SNV IS CHOSEN AS
  #      CAUSAL FOR ALL FAMILIES IN STUDY.  DOES NOT CONSIDER ALLELE FREQUENCY.
  if (is.null(SNV_map$is_CRV)) {
    SNV_map$is_CRV = FALSE
    SNV_map$is_CRV[sample(1:nrow(SNV_map), size = 1)] = TRUE
    warning("The variable is_CRV is missing from SNV_map.",
            "\n ...Randomly selecting one SNV to be the causal rare variant for all pedigrees")
  }

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
  #markers not carried by any member of our study.
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
