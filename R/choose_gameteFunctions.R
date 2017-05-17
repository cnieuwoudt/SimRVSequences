#' Simulate inheritance of parental gamete to offspring
#'
#' Simulate inheritance of parental gamete to offspring based on rare variant statuses of parent and offspring.
#'
#' Here we use the RV statuses of the parent and offspring to determine which of the parental gametes are appropriate options for transmission.  Upon reducing the sample space appropriately we choose from the remaining options with equal probability.
#' \enumerate{
#' \item If the parent \strong{is not} a carrier of the rare variant then we choose any of the four gametes with equal probability since the offpring could not have inherited the rare variant from this parent.
#' \item If the parent \strong{is} a carrier of the rare variant and
#' \itemize{
#' \item the offspring \strong{is not} a carrier of the rare variant, we choose with equal probability from the two gametes that \strong{do not} contain the rare variant.
#' \item the offspring \strong{is} a carrier of the rare variant, we choose with equal probability from the two gametes that \strong{do}contain the rare variant.
#' }
#'}
#'
#' @inheritParams sim_gameteFormation
#' @param RV_locus Numeric list of length 2. A list containing (1) the chromosome upon which the rare variant resides (2) the position in cM where the rare variant resides.
#' @param parent_RValleles Numeric list of length 2. The paternal and maternal alleles at the disease locus (1 = RV inherited, 0 otherwise)
#' @param offspring_RVstatus Numeric. 1 if offspring inherits the RV from parent, 0 otherwise.
#'
#' @return A list containing (1) a list of inherited haplotype codings (2) the chiasmata locations
#' @export
#'
#' @examples
#' library(SimRVPedigree)
#' #Read in age-specific hazards
#' data(AgeSpecific_Hazards)
#'
#' # Simulate pedigree ascertained for multiple affected individuals
#' set.seed(13)
#' ex_RVped <- sim_RVped(onset_hazard = AgeSpecific_Hazards[, 1],
#'                       death_hazard = AgeSpecific_Hazards[, c(2,3)],
#'                       part = seq(0, 100, by = 1),
#'                       RR = 15, FamID = 1,
#'                       founder_byears = c(1900, 1910),
#'                       ascertain_span = c(1900, 2015),
#'                       num_affected = 2, stop_year = 2015,
#'                       recall_probs = c(1, 0.75, 0.5))[[2]]
#' library(kinship2)
#' ex_ped <- with(ex_RVped, pedigree(id = ID,
#'                                   dadid = dad_id,
#'                                   momid = mom_id,
#'                                   sex = gender + 1,
#'                                   affected = cbind(Affected = affected,
#'                                                    Proband = proband,
#'                                                    RV_status = DA1 + DA2),
#'                                  status = ifelse(is.na(death_year), 0, 1)))
#' plot(ex_ped)
#'
#' # say we are simulating gamete inheritance from the parent
#' # with ID 1 to the offspring with ID 3
#' my_chrom_map = data.frame(chrom     = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
#'                                         13, 14, 15, 16, 17, 18, 19, 20, 21, 22),
#'                           start_pos = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'                                         0, 0, 0, 0, 0, 0, 0, 0, 0),
#'                           end_pos   = c(250, 200, 190, 180, 170, 160, 150, 140, 130, 120,
#'                                         250, 200, 190, 180, 170, 160, 150, 140, 130, 120,
#'                                         100, 100))
#' my_chrom_map$center = round(my_chrom_map$end_pos*0.2)
#'
#' my_chrom_map = data.frame(chrom     = c(1),
#'                           start_pos = c(0),
#'                           end_pos   = c(250),
#'                           center    = c(50))
#'
#' sim_gameteInheritance(RV_locus = c(1, 150),
#'                       parent_RValleles = c(0, 1),
#'                       offspring_RVstatus = c(0),
#'                       chrom_map = my_chrom_map,
#'                       allele_IDs = c(2, 3))
#'
#' my_parent_RValleles = ex_RVped[which(ex_RVped$ID == 1),
#'                                which(colnames(ex_RVped) %in% c("DA1", "DA2"))]
#'
#' my_offspring_RVstatus = sum(ex_RVped[which(ex_RVped$ID == 3),
#'                                      which(colnames(ex_RVped) %in%
#'                                      c("DA1", "DA2"))])
#'
#' my_allele_IDs = c(2, 3)
#'
#' sim_gameteInheritance(RV_locus = c(1, 150),
#'                       parent_RValleles = my_parent_RValleles,
#'                       offspring_RVstatus = my_offspring_RVstatus,
#'                       chrom_map = my_chrom_map,
#'                       allele_IDs = c(2, 3))
#'
#' dad_combos <- ex_RVped[which(!is.na(ex_RVped$dad_id)), c(2, 4)]
#' mom_combos <- ex_RVped[which(!is.na(ex_RVped$mom_id)), c(2, 5)]
#'
#' colnames(dad_combos)[2] <- colnames(mom_combos)[2] <- "par_ID"
#' po_combos <- rbind(dad_combos, mom_combos)
#'
#' system.time( for (i in 1:nrow(po_combos)) {
#'  sim_gameteInheritance(RV_locus = my_RV_locus,
#'                       parent_RValleles = ex_RVped[which(ex_RVped$ID == po_combos[i, 2]), c(7, 8)],
#'                       offspring_RVstatus = sum(ex_RVped[which(ex_RVped$ID == po_combos[i, 1]), c(7, 8)]),
#'                       chrom_map = my_chrom_map,
#'                       allele_IDs = c(2, 3))
#' })
#'
sim_gameteInheritance <- function(RV_locus, parent_RValleles,
                                  offspring_RVstatus,
                                  chrom_map, allele_IDs,
                                  burn_in = 1000, gamma_params = c(2.63, 2.63/0.5)) {

  parental_gametes <- sim_gameteFormation(chrom_map, allele_IDs, burn_in, gamma_params)

  # check to see if parent is a carrier of the RV, if not
  # we can choose between the 4 gametes with equal probability
  # since we know that the offspring could not have inherited the RV
  # from this parent
  if (sum(parent_RValleles) == 0) {
    inherited_Ggrp <- sample(1, x = c("A", "B", "C", "D"))
  } else {
    # determine which of the parental gametes are acceptable choices for
    # transmission to offspring given the offspring's RV status.
    # NOTE: The offspring can only inherit 1 copy of the RV so we don't have to
    # worry if the RV is maternally/paternally inherited by an offspring, we
    # only have to worry about maternal/paternal inheritance in the parent

    # get the list location corresponding to the chromosome that the RV is located on
    RV_chromLoc <- which(chrom_map[, 1] == RV_locus[1, 1])

    # Store the halpotype df for the RV chromosome
    RV_chromHaps <- parental_gametes[[1]][[RV_chromLoc]]

    #store the chiasmata locations for the RV chromosome
    RV_chromChias <- parental_gametes[[2]][[RV_chromLoc]]

    #Identify the parental allele ID that corresponds to the allele
    #inherited by offspring at the RV locus
    RV_allele <- allele_IDs[which(parent_RValleles == offspring_RVstatus)]

    #determine which halpotype interval the RV lies in
    if(length(RV_chromChias) == 0){
      col_loc <- 1
    } else if (RV_chromChias[length(RV_chromChias)] <  RV_locus[1, 2]){
      col_loc <- length(RV_chromChias) + 1
    } else{
      col_loc <- min(which(((RV_chromChias - RV_locus[1, 2]) > 0) == TRUE))
    }

    #determine which gamete groups contain the appropriate
    # RV, given the offsprings RV status
    inherited_Ggrp <- sample(1, x = RV_chromHaps[which(RV_chromHaps[, col_loc] == RV_allele),
                                                 ncol(RV_chromHaps)])
  }

  inherited_haplotypes <- lapply(c(1:nrow(chrom_map)), function(x){
    parental_gametes[[1]][[x]][which(parental_gametes[[1]][[x]][, ncol(parental_gametes[[1]][[x]])] == inherited_Ggrp),
                               -ncol(parental_gametes[[1]][[x]])]
    })

  return(list(haplotypes = inherited_haplotypes, cross_locations = parental_gametes[[2]]))

}
