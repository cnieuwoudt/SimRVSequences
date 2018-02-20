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
#' set.seed(1960)
#' ex_RVped <- sim_RVped(hazard_rates = hazard(hazardDF = AgeSpecific_Hazards),
#'                       GRR = 10, RVfounder = TRUE, FamID = 1,
#'                       founder_byears = c(1900, 1920),
#'                       ascertain_span = c(1995, 2015),
#'                       num_affected = 2, stop_year = 2017,
#'                       recall_probs = c(1, 1, 0))[[2]]
#' plot(ex_RVped)
#'
#' # say we are simulating gamete inheritance from the parent
#' # with ID 1 to the offspring with ID 4
#'
#' data(mark_map)
#' data(hg_chrom)
#' head(mark_map)
#' head(hg_chrom)
#'
#' my_chrom_map = hg_chrom[17, ]
#' my_chrom_map$start_pos = convert_BP_to_cM(my_chrom_map$start_pos)
#' my_chrom_map$end_pos = convert_BP_to_cM(my_chrom_map$end_pos)
#' my_chrom_map$center = convert_BP_to_cM(my_chrom_map$center)
#'
#' sim_gameteInheritance(RV_locus = mark_map[1, 1:2],
#'                       parent_RValleles = c(0, 1),
#'                       offspring_RVstatus = c(0),
#'                       chrom_map = my_chrom_map,
#'                       allele_IDs = c(2, 3))
#'
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

    # get the list location corresponding of the chromosome that the RV is on.
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
