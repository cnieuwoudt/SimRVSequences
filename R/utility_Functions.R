#' Assign gamete order
#'
#' For internal use.
#'
#' @return A list corresponding to assigned gamete order.
#' @export
#'
#' @examples
#' gam_order()
#'
gam_order <- function(){

  combos <- matrix(c(rep(c("A", "B"), each = 2), rep(c("C", "D"), 2),
                     rep(c("B", "A"), each = 2), rep(c("D", "C"), 2),
                     rep(c("C", "D"), 2), rep(c("A", "B"), each = 2),
                     rep(c("D", "C"), 2), rep(c("B", "A"), each = 2)),
                   nrow = 8)

  return(combos[sample(1, x = c(1:8)), ])
}


#' Determine number of chiasmata before centromere
#'
#' For internal use.
#'
#' @param chiasmata_pos The simulated chiasmata postions
#' @param center_loc The centromere locations
#'
#' @return The number of chiasmata before the centromere
#' @export
#'
#' @examples
#' ex_chias <- sim_chiasmataPositions(chrom_map = data.frame(start = 47, stop = 198))
#' ex_chias
#' chias_count_BC(ex_chias, 78)
chias_count_BC <- function(chiasmata_pos, center_loc){
  ifelse(length(which(chiasmata_pos < center_loc)) > 0,
         max(which(chiasmata_pos < center_loc)),
         0)
}

#' Get parent and offspring information from a pedigree
#'
#' For internal use.
#'
#' @param ped_file Data.frame. The pedigree file, must have same format as pedigree simulated with \code{sim_RVped}
#' @param offspring_ID Numeric. The ID of the offspring.
#' @param parent_ID Numeric. The ID of the parent.
#'
#'
#' @return A list containing the parent's paternal and maternal alleles at the disease locus, and the RV status of the offspring
#' @export
#' @importFrom reshape2 melt
#'
get_parOffInfo <- function(ped_file){

  mdata <- melt(ped_file[which(!is.na(ped_file$dad_id)),
                         which(colnames(ped_file) %in%
                                 c("ID", "dad_id", "mom_id", "Gen"))],
                id = c("ID", "Gen"))
  mdata <- mdata[, ]
  colnames(mdata) = c("offspring_ID", "Gen", "parent", "parent_ID")


  mdata$Off_RVstatus <- apply(as.data.frame(mdata[, 1]), 1, function(x){
    sum(ped_file[which(ped_file$ID == x), which(colnames(ped_file) %in% c("DA1", "DA2"))])
  })

  mdata$Par_DA1 <- apply(as.data.frame(mdata[, 4]), 1, function(x){
    sum(ped_file[which(ped_file$ID == x), which(colnames(ped_file) == "DA1")])
  })

  mdata$Par_DA2 <- apply(as.data.frame(mdata[, 4]), 1, function(x){
    sum(ped_file[which(ped_file$ID == x), which(colnames(ped_file) == "DA2")])
  })

  return(mdata)
}


#' Reduce chiasmata vector to crossovers that transmitted gamete participates in based on the allele vector.
#'
#' For internal use.
#'
#' @param gamete_haplo Numeric vector. The inherited haplotype.
#' @param chias_locations  Numeric vector.  Chiasmata locations.
#'
#' @return The locations of crossovers
#' @export
#'
#' @examples
#' haplo_vec <- sample(x = c(2, 3), size = 10, replace = T)
#' chias_vec <- cumsum(rgamma(9, shape = 2.63, rate = 2*2.63))
#' reduce_to_events(gamete_haplo = haplo_vec,
#'                  chias_locations = chias_vec)
#'
reduce_to_events <- function(gamete_haplo, chias_locations){
  if(sum(gamete_haplo == gamete_haplo[1]) == length(gamete_haplo)){
    cross_loc = numeric(0)
  } else {
    keep_ind <- c()
    i = 1

    #find position of last match to first element of gamete_haplo
    keep_ind[i] <- Position(function(x){x != gamete_haplo[1]}, gamete_haplo) - 1
    d = keep_ind[i]

    while (d < length(gamete_haplo)) {
      rhap <- gamete_haplo[c((keep_ind[i] + 1) : length(gamete_haplo))]
      if ( sum(rhap == rhap[1]) == length(rhap) | length(rhap) == 1 ) {
        #if the remaining items to check have length 1 or are
        #all the same element we are done so break while loop
        d <- length(gamete_haplo) + 1
      } else {
        keep_ind[i + 1] <- keep_ind[i] +
          (Position(function(x){x != rhap[1]}, rhap) - 1)

        d <- keep_ind[i + 1]
        i = i + 1
      }
    }
    cross_loc <- chias_locations[keep_ind]
  }
  return(cross_loc)
}

#' Determine if a number os odd
#'
#' @param x Numeric.
#'
#' @return Boolean
#' @export
#'
is_odd <- function(x) {x %% 2 != 0}

#' Remove unaffected relatives
#'
#' Remove unaffected relatives
#'
#' @param ped_file data.frame. A pedigree.
#'
#' @return \code{retA_ped} A pedigree containing only affected members, obligate carriers, and founders.
#' @export
#'
#' @references OUR MANUSCRIPT
#' @references Thompson, E. (2000). \emph{Statistical Inference from Genetic Data on Pedigrees.} NSF-CBMS Regional Conference Series in Probability and Statistics, 6, I-169. Retrieved from http://www.jstor.org.proxy.lib.sfu.ca/stable/4153187
#'
#' @examples
#' library(SimRVPedigree)
#' #Read in example pedigrees
#' data(EgPeds)
#'
#' library(kinship2)
#' #plot full pedigree
#' ex_pedigree <- pedigree(id = EgPeds$ID,
#'                         dadid = EgPeds$dad_id,
#'                         momid = EgPeds$mom_id,
#'                         sex = (EgPeds$gender + 1),
#'                         affected = EgPeds$affected,
#'                         famid = EgPeds$FamID)
#' plot(ex_pedigree['1'])
#'
#' #reduce to affected only pedigree
#' Apeds = affected_onlyPed(EgPeds[which(EgPeds$FamID == 1), ])
#' red_ped <- pedigree(id = Apeds$ID,
#'                     dadid = Apeds$dad_id,
#'                     momid = Apeds$mom_id,
#'                     sex = (Apeds$gender + 1),
#'                     affected = Apeds$affected,
#'                     famid = Apeds$FamID)
#' plot(red_ped['1'])
#'
affected_onlyPed = function(ped_file){

  #create new ped file with affecteds only
  retA_ped <- ped_file[which(ped_file$affected == 1), ]

  if (nrow(retA_ped) == 0) {
    warning("No affecteds to assign affected generation")
    return(retA_ped)
  } else {
    d <- 0
    while (d == 0) {
      #find the dad IDs that are required but have been removed
      miss_dad  <- !is.element(retA_ped$dad_id,
                               retA_ped$ID[which(retA_ped$gender == 0)])
      readd_dad <- retA_ped$dad_id[miss_dad]
      readd_dad <- unique(readd_dad[!is.na(readd_dad)])

      #find the mom IDs that are required but have been removed
      miss_mom  <- !is.element(retA_ped$mom_id,
                               retA_ped$ID[which(retA_ped$gender == 1)])
      readd_mom <- retA_ped$mom_id[miss_mom]
      readd_mom <- unique(readd_mom[!is.na(readd_mom)])

      #check to see if we need to readd anyone
      if (length(c(readd_dad, readd_mom)) == 0) {
        d <- 1
      } else {
        #Now pull the rows containing the required parents
        # from the original ped_file
        readd <- ped_file[which(ped_file$ID %in% c(readd_dad, readd_mom)), ]

        #combine with affected ped file
        retA_ped <- rbind(retA_ped, readd)
      }
    }
  }
  return(retA_ped)
}

#' Set Founder Genotypes and Return Reduced DF of Genotypes
#'
#'
#' @return list of (1) family genotypes and (2) a reduced pool of genos
#' @export
#'
get_FGenos <- function(founder_ids, RV_founder, FamID, founder_genotypes, FamRV) {
  #choose row corresponding to seq_data for RV founder
  RV_founder_row <- sample(size = 1,
                           x = which(founder_genotypes[ , which(colnames(founder_genotypes) == FamRV)] == 1))

  #choose row corresponding to seq_data for other founders
  NOTRV_founder_row <- sample(size = 2*length(founder_ids) + 1,
                           x = which(founder_genotypes[ , which(colnames(founder_genotypes) == FamRV)] == 0),
                           replace = FALSE)

  #combine into 1 DF
  fam_genos <- rbind(founder_genotypes[RV_founder_row, ],
                     founder_genotypes[NOTRV_founder_row, ])

  #update ID variable
  fam_genos$ID <- c(rep(RV_founder, 2), rep(founder_ids, each = 2))

  #ramdomly permute RV founders rows so that half of the time the RV is a maternally inherited and the other half of the time it is paternally inherited.
  fam_genos[c(1,2), ] <- fam_genos[sample(x = c(1, 2), size = 2, replace = F), ]

  #reduce genotypes DF so that these rows cannot be chosen again
  red_genotypes <- founder_genotypes[-c(RV_founder_row, NOTRV_founder_row), ]

  my_return <- list(fam_genos, red_genotypes)
  return(my_return)
}
