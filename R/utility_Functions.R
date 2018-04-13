#' Assign gamete order
#'
#' For internal use.
#'
#' @return A list corresponding to assigned gamete order.
#' @keywords internal
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
#' @keywords internal
chias_count_BC <- function(chiasmata_pos, center_loc){
  ifelse(length(which(chiasmata_pos < center_loc)) > 0,
         max(which(chiasmata_pos < center_loc)),
         0)
}

#' Get parent and offspring information from a pedigree
#'
#' \strong{For internal use.}
#'
#' @param ped_file Data.frame. The pedigree file, must have same format as pedigree simulated with \code{sim_RVped}
#'
#' @return A list containing the parent's paternal and maternal alleles at the disease locus, and the RV status of the offspring
#' @keywords internal
#' @importFrom reshape2 melt
get_parOffInfo <- function(ped_file){

  mdata <- melt(ped_file[which(!is.na(ped_file$dadID)),
                         which(colnames(ped_file) %in%
                                 c("ID", "dadID", "momID", "Gen"))],
                id = c("ID", "Gen"))
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
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' haplo_vec <- sample(x = c(2, 3), size = 10, replace = TRUE)
#' chias_vec <- cumsum(rgamma(9, shape = 2.63, rate = 2*2.63))
#' reduce_to_events(gamete_haplo = haplo_vec,
#'                  chias_locations = chias_vec)
#'}
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

#' Determine if input is an odd number
#'
#' @param x Numeric.
#'
#' @return Boolean
#' @keywords internal
#'
is_odd <- function(x) {x %% 2 != 0}


#' Determine if input is an integer
#'
#' @param x Numeric.
#'
#' @return Boolean
#' @keywords internal
#'
is_int <- function(x) {x %% 1 == 0}

#' Remove unaffected relatives
#'
#' Remove unaffected relatives
#'
#' @param ped_file data.frame. A pedigree.
#'
#' @return \code{retA_ped} A pedigree containing only affected members, obligate carriers, and founders.
#' @export
#'
#' @examples
#' library(SimRVPedigree)
#' #Read in example pedigrees and create ped object
#' data(EgPeds)
#' ex_peds <- new.ped(EgPeds)
#'
#' #plot full pedigree
#' plot(ex_peds[which(ex_peds$FamID == 1), ])
#'
#' #reduce to affected only pedigree
#' Apeds = affected_onlyPed(ex_peds[which(ex_peds$FamID == 1), ])
#' plot(Apeds)
#'
affected_onlyPed = function(ped_file){

  #create new ped file with affecteds only
  retA_ped <- ped_file[ped_file$affected, ]

  if (nrow(retA_ped) == 0) {
    warning(paste0("Disease-affected relatives are not present in pedigree with FamID ",
                   sep = "", ped_file$FamID[1]))
    return(retA_ped)
  } else {
    d <- 0
    while (d == 0) {
      #find the dad IDs that are required but have been removed
      miss_dad  <- !is.element(retA_ped$dadID,
                               retA_ped$ID[which(retA_ped$sex == 0)])
      readd_dad <- retA_ped$dadID[miss_dad]
      readd_dad <- unique(readd_dad[!is.na(readd_dad)])

      #find the mom IDs that are required but have been removed
      miss_mom  <- !is.element(retA_ped$momID,
                               retA_ped$ID[which(retA_ped$sex == 1)])
      readd_mom <- retA_ped$momID[miss_mom]
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

#' Convert from basepairs to centimorgan
#'
#' Convert from basepairs to centimorgan
#'
#' @param pos_BP Numeric.  The position in basepairs.
#'
#' @return pos_CM The postion in centiMorgans
#' @export
#'
convert_BP_to_cM <- function(pos_BP){ pos_BP/1000000 }

#' Convert from centiMorgan to basepairs
#'
#' Convert from centiMorgan to basepairs
#'
#' @param pos_CM Numeric.  The position in centiMorgan.
#'
#' @return pos_BP The postion in basepairs
#' @export
#'
convert_CM_to_BP <- function(pos_CM){ pos_CM*1000000 }
