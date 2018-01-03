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

#' Determine if input is an odd number
#'
#' @param x Numeric.
#'
#' @return Boolean
#' @export
#'
is_odd <- function(x) {x %% 2 != 0}


#' Determine if input is an integer
#'
#' @param x Numeric.
#'
#' @return Boolean
#' @export
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
    warning(paste0("No disease-affected relative present in pedigree with FamID ",
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

#' Set Founder Genotypes and Return Reduced DF of Genotypes
#'
#' NEVER USED??? SHOULD PROBABLY BE REMOVED
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

#' Estimate haplotype distribution
#'
#' Estimate haplotype distribution
#'
#' @inheritParams sim_RVstudy
#' @param pop_haplos Data.frame.  A data frame containing the population haplotypes
#'
#' @return haplo_dist The haplotype distribution
#' @export
#' @importFrom plyr count
#'
#' @examples
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
#' founder_seq <- as.data.frame(matrix(sample(2000*nrow(mark_map),
#'                                     x = c(0, 1),
#'                                     replace = T,
#'                                     prob = c(0.95, 0.05)),
#'                              nrow = 2*1000))
#' colnames(founder_seq) = as.character(mark_map$marker)
#'
#' hdist = estimate_haploDist(founder_seq, mark_map)
#' hdist[[1]]
#' hdist[[2]]
#'
#'
estimate_haploDist <- function(pop_haplos, marker_map){
  haplo_dist <- list()
  for (i in 1:length(unique(marker_map$chrom))) {
    haplo_dist[[i]] <- count(pop_haplos[, c(which(marker_map$chrom == unique(marker_map$chrom)[i]))])
    haplo_dist[[i]]$freq <- haplo_dist[[i]]$freq/nrow(pop_haplos)
    colnames(haplo_dist[[i]]) <- c(marker_map$marker[marker_map$chrom == unique(marker_map$chrom)[i]], "prob")
  }

  return(haplo_dist)
}


#' Condition Haplotype distribution
#'
#' @param Chaplo_dist haplotype distribution for a single chromosome?
#' @param RV_marker familial RV marker
#' @param RV_status 0 or 1, 1 if RV is inherited
#'
#' @return The conditioned haplotype distribution.
#' @export
#'
condition_haploDist <- function(Chaplo_dist, RV_marker, RV_status){

  #determine which rows of Chaplo_dist are appropriate given rv status
  keep_rows <- which(Chaplo_dist[, which(colnames(Chaplo_dist) == RV_marker)] == RV_status)

  cond_haploDist <- Chaplo_dist[keep_rows, ]
  cond_haploDist$prob <- cond_haploDist$prob/sum(cond_haploDist$prob)

  return(cond_haploDist)
}


#' Draw Founder Genotypes from Haplotype Distribution Given Familial Risk Variant
#'
#'
#' @return list of familial founder genotypes
#' @export
#'
sim_FGenos <- function(founder_ids, RV_founder, FamID,
                       haplotype_dist, FamRV, marker_map) {

  #store chromosomal postion (i.e. list position in haplotype_dist) of risk variant
  RV_chrom_pos <- which(unique(marker_map$chrom) == marker_map$chrom[marker_map$marker == FamRV])

  fam_genos <- list()
  #NOTE: If we want to change from SNVs to another variant
  #type we will need to change the functionality of RV_status below

  for (k in 1:length(unique(marker_map$chrom))) {
    if (k == RV_chrom_pos) {
      #Since the risk variant is located on this chromosome, we must condition the
      #haplotype distribution on RV status before drawing haplotype
      RVdist <- condition_haploDist(haplotype_dist[[k]],
                                    RV_marker = FamRV,
                                    RV_status = 1)

      RV_founder_dat = RVdist[sample(x = c(1:nrow(RVdist)),
                                     size = 1,
                                     prob = RVdist$prob),
                              -ncol(RVdist)]

      NRVdist <- condition_haploDist(haplotype_dist[[k]],
                                     RV_marker = FamRV,
                                     RV_status = 0)

      NRV_founder_dat <- NRVdist[sample(x = c(1:nrow(NRVdist)),
                                        size = (2*length(founder_ids) + 1),
                                        replace = TRUE,
                                        prob = NRVdist$prob),
                                 -ncol(NRVdist)]

      fam_genos[[k]] <- rbind(RV_founder_dat, NRV_founder_dat)
    } else {

      #since the familial RV does not reside on the kth chromosome we
      #sample founder haplotypes according to from the haplotype distribution

      fam_genos[[k]] = haplotype_dist[[k]][sample(x = c(1:nrow(haplotype_dist[[k]])),
                                                  size = 2*(length(founder_ids) + 1),
                                                  replace = TRUE,
                                                  prob = haplotype_dist[[k]]$prob),
                                           -ncol(haplotype_dist[[k]])]
    }
  }

  founder_genos <- do.call("cbind", fam_genos)
  founder_genos$ID <- rep(c(RV_founder, founder_ids), each = 2)
  rownames(founder_genos) <- NULL

  #ramdomly permute RV founders rows so that the RV is not always paternally inherited.
  founder_genos[c(1,2), ] <- founder_genos[sample(x = c(1, 2), size = 2, replace = F), ]

  return(founder_genos)
}
