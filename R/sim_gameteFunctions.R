#' Simulate chiasmata positions
#'
#' Simulate chiasmata positions among a group of four chromatids.
#'
#'
#' @param chrom_map Data.frame with 1 row and 2 columns. The two columns represent the start and stop positions (in cM) over which to simulate recombination.
#' @param gamma_params Numeric list of length 2. The respective shape and scale parameters gamma distribution used to simulate distance between chiasmata, default from Vorrips.
#'
#' @return A list containing chiasmata postions and the number of chiasmata that occur before the centromere.
#' @export
#'
#' @importFrom stats rexp
#' @importFrom stats rgamma
#'
#' @examples
#'
#' sim_chiasmataPositions(chrom_map = data.frame(start = 47, stop = 198, center = 100))
#'
#' set.seed(1)
#' system.time(for(i in 1:10000){
#' sim_chiasmataPositions(chrom_map = data.frame(start = 0, stop = 200, center = 50))
#' })
#'
#'
sim_chiasmataPositions <- function(chrom_map,
                                   mean_dist = 0.5,
                                   gamma_params = c(2.63, 2.63/0.5)){

  # generate first chiasmata position using burn in process
  # suggested in Voorrips 2012, see distToFirstChiasma in Multivalent.java
  # and ranExp, ranGamma, and ranDistInterference in Tools.java
  #Muliply by 10 as suggested by Voorrips, and then by 1000 to convert to cM
  burnDist <- 1000*mean_dist
  try_pos <- rexp(1, rate = 1/mean_dist)*100

  while(length(try_pos[which(try_pos > burnDist)]) == 0){
    try_pos <- c(try_pos, try_pos[length(try_pos)] +
                   rexp(10, rate = 1/mean_dist)*100
                   #cumsum(rgamma(5 + (burnDist*gamma_params[2])/(gamma_params[1]*100),
                    #             shape = gamma_params[1],
                    #             rate = gamma_params[2])*100)
                 )
  }

  chiasmata_pos <- try_pos[min(which(try_pos > burnDist))] - burnDist

  #simulate chiasmata along the length of the chromosome
  while (length(chiasmata_pos[which(chiasmata_pos >= chrom_map[1, 2] )]) == 0){
    chiasmata_pos <- c(chiasmata_pos,
                       chiasmata_pos[length(chiasmata_pos)] +
                         rexp(10, rate = 1/mean_dist)*100
                         #cumsum(rgamma(max(5, 5*round((diff(as.numeric(chrom_map)))/50)),
                          #             shape = gamma_params[1],
                          #             rate = gamma_params[2])*100)
                       )
  }

  chiasmata_pos <- chiasmata_pos[which(chiasmata_pos < chrom_map[1, 2] )]

  return(chiasmata_pos)

}



#' Simulate recombination along a bundle of four chromatids.
#'
#' Simulate recombination along a bundle of four chromatids.
#'
#' Given the possible chiasmata positions returned from \code{sim_chiasmataPostions}, we want to randomly select two non-sister chromatids to participate in each recombination event.  We assume no chromatid interference so that the non-sister chromatids participating in a crossover event are independent of those chosen in previous crossover events.
#'
#' @param num_chiasmata Numeric. The number of chiasmata to simulate among the chromatid bundle.
#' @param allele_IDs Numeric list of length 2. The identification numbers for the respective paternal and maternal alleles of the individual for whom we wish to simulate recombination.
#'
#' @return haploid_mat. A matrix with rows representing recombined haplotypes.
#' @export
#'
#' @seealso \code{\link{sim_chiasmataPostions}}
#'
#' @examples
#' my_chrom_map <- data.frame(start = 0, stop = 200, center = 40)
#' sim_chias_pos <- sim_chiasmataPositions(my_chrom_map)
#' chias_count_BC(sim_chias_pos, 40)
#' sim_haploidFormation(num_chiasmata = length(sim_chias_pos),
#'                      before_center = chias_count_BC(sim_chias_pos, my_chrom_map[1,3]),
#'                      allele_IDs = c(2, 3))
#'
#' system.time(for (i in 1:10000) {
#' sim_chias_pos <- sim_chiasmataPositions(my_chrom_map)
#' sim_haploidFormation(num_chiasmata = length(sim_chias_pos),
#'                      before_center = chias_count_BC(sim_chias_pos, my_chrom_map[1,3]),
#'                      allele_IDs = c(2, 3))
#' })
#'
sim_haploidFormation <- function(num_chiasmata,
                                 before_center,
                                 allele_IDs) {

  #each column in haploid_mat represents the alleles on either side of a chiasmata
  haploid_mat <- matrix(rep(allele_IDs, each = 2*(num_chiasmata + 1)),
                        nrow = 4, byrow = T)

  #choose non-sister chromatids to participate in each chiasmata event
  choose_gam = matrix(c(sample(c(1, 2), num_chiasmata, replace = T),
                        sample(c(3, 4), num_chiasmata, replace = T)), ncol = 2)

  #swap allele sequences at sucessive chiasmata locations from left to right
  if(num_chiasmata > 0){
    for (i in 1:num_chiasmata) {
      haploid_mat[choose_gam[i, ], c(1: i)] <- haploid_mat[rev(choose_gam[i, ]), c(1: i)]
    }
  }

  #reorder the rows so that sister chromatids are together at the centromeres
  haploid_mat <- as.data.frame(haploid_mat[order(haploid_mat[, (before_center + 1)]), ])
  haploid_mat$gamete_grp <- gam_order()

  return(haploid_mat)

}


#' Simulate gamete formation
#'
#' Simulate gamete formation
#'
#' @inheritParams sim_haploidFormation
#'
#' @param chrom_map Data.frame.  A data.frame consisting of three columns: column 1 contains the chromosome numbers, column 2 start postion of chromosome (in cM), column 3 end position of chromosome (in cM).
#'
#' @return chrom_haps and gamete_group
#' @export
#'
#' @examples
#' my_chrom_map = data.frame(chrom = c(1, 2, 3),
#'                           start_pos = c(0, 0, 0),
#'                           end_pos = c(250, 120, 75),
#'                           center = c(50, 24, 15))
#'
#' sim_gameteFormation(my_chrom_map, allele_IDs = c(0, 1))
#'
#' system.time(for (i in 1:10000) {
#'   sim_gameteFormation(my_chrom_map, c(0, 1))
#' })
sim_gameteFormation <- function(chrom_map, allele_IDs) {

  chrom_chiasmataPos <- lapply(c(1:nrow(chrom_map)),
                               function(x){
                                 sim_chiasmataPositions(chrom_map = chrom_map[x, -1])
                                 })

  BC_count <- lapply(c(1:nrow(chrom_map)),
                     function(x){
                       chias_count_BC(chrom_chiasmataPos[[x]],
                                      center_loc = chrom_map[x, 4])
                     })

  chrom_haps <- lapply(c(1:nrow(chrom_map)),
                       function(x){
                         sim_haploidFormation(num_chiasmata = length(chrom_chiasmataPos[[x]]),
                                              before_center = BC_count[[x]],
                                              allele_IDs)
                         })

  fun_return <- list(chrom_haps, chrom_chiasmataPos)

  return(fun_return)

}

