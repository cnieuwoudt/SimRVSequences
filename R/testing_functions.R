#' Simulate chiasmata positions
#'
#' Simulate chiasmata positions
#'
#' Issues: from what exponential distribution should the original start pos come from? Currently using exponential with mean 0.5
#'
#' Also, need to verify that parmaters chosen for gamma distribution result in an average distance of 1 crossover per M.
#'
#' @param chrom_map Data.frame with 1 row and 2 columns. The two columns represent the start and stop positions (in cM) over which to simulate recombination.
#' @param mean_dist Numeric. The mean of the expoential distribution.
#'
#' @return A list containing chiasmata postions.
#' @export
#'
#' @examples
#'
#' sim_chiasmataPositionsNI(chrom_map = data.frame(start = 47, stop = 198, center = 100))
#'
#' set.seed(1)
#' system.time(for(i in 1:10000){
#' sim_chiasmataPositionsNI(chrom_map = data.frame(start = 0, stop = 200, center = 50))
#' })
#'
#'
sim_chiasmataPositionsNI <- function(chrom_map,
                                     mean_dist = 0.5){

  chiasmata_pos <- rexp(1, rate = 1/mean_dist)*100

  #simulate chiasmata along the length of the chromosome
  while (length(chiasmata_pos[which(chiasmata_pos >= chrom_map[1, 2] )]) == 0){
    chiasmata_pos <- c(chiasmata_pos,
                       chiasmata_pos[length(chiasmata_pos)] +
                         cumsum(rexp(max(5, 5*round((diff(as.numeric(chrom_map)))/50)),
                                     rate = 1/mean_dist)*100))
  }

  chiasmata_pos <- chiasmata_pos[which(chiasmata_pos < chrom_map[1, 2] )]

  return(chiasmata_pos)

}
