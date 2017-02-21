
#' Simulate chiasmata positions
#'
#' Simulate chiasmata positions
#'
#' Issues: from what exponential distribution should the original start pos come from? Currently using exponential with mean 0.5
#'
#' @param chrom_length Numeric. The chromosome length, in M, over which to simulate chiasmata.
#' @param gamma_params Numeric list of lenght 2. The respective shape and scale parameters gamma distribution used to simulate distance between chiasmata, default from Vorrips.
#'
#' @return A vector of chiasmata postions.
#' @export
#'
#' @examples
#'
#' sim_chiasmataPositions(chrom_length = 25)
#'
#' set.seed(1)
#' system.time(for(i in 1:10000){
#' sim_chiasmataPositions(chrom_length = 500)
#' })
#'
sim_chiasmataPositions <- function(chrom_length,
                                   gamma_params = c(2.63, 2.63/2)){
  #generate random starting point before chrom_startPos
  current_pos <- 0 - rexp(1, 0.5)

  while(current_pos < 0) {
    first_chias <- rgamma(1, shape = gamma_params[1],
                          scale = gamma_params[2])
    if(current_pos + first_chias > 0 ){
      chiasmata_pos <- current_pos <- current_pos + first_chias
    }
  }

  chiasmata_pos <- c(chiasmata_pos,
                     chiasmata_pos + cumsum(rgamma(round(chrom_length),
                                                   shape = gamma_params[1],
                                                   scale = gamma_params[2])))

  #make sure we simulated enough chiasmata along the chromosome
  end_reached = FALSE
  while ( end_reached == FALSE){
    if ( max(chiasmata_pos) >= chrom_length ) {
      end_reached = TRUE
      chiasmata_pos <- chiasmata_pos[which(chiasmata_pos < chrom_length)]
    } else {
      chiasmata_pos <- c(chiasmata_pos,
                         chiasmata_pos + cumsum(rgamma(round(chrom_length)/2,
                                                       shape = gamma_params[1],
                                                       scale = gamma_params[2])))
    }

  }

  return(chiasmata_pos)

}
