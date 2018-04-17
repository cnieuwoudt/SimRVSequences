#' Create chromosome map from marker map
#'
#' @inheritParams sim_RVstudy
#'
#' @return a dataframe catalouging the start and stop positions, in base pairs, for each chromosome.  We use this information to determine what regions to simulate recombination over.
#' @keywords internal
create_chrom_map <- function(marker_map){
  chrom_map <- do.call(rbind, lapply(sort(unique(marker_map$chrom)), function(x){
    c(x, range(marker_map$position[marker_map$chrom == x]))
  }))

  chrom_map <- as.data.frame(chrom_map)
  colnames(chrom_map) = c("chrom", "start_pos", "end_pos")
  return(chrom_map)
}
