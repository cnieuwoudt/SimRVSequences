#' Create chromosome map from marker map
#'
#' @inheritParams sim_RVstudy
#'
#' @return a dataframe catalouging the start and stop positions, in base pairs, for each chromosome.  We use this information to determine what regions to simulate recombination over.
#' @export
create_chrom_map <- function(marker_map){
  chrom_map <- do.call(rbind, lapply(sort(unique(marker_map$chrom)), function(x){
    c(x, range(marker_map$position[marker_map$chrom == x]))
  }))

  chrom_map <- as.data.frame(chrom_map)
  colnames(chrom_map) = c("chrom", "start_pos", "end_pos")
  return(chrom_map)
}


#' Check marker_map for possible issues
#'
#' @inheritParams sim_RVstudy
#' @keywords internal
#'
check_marker_map <- function(marker_map){
  #check to see if marker_map contains the column information we expect
  # and check to see if we have any missing values.

  ## Check colID variable
  if(!c("colID") %in% colnames(marker_map)){
    stop("The variable 'colID' is missing from marker_map.")
  }
  if(any(is.na(marker_map$colID))){
    stop("The variable 'colID' in the marker_map dataset contains missing values.")
  }

  ## Check chrom variable
  if(!c("chrom") %in% colnames(marker_map)){
    stop("The variable 'chrom' is missing from marker_map.")
  }
  if(any(is.na(marker_map$chrom))){
    stop("The variable 'chrom' in the marker_map dataset contains missing values.")
  }

  ## Check position variable
  if(!c("position") %in% colnames(marker_map)){
    stop("The variable 'position' is missing from marker_map.")
  }
  if(any(is.na(marker_map$position))){
    stop("The variable 'position' in the marker_map dataset contains missing values.")
  }

  ## Check marker variable
  if(!c("marker") %in% colnames(marker_map)){
    stop("The variable 'marker' is missing from marker_map.")
  }
  if(any(is.na(marker_map$marker))){
    stop("The variable 'marker' in the marker_map dataset contains missing values.")
  }

}
