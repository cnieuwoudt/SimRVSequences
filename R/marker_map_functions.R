#' Create chromosome map from marker map
#'
#' @inheritParams sim_RVstudy
#'
#' @return a dataframe catalouging the start and stop positions, in base pairs, for each chromosome.  We use this information to determine what regions to simulate recombination over.
#' @export
create_chrom_map <- function(SNV_map){
  chrom_map <- do.call(rbind, lapply(sort(unique(SNV_map$chrom)), function(x){
    c(x, range(SNV_map$position[SNV_map$chrom == x]))
  }))

  chrom_map <- as.data.frame(chrom_map)
  colnames(chrom_map) = c("chrom", "start_pos", "end_pos")
  return(chrom_map)
}


#' Check SNV_map for possible issues
#'
#' @inheritParams sim_RVstudy
#' @keywords internal
#'
check_SNV_map <- function(SNV_map){
  #check to see if SNV_map contains the column information we expect
  # and check to see if we have any missing values.

  ## Check colID variable
  if(!c("colID") %in% colnames(SNV_map)){
    stop("The variable 'colID' is missing from SNV_map.")
  }
  if(any(is.na(SNV_map$colID))){
    stop("The variable 'colID' in the SNV_map dataset contains missing values.")
  }

  ## Check chrom variable
  if(!c("chrom") %in% colnames(SNV_map)){
    stop("The variable 'chrom' is missing from SNV_map.")
  }
  if(any(is.na(SNV_map$chrom))){
    stop("The variable 'chrom' in the SNV_map dataset contains missing values.")
  }

  ## Check position variable
  if(!c("position") %in% colnames(SNV_map)){
    stop("The variable 'position' is missing from SNV_map.")
  }
  if(any(is.na(SNV_map$position))){
    stop("The variable 'position' in the SNV_map dataset contains missing values.")
  }

  ## Check marker variable
  if(!c("marker") %in% colnames(SNV_map)){
    stop("The variable 'marker' is missing from SNV_map.")
  }
  if(any(is.na(SNV_map$marker))){
    stop("The variable 'marker' in the SNV_map dataset contains missing values.")
  }

}
