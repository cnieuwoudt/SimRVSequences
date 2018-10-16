#' Check SNV_map for possible issues
#'
#' INTENDED FOR INTERNAL USE ONLY
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

#' Checks individual pedigrees for proper format.
#'
#' @param ped_file data.frame The pedigree.
#'
#' @keywords internal
#'
check_ped <- function(ped_file){

  moms <- unique(ped_file$momID[!is.na(ped_file$momID)])
  dads <- unique(ped_file$dadID[!is.na(ped_file$dadID)])

  if (any(ped_file$sex[which(ped_file$ID %in% moms)] != 1) |
      any(ped_file$sex[which(ped_file$ID %in% dads)] != 0)){

    wrong_sex <- c(ped_file$ID[which(ped_file$sex[which(ped_file$ID %in% dads)] != 0)],
                   ped_file$ID[which(ped_file$sex[which(ped_file$ID %in% moms)] != 1)])

    stop(paste0('Sex improperly specifed ID: ', sep = '', wrong_sex, '. \nPlease ensure that for males: sex = 0; and for females: sex = 1.'))
  }

  if (any(!moms %in% ped_file$ID) | any(!dads %in% ped_file$ID)) {

    wrong_par <- c(ped_file$ID[which(ped_file$momID == moms[which(!moms %in% ped_file$ID)])],
                   ped_file$ID[which(ped_file$dadID == dads[which(!dads %in% ped_file$ID)])])

    stop(paste0('ID: ', sep = '', wrong_par, '.  Non-founders must have a mother and a father. Founders have neither.'))
  }

  if (any(!is.na(ped_file$momID[is.na(ped_file$dadID)])) |
      any(!is.na(ped_file$dadID[is.na(ped_file$momID)]))) {
    stop("Non-founders must have both a mother and a father, while founders have neither.")
  }
}


#' Checks ped_files for expected info and format.
#'
#' @inheritParams sim_RVstudy
#'
#' @keywords internal
#'
check_peds <- function(ped_files){

  if (!"FamID" %in% colnames(ped_files) |
      !"ID" %in% colnames(ped_files) |
      !"dadID" %in% colnames(ped_files) |
      !"momID" %in% colnames(ped_files) |
      !"sex" %in% colnames(ped_files) |
      !"affected" %in% colnames(ped_files) |
      !"DA1" %in% colnames(ped_files) |
      !"DA2" %in% colnames(ped_files)) {
    stop('ped_files must contain the following variables: FamID, ID, dadID, momID, sex, affected, DA1, DA2')
  }

  if (any(is.na(ped_files$ID))) {
    stop('ID contains missing values.  Please ensure all individuals have a valid ID.')
  }

  if (!is.logical(ped_files$affected)) {
    stop('In ped_files: expecting affected to be a logical variable.')
  }

}
