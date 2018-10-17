#' Check SNV_map for possible issues
#'
#' \strong{For internal use.}
#'
#' @inheritParams sim_RVstudy
#' @export
#'
check_SNV_map <- function(SNV_map){
  #check to see if SNV_map contains the column information we expect
  # and check to see if we have any missing values.

  ## Check colID variable
  if (!c("colID") %in% colnames(SNV_map)) {
    stop("The variable 'colID' is missing from SNV_map.")
  }
  if (any(is.na(SNV_map$colID))) {
    stop("The variable 'colID' in the SNV_map dataset contains missing values.")
  }

  ## Check chrom variable
  if (!c("chrom") %in% colnames(SNV_map)) {
    stop("The variable 'chrom' is missing from SNV_map.")
  }
  if (any(is.na(SNV_map$chrom))) {
    stop("The variable 'chrom' in the SNV_map dataset contains missing values.")
  }

  ## Check position variable
  if (!c("position") %in% colnames(SNV_map)) {
    stop("The variable 'position' is missing from SNV_map.")
  }

  if (any(is.na(SNV_map$position))) {
    stop("The variable 'position' in the SNV_map dataset contains missing values.")
  }

  ## Check marker variable
  if (!c("marker") %in% colnames(SNV_map)) {
    stop("The variable 'marker' is missing from SNV_map.")
  }
  if (any(is.na(SNV_map$marker))) {
    stop("The variable 'marker' in the SNV_map dataset contains missing values.")
  }

  #when is_CRV is specified, check to see that it is TRUE for
  #at least one marker, and that there are haplotypes that carry
  #the SNV
  if (!is.null(SNV_map$is_CRV)) {
    if (sum(SNV_map$is_CRV) == 0) {
      stop("In SNV_map: is_CRV is FALSE for all markers.")
    }
  }

}

#' Checks individual pedigrees for proper format.
#'
#' \strong{For internal use.} Checks individual pedigrees for formatting (i.e. mom/dad properly specified, etc.)
#'
#' @param ped_file data.frame The pedigree.
#'
#' @export
#'
check_ped <- function(ped_file){

  #gather all mom and dad IDs for non-founders
  moms <- unique(ped_file$momID[!is.na(ped_file$momID)])
  dads <- unique(ped_file$dadID[!is.na(ped_file$dadID)])

  #check to see if moms are female and dads are male
  if (any(ped_file$sex[which(ped_file$ID %in% moms)] != 1) |
      any(ped_file$sex[which(ped_file$ID %in% dads)] != 0)){

    wrong_sex <- c(ped_file$ID[which(ped_file$sex[which(ped_file$ID %in% dads)] != 0)],
                   ped_file$ID[which(ped_file$sex[which(ped_file$ID %in% moms)] != 1)])

    stop(paste0('Sex improperly specifed ID: ', sep = '', wrong_sex,
                '. \n Please ensure that for males: sex = 0; and for females: sex = 1.'))
  }

  #check to see that the moms and dads are actually included in the pedigree
  #that is check to see that the IDs of moms and dads are properly specified
  if (any(!moms %in% ped_file$ID) | any(!dads %in% ped_file$ID)) {

    wrong_par <- c(ped_file$ID[which(ped_file$momID == moms[which(!moms %in% ped_file$ID)])],
                   ped_file$ID[which(ped_file$dadID == dads[which(!dads %in% ped_file$ID)])])

    stop(paste0('ID: ', sep = '', wrong_par,
                '.  Non-founders must have a mother and a father. Founders have neither.'))
  }

  #check to see that both parents are missing for founders
  if (any(!is.na(ped_file$momID[is.na(ped_file$dadID)])) |
      any(!is.na(ped_file$dadID[is.na(ped_file$momID)]))) {
    stop("Non-founders must have both a mother and a father, while founders have neither.")
  }


  #check to see that when founders do not introduce a cRV at the familial
  #disease locus that offspring do not carry it.
  #Essentially, this is a check for de novo mutations between founders
  #and non-founders. This will NOT catch a

  #dadIDs of the non-founders who inherited a cRV from dad
  inhrt_fromDad <- unique(ped_file$dadID[ped_file$DA1 == 1 & !is.na(ped_file$dadID)])

  #momIDs of the non-founders who inherited a cRV from mom
  inhrt_fromMom <- unique(ped_file$momID[ped_file$DA2 == 1 & !is.na(ped_file$dadID)])

  if (length(inhrt_fromDad) > 0) {
    #count the number of cRVs held by each dad from whom a non-founder inherited a cRV
    dadRVcounts <- sapply(inhrt_fromDad, function(x){
      sum(ped_file[ped_file$ID == x, c("DA1", "DA2")])
    })

    if (any(dadRVcounts == rep(0, length(inhrt_fromDad)))) {
      stop("\n Detecting de novo mutation. \n Please check that variable DA1, which represents the \n paternally inherited allele, is properly specified in ped_files.")
    }
  }

  if (length(inhrt_fromMom) > 0) {
    #count the number of cRVs held by each dad from whom a non-founder inherited a cRV
    momRVcounts <- sapply(inhrt_fromMom, function(x){
      sum(ped_file[ped_file$ID == x, c("DA1", "DA2")])
    })

    if (any(momRVcounts == rep(0, length(inhrt_fromMom)))) {
      stop("\n Detecting de novo mutation. \n Please check that variable DA2, which represents the \n maternally inherited allele, is properly specified in ped_files.")
    }
  }

  #check to make sure that only 1 founder introduced the causal rare variant
  if (sum(ped_file[is.na(ped_file$dadID), c("DA1", "DA2")]) > 1) {
    stop("Assumption violated.",
         "\n Reformat ped_files so that only one cRV is introduced per pedigree.")
  }
}


#' Checks ped_files for expected info and format.
#'
#' \strong{For internal use.}
#'
#' @inheritParams sim_RVstudy
#'
#' @export
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
