#' Construct offspring sequence from parental allele vector
#'
#' For internal use.
#'
#' @param parental_genotypes The parental genotype sequence information.
#' @param Cmarker_map Data.frame. Must contain three columns with: column 1: marker names, must be listed in the same order as in the founder genotype file, column 2: the chromosomal position of the marker, column 3: the position of the marker in cM.
#' @param inherited_haplotype The inherited haplotype sequence.
#' @param chiasmata_locations The chiasmata locations.
#' @param REDchrom_map Data.frame.  The chromosome map, reduced to the chromosome in question.
#'
#' @return offspring_seq
#' @export
#'
#'
#'
#'
reconstruct_fromHaplotype <- function(parental_genotypes,
                                      Cmarker_map,
                                      inherited_haplotype,
                                      chiasmata_locations,
                                      REDchrom_map){

  if (length(chiasmata_locations) > 0){

    #determine which inherited haplotype participated in chiasmata
    switch_alleles_loc <- c(REDchrom_map$start_pos,
                            reduce_to_events(as.numeric(inherited_haplotype), chiasmata_locations),
                            REDchrom_map$end_pos + 1) #1 added here in case of marker at end of chromosome

    #determine first allele ID in inherited_haplotype
    #set the offspring's haplotype sequence to the other (original) parental haplotype
    offspring_seq <- parental_genotypes[ifelse(inherited_haplotype[1, 1] == 1, 2, 1), ]

    #store the first allele ID in the haplotype sequence
    switch_alle <- inherited_haplotype[1, 1]

    #cycle through all swaps
    for(i in 1:(length(switch_alleles_loc) %/% 2)){
      start_switch <- switch_alleles_loc[2*i - 1]
      end_switch <- switch_alleles_loc[2*i]

      #switch alleles between chiasmata
      offspring_seq[which(Cmarker_map$position >= start_switch & Cmarker_map$position < end_switch)] <-
        parental_genotypes[switch_alle, which(Cmarker_map$position >= start_switch & Cmarker_map$position < end_switch)]
    }
  } else {
    offspring_seq <- parental_genotypes[inherited_haplotype[1], ]
  }

  return(offspring_seq)
}

#' Simulate sequence data for a pedigree
#'
#' @inheritParams sim_gameteInheritance
#' @param ped_file Data frame. Must match format of pedigree simulated by sim_RVped
#' @param marker_map Dataframe. Must contain three columns with: column 1: marker names, must be listed in the same order as in the founder genotype file, column 2: the chromosomal position of the marker, column 3: the position of the marker in cM.
#' @param RV_marker character. The marker name of the RV locus.
#' @param founder_genos Dataframe.  A dataframe with rows corresponding to founders, and colums corresponding to markers.  Markers must be listed in same order as \code{marker_map}.
#'
#' @return offspring_sequences
#' @export
#'
#' @examples
#' #FIND SHORT WORKING EXAMPLE
#'
sim_RVseq <- function(ped_file, founder_genos,
                      marker_map, chrom_map, RV_marker,
                      burn_in = 1000, gamma_params = c(2.63, 2.63/0.5)){

  #Get parent/offspring information
  #i.e. for each offspring find RV_status,
  #parent IDs, and parent alleles at RV locus
  PO_info <- get_parOffInfo(ped_file)
  PO_info <- PO_info[order(PO_info$Gen, PO_info$offspring_ID),]

  ped_genos <- founder_genos[[1]]
  ped_geno_IDs <- founder_genos[[2]]


  #for each offspring simulate transmission of parental data
  for (i in 1:nrow(PO_info)) {
    #determine the chromosome number and location of the familial RV locus
    #then store as a dataframe with chrom in the first column
    RVL <- marker_map[which(marker_map$marker == RV_marker),
                      which(colnames(marker_map) %in% c("chrom", "position"))]

    if(colnames(RVL[1]) != "chrom"){
      RVL <- RVL[, c(2, 1)]
    }

    #simulate recombination events for this parent offspring pair
    loop_gams <- sim_gameteInheritance(RV_locus = RVL,
                                       parent_RValleles = PO_info[i, c(6, 7)],
                                       offspring_RVstatus = PO_info[i, 5],
                                       chrom_map,
                                       allele_IDs = c(1, 2),
                                       burn_in, gamma_params)

    #construct offspring's inherited material from this parent
    loop_seq <- lapply(c(1:nrow(chrom_map)),
                       function(x){
                         reconstruct_fromHaplotype(parental_genotypes =
                                                     ped_genos[which(ped_geno_IDs == PO_info[i, 4]),
                                                                        which(marker_map$chrom == chrom_map$chrom[x])],
                                                   Cmarker_map = marker_map[which(marker_map$chrom == chrom_map$chrom[x]),],
                                                   inherited_haplotype = loop_gams$haplotypes[[x]],
                                                   chiasmata_locations = loop_gams$cross_locations[[x]],
                                                   REDchrom_map = chrom_map[x, ])
                       })

    #append ID for this haplotype to the list of IDs
    ped_geno_IDs <- c(ped_geno_IDs, PO_info[i, 1])

    #append this haplotype to the other familial haplotypes
    ped_genos <- rbind(ped_genos, unlist(loop_seq))
  }

  #create a data.frame to store identifying info
  geno_map <- data.frame(FamID = rep(ped_file$FamID[1], length(ped_geno_IDs)),
                         ID = ped_geno_IDs,
                         affected =  rep(FALSE, length(ped_geno_IDs)))
  #identify affected individuals
  geno_map$affected[geno_map$ID %in% ped_file$ID[ped_file$affected]] <- TRUE

  #Return the genomes matrix and a data.frame continating identifying
  #information for the of IDs to identify the
  #family member to whom
  return(list(ped_genos = ped_genos, geno_map = geno_map))
}
