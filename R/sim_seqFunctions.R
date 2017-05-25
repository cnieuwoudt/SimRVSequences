#' Construct offspring sequence from parental allele vector
#'
#' For internal use.
#'
#' May need to reasses this code, particularly at line 27.
#' COMMENT THIS CODE ASAP
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

    #determine which chiasmata inherited haplotype participated in
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
      offspring_seq[, which(Cmarker_map$position >= start_switch & Cmarker_map$position < end_switch)] <-
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
#' library(SimRVPedigree)
#' library(kinship2)
#' #Read in age-specific hazards
#' data(EgPeds)
#'
#' ex_RVped <- EgPeds[which(EgPeds$FamID == 4), ]
#'
#' # Define pedigree object for trimmed pedigree, i.e, pedigree with
#' #proband selected and relatives trimmed
#' TrimRVped <- pedigree(id = ex_RVped$ID,
#'                       dadid = ex_RVped$dad_id,
#'                       momid = ex_RVped$mom_id,
#'                       sex = (ex_RVped$gender + 1),
#'                       affected = cbind(Affected = ex_RVped$affected,
#'                                        Proband = ex_RVped$proband,
#'                                        RV_status = ex_RVped$DA1 +
#'                                                    ex_RVped$DA2),
#'                       famid = ex_RVped$FamID)['4']
#'
#'  plot(TrimRVped)
#'  pedigree.legend(TrimRVped, location = "topleft",  radius = 0.25)
#'
#' my_chrom_map = data.frame(chrom     = c(1),
#'                           start_pos = c(0),
#'                           end_pos   = c(250),
#'                           center = c(50))
#' my_chrom_map
#' my_RV_marker <- "1_50"
#'
#' link_map <- data.frame(chromosome = c(1, 1, 1, 1, 1, 1, 1),
#'                        position = c(20, 50, 100, 125, 175, 200, 250))
#' link_map$marker <- paste0(link_map$chromosome, sep = "_", link_map$position)
#' link_map <- link_map[, c(3, 1, 2)]
#' link_map
#'
#'
#' founder_seq2 <- matrix(rep(letters[1:(2*nrow(link_map))], length(which(is.na(ex_RVped$dad_id)))),
#'                        nrow = 2*length(which(is.na(ex_RVped$dad_id))),
#'                        byrow = T)
#' colnames(founder_seq2) = as.character(link_map$marker)
#'
#' founder_seq2[1, which(colnames(founder_seq2) == my_RV_marker)] <- 'X'
#' founder_seq2 <- as.data.frame(founder_seq2)
#' founder_seq2$ID = rep(ex_RVped$ID[which(is.na(ex_RVped$dad_id))], each = 2)
#' founder_seq2
#'
#' set.seed(6)
#' ped_seq <- sim_RVseq(ped_file = ex_RVped,
#'                      founder_genos = founder_seq2,
#'                      marker_map = link_map,
#'                      chrom_map = my_chrom_map,
#'                      RV_marker = my_RV_marker)
#' ped_seq
#'
#' set.seed(6)
#' system.time(sim_RVseq(ped_file = ex_RVped,
#'                       founder_genos = founder_seq2,
#'                       marker_map = link_map,
#'                       chrom_map = my_chrom_map,
#'                       RV_marker = my_RV_marker))
#'
sim_RVseq <- function(ped_file, founder_genos,
                      marker_map, chrom_map, RV_marker,
                      burn_in = 1000, gamma_params = c(2.63, 2.63/0.5)){

  #Get parent/offspring information
  #i.e. for each offspring find RV_status,
  #parent IDs, and parent alleles at RV locus
  PO_info <- get_parOffInfo(ped_file)
  PO_info <- PO_info[order(PO_info$Gen, PO_info$offspring_ID),]

  ped_genos <- founder_genos

  #for each offspring simulate transmission of parental data
  for (i in 1:nrow(PO_info)) {
    loop_gams <- sim_gameteInheritance(RV_locus = marker_map[which(marker_map$marker == RV_marker), c(2:3)],
                                       parent_RValleles = PO_info[i, c(6, 7)],
                                       offspring_RVstatus = PO_info[i, 5],
                                       chrom_map,
                                       allele_IDs = c(1, 2),
                                       burn_in, gamma_params)

    loop_seq <- lapply(c(1:nrow(chrom_map)),
                       function(x){
                         reconstruct_fromHaplotype(parental_genotypes =
                                                     ped_genos[which(ped_genos$ID == PO_info[i, 4]),
                                                                        which(marker_map$chromosome == chrom_map$chrom[x])],
                                                   Cmarker_map = marker_map[which(marker_map$chromosome == chrom_map$chrom[x]),],
                                                   inherited_haplotype = loop_gams$haplotypes[[x]],
                                                   chiasmata_locations = loop_gams$cross_locations[[x]],
                                                   REDchrom_map = chrom_map[x, ])
                       })
    loop_seq[[nrow(chrom_map) + 1]] = data.frame(ID = PO_info[i, 1])

    ped_genos <- rbind(ped_genos, do.call("cbind", loop_seq))
  }

  #ped_genos <- ped_genos[order(ped_genos$ID),]
  ped_genos$FamID <- ped_file$FamID[1]
  ped_genos$FamRV <- RV_marker
  #rownames(ped_genos) = NULL

  return(ped_genos)
}
