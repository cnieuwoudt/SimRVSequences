#' Create recombination map
#'
#' Create a recombination map for use with SLiM 2.0 (Messer)
#'
#' @param exon_df
#'
#' @return rc_map A recombination map that may be used in conjunction with SLiM 2.0 (cite Messer), provided that the end position is shifted forward by one position.  See example.
#' @export
#'
#' @examples
#' data(hg_exons)
#' smap <- create_slimMap(hg_exons)
#' head(smap)
#'
#' NOTE: The first posistion in an eidos array begins at zero,
#' hence, to use with Slim, users must make the following change
#' smap$endPos <- smap$endPos - 1
#' head(smap)
#'
create_slimMap <- function(exon_df){
  #split into dataset for each chromosome
  bychr <- lapply(sort(unique(exon_df$chrom)), function(x){
    subset(exon_df, chrom == x)
  })

  #compile data on introns
  int_dist <- lapply(bychr, function(x){
    data.frame(chrom = x$chrom,
               dist = c((x$exonStart[1] - 1), (x$exonStart[-1] - x$exonEnd[-nrow(x)] - 1)),
               no = c(1:nrow(x)),
               recRate = 1E-8*c(5E+7, (x$exonStart[-1] - x$exonEnd[-nrow(x)] - 1)),
               mutRate = rep(0, nrow(x)),
               type = rep("intron", nrow(x)),
               simDist = rep(1, nrow(x)))
  })
  #set the recombination rate to zero for the first chromosome
  #Note: the first recombination rate is 0.5 for all other
  #chromosomes so that successive chromosomes are unlinked
  int_dist[[1]]$recRate[1] = 0

  #compile data on exons
  ex_dist <- lapply(bychr, function(x){
    data.frame(chrom = x$chrom,
               dist = c(x$exonEnd - x$exonStart + 1),
               no = c(1:nrow(x)),
               recRate = rep(1E-8, nrow(x)),
               mutRate = rep(1E-8, nrow(x)),
               type = rep("exon", nrow(x)),
               simDist = c(x$exonEnd - x$exonStart + 1))
  })

  #combine intron and exon data
  ie_dist <- lapply(1:length(int_dist), function(x){
    rbind(int_dist[[x]], ex_dist[[x]])
  })

  #order by type and number so that they appear in proper order
  ie_dist <- lapply(ie_dist, function(x){
    x[order(x$no, x$type), ]
  })

  #recombine datasets and return
  rc_map <- do.call(rbind, ie_dist)
  rc_map$endPos <- cumsum(rc_map$simDist)
  row.names(rc_map) = NULL
  return(rc_map)
}


#' Create list of haplotype matrices
#'
#' @param markermap The marker map returned by remao
#' @param genomat The full sparse matrix of genotypes
#'
#' @return A list of haplotype matrices
#' @export
create_haplotypeSet <- function(markMap, genoMat){
  lapply(sort(unique(markMap$chrom)), function(x){
    genoMat[, which(markMap$chrom == x)]
  })
}


#' Re-map slim mutations
#'
#' @param mutationDF The Mutation data frame returned by reMap_mutations
#' @param genoMat The sparse genotypes matrix returned by reMap_mutations
#' @param recombMap The recombination map provided to slim
#'
#' @return A re-mapped mutation data frame
#' @export
reMap_mutations <- function(mutationDF, genoMat, recombMap){
  bychr <- lapply(sort(unique(recombMap$chrom)), function(x){
    subset(recombMap, chrom == x)
  })

  #subset by introns for each chromosome
  bychr_int <- lapply(bychr, function(x){
    subset(x, type == "intron")
  })

  #get chrom starts and stops
  chr_start <- unlist(lapply(bychr, function(x){ x$endPos[1] }))
  chr_end <- unlist(lapply(bychr, function(x){ x$endPos[nrow(x)] }))

  #determine which chromosome each mutation falls on
  mutationDF$chrom <- cut(mutationDF$position, breaks = c(1, chr_end), labels = FALSE)
  #renumber for chromosomes included in recombMap
  mutationDF$chrom <- sort(unique(recombMap$chrom))[mutationDF$chrom]


  #create separate data frames for each chromosome
  mut_by_chrom <- lapply(sort(unique(mutationDF$chrom)), function(x){
    subset(mutationDF, chrom == x)
  })


  for(i in 1:length(mut_by_chrom)){
    #shift mutations and map so that the first intron starts at
    #position 1 in each chromosome
    mut_by_chrom[[i]]$position <- mut_by_chrom[[i]]$position - (chr_start - 1)[i]
    bychr[[i]]$newEndPos <- bychr[[i]]$endPos - (chr_start - 1)[i]

    #determine which exon each mutation falls in
    mut_by_chrom[[i]]$ex_num <- (cut(mut_by_chrom[[i]]$position,
                                     breaks = bychr[[i]]$newEndPos,
                                     labels = FALSE) + 1)/2

    #get cumulative intron distance
    bychr_int[[i]]$cumDist <- cumsum(bychr_int[[i]]$dist)

    mut_by_chrom[[i]]$position <- mut_by_chrom[[i]]$position +
      bychr_int[[i]]$cumDist[mut_by_chrom[[i]]$ex_num] -
      mut_by_chrom[[i]]$ex_num

    mut_by_chrom[[i]]$colID <- c(1:nrow(mut_by_chrom[[i]]))

  }

  markMap <- do.call(rbind, mut_by_chrom)

  genoList <- create_haplotypeSet(markMap, genoMat)
  return(list(markMap = markMap, genoList = genoList))
}

