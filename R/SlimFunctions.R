#' Create recombination map
#'
#' Create a recombination map for use with SLiM 2.0 (Messer)
#'
#' In addition to allowing users to specify recombination hotspots, the recombination map provided to SLiM can be used to simulate mutations over unlinked regions (i.e. in different chromosomes) or in linked but non-contiguous regions (i.e in exon-only data).  The \code{\link{create_slimMap}} function may be used to create a recombination map to exon-only data with SLiM.
#'
#' The argument \code{exon_df} must contain the following variables:
#' \tabular{lll}{
#' \strong{name} \tab \strong{type} \tab \strong{description} \cr
#' \code{chrom} \tab numeric \tab chromosome identification number\cr
#' \code{exonStart} \tab numeric \tab the exon's starting position, in bp\cr
#' \code{exonStop} \tab numeric \tab the exon's ending position, in bp\cr
#' }
#'
#'
#' The returned data frame, \code{rc_map}, contains the information required by slim as well as information that is used to remap mutations to their correct positions after simulation.  The variables \code{recRate}, \code{mutRate}, and \code{endPos} are the only variables that are required by SLiM.
#'
#' @param exon_df Data frame. A data frame cataloguing exon data. See details.
#' @param mutation_rate Numeric.  The per site per generation mutation rate, by default \code{mutation_rate} = 1E-8.
#' @param recomb_rate Numeric.  The per site per generation mutation rate, by default \code{recomb_rate} = 1E-8.
#'
#' @return rc_map A recombination map that may be used in conjunction with SLiM 2.0 (cite Messer), provided that the end position is shifted forward by one position.  See details and example.
#' @export
#'
#' @examples
#' data(hg_exons)
#' smap <- create_slimMap(hg_exons)
#' head(smap)
#'
#' # NOTE: The first position in an eidos array begins at zero, not one.
#' # Therefore, users must shift endPos up one unit before supplying to Slim
#'
#' smap$endPos <- smap$endPos - 1
#' head(smap)
#'
create_slimMap <- function(exon_df, mutation_rate = 1E-8, recomb_rate = 1E-8){
  #split into dataset for each chromosome
  bychr <- lapply(sort(unique(exon_df$chrom)), function(x){
    subset(exon_df, chrom == x)
  })

  #compile data on introns
  int_dist <- lapply(bychr, function(x){
    data.frame(chrom = x$chrom,
               dist = c((x$exonStart[1] - 1), (x$exonStart[-1] - x$exonEnd[-nrow(x)] - 1)),
               no = c(1:nrow(x)),
               recRate = recomb_rate*c(0.5/recomb_rate, (x$exonStart[-1] - x$exonEnd[-nrow(x)] - 1)),
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
               recRate = rep(recomb_rate, nrow(x)),
               mutRate = rep(mutation_rate, nrow(x)),
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
#' Intended for internal use
#'
#' @param mutationDF The Mutation data frame returned by reMap_mutations
#' @param recomb_map The recombination map provided to slim
#'
#' @return A re-mapped mutation data frame
#' @export
reMap_mutations <- function(mutationDF, recomb_map){
  #split into data for different chromosomes, because it
  #makes myhead hurt to think this as 1 chomosome
  bychr <- lapply(sort(unique(recomb_map$chrom)), function(x){
    subset(recomb_map, chrom == x)
  })

  #subset by introns for each chromosome
  bychr_int <- lapply(bychr, function(x){
    subset(x, type == "intron")
  })

  #get chrom starts and stops
  chr_start <- unlist(lapply(bychr, function(x){ x$endPos[1] }))
  chr_end <- unlist(lapply(bychr, function(x){ x$endPos[nrow(x)] }))

  #determine which chromosome each mutation falls on
  mutationDF$chrom <- cut(mutationDF$position,
                          breaks = c(1, chr_end),
                          labels = FALSE)

  #renumber for chromosomes included in recomb_map
  mutationDF$chrom <- sort(unique(recomb_map$chrom))[mutationDF$chrom]


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

  mut_dat <- do.call(rbind, mut_by_chrom)
  return(mut_dat[, c(1:4)])
}

#' Read SLiM 2.0 Output
#'
#' Extract SNP data from SLiM output.
#'
#' The \code{read_slim} function is used to extract SNP (single nucleotide polymorphism) data from the text file produced by SLiM's outputFull() method. Currently, we do not support output in MS or VCF format (i.e. output produced by SLiM's outputMSSample() or outputVCFSample() methods).
#'
#' In addition to allowing users to specify recombination hotspots, the recombination map provided to SLiM can be used to simulate mutations over unlinked regions (i.e. in different chromosomes) or in linked but non-contiguous regions (i.e in exon-only data).  The \code{\link{create_slimMap}} function may be used to create a recombination map to exon-only data with SLiM.  In this case, simply supply the data frame returned by \code{create_slimMap} to the \code{recomb_map} argument of \code{read_slim} to map mutations to their actual locations and/or their appropriate chromosomes.  If \code{create_slimMap} was not used to create the recombination map used with SLiM 2.0, users must take care to ensure that \code{recomb_map} is of the same form as the output produced by \code{\link{create_slimMap}}.
#'
#' The first item returned by \code{read_slim} is a data frame named \code{Mutations}, which catalouges the SNP IDs, genomic positions, chromosomal positions, and allele frequency.  The variable \code{colID} references the SNP's ID which is also its column position in the sparse genotype matrix. The variable \code{position} is the genomic position of the SNP (in bp), and \code{afreq} is the allele frequency of the SNP.  If \code{recomb_map} is provided and contains multiple chromosomes, \code{chrom} identifies the chromosome that the mutation resides on.  When \code{recomb_map} is not provided, we assume that all mutations reside on the first chromosome so that \code{chrom} is 1 for every mutation.
#'
#' The second item returned is a sparse matrix named \code{Genomes}.  This matrix contains two rows for each individual in the (haploid) population. That is, each row is one of the haplotypes for a single individual in the population.  For example, the first individual's maternally and paternally inherited haplotypes are stored in rows one and two, respectively.  In general, the \eqn{i^{th}} individual's maternally inherited haplotype is stored in row \eqn{2i-1} and the paternally inherited haplotype is stored in row \eqn{2i}. (ALTERNATIVE DESC: The third and fourth rows contain haplotypes for the second individual, fifth and sixth rows contain haplotypes for the second individual, and so on. )
#'
#'
#' NOTE TO SELF:
#' \itemize{
#'  \item For file extension internal check use: importFrom tools file_ext
#'  \item If possible, add test: two rows for each person in genotypes matrix
#' }
#'
#' @param file_path character.  The file path of the SLiM 2.0 output file.
#' @param keep_maf numeric. The largest allele frequency for retained SNPs.  All variants with allele frequency greater than \code{keep_maf} will be removed.  Please note, removing common variants is recommended for large datasets due to the limitations of data allocation in \code{R}.
#' @param recomb_map data frame. The recombination map provided to SLiM 2.0, see details.
#'
#' @return  A list containing:
#' @return \item{\code{Mutations} }{A dataframe containing SNP information, see details..}
#' @return \item{\code{Genomes} }{A sparse matrix of haplotypes, see details.}
#' @importFrom Matrix sparseMatrix
#' @export
#'
#' @examples
#' \dontrun{
#' sout = read_slim(file_path = "C:/Data/Slim/SlimFullFix_out.txt", keep_maf = 0.01, recomb_map = create_slimMap(hg_exons))
#' sout = read_slim(file_path = "C:/Data/Slim/SLiMtest_output.txt", keep_maf = 0.01)
#' }
read_slim <- function(file_path, keep_maf = 0.01, recomb_map = NULL){
  print("Reading Slim File")
  exDat = readLines(file_path)

  #The default output for Slim  is a .txt file with the following headings:
  # - OUT: Contains number of generations, type (i.e. "A" for autosome), and
  #   the file name
  #
  # - Populations:
  #   next line is pop description, "p1", pop size, and
  #   type (i.e. "H" hermaphroditic)
  #
  # - Mutations:
  #   each mutation gets a separate line with (1) temp id, (2) permanent id,
  #   (3) mutation type,  (4) base position, (5) selection coefficient, 0 for
  #   neutral model, (6) dominance coefficient (0.5 for neutral), (7) id of
  #   sub-population in which the mutation arose, (8) the generation in which
  #   the mutation arose, and (9) the prevalence of the mutation.
  #
  # - Individuals:
  #   lists which genomes belong to which individual, highly redundant as this
  #   information is also listed at the beginning of each genome.
  #
  # - Genomes:
  #   each line lists one of the genomes for an individual, i.e. individual 1's
  #   genomes are contained in lines 1 & 2, individual 2's genomes are contained
  #   in lines 3 & 4, etc. Each line begins with the genome id, followed by the
  #   type (i.e. "A" for autosome), and a list of mutations.  The mutations are
  #   identified by the temporary ID (tempID) specified in the Mutations output.

  #find heading locations
  PopHead <- which(exDat == "Populations:")
  MutHead <- which(exDat == "Mutations:")
  IndHead <- which(exDat == "Individuals:")
  GenHead <- which(exDat == "Genomes:")

  popCount <- as.numeric(unlist(strsplit(exDat[PopHead+1], split = " "))[2])
  #-----------#
  # Mutations #
  #-----------#
  print("Creating Mutations Dataset")
  #create mutation dataset from slim's Mutation output
  #only retaining the tempID, position, and prevalence of each mutation
  MutData <- do.call(rbind,
                     lapply((MutHead + 1):(IndHead - 1), function(x){
                       as.numeric(unlist(strsplit(exDat[x], split = " "))[c(1, 4, 9)])
                     })
  )

  MutData <- as.data.frame(MutData)
  colnames(MutData) <- c("tempID", "position", "prevalence")

  #add 1 to temp ID so that we can easily associate mutations
  #to columns by default slim's first tempID is 0, not 1.
  MutData$tempID <- MutData$tempID + 1
  #First position in slim is 0, not 1
  MutData$position <- MutData$position + 1

  #calculate the minor allele frequency
  MutData$afreq <- MutData$prevalence/(2*popCount)

  #order Mutation dataset by tempID, so that (later) we can order
  #the mutations in each haplotypes by increasing genomic position
  MutData <- MutData[order(MutData$tempID), ]
  MutData$colID <- cumsum(MutData$afreq <= keep_maf)*(MutData$afreq <= keep_maf)

  #create dataframe of rare mutations only
  RareMutData <- subset(MutData, colID > 0)

  #-----------#
  # Genotypes #
  #-----------#
  print("Creating Sparse Genotypes Matrix")
  #determine future row and column position of each mutation listed in genomes
  #row will correspond to person, column will correspond to the tempID of the
  #mutation
  jpos <- lapply(1:(2*popCount), function(x){
    extract_tempIDs(mutString = exDat[GenHead + x],
                    rarePos = MutData$colID)
  })

  ipos <- lapply(1:length(jpos), function(x){
    rep(x, length(jpos[[x]]))
  })

  #create sparse matrix containing mutations(columns) for each individual(row)
  GenoData <- sparseMatrix(i = unlist(ipos),
                           j = unlist(jpos),
                           x = rep(1, length(unlist(jpos))),
                           dims = c(2*popCount, nrow(RareMutData)))

  #order by genomic postion of rare mutation
  GenoData <- GenoData[, order(RareMutData$position)]
  RareMutData <- RareMutData[order(RareMutData$position),]

  #Re-format tempID so that it corresponds to the column
  # (in GenoData) that the mutation is stored in
  RareMutData$colID <- 1:nrow(RareMutData)
  RareMutData <- RareMutData[, c(5, 2, 4)]

  if (!is.null(recomb_map)) {
    print("Remapping Mutations")
    RareMutData <- reMap_mutations(mutationDF = RareMutData,
                                   recomb_map)
  } else {
    RareMutData$chrom <- 1
  }

  row.names(RareMutData) = NULL

  return(list(Mutations = RareMutData[, c(1, 4, 2, 3)], Genomes = GenoData))
}



#' Determine i and j positions of mutations for sparse matrix
#'
#' @param mutString character. String containing mutations
#' @param indPos numeric. row number of individual
#' @param rarePos numeric vector. position of variations
#'
#' @return data.frame with x and y positions of mutations
#' @keywords internal
extract_tempIDs <- function(mutString, rarePos){
  tids <- as.numeric(unlist(strsplit(mutString, split = " "))[-c(1:2)]) + 1
  #Subset colID by tids (tempID) and retain the colIDs that are non-zero,
  #i.e. the rare varaints
  rarePos[tids][rarePos[tids] > 0]
}
