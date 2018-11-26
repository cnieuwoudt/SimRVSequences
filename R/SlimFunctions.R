#' Create recombination map
#'
#' Create a recombination map for use with SLiM (Haller and Messer).
#'
#' One of SLiM's many features is that it can simulate recombination hotspots using a user-specified recombination map.  This recombination map may be utilized to simulate mutations over unlinked regions (i.e. in different chromosomes) or in linked but non-contiguous regions (i.e in exon-only data).  The \code{create_slimMap} function may be used to generate the recombination map required by SLiM to simulate exon-only SNV data.
#'
#' We expect that \code{exon_df} does not contain any overlapping segments.  Prior to supplying the exon data to \code{create_slimMap} users must combine overlapping exons into a single observation.  The \code{\link{combine_exons}} function may be used to accomplish this task.
#'
#' The argument \code{exon_df} must contain the following variables:
#' \tabular{lll}{
#' \strong{name} \tab \strong{type} \tab \strong{description} \cr
#' \code{chrom} \tab numeric \tab chromosome identification number\cr
#' \code{exonStart} \tab numeric \tab the position of the first base pair in the exon\cr
#' \code{exonStop} \tab numeric \tab the position of the last base pair in the exon\cr
#' }
#'
#' The data frame returned by \code{create_slimMap} contains variables required by SLiM to simulate exon-only data.  Additionally, the returned data frame also includes variables that are required to re-map mutations to their correct positions when importing SLiM data to \code{R}.  The variables contained in the returned data frame are described as follows.
#' \describe{
#' \item{\code{chrom}}{The chromosome number.}
#' \item{\code{segLength}}{The length of the segment in base pairs.  We assume that segments contain the positions listed in \code{exonStart} and \code{exonEnd}.  Therefore, for a combined exon segment, \code{segLength} is calculated as \code{exonEnd - exonStart + 1}.}
#' \item{\code{recRate}}{The per-site per-generation recombination rate.  Following Harris, segments between exons on the same chromosome are simulated as a single base pair with \code{rec_rate} equal to recombination rate multiplied by the number of base pairs in the segment.  For each chromosome, a single site is created between the last exon on the previous chromosome and the first exon of the current chromosome.  This site will have recombination rate 0.5 to accommodate unlinked chromosomes.}
#' \item{\code{mutRate}}{The per-site per-generation mutation rate.  Since we are interested in exon-only data, the mutation rate outside exons is set to zero.}
#' \item{\code{exon}}{A logical variable that is \code{TRUE} if the segment is an exon and \code{FALSE} otherwise.}
#' \item{\code{simDist}}{The simulated exon length, in base pairs. When \code{exon = TRUE}, \code{simDist = segLength}; however, when \code{exon = FALSE}, \code{simDist = 1} since segments between exons on the same chromosome are simulated as a single base pair.}
#' \item{\code{endPos}}{The simulated end position, in base pairs, of the segment.}
#' }
#'
#' Only three of the variables returned by \code{create_slimMap} are required by SLiM to simulate exon-only data: \code{recRate}, \code{mutRate}, and \code{endPos}.  The other variables seen in the output above are used by our \code{\link{read_slim}} function to re-map mutations to their correct positions when importing SLiM data to \code{R}.
#'
#' Please note SLiM is written in a scripting language called Eidos. Unlike an \code{R} array, the first position in an Eidos array is 0.  Therefore, we must shift the variable \code{endPos} forward 1 unit before supplying this data to SLiM. See example.
#'
#' @param exon_df Data frame. A data frame that contains the positions of each exon to simulate.  This data frame must contain the variables \code{chrom}, \code{exonStart}, and \code{exonEnd}.  See details.
#' @param mutation_rate Numeric.  The per-site per-generation mutation rate, assumed to be constant across the genome. By default, \code{mutation_rate= 1E-8}, as in Harris.
#' @param recomb_rate Numeric.  The per-site per-generation mutation rate, assumed to be constant across the genome. By default, \code{mutation_rate= 1E-8}, as in Harris.
#'
#' @return A recombination map that may be used in conjunction with SLiM (Haller and Messer 2017).  See details and example.
#'
#' @references Benjamin Haller and Phillip W. Messer (2017). \emph{Slim 2: Flexible, interactive forward genetic simulations}. Molecular Biology and Evolution; 34(1), pp. 230-240.
#' @references Kelly Harris and Rasmus Nielsen (2016). \emph{The genetic cost of neanderthal introgression}. Genetics, 203(2): pp. 881-891.
#'
#' @export
#'
#' @seealso \code{\link{combine_exons}}
#'
#' @examples
#' #load hg_exons data
#' data(hg_exons)
#'
#' #since the exons in hg_exons have already been combined into
#' #overlapping exons, we supply hg_exons to create_slimMap
#' slimMap <- create_slimMap(hg_exons)
#' head(slimMap)
#'
#' # restrict output to the variables required by SLiM
#' slimMap <- slimMap[, c("recRate", "mutRate", "endPos")]
#'
#' # shift endPos up by one unit
#' slimMap$endPos <- slimMap$endPos - 1
#'
#' # print first four rows of slimMap
#' head(slimMap, n = 4)
#'
create_slimMap <- function(exon_df, mutation_rate = 1E-8, recomb_rate = 1E-8){
  #split into dataset for each chromosome
  bychr <- lapply(sort(unique(exon_df$chrom)), function(x){
    exon_df[exon_df$chrom == x, ]
  })

  #compile data on introns
  int_dist <- lapply(bychr, function(x){
    data.frame(chrom = x$chrom,
               segLength = c((x$exonStart[1] - 1), (x$exonStart[-1] - x$exonEnd[-nrow(x)] - 1)),
               no = c(1:nrow(x)),
               recRate = recomb_rate*c(0.5/recomb_rate, (x$exonStart[-1] - x$exonEnd[-nrow(x)] - 1)),
               mutRate = rep(0, nrow(x)),
               type = rep("intron", nrow(x)),
               exon = rep(FALSE, nrow(x)),
               simDist = rep(1, nrow(x)))
  })
  #set the recombination rate to zero for the first chromosome
  #Note: the first recombination rate is 0.5 for all other
  #chromosomes so that successive chromosomes are unlinked
  int_dist[[1]]$recRate[1] = 0

  #compile data on exons
  ex_dist <- lapply(bychr, function(x){
    data.frame(chrom = x$chrom,
               segLength = c(x$exonEnd - x$exonStart + 1),
               no = c(1:nrow(x)),
               recRate = rep(recomb_rate, nrow(x)),
               mutRate = rep(mutation_rate, nrow(x)),
               type = rep("exon", nrow(x)),
               exon = rep(TRUE, nrow(x)),
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
  return(rc_map[, -c(3, 6)])
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
    recomb_map[recomb_map$chrom == x, ]
  })

  #subset by introns for each chromosome
  bychr_int <- lapply(bychr, function(x){
    x[!x$exon, ]
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
    mutationDF[mutationDF$chrom == x, ]
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
    bychr_int[[i]]$cumDist <- cumsum(bychr_int[[i]]$segLength)

    mut_by_chrom[[i]]$position <- mut_by_chrom[[i]]$position +
      bychr_int[[i]]$cumDist[mut_by_chrom[[i]]$ex_num] -
      mut_by_chrom[[i]]$ex_num

  }

  mut_dat <- do.call(rbind, mut_by_chrom)
  return(mut_dat[, c(1:4)])
}

#' Import SLiM data to R
#'
#'To import SLiM data into \code{R}, we provide the \code{read_slim} function, which has been tested for SLiM versions 2.0-3.1.  Presently, the \code{read_slim} function is only appropriate for single-nucleotide variant (SNV) data produced by SLiM's outputFull() method.  We do not support output in MS or VCF data format, i.e. produced by outputVCFsample() or outputMSSample() in SLiM.
#'
#' In addition to reducing the size of the data, the argument \code{keep_maf} has practicable applicability.  In family-based studies, common SNVs are generally filtered out prior to analysis.  Users who intend to study common variants in addition to rare variants may need to run chromosome specific analyses to allow for allocation of large data sets in \code{R}.
#'
#' The argument \code{recomb_map} is used to remap mutations to their actual locations and chromosomes.  This is necessary when data has been simulated over non-contiguous regions such as exon-only data.  If \code{\link{create_slimMap}} was used to create the recombination map for SLiM, simply supply the output of \code{create_slimMap} to \code{recomb_map}.  If \code{recomb_map} is not provided we assume that the SNV data has been simulated over a contiguous segment starting with the first base pair on chromosome 1.
#'
#' The data frame \code{pathway_df} allows users to identify SNVs located within a pathway of interest.  When supplied, we expect that \code{pathwayDF} does not contain any overlapping segments.  \emph{All overlapping exons in \code{pathway_df} MUST be combined into a single observation.  Users may combine overlapping exons with the \code{\link{combine_exons}} function.}
#'
#' The \code{read_slim} function returns a list containing two items:
#' \enumerate{
#' \item \code{Haplotypes} A sparse matrix of class dgCMatrix. The columns in {Haplotypes} represent distinct SNVs, while the rows repesent individual haplotypes. We note that this matrix contains two rows of data for each diploid individual in the population: one row for the maternally ihnherited haplotype and the other for the paternally inherited haplotype.
#' \item \code{Mutations} A data frame cataloging SNVs in \code{Haplotypes}. The variables in the \code{Mutations} data set are described as follows:
#' \tabular{ll}{
#' \code{colID} \tab Associates the rows, i.e. SNVs, in \code{Mutations} to the columns of \code{Haplotypes}. \cr
#' \code{chrom} \tab The chromosome that the SNV resides on. \cr
#' \code{position} \tab The position of the SNV in base pairs. \cr
#' \code{afreq} \tab The derived allele frequency of the SNV. \cr
#' \code{marker} \tab A unique character identifier for the SNV.\cr
#' \code{pathwaySNV} \tab Identifies SNVs located within the pathway of interest as \code{TRUE}. Note that this variable is omitted when users do not supply \code{pathway_df} to \code{read_slim}. \cr
#' }}
#'
#'
#' @param file_path character.  The file path of the .txt output file created by the outputFull() method in SLiM.
#' @param keep_maf numeric. The largest allele frequency for retained SNVs, by default \code{keep_maf = 0.01}.  All variants with allele frequency greater than \code{keep_maf} will be removed. Please note, removing common variants is recommended for large data sets due to the limitations of data allocation in R. See details.
#' @param recomb_map data frame. (Optional) A recombination map of the same format as the data frame returned by \code{\link{create_slimMap}}. See details.
#' @param pathway_df data frame. (Optional) A data frame that contains the positions for each exon in a pathway of interest.  See details.
#'
#' @return  A list containing:
#' @return \item{\code{Haplotypes} }{A sparse matrix of haplotypes. See details.}
#' @return \item{\code{Mutations}}{A data frame cataloging SNVs in \code{Haplotypes}. See details.}
#' @importFrom Matrix sparseMatrix
#' @export
#'
#' @references Haller, B., Messer, P. W. (2017). \emph{Slim 2: Flexible, interactive forward genetic simulations}. Molecular Biology and Evolution; 34(1), pp. 230-240.
#' @references Douglas Bates and Martin Maechler (2018). \strong{Matrix: Sparse and Dense Matrix Classes and Methods}.
#' \emph{R package version 1.2-14}. https://CRAN.R-project.org/package=Matrix
#'
#' @seealso \code{\link{create_slimMap}}, \code{\link{combine_exons}}, \code{\link{dgCMatrix-class}}
#'
#' @examples
#'
#' # If create_slimMap was used to create the recombination
#' # map for SLiM from the hg_exons data set, and if the .txt file
#' # produced by SLiM's outputFull() method is saved as "slimOut.txt"
#' # in the current working directory we import "slimOut.txt" to R
#' # using the following command.
#'
#' \dontrun{
#' s_out <- read_slim(file_path  = "slimOut.txt",
#'                    recomb_map = create_slimMap(hg_exons))
#' }
#'
read_slim <- function(file_path, keep_maf = 0.01,
                      recomb_map = NULL,
                      pathway_df = NULL){
  #NOTE: Time to read file ~19 secs
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
  #   lists which genomes belong to which individual, redundant since this
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

  popCount <- as.numeric(unlist(strsplit(exDat[PopHead + 1], split = " "))[2])

  #-----------#
  # Mutations #
  #-----------#
  #NOTE: time to create Mutations data set ~4 secs
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

  #calculate the derived allele frequency
  MutData$afreq <- MutData$prevalence/(2*popCount)

  #order Mutation dataset by tempID, so that (later) we can order
  #the mutations in each haplotypes by increasing genomic position
  MutData <- MutData[order(MutData$tempID), ]

  #Identify the future sparseMatrix column ID of variants with <= keep_maf.
  #Variants with large maf are assigned column ID 0, i.e. thrown away.
  #Variants with sufficiently small maf are assigned increasing colIDs.
  #These are like new tempIDs, they do not reflect genomic position, yet.
  MutData$colID <- cumsum(MutData$afreq <= keep_maf)*(MutData$afreq <= keep_maf)

  #Using the identified colID, create dataframe of rare mutations only
  RareMutData <- MutData[MutData$colID > 0, ]

  #-----------#
  # Genotypes #
  #-----------#
  #NOTE: jpos is the computationally expensive task (uses strsplit)
  #this chuck takes ~2.2 minutes to run (for complete exon-only data)
  #This is an improvement from the old time: ~ 6 mins
  print("Creating Sparse Haplotypes Matrix")
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
  RareMutData <- RareMutData[, c(5, 2, 3, 4)]

  #-----------------------------#
  # Re-code Identical Mutations #
  #-----------------------------#
  # Assuming that all mutations are of the same type, different mutations at
  # the same site are actually identical mutations from different lineages.
  # For simplicity, we recode these mutations so that they are only cataloged
  # once.
  # TALK TO JINKO BEFORE YOU DO THIS
  if (any(duplicated(RareMutData$position))) {
    print("Recoding Identical Mutations")
    com_id_muts <- combine_identicalmutations(mutmap = RareMutData,
                                              hapmat = GenoData,
                                              pCount = popCount,
                                              keep_maf)
    RareMutData <- com_id_muts[[1]]
    GenoData <- com_id_muts[[2]]
  } else {
    RareMutData <- RareMutData[, -3]
  }


  #------------------#
  # Re-map Mutations #
  #------------------#
  if (!is.null(recomb_map)) {
    print("Remapping Mutations")
    RareMutData <- reMap_mutations(mutationDF = RareMutData,
                                   recomb_map)
  } else {
    RareMutData$chrom <- 1
  }

  row.names(RareMutData) = NULL

  RareMutData$marker <- paste0(RareMutData$chrom, sep = "_", RareMutData$position)

  #reduce RareMutData, yet again, to the columns we actually need
  #really should clean this up soon
  RareMutData <- RareMutData[, c(1, 4, 2, 3, 5)]


  #----------------------#
  # Identify Pathway RVs #
  #----------------------#
  #if pathway data has been supplied, identify pathway SNVs
  if (!is.null(pathway_df)) {
    print("Identifying Pathway SNVs")
    RareMutData <- identify_pathwaySNVs(markerDF = RareMutData,
                                        pathwayDF = pathway_df)
  }


  return(list(Haplotypes = GenoData, Mutations = RareMutData))
}



#' Determine i and j positions of mutations for sparse matrix
#'
#' @param mutString character. String containing mutations
#' @param rarePos numeric vector. position of variations
#'
#' @return data.frame with x and y positions of mutations
#' @export
extract_tempIDs <- function(mutString, rarePos){
  tids <- as.numeric(strsplit(mutString, split = " ", fixed = TRUE)[[1]][-c(1:2)]) + 1
  #Subset colID by tids (tempID) and retain the colIDs that are non-zero,
  #i.e. the rare varaints
  rarePos[tids][rarePos[tids] > 0]
}


#' Combine identical mutations
#'
#' Assuming that all mutations are of the same type (in SLIM simulation), different mutations at the same site are actually identical mutations from different lineages. This function re-code these mutations so that they are only cataloged once.
#'
#' @param mutmap data.frame The SNV_map with identical mutations
#' @param hapmat sparseMatrix The sparseMatrix of haplotypes
#' @param pCount the population count
#' @param keep_maf numeric. The largest allele frequency for retained SNVs, by default \code{keep_maf = 0.01}.  All variants with allele frequency greater than \code{keep_maf} will be removed.
#' @importFrom Matrix rowSums
#'
#' @return  A list containing:
#' @return \item{\code{hapmat} }{A sparse matrix of haplotypes. See details.}
#' @return \item{\code{mutmap}}{A data frame cataloging SNVs in \code{hapmap}.}
#' @export
#'
combine_identicalmutations <- function(mutmap, hapmat, pCount, keep_maf){

  # identify the positions at which identical mutations
  # from different lineages exist.
  im_pos = unique(mutmap$position[duplicated(mutmap$position)])

  # find the column locations for the identical mutations
  col_loc <- lapply(im_pos, function(x){
    mutmap$colID[which(mutmap$position == x)]
  })

  #find the combined SNV data
  comb_mut <- lapply(col_loc, function(x){
    rowSums(hapmat[, x])
  })

  keep_SNVcol <- unlist(lapply(col_loc, function(x){
    x[1]
  }))

  #replace with combined haplotype data
  hapmat[, keep_SNVcol] <- do.call(cbind, comb_mut)

  #determine the superfluous columns and remove them, since we have
  #already accounted for them by combining the columns for this SNV
  remove_cols <- unlist(lapply(col_loc, function(x){
    x[-1]
  }))

  #remove the columns of the identical SNVs
  hapmat <- hapmat[, -remove_cols]

  #combine the SNVs in mutmap, i.e. calculate prevalence and derived allele frequency
  #and relable column ID
  for (i in 1:length(col_loc)) {
    mutmap$prevalence[col_loc[[i]][1]] <- sum(mutmap$prevalence[col_loc[[i]]])
    mutmap$afreq[col_loc[[i]][1]] <- mutmap$prevalence[col_loc[[i]][1]]/(2*pCount)
  }

  mutmap <- mutmap[-remove_cols, ]
  mutmap$colID <- 1:nrow(mutmap)

  if (any(mutmap$afreq > keep_maf)) {
    keep_cols <- which(mutmap$afreq <= keep_maf)
    hapmat <- hapmat[, keep_cols]
    mutmap <- mutmap[keep_cols, ]
    mutmap$colID <- 1:nrow(mutmap)
  }

  return(list(mutmap[, c("colID", "position", "afreq")], hapmat))
}
