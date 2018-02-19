#' Read arlequin project file
#'
#' Read arlequin project file with DataType SNP only, and extract the SNPs as a dataframe well as their postions.
#'
#' @param file character.  The name of the arlequin project file, note that the file name must include the .arp file extension.
#'
#' @return  A list containing:
#' @return \item{\code{SNPs} }{A dataframe of SNPs, with columns for each loci in Postions and rows for individual haplotypes.}
#' @return \item{\code{Positions} }{A numeric list of the postions, in bp, of the SNPs contained in SNPs}
#' @export
#'
#' @examples
#' \dontrun{
#' setwd("C:/Users/cnieuwoudt/Documents/fsc26_win64/1PopDNAtest")
#' exDat = read.arp("1PopDNAtest_1_1.arp")
#' head(exDat[[1]])
#' head(exDat[[2]])
#' }
#'
read.arp <- function(file){
  exDat = readLines(file)

  #---------------#
  # SNP Positions #
  #---------------#
  #Assumes that postions are listed at line 20 - will need a way to check this
  #extract postions of SNPs and store as numeric vector
  Positions = unlist(strsplit(exDat[20], split = ", "))

  # Positions[1]
  # class(Positions[1])
  #remove # before first position
  Positions[1] <- paste0(unlist(strsplit(Positions[1], split = ""))[-1],
                         collapse = "")
  #convert to numeric variable
  Positions <- as.numeric(Positions)
  # head(Positions)

  #------#
  # SNPs #
  #------#
  #Store positons of first and last haplotype rows
  #Here I an using the structure of the arlequin file
  #to determine these positions
  firstHaplo <- which(exDat == "\t\tSampleData= {") + 1
  lastHaplo <- which(exDat == "}") - 2

  SNPs <- do.call(rbind, lapply(exDat[firstHaplo:lastHaplo], function(x){
    as.numeric(unlist(strsplit(unlist(strsplit(x, split = " "))[2], split = "")))
  }))

  return(list(SNPs = SNPs, Positions = Positions))
}
