#' Returns the row and columns locations of mutations for a person
#'
#' @param person_index The column location of the person
#' @param person_genos The genotype vector for the person
#'
#' @return A list containing two vectors, the first contains the row locations of mutations, the second contains the column locations of mutations.
#' @keywords internal
#'
get_SMindex_by_person <- function(person_index, person_genos){
  #j_1 stores the positions where the
  #individual with person_index carries an SNV
  #on their 1st haplotype
  #
  #These will become the non-zero column entries
  #in our sparseMatrix (i.e. argument j to sparseMatrix)
  j_1 <- which(person_genos %in% c("1|0", "1|1"))

  #j_2 stores the positions where the
  #individual with person_index carries an SNV
  #on their 2nd haplotype
  #
  #These will also become non-zero column entries
  #in our sparseMatrix (i.e. argument j to sparseMatrix)
  j_2 <- which(person_genos %in% c("0|1", "1|1"))


  # create the vector of row positions for person with person_index,
  # which will be supplied to sparseMatrix
  i_pos <- c(rep((2*person_index - 1), length(j_1)),
             rep((2*person_index), length(j_2)))

  return(list(i_pos = i_pos,
              j_pos = c(j_1, j_2)))

}


#' Convert genotypes to sparse Matrix format.
#'
#' @param genotypes The genotypes returned by read.vcfR
#'
#' @return A sparseMatrix
#' @export
#'
genos2sparseMatrix <- function(genotypes){
  # Get the row and column location of each mutation
  # person-by-person (i.e. row-by-row) from the genotypes
  # matrix returned by read.vcfR
  #
  # NOTE: index_by_person, below, is a list of lists
  # the first item in the first list are the row postions of SNVs for the first person
  # and the second item in the first list are the column postions of SNVs for the first person
  index_by_person <- lapply(1:ncol(genotypes), function(x){
    get_SMindex_by_person(x, genotypes[, x])
  })

  #input the row and column data into the sparseMatrix
  SM_format <- sparseMatrix(i = unlist(lapply(index_by_person, `[[`, 1)),
                            j = unlist(lapply(index_by_person, `[[`, 2)),
                            x = rep(1, length(unlist(lapply(index_by_person, `[[`, 1)))))

  row.names(SM_format) = rep(colnames(genotypes), each = 2)

  return(SM_format)
}
