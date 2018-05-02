library(testthat)
context("reconstruct_fromHaplotype")

data(EXmut)

#Chromosome info we will need
#store the mutation data for chromosome 1
C1mut = EXmut[which(EXmut$chrom == 1), ]

#reformat so that mutations occur at consecutive positions
#this helps with testing.
C1mut$position = seq(1:nrow(C1mut))

#this will be the reduced chromosome map for chr 1
C1map = data.frame(chrom = 1,
                   start_pos = min(C1mut$position),
                   end_pos = max(C1mut$position))

test_that("reconstruct_fromHaplotype returns an identical sequence when no crossovers occur", {

  #initialize event_loc to an empty list, which is what we
  #observe when no crossovers are simulated
  event_loc = c()

  #create inherited haplotype matrix,
  #this will contain 1 row amd one column
  #the only entry will be the haplotype to inherit
  inherit_hap = matrix(sample(1:2, size = 1),
                       nrow = 1)

  #matrix with 1 stored for every mutation from paternal haplotype
  #and 2 for every mutation from maternal haplotype
  parent_hap = as.data.frame(matrix(c(rep(1, nrow(C1mut)), rep(2, nrow(C1mut))),
                                    nrow = 2, byrow = TRUE))


  #reconstruct the offspring sequence given the locations of crossovers
  #and the participating haplotype sequences.
  inherited_genomic_seq <- reconstruct_fromHaplotype(parental_genotypes = parent_hap,
                                                     Cmarker_map = C1mut,
                                                     inherited_haplotype = inherit_hap,
                                                     chiasmata_locations = event_loc,
                                                     REDchrom_map = C1map)
  #since no crossovers occur every position should be equal
  #to the value stored in inherit_hap
  expect_true(all(inherited_genomic_seq == inherit_hap))
})


test_that("reconstruct_fromHaplotype contains the correct number of swaps", {

  #sample crossover locations over the area of interest in this chromosome
  event_loc = sort(runif(sample(1:6, size = 1), min(C1mut$position), max(C1mut$position)))


  #sample gametes for each crossover.
  #Note that some crossover events will be trivial if this gamete did
  #not pariticipate in the crossover.
  inherit_hap = matrix(sample(1:2,
                              size = (length(event_loc) + 1),
                              replace = TRUE),
                       nrow = 1)

  #matrix with 1 stored for every mutation from paternal haplotype
  #and 2 for every mutation from maternal haplotype
  parent_hap = as.data.frame(matrix(c(rep(1, nrow(C1mut)), rep(2, nrow(C1mut))),
                                    nrow = 2, byrow = TRUE))

  inherited_genomic_seq <- reconstruct_fromHaplotype(parental_genotypes = parent_hap,
                                                     Cmarker_map = C1mut,
                                                     inherited_haplotype = inherit_hap,
                                                     chiasmata_locations = event_loc,
                                                     REDchrom_map = C1map)


  inherit_hap
  event_loc
  inherited_genomic_seq
  rle(as.numeric(inherit_hap))$values
  rle(as.numeric(inherited_genomic_seq))$values


  expect_equal(rle(as.numeric(inherit_hap))$values, rle(as.numeric(inherited_genomic_seq))$values)

})
