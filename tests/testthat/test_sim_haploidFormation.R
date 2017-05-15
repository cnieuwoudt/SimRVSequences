library(testthat)
context("sim_hapliodFormation")
test_that("sim_hapliodFormation has ncols = length(chiasmata positions) + 2", {

  my_chrom_map <- data.frame(start = 0, stop = 250, center = 50)
  sim_chias_pos <- sim_chiasmataPositions(my_chrom_map)
  chias_count_BC(sim_chias_pos, my_chrom_map[1,3])
  my_haps <- sim_haploidFormation(num_chiasmata = length(sim_chias_pos),
                                  before_center = chias_count_BC(sim_chias_pos,
                                                                 my_chrom_map[1,3]),
                                  allele_IDs = c(2, 3))
  expect_equal(ncol(my_haps) - 2, length(sim_chias_pos))
})

test_that("sim_hapliodFormation always returns a unique gamete group for each recombined haploid", {

  my_chrom_map <- data.frame(start = 0, stop = 250, center = 50)
  sim_chias_pos <- sim_chiasmataPositions(my_chrom_map)
  chias_count_BC(sim_chias_pos, my_chrom_map[1,3])
  my_haps <- sim_haploidFormation(num_chiasmata = length(sim_chias_pos),
                                  before_center = chias_count_BC(sim_chias_pos,
                                                                 my_chrom_map[1,3]),
                                  allele_IDs = c(2, 3))

  expect_equal(sort(my_haps$gamete_grp), LETTERS[1:4])
})
