library(testthat)
context("sim_RVstudy")

data("study_peds")
data("EXhaps")
data("EXmuts")

#SMALL TOY DATASETS FOR TESTING
toy_haps <- sparseMatrix(i = seq(1:10), j = seq(1:10), x = rep(1, 10))
toy_muts <- data.frame(colID = seq(1:10),
                       chrom = rep(1, 10),
                       position = round(seq(1001, 2000001, length.out = 10)*1000),
                       pathwaySNV = sample(x = c(TRUE, FALSE), size = 10,
                                           replace = TRUE, prob = c(0.2, 0.8)),
                       is_CRV = sample(c(rep(FALSE, 9), TRUE), size = 10))

toy_muts$marker = paste0(toy_muts$chrom, sep = "_", toy_muts$position)


test_that("rows of haplo_map are equal to rows ped_haplos", {

  study_seq <- sim_RVstudy(ped_files = study_peds,
                           SNV_map = toy_muts,
                           haplos = toy_haps,
                           remove_wild = FALSE,
                           affected_only = TRUE)

  expect_equal(nrow(study_seq$ped_haplos), nrow(study_seq$haplo_map))
})


test_that("warning if is_CRV is missing from SNV_map", {

  expect_warning(sim_RVstudy(ped_files = study_peds,
                             SNV_map = toy_muts[, -5],
                             haplos = toy_haps,
                             remove_wild = FALSE,
                             affected_only = TRUE))
})

