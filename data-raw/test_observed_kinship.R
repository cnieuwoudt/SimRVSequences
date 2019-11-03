# NOTE:
# to run this script you must temporarily export sim_seq

#create a pedigree with only two first cousins
D2_ped = data.frame(FamID = rep(1, 8),
                    ID = 1:8,
                    sex = c(0, 1, 0, 1, 0, 1, 0, 0),
                    dadID = c(NA, NA, NA, 1, NA, 1, 3, 5),
                    momID = c(NA, NA, NA, 2, NA, 2, 4, 6),
                    affected = c(F, F, F, F, F, F, T, T),
                    Gen = c(1, 1, 2, 2, 2, 2, 3, 3),
                    DA1 = rep(0, 8),
                    DA2 = rep(0, 8))
plot(new.ped(D2_ped))


#create some dummy founder data to drop down the pedigree
my_founder_genos = list()
my_founder_genos[[1]] = matrix(c(rep(1, 1),
                                 rep(2, 1),
                                 rep(3, 1),
                                 rep(4, 1),
                                 rep(5, 1),
                                 rep(6, 1),
                                 rep(7, 1),
                                 rep(8, 1)),
                               nrow = 8, ncol = 1,
                               byrow = TRUE)
my_founder_genos[[2]] = c(1, 1, 2, 2, 3, 3, 5, 5)

#initialize chrom map (simulatiing over approximately 250 cM)
my_chrom_map = data.frame(chrom = 1, start_pos = 1, end_pos = 250)

#store all combinations of individuals to compare
sim_set <- data.frame(relationship = c("self", "parent-offspring", "full-sibs", "avuncular", "grandparent-grandchild", "first-cousins"),
                      ID1 = c(6, 4, 4, 6, 1, 7),
                      ID2 = c(6, 7, 6, 7, 8, 8),
                      expected_kin = c(0.5, 0.25, 0.25, 0.125, 0.125, 0.0625),
                      observed = rep(NA, 6))

sim_kin = list()

nreps = 10000
set.seed(23957)
for(k in 1:nrow(sim_set)){
  #initialize list to hold results
  sim_kin[[k]] = rep(NA, nreps)
  pb <- txtProgressBar(min = 0, max = nreps, style = 3)
  for(i in 1:nreps){
    #randomly sample the location of the SNV to drop down the pedigree
    my_SNV_map = data.frame(colID = 1,
                            chrom = 1,
                            position = runif(1, 1, 250)) #ampling position in cM
    my_SNV_map$marker = paste0(my_SNV_map$chrom, "_", my_SNV_map$position)

    #sample the sequence data for the pedigree
    fam_seq = sim_seq(ped_file = D2_ped, founder_genos = my_founder_genos,
                      SNV_map = my_SNV_map, chrom_map = my_chrom_map, RV_marker = '1_1',
                      burn_in = 1000, gamma_params = c(2.63, 2.63/0.5))

    #pull the alleles for each pair of individuals
    ID1_alleles <- fam_seq$ped_genos[fam_seq$geno_map$ID == sim_set$ID1[k], ]
    ID2_alleles <- fam_seq$ped_genos[fam_seq$geno_map$ID == sim_set$ID2[k], ]

    #Randomly sample 1 allele from each individual and record 1 if
    #the sampled alleles are the same
    sim_kin[[k]][i] = (sample(ID1_alleles, size = 1) == sample(ID2_alleles, size = 1))*1

    setTxtProgressBar(pb, i)
  }
  close(pb)
}

sim_kin_df = do.call(cbind, sim_kin)
save(sim_kin_df, file = "C:/Data/SimSeqValidation/sim_kin_df.rdata", compress='xz')
load(file = "C:/Data/SimSeqValidation/sim_kin_df.rdata")

sim_set$observed <- apply(sim_kin_df, 2, mean)
sim_set$std <- apply(sim_kin_df, 2, function(x){mean(x)*(1-mean(x))/length(x)})


sim_set$p.val <- sapply(1:nrow(sim_set), function(x){
  binom.test(x = sum(sim_kin_df[, x]),
             p = sim_set$expected_kin[x],
             n = length(sim_kin_df[, x]), alternative = "two.sided")$p.value
})

sim_set
