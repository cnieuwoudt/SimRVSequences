# NOTE:
# must temporarily export sim_gameteInheritacne and reduce_to_evnets

simRecomb_test <- function(num_reps, chr_start, chr_stop, chr_center,
                           int_widths, int_start_pos,
                           conditional_gene_drop = FALSE,
                           burn_in = 1000, gamma_params = c(2.63, 2.63/0.5)){


  #sim_gameteInheritance simulates recombination and inheritance for every
  #chrom in chrom_map
  if (conditional_gene_drop){
    #note have to simuate each chromosome separately,
    #which is why you need the outer lappy
    loop_gams <- lapply(1:num_reps, function(x){
      sim_gameteInheritance(RV_locus = data.frame(chr = 1,
                                                  pos = sample(1:chr_stop, size = 1)), #simulating randomly in each trial
                            parent_RValleles = sample(x = c(0, 0, 1), size = 2), #setting for conditional gene-drop
                            offspring_RVstatus = sample(x = c(0, 1), size = 1),#again, setting for conditional gene-drop
                            chrom_map = data.frame(chrom     = 1,
                                                   start_pos = chr_start,
                                                   end_pos   = chr_stop,
                                                   center    = chr_center),
                            allele_IDs = c(2, 3),
                            burn_in, gamma_params)
    })

  } else {
    #note have to simuate each chromosome separately,
    #which is why you need the outer lappy
    loop_gams <- lapply(1:num_reps, function(x){
      sim_gameteInheritance(RV_locus = data.frame(chr = 1,
                                                  pos = round((chr_stop-chr_start)/2)),  #it doesn't matter what we set this to here
                            parent_RValleles = c(0, 0), #setting for conditional gene-drop
                            offspring_RVstatus = c(0), #again, setting for conditional gene-drop
                            chrom_map = data.frame(chrom     = 1,
                                                   start_pos = chr_start,
                                                   end_pos   = chr_stop,
                                                   center    = chr_center),
                            allele_IDs = c(2, 3),
                            burn_in, gamma_params)
    })

  }



  #The
  loop_events <- lapply(1:num_reps, function(x){
    reduce_to_events(gamete_haplo = as.numeric(loop_gams[[x]]$haplotypes[[1]]),
                     chias_locations = as.numeric(loop_gams[[x]]$cross_locations[[1]]))
  })


  #count the number of crossovers inside the interval of interest.
  #NOTE: since we are only supplying the interval of interest to cut,
  #all crossovers not inside the interval will be assigned NA.
  nChias_inInt <- sapply(1:num_reps, function(x){
    sum(!is.na(cut(x = loop_events[[x]], breaks = c(int_start_pos, int_start_pos+int_widths), include.lowest = TRUE)))
  })


  return(nChias_inInt)
}




#lets start with a simple run
#compare expected to observed recombination frequency over the entire length of the chromosome

chrom_lengths <- c(50, 100, 150, 200)
int_widths <- c(1, 5, 10, 50, 100, 150, 200)
int_loc <- c("top", "center", "bottom")

test_vals = expand.grid(c_length = chrom_lengths, int_width = int_widths, int_loc = int_loc)
head(test_vals)
test_vals = test_vals[test_vals$int_width <= test_vals$c_length, ]
test_vals$chr_starts = rep(0, nrow(test_vals))
test_vals$chr_stops = test_vals$c_length
test_vals$chr_centers = test_vals$c_length*0.2
test_vals$int_start_pos = ifelse(test_vals$int_loc == "top", 0,
                                 ifelse(test_vals$int_loc == "center",
                                        (test_vals$c_length/2)-(test_vals$int_width/2),
                                        test_vals$c_length - test_vals$int_width))
test_vals$expected = 0.5*tanh(2*(test_vals$int_width)/100)
row.names(test_vals) = NULL
head(test_vals, n = 10)

#---------------------------------#
# With Recombinaiton interference #
# regular gene-drop model         #
#---------------------------------#
nreps = 100000
test_CI_gda <- list()
pb <- txtProgressBar(min = 0, max = nrow(test_vals), style = 3)
set.seed(9090)
for(i in 1:nrow(test_vals)){
  test_CI_gda[[i]] <- simRecomb_test(num_reps = nreps,
                                  chr_start = test_vals$chr_starts[i],
                                  chr_stop = test_vals$chr_stops[i],
                                  chr_center = test_vals$chr_centers[i],
                                  int_widths = test_vals$int_width[i],
                                  int_start_pos = test_vals$int_start_pos[i])
  setTxtProgressBar(pb, i)
}

close(pb)

#combine results
test_vals_CI_gda = test_vals
test_vals_CI_gda$RC_observed = sapply(test_CI_gda, function(x){
  #x %% 2 will returns 1 if number of recombiations are odd
  #this function returns the proportion of odd number crossovers in interval
  mean(x %% 2)
})
test_vals_CI_gda$RC_count = sapply(test_CI_gda, function(x){
  #x %% 2 will returns 1 if number of recombiations are odd
  #this function returns the proportion of odd number crossovers in interval
  sum(x %% 2)
})
test_vals_CI_gda$pvalue_binom = sapply(1:length(test_CI_gda), function(y){
  binom.test(x = test_vals_CI_gda$RC_count[y],
             p = test_vals_CI_gda$expected[y],
             n = length(test_CI_gda[[y]]), alternative = "two.sided")$p.value
})

test_vals_CI_gda$ave_numRC = sapply(test_CI_gda, mean)



test_vals_CI_gda

#gda = gene drop algorithm
test_recomb_gda = test_vals_CI_gda[order(test_vals_CI_gda$c_length, test_vals_CI_gda$int_loc, test_vals_CI_gda$int_width, decreasing = FALSE), ]
row.names(test_recomb_gda) = NULL

save(test_recomb_gda, file = "C:/Data/SimSeqValidation/test_recomb_gda.rdata", compress='xz')
load(file = "C:/Data/SimSeqValidation/test_recomb_gda.rdata")

#---------------------------------#
# With Recombinaiton interference #
# conditional gene-drop model     #
#---------------------------------#
nreps = 100000
test_CI_cgda <- list()

set.seed(411711) #round 1
set.seed(345897) #round 2
pb <- txtProgressBar(min = 0, max = nrow(test_vals), style = 3)
for(i in 1:nrow(test_vals)){
  test_CI_cgda[[i]] <- simRecomb_test(num_reps = nreps,
                                       chr_start = test_vals$chr_starts[i],
                                       chr_stop = test_vals$chr_stops[i],
                                       chr_center = test_vals$chr_centers[i],
                                       int_widths = test_vals$int_width[i],
                                       int_start_pos = test_vals$int_start_pos[i],
                                       conditional_gene_drop = TRUE)
  setTxtProgressBar(pb, i)
}
close(pb)

#combine results
test_vals_CI_cgda = test_vals
test_vals_CI_cgda$RC_observed = sapply(test_CI_cgda, function(x){
  #x %% 2 will returns 1 if number of recombiations are odd
  #this function returns the proportion of odd number crossovers in interval
  mean(x %% 2)
})

test_vals_CI_cgda$RC_count = sapply(test_CI_cgda, function(x){
  #x %% 2 will returns 1 if number of recombiations are odd
  #this function returns the proportion of odd number crossovers in interval
  sum(x %% 2)
})

test_vals_CI_cgda$pvalue_binom = sapply(1:length(test_CI_cgda), function(y){
  binom.test(x = sum(test_CI_cgda[[y]] %% 2),
             p = test_vals_CI_cgda$expected[y],
             n = length(test_CI_cgda[[y]]), alternative = "two.sided")$p.value
})

test_vals_CI_cgda$ave_numRC = sapply(test_CI_cgda, mean)


#cgda = conditional gene drop algorithm
test_recomb_cgda = test_vals_CI_cgda[order(test_vals_CI_cgda$c_length, test_vals_CI_cgda$int_loc, test_vals_CI_cgda$int_width, decreasing = FALSE), ]
row.names(test_recomb_cgda) = NULL
test_recomb_cgda

save(test_recomb_cgda, file = "C:/Data/SimSeqValidation/test_recomb_cgda.rdata", compress='xz')

#for round 2, save as test_recomb_cgda2
test_recomb_cgda2 = test_recomb_cgda
save(test_recomb_cgda2, file = "C:/Data/SimSeqValidation/test_recomb_cgda2.rdata", compress='xz')

load(file = "C:/Data/SimSeqValidation/test_recomb_cgda.rdata")
