#vignette timing pedigrees and SNV data
#simulate data for a sample of pedigrees
# DO NOT RUN #
# instead load simulated pedigrees
library(SimRVPedigree)

# #Create hazard object from AgeSpecific_Hazards data
# data(AgeSpecific_Hazards)
# my_HR = hazard(AgeSpecific_Hazards)
#
# # load libraries needed to simulate pedigrees in parallel.
# library(doParallel)
# library(doRNG)
#
# npeds <- 100    #set the number of pedigrees to generate
#
# cl <- makeCluster(8)   # create cluster
# registerDoParallel(cl) # register cluster
#
#
# #simulate a sample of five pedigrees using foreach
# time_peds = foreach(i = seq(npeds), .combine = rbind,
#                     .packages = c("SimRVPedigree"),
#                     .options.RNG = 844090518
# ) %dorng% {
#   # Simulate pedigrees ascertained for at least three disese-affected individuals,
#   # according to the age-specific hazard rates in the `AgeSpecific_Hazards` data
#   # set, ascertained from 1980 to 2018, with seed-founder birth year spanning
#   # from 1900 to 1920, stop year set to 2018, and with genetic relative-risk 50.
#   sim_RVped(hazard_rates = my_HR,
#             GRR = 50, FamID = i,
#             RVfounder = TRUE,
#             founder_byears = c(1900, 1920),
#             ascertain_span = c(1980, 2018),
#             stop_year = 2018,
#             recall_probs = c(1, 0.5, 0),
#             num_affected = 3)[[2]]}
#
# stopCluster(cl) #shut down cluster
#
#
# #save simulated pedigrees
# save(time_peds, file="data-raw/time_peds.rdata", compress='xz')

#load pedigrees for timing
load(file="data-raw/time_peds.rdata")
class(time_peds)


library(SimRVSequences)
data(hg_exons)
data(hg_apopPath)
data(study_peds)

#create mutation and genome set from slim output
a = Sys.time()
sout = read_slim(file_path = "C:/Data/Slim/SlimFINALout.txt",
                 keep_maf = 0.01,
                 recomb_map = create_slimMap(hg_exons),
                 pathway_df = hg_apopPath)
b = Sys.time()
difftime(b, a, units = "mins")

class(sout)

#choose pool of causal rare variants
dim(sout$Mutations)
table(sout$Mutations$afreq[sout$Mutations$pathwaySNV])

set.seed(411)
crv_rows <- sample(which(sout$Mutations$pathwaySNV & sout$Mutations$afreq == 5e-05),
                   size = 20, replace = FALSE)
sout$Mutations$is_CRV <- FALSE
sout$Mutations$is_CRV[crv_rows] <- TRUE

sum(sout$Mutations$afreq[sout$Mutations$is_CRV])

