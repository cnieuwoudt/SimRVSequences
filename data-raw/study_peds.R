# load the SimRVPedigree library
library(SimRVPedigree)

#Create hazard object from AgeSpecific_Hazards data
data(AgeSpecific_Hazards)
my_HR = hazard(AgeSpecific_Hazards)

library(doParallel)
library(doRNG)

npeds <- 5    #set the number of pedigrees to generate

cl <- makeCluster(5)   # create cluster
registerDoParallel(cl) # register cluster


start_time <- Sys.time()
study_peds = foreach(i = seq(npeds), .combine = rbind,
                  .packages = c("kinship2", "SimRVPedigree"),
                  .options.RNG = 844090518
) %dorng% {
  # Simulate pedigrees ascertained for at least three disese-affected individuals,
  # according to the age-specific hazard rates in the `AgeSpecific_Hazards` data
  # set, ascertained from 1980 to 2018, with seed-founder birth year spanning
  # from 1900 to 1920, stop year set to 2018, and with genetic relative-risk 50.
  sim_RVped(hazard_rates = my_HR,
            GRR = 50, FamID = i,
            RVfounder = TRUE,
            founder_byears = c(1900, 1920),
            ascertain_span = c(1980, 2018),
            stop_year = 2018,
            recall_probs = c(1, 0.5, 0),
            num_affected = 3)[[2]]}
stop_time <- Sys.time()
stopCluster(cl)

difftime(stop_time, start_time, units = "mins")

# par(mfrow = c(2, 3))
# plot(study_peds[study_peds$FamID == 1,])
# plot(study_peds[study_peds$FamID == 2,])
# plot(study_peds[study_peds$FamID == 3,])
# plot(study_peds[study_peds$FamID == 4,])
# plot(study_peds[study_peds$FamID == 5,])

save(study_peds, file="data/study_peds.rdata", compress='xz')
