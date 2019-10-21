#simulate data for a sample of pedigrees
# DO NOT RUN #
# instead load simulated pedigrees
library(SimRVPedigree)

#Create hazard object from AgeSpecific_Hazards data
data(AgeSpecific_Hazards)
my_HR = hazard(AgeSpecific_Hazards)

# load libraries needed to simulate pedigrees in parallel.
library(doParallel)
library(doRNG)

npeds <- 100    #set the number of pedigrees to generate

cl <- makeCluster(8)   # create cluster
registerDoParallel(cl) # register cluster


#simulate a sample of five pedigrees using foreach
s_peds = foreach(i = seq(npeds), .combine = rbind,
                 .packages = c("SimRVPedigree"),
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

stopCluster(cl) #shut down cluster


# save(time_peds, file="data-raw/time_peds.rdata", compress='xz')
# save(time_peds, file="C:/Users/cnieuwoudt/Google Drive/SeqSim/R code/time_peds.rdata", compress='xz')


#import pedigrees
time_peds <- read.csv("C:/Users/cnieuwoudt/Google Drive/SeqSim/R code/time_peds.csv")


#load pedigrees for timing
load(file="data-raw/time_peds.rdata")
class(time_peds)

#-------------------------#
# Timing table - Vignette #
#-------------------------#
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

# sout$Mutations$is_CRV <- FALSE
# sout$Mutations$is_CRV[sout$Mutations$marker == "1_170552964"] <- TRUE


# create different sized pedigree samples to evaluate timing #
#modify FamIDs so that we can combine with original 100
time_peds_dummy200 <- time_peds
time_peds_dummy200$FamID <- time_peds_dummy200$FamID + 100
time_peds_dummy200 <- rbind(time_peds, time_peds_dummy200)

#data frame of times to populate for each set of settings
time_dat <- data.frame(time_ped10 = rep(NA, 10),
                       time_ped50 = rep(NA, 10),
                       time_ped100 = rep(NA, 10),
                       time_ped150 = rep(NA, 10),
                       time_ped200 = rep(NA, 10))

simStudy_settings <- data.frame(affOnly = c(TRUE, TRUE, FALSE, FALSE),
                                remWild = c(TRUE, FALSE, TRUE, FALSE))

num_peds = c(10, 50, 100, 150, 200)

#one entry in time_res for each set of settings in simStudy_settings
#the jth column in time_dat is for the jth setting in num_peds
#the ith row in time_dat is for the ith run (of 10)
time_res = list()
time_res[[1]] = time_dat
time_res[[2]] = time_dat
time_res[[3]] = time_dat
time_res[[4]] = time_dat

set.seed(2304789)
for(i in 1:length(num_peds)){
  for(k in 1:10){
    #sample appropriate number of pedigrees from sample of 200 pedigrees
    tpeds = time_peds_dummy200[time_peds_dummy200$FamID %in%
                                 sample(time_peds_dummy200$FamID, size = num_peds[[i]]), ]

    #time sim_RVstudy for each set of settings
    for(j in 1:nrow(simStudy_settings)){
      a = Sys.time()
      study_seq <- sim_RVstudy(ped_files = tpeds,
                               SNV_data = sout,
                               affected_only = simStudy_settings$affOnly[[j]],
                               remove_wild = simStudy_settings$remWild[[j]])
      b = Sys.time()

      time_res[[j]][k, i] = difftime(b, a, units = "mins")

    }
  }
}

time_res1 = time_res[[1]]
time_res2 = time_res[[2]]
time_res3 = time_res[[3]]
time_res4 = time_res[[4]]

# save(time_res1, file="C:/Users/cnieuwoudt/Google Drive/SeqSim/R code/time_res1.rdata", compress='xz')
# save(time_res2, file="C:/Users/cnieuwoudt/Google Drive/SeqSim/R code/time_res2.rdata", compress='xz')
# save(time_res3, file="C:/Users/cnieuwoudt/Google Drive/SeqSim/R code/time_res3.rdata", compress='xz')
# save(time_res4, file="C:/Users/cnieuwoudt/Google Drive/SeqSim/R code/time_res4.rdata", compress='xz')


load(file="C:/Users/cnieuwoudt/Google Drive/SeqSim/R code/time_res1.rdata", compress='xz')
load(file="C:/Users/cnieuwoudt/Google Drive/SeqSim/R code/time_res2.rdata", compress='xz')
load(file="C:/Users/cnieuwoudt/Google Drive/SeqSim/R code/time_res3.rdata", compress='xz')
load(file="C:/Users/cnieuwoudt/Google Drive/SeqSim/R code/time_res4.rdata", compress='xz')


my.cols = c(SFUred = rgb(red = 166/225, green = 25/225, blue = 46/225), #SFUred
            SFUblack = rgb(red = 84/225, green = 88/225, blue = 90/225), #SFUblack
            SFUblue = rgb(red = 0/225, green = 112/225, blue = 150/225), #SFUblue
            SFUgold = rgb(red = 193/225, green = 160/225, blue = 30/225),#SFUgold
            SFUgreen = rgb(red = 112/225, green = 140/225, blue = 50/225),
            SFUteal = rgb(red = 0, green = 128/225, blue = 128/225),
            SFUmidnight = rgb(red = 25/225, green = 25/225, blue = 112/225)) #SFUgreen

win.graph(h = 5, w = 6)
plot(x = num_peds, y = apply(time_res1, 2, mean),
     xlab = "Number of Pedigrees",
     xaxt = "n",
     ylab = "Time (in minutes)",
     main = "Timing by Number of Pedigree and Settings",
     col = my.cols[[1]],
     pch = 19,
     lwd = 1,
     type = "b",
     ylim = c(0, 7))
#create legend to detail population
legend("topleft", #title = "settings",
       legend = c("affected_only = TRUE, remove_wild = TRUE",
                  "affected_only = TRUE, remove_wild = FALSE",
                  "affected_only = FALSE, remove_wild = TRUE",
                  "affected_only = FALSE, remove_wild = FALSE"),
       col = c(my.cols[1], my.cols[3], my.cols[6], my.cols[7]), lwd = 2)

axis(side = 1, at = num_peds,
     labels = as.character(num_peds),
     cex = 1.2, line = 0, lwd = 1)

points(x = num_peds, y = apply(time_res2, 2, mean),
       col = my.cols[[3]],
       pch = 15,
     type = "b")
points(x = num_peds, y = apply(time_res3, 2, mean),
       col = my.cols[[6]],
       pch = 17,
       type = "b")
points(x = num_peds, y = apply(time_res4, 2, mean),
       col = my.cols[[7]],
       pch = 18,
       type = "b")

