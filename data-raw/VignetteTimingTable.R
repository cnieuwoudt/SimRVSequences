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

#----------------------#
#   Time sim_RVstudy   #
# affected_only = TRUE #
# remove_wild = TRUE   #
#----------------------#
timing1 <- rep(NA, 5)
for(i in 1:5){
  a = Sys.time()
  study_seq <- sim_RVstudy(ped_files = time_peds,
                           SNV_data = sout,
                           affected_only = TRUE,
                           remove_wild = TRUE)
  b = Sys.time()

  timing1[i] = difftime(b, a, units = "mins")
}

# study_seq$SNV_map[duplicated(study_seq$SNV_map$marker), ]
# nrow(study_seq$SNV_map[duplicated(study_seq$SNV_map$marker), ])
# study_seq$SNV_map[study_seq$SNV_map$marker == "1_170552964", ]
# study_seq$SNV_map[study_seq$SNV_map$marker %in% study_seq$SNV_map$marker[duplicated(study_seq$SNV_map$marker)], ]


mean(timing1)
#0.976205537

#----------------------#
#   Time sim_RVstudy   #
# affected_only = TRUE #
# remove_wild = FALSE  #
#----------------------#
timing2 <- rep(NA, 5)
for(i in 1:5){
  a = Sys.time()
  study_seq <- sim_RVstudy(ped_files = time_peds,
                           SNV_map = sout$Mutations,
                           haplos = sout$Haplotypes,
                           affected_only = TRUE,
                           remove_wild = FALSE)
  b = Sys.time()

  timing2[i] = difftime(b, a, units = "mins")
}
mean(timing2)
#1.26697275

#----------------------#
#   Time sim_RVstudy    #
# affected_only = FALSE #
# remove_wild = TRUE    #
#----------------------#
timing3 <- rep(NA, 5)
for(i in 1:5){
  a = Sys.time()
  study_seq <- sim_RVstudy(ped_files = time_peds,
                           SNV_map = sout$Mutations,
                           haplos = sout$Haplotypes,
                           affected_only = FALSE,
                           remove_wild = TRUE)
  b = Sys.time()

  timing3[i] = difftime(b, a, units = "mins")
}
mean(timing3)
#2.0315283

#----------------------#
#   Time sim_RVstudy    #
# affected_only = FALSE #
# remove_wild = FALSE   #
#----------------------#
timing4 <- rep(NA, 5)
for(i in 1:5){
  a = Sys.time()
  study_seq <- sim_RVstudy(ped_files = time_peds,
                           SNV_map = sout$Mutations,
                           haplos = sout$Haplotypes,
                           affected_only = FALSE,
                           remove_wild = FALSE)
  b = Sys.time()

  timing4[i] = difftime(b, a, units = "mins")
}

mean(timing4)
#3.05119482



#--------------------#
# pedigree reduction #
#--------------------#
p1 = new.ped(time_peds)
p1inf <- summary(p1)$family_info

study_seq <- sim_RVstudy(ped_files = time_peds,
                         SNV_map = sout$Mutations,
                         haplos = sout$Haplotypes,
                         affected_only = TRUE,
                         remove_wild = TRUE)

summary(study_seq)
class(study_seq$ped_files)
dim(study_seq$ped_files)
dim(time_peds)

View(study_seq$ped_files)

p2inf <- summary(new.ped(study_seq$ped_files))$family_info
head(p2inf)

head(p1inf)

#total relatives - no reduction
sum(p1inf$totalRelatives)
sum(p2inf$totalRelatives)

#ave fam reduction
mean(p1inf$totalRelatives - p2inf$totalRelatives)


#------------------#
# orginal SNV data #
#------------------#
exDat = readLines("C:/Data/Slim/SlimFINALout.txt")
PopHead <- which(exDat == "Populations:")
MutHead <- which(exDat == "Mutations:")
IndHead <- which(exDat == "Individuals:")
GenHead <- which(exDat == "Genomes:")

length(exDat)
#each row in the Mutations section represents an individual SNV, create the original
#mutations matrix
MutData <- do.call(rbind,
                   lapply((MutHead + 1):(IndHead - 1), function(x){
                     as.numeric(unlist(strsplit(exDat[x], split = " "))[c(1, 4, 9)])
                   })
)

MutData <- as.data.frame(MutData)
colnames(MutData) <- c("tempID", "position", "prevalence")
nrow(MutData)

#store the mutation types
MutData$type <- sapply((MutHead + 1):(IndHead - 1), function(x){
  unlist(strsplit(exDat[x], split = " ", fixed = TRUE))[3]})
head(MutData)

MutData$afreq = MutData$prevalence/20000
m2 = MutData
remove(MutData)



#------------#
# SNV limits #
#------------#
#create ped object
#time_peds = new.ped(time_peds)
#tested keep_maf: 0.01, 0.02, 0.05, 0.07, 0.1,

#too big 0.5, 0.25, 0.15

library(SimRVSequences)
data(hg_exons)
data(hg_apopPath)
data(study_peds)

#create mutation and genome set from slim output
a = Sys.time()
sout = read_slim(file_path = "C:/Data/Slim/SlimFINALout.txt",
                 keep_maf = 0.15,
                 recomb_map = create_slimMap(hg_exons),
                 pathway_df = hg_apopPath,
                 single_SNVtype = FALSE)
b = Sys.time()
difftime(b, a, units = "mins")
nrow(sout[[2]])

#----------------------#
# Single Pedigree Time #
#      sim_RVstudy     #
# affected_only = TRUE #
# remove_wild = TRUE   #
#----------------------#
time_single <- rep(NA, 100)
for(i in 1:100){
  a = Sys.time()
  study_seq <- sim_RVstudy(ped_files = time_peds[time_peds$FamID == i, ],
                           SNV_map = sout$Mutations,
                           haplos = sout$Haplotypes,
                           affected_only = TRUE,
                           remove_wild = TRUE)
  b = Sys.time()

  time_single[i] = difftime(b, a, units = "secs")
}
mean(time_single)
plot(time_single)
