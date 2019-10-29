#load pedigrees and SNV data
source("data-raw/VignetteTimingData.R")


#------------------------#
# Single Pedigree Timing #
#------------------------#
#time sv_RVstudy for the various timing settings
#data frame of times to populate for each set of settings
time_dat_1P <- data.frame(time_set1 = rep(NA, 100),
                       time_set2 = rep(NA, 100),
                       time_set3 = rep(NA, 100),
                       time_set4 = rep(NA, 100))

simStudy_settings <- data.frame(affOnly = c(TRUE, TRUE, FALSE, FALSE),
                                remWild = c(TRUE, FALSE, TRUE, FALSE))


for(i in 1:100){ #loop over the pedigrees
  #time sim_RVstudy for each set of settings
  for(j in 1:nrow(simStudy_settings)){
    a = Sys.time()
    study_seq <- sim_RVstudy(ped_files = time_peds[time_peds$FamID == i, ],
                             SNV_data = sout,
                             affected_only = simStudy_settings$affOnly[[j]],
                             remove_wild = simStudy_settings$remWild[[j]])
    b = Sys.time()

    time_dat_1P[i, j] = difftime(b, a, units = "mins")
  }
}


#---------------------#
# 100 Pedigree Timing #
#---------------------#
#time sv_RVstudy for the various timing settings
#data frame of times to populate for each set of settings
time_dat_100P <- data.frame(time_set1 = rep(NA, 1),
                          time_set2 = rep(NA, 1),
                          time_set3 = rep(NA, 1),
                          time_set4 = rep(NA, 1))

simStudy_settings <- data.frame(affOnly = c(TRUE, TRUE, FALSE, FALSE),
                                remWild = c(TRUE, FALSE, TRUE, FALSE))

for(j in 1:nrow(simStudy_settings)){
  a = Sys.time()
  study_seq <- sim_RVstudy(ped_files = time_peds,
                           SNV_data = sout,
                           affected_only = simStudy_settings$affOnly[[j]],
                           remove_wild = simStudy_settings$remWild[[j]])
  b = Sys.time()

  time_dat_100P[1, j] = difftime(b, a, units = "mins")
}

time_dat_100P

#-------------------------------------#
# pedigree reduction (100 ped sample) #
#-------------------------------------#
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
