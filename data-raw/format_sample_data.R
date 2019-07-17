file_path <- 'https://raw.githubusercontent.com/cnieuwoudt/1000-Genomes-Exon-Data/master/Vignette%20Data/SampleInfo1.csv'
S1 <- read.csv(file_path)


#file_path = paste0("C:/Users/Christina/Documents/GitHub/1000-Genomes-Exon-Data/Exon Data/exons_chr", 21, ".vcf.gz", sep = "")

file_path = "C:/Users/cnieuwoudt/Google Drive/SeqSim/Group Meetings/exons_chr21.vcf.gz"

library(vcfR)
#import the vcf file using read.vcfR from the vcfR package
vcf <- read.vcfR(file_path)


str(vcf)
vcf@gt[1:3, 1:10]

vcf_sampleIDs = colnames(vcf@gt[, -1])

#Reduce the Samples data to include only the individuals who are in
#the vcf file
S_red <- S1[which(S1$Sample %in% vcf_sampleIDs), ]
head(S_red)
S_red <- S_red[order(S_red$Sample, S_red$Family.ID), ]



#--------------------#
# Reduce to Families #
#--------------------#
#List of repeated Family.IDs
#subset S_red further to include only non-duplicated FamIDs
NU_Fids <- unique(S_red$Family.ID[duplicated(S_red$Family.ID)])

#S_Fams holds the individuals who belong to a family
S_Fams <- S_red[S_red$Family.ID %in% NU_Fids, ]
S_Fams <- S_Fams[order(S_Fams$Sample), ]
View(S_Fams)

#Unfortunately, Family.ID is not used as expected and therefore is not
#useful in this context.  Instead, we will manually check each relationship
#variable to see which individuals need to be assessed for removal.


#-------------------------#
# Unexpected Parent Child #
#-------------------------#
# Sample IDs of unexpected parent/child
unique(S_red$Unexpected.Parent.Child)

# determine rows that hold individuals who have been listed as an unexpected parent/child
S_red[which(S_red$Sample %in% unique(S_red$Unexpected.Parent.Child)[-1]), ]
# None of these individuals are in the sample.
# No one to remove, it appears these individuals are not included in the sample.

#---------------#
# Non Paternity #
#---------------#
# Sample IDs of unexpected parent/child
unique(S_red$Non.Paternity)

# determine rows that hold individuals who have been listed as a non-parent
S_red[which(S_red$Sample %in% unique(S_red$Non.Paternity)[-1]), ]
# None of these individuals are in the sample.
# No one to remove, it appears these individuals are not included in the sample.


#----------#
# Siblings #
#----------#
# Sample IDs of siblings
SibIDs <- unlist(strsplit(unique(S_red$Siblings), split = ", "))

# determine rows that hold individuals who have been listed as a sibling
S_red[which(S_red$Sample %in% SibIDs), ]
# None of these individuals are in the sample.
# No one to remove, it appears these individuals are not included in the sample.

#--------------#
# Grandparents #
#--------------#
# Sample IDs of grandparents
unique(S_red$Grandparents)

# determine rows that hold individuals who have been listed as a grandparent
S_red[which(S_red$Sample %in% unique(S_red$Grandparents)[-1]), ]
# None of these individuals are in the sample.
# No one to remove, it appears these individuals are not included in the sample.



#-----------#
# Avuncular #
#-----------#
# Sample IDs of avualcular relations
unique(S_red$Avuncular)

AvuncIDs <- c("HG01442", "NA06997", "NA19192", "NA19469",
              "NA19672", "NA19660", "NA19714", "NA20337",
              "NA20363", "NA20907")

# determine rows that hold individuals who have been listed as avuncular
S_red[which(S_red$Sample %in% AvuncIDs), ]
# None of these individuals are in the sample.
# No one to remove, it appears these individuals are not included in the sample.


#---------------#
# Half Siblings #
#---------------#
# Sample IDs of half siblings
unique(S_red$Half.Siblings)

# determine rows that hold individuals who have been listed as a half sib
S_red[which(S_red$Sample %in% unique(S_red$Half.Siblings)[-1]), ]
# None of these individuals are in the sample.
# No one to remove, it appears these individuals are not included in the sample.


#----------------------#
# Unknown Second Order #
#----------------------#
# Sample IDs of individuals with the specified relationship
unique(S_red$Unknown.Second.Order)

SecOrd_ids <- c(unique(S_red$Unknown.Second.Order)[-c(1, 3, 4, 13, 15, 22)],
                #unlist IDs separated by commas
                unlist(strsplit(unique(S_red$Unknown.Second.Order)[c(3, 4, 22)], split = ", ")),
                #manually insert the last two
                "NA07031", "NA19101")


# determine rows that hold individuals with the specified relationship
S_red[which(S_red$Sample %in% SecOrd_ids), ]
# None of these individuals are in the sample.
# No one to remove, it appears these individuals are not included in the sample.


#-------------#
# Third Order #
#-------------#
# Sample IDs of individuals with the specified relationship
unique(S_red$Third.Order)

ThiOrd_ids <- c(unique(S_red$Third.Order)[-c(1, 3, 4, 18, 23, 24, 34, 28, 29, 30)],
                #unlist IDs separated by commas
                unlist(strsplit(unique(S_red$Third.Order)[c(3, 4, 18, 23, 24, 34)], split = ", ")),
                #manually insert the last two
                "NA12264", "NA19178", "NA06993")


# determine rows that hold individuals with the specified relationship
S_red[which(S_red$Sample %in% ThiOrd_ids), ]

#
#  Next, we will need to sift through these individuals on
#  a case-by-case basis to determine who will be to kept.
#


# To get an idea of what we are working with lets subset the data
# to contain the individuals with a specified third order relationship.
# Also, I'm removing the following variables as they are missing for all individuals:
#   Non.Paternity, Siblings, Grandparents, Avuncular, Half.Siblings, Unknown.Second.Order
SamDat_TO <- S_red[which(S_red$Third.Order != "" | S_red$Sample %in% ThiOrd_ids),
                   -which(colnames(S_red) %in% c("Non.Paternity", "Siblings",
                                                 "Grandparents", "Avuncular",
                                                 "Half.Siblings", "Unknown.Second.Order"))]


SamDat_TO <- SamDat_TO[order(SamDat_TO$Population, SamDat_TO$Sample, SamDat_TO$Family.ID), ]
row.names(SamDat_TO) = NULL

View(SamDat_TO)
SamDat_TO$FamID = NA
unique(SamDat_TO$Population)

#-----#
# CEU #
#-----#
SamDat_TO[SamDat_TO$Population == "CEU", c(1, 8)]
SamDat_TO$FamID[1:2] = 1:2

#-----#
# ESN #
#-----#
SamDat_TO[SamDat_TO$Population == "ESN", c(1, 8)]

SamDat_TO$FamID[3:4] = 3
SamDat_TO$FamID[c(5,9)] = 4
SamDat_TO$FamID[6:8] = 5

#-----#
# GIH #
#-----#
SamDat_TO[SamDat_TO$Population == "GIH", c(1, 8)]
SamDat_TO[SamDat_TO$Sample == "NA21134", ]
SamDat_TO$FamID[10] = 6


#-----#
# GWD #
#-----#
SamDat_TO[SamDat_TO$Population == "GWD", c(1, 8)]
SamDat_TO$FamID[SamDat_TO$Population == "GWD"] = 7

#-----#
# LWK #
#-----#
SamDat_TO[SamDat_TO$Population == "LWK", c(1, 8)]
SamDat_TO[SamDat_TO$Sample == "NA19039", ]
SamDat_TO$FamID[c(15, 19)] = 8
SamDat_TO$FamID[c(16, 20)] = 9
SamDat_TO$FamID[c(17, 18)] = 10

#-----#
# MSL #
#-----#
SamDat_TO[SamDat_TO$Population == "MSL", c(1, 8)]
SamDat_TO[SamDat_TO$Sample == "HG03383", ]
SamDat_TO$FamID[c(21, 22)] = 11
SamDat_TO$FamID[c(23)] = 12
SamDat_TO$FamID[c(24, 25, 26)] = 13
SamDat_TO$FamID[c(27:31)] = 14
SamDat_TO$FamID[c(32:33)] = 15

#-----#
# PJL #
#-----#
SamDat_TO[SamDat_TO$Population == "PJL", c(1, 8)]
SamDat_TO$FamID[c(34:35)] = 16
SamDat_TO$FamID[c(36, 41)] = 17
SamDat_TO$FamID[c(37:38)] = 18
SamDat_TO$FamID[c(39:40)] = 19

#-----#
# YRI #
#-----#
SamDat_TO[SamDat_TO$Population == "YRI", c(1, 8)]
S_red[S_red$Sample == "NA19178",]
SamDat_TO$FamID[c(42)] = 20

SamDat_TO = SamDat_TO[order(SamDat_TO$FamID, SamDat_TO$Sample), ]


Samples_TO <- rep(NA, 20)
set.seed(2934756)
for(i in 1:20){
  if(length(SamDat_TO$Sample[SamDat_TO$FamID == i]) == 1){
    Samples_TO[i] = SamDat_TO$Sample[SamDat_TO$FamID == i]
  } else {
    Samples_TO[i] = sample(x = SamDat_TO$Sample[SamDat_TO$FamID == i], size = 1)
  }
}

#All of the individuals who do not have any relatives in the sample
NotTO_IDs <- S_red$Sample[S_red$Third.Order == "" & !(S_red$Sample %in% ThiOrd_ids)]


#combine the individuals who did not have any relatives with the
#thrid order relatives we just sampled
keep.samples = c(NotTO_IDs, Samples_TO)
SampleData <- S_red[S_red$Sample %in% keep.samples, ]
View(SampleData)


#-------------#
# Third Order #
#-------------#
# check to make sure we don't have any pesky third order cousins included
# Sample IDs of individuals with the specified relationship
unique(SampleData$Third.Order)

ThiOrd_ids <- c(unique(SampleData$Third.Order)[-c(1, 2, 19, 15, 16, 17)],
                #unlist IDs separated by commas
                unlist(strsplit(unique(SampleData$Third.Order)[c(2, 19)], split = ", ")),
                #manually insert the last three
                "NA06993", "NA12264", "NA19178")


# determine rows that hold individuals with the specified relationship
SampleData[which(SampleData$Sample %in% ThiOrd_ids), ]
# YAY!!!!!!

#reduce to columns we are interested in
SampleData <- SampleData[, c(1, 3, 4, 5)]

write.table(SampleData, "C:/Users/cnieuwoudt/Documents/GitHub/1000-Genomes-Exon-Data/Formatted-SNVdata/SampleData.csv",
            col.names = T, row.names = F, sep = ",")

