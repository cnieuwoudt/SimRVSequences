library(SimRVSequences)
library(Matrix)
library(vcfR)

#import Sample Data
#These are the individuals we have sampled who are unrelated
SampleData <- read.csv("~/GitHub/1000-Genomes-Exon-Data/Formatted-SNVdata/SampleData.csv", stringsAsFactors=FALSE)

sampLength <- c()
SampID_order <- list()
dimHaplotypes <- list()

#set the chromosome number of the exon data to import
chrom_num = 22


file_path = paste0("C:/Users/cnieuwoudt/Documents/GitHub/1000-Genomes-Exon-Data/Exon Data/exons_chr",
                   chrom_num, ".vcf.gz", sep = "")

#import the vcf file using read.vcfR from the vcfR package
vcf <- read.vcfR(file_path)

#-------------------------------#
# extract and format Haplotypes #
#-------------------------------#
#store the genotype data
#NOTE: removing the first column since it does not hold genotype data
genotypes <- vcf@gt[, -1]

#reduce the genotypes to contain only the individuals we are retaining
#NOTE: people are columns.
genotypes <- genotypes[, which(colnames(genotypes) %in% SampleData$Sample)]

sampLength[[chrom_num]] = dim(genotypes)[2]
SampID_order[[chrom_num]] = unique(colnames(genotypes))

#convert genotypes to sparseMatrix format
Haplotypes <- genos2sparseMatrix(genotypes)

#------------------------------#
# extract and format Mutations #
#------------------------------#

#extract and store the mutation data using vcfR2tidy
#NOTE: can use info_types to choose which INFO variables are included
#NOTE: specifying info_only = TRUE since we do not need to re-process
#the genotype data (this was accomplished when formatting "Haplotypes")
muts <- vcfR2tidy(vcf, info_only = TRUE,
                  info_fields = c("AF", "AC", "NS", "AN",
                                  "EAS_AF", "EUR_AF", "AFR_AF", "AMR_AF", "SAS_AF",
                                  "DP"),
                  info_types = c(AF = "n", AC = "i", NS = "i", AN = "i",
                                 EAS_AF = "n", EUR_AF = "n", AFR_AF = "n", AMR_AF = "n", SAS_AF = "n",
                                 DP = "i"))


#muts$meta contains variable descriptions for the INFO variables
#muts$meta


muts$fix$AF <- muts$fix$AC/muts$fix$AN
#head(muts$fix)

#store as dataframe, remove columns that contain only NA values,
#and create colID variable to identify the column positon of the mutation
#RECALL: mutations are columns and individuals are rows in the sparseMatrix "Haplotypes"
Mutations <- as.data.frame(muts$fix)
Mutations <- cbind(seq(1:dim(Haplotypes)[2]),
                   Mutations[, sapply(1:ncol(Mutations), function(x){any(!is.na(Mutations[, x]))})])
#head(Mutations)

#rename columns for consistency with read_slim output
#According to muts$meta the variable "AF" is the estimated allele freq,
#i.e. the variable we previously named "afreq" in the Mutations output
#returned by read_slim.
colnames(Mutations)[c(1:3, 7)] = c("colID", "chrom", "position", "afreq")
#head(Mutations)


#store chromosome number as integer
Mutations$chrom <- as.integer(Mutations$chrom)

# Note that some of the SNVs have an allele count of zero, i.e. AC = 0.
# Meaning, no individual in the sample carries the SNV.  Also, note that
# since we removed 22 third order relatives, any SNVs that were carried by
# only these individuals will no longer be present in the haplotype data, i.e.
# the column sums for such SNVs will be zero.
#which(Mutations$AC == 0)

library(Matrix)
#which(colSums(Haplotypes) == 0)

# To save space, it makes sense to remove these SNVs from the data.
# remove SNVs that are not carried by any member of the sample
remove_cols <- which(colSums(Haplotypes) == 0)
Haplotypes <- Haplotypes[, -remove_cols]
Mutations <- Mutations[-remove_cols, ]
Mutations$colID = seq(1:nrow(Mutations))
#TODO: get rid of the variable colID. Really, what is the point of this
#variable? SNVs should be in same order in both haplotypes and
#mutations anyway....

dimHaplotypes[[chrom_num]] = dim(Haplotypes)
#Mutations <- Mutations[, c(1:5, 7:16)]

#Assign object name based on chromosoem number
object_name <- paste0("SNVdata_chrom", chrom_num, sep = "")

#assign name based on chromosome number
assign(object_name,
       SNVdata(Haplotypes = Haplotypes, Mutations = Mutations))

#store the formatted SNVdata object in the "Formatted SNVdata" folder
save(SNVdata_chrom22,
     file = paste0("C:/Users/cnieuwoudt/Documents/GitHub/1000-Genomes-Exon-Data/Formatted-SNVdata/SNVdata_chrom", chrom_num, ".rda", sep = ""))

remove(SNVdata_chrom22)
remove(genotypes)
remove(Haplotypes)
remove(Mutations)
remove(muts)
remove(vcf)
remove(remove_cols)

#----------------------------------#
# After all chromomsomes processed #
#----------------------------------#
# test import_SNVdata function
ex_data <- import_SNVdata(1:22, pathway_df = hg_apopPath)

# Quick quality checks
all(sapply(SampID_order, identical, SampID_order[[1]]))
all(sapply(dimHaplotypes, function(x){x[[1]]}) == 5052)
all(sampLength == 2526)

#number of retained SNVs in for each chromosome
sapply(dimHaplotypes, function(x){x[[2]]})
#NOTE: follow up with Wendy about chromosome 19
