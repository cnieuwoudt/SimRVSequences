library(SimRVSequences)
library(Matrix)
library(vcfR)

#import Sample Data
#These are the individuals we have sampled who are unrelated
SampleData <- read.csv("~/GitHub/1000-Genomes-Exon-Data/Formatted-SNVdata/SampleData.csv", stringsAsFactors=FALSE)

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
dim(genotypes)

#convert genotypes to sparseMatrix format
Haplotypes <- genos2sparseMatrix(genotypes)
dim(Haplotypes)
Haplotypes[1:10, 1:10]

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
                                  "VT", "EX_TARGET", "DP"),
                  info_types = c(AF = "n", AC = "i", NS = "i", AN = "i",
                                 EAS_AF = "n", EUR_AF = "n", AFR_AF = "n", AMR_AF = "n", SAS_AF = "n",
                                 DP = "i"))


#muts$meta contains variable descriptions for the INFO variables
muts$meta


#store as dataframe, remove columns that contain only NA values,
#and create colID variable to identify the column positon of the mutation
#RECALL: mutations are columns and individuals are rows in the sparseMatrix "Haplotypes"
Mutations <- as.data.frame(muts$fix)
Mutations <- cbind(seq(1:dim(Haplotypes)[2]),
                   Mutations[, sapply(1:ncol(Mutations), function(x){any(!is.na(Mutations[, x]))})])
head(Mutations)

#remove unnecessary columns ???, i.e. FILTER, EX_TARGET, and VT
#From metadata:
#EX_TARGET indicates whether a variant is within the exon pull down target boundaries
#VT indicates what type of variant the line represents
unique(Mutations$VT)
unique(Mutations$FILTER)


#rename columns for consistency with read_slim output
#According to muts$meta the variable "AF" is the estimated allele freq,
#i.e. the variable we previously named "afreq" in the Mutations output
#returned by read_slim.
colnames(Mutations)[c(1:3, 7)] = c("colID", "chrom", "position", "afreq")
head(Mutations)

#Chrom must be a numeric variable.
#Convert chrom to numeric
Mutations$chrom <- as.numeric(Mutations$chrom)

#How rare are the SNVs in the exported exon data?
nrow(Mutations)
length(which(Mutations$afreq <= 0.01))

#from the following output rouhgly 93% of the SNVs are rare (AF <= 0.01) for chr22
#Did we filter for rarity? Is it normal for so many SNVs to be rare?
length(which(Mutations$afreq <= 0.01))/nrow(Mutations)

#Do we have data for every person, i.e. no missing values for any mutation
unique(Mutations$NS)  #GOOD! no missing values
unique(Mutations$AN)  #GOOD! no missing values

head(Mutations)
Mutations <- Mutations[, c(1:5, 7:15, 17)]

#Assign object name based on chromosoem number
object_name <- paste0("SNVdata_chrom", chrom_num, sep = "")

#assign name based on chromosome number
assign(object_name,
       SNVdata(Haplotypes = Haplotypes, Mutations = Mutations))

#store the formatted SNVdata object in the "Formatted SNVdata" folder
save(SNVdata_chrom22,
     file = paste0("C:/Users/cnieuwoudt/Documents/GitHub/1000-Genomes-Exon-Data/Formatted-SNVdata/SNVdata_chrom", chrom_num, ".rda", sep = ""))

# #store the formatted SNVdata object in the "Formatted SNVdata" folder
# save(vcf_chrom21,
#      file = paste0("C:/Users/cnieuwoudt/Documents/GitHub/1000-Genomes-Exon-Data/Formatted-SNVdata/vcf_chrom", chrom_num, ".rda", sep = ""))

# remove(SNVdata_chrom9)
#
# #load the data from my computer, with a new name, to see if format was preserved
# load(file = paste0("C:/Users/cnieuwoudt/Documents/GitHub/1000-Genomes-Exon-Data/Formatted SNVdata/SNVdata_chrom", chrom_num, ".rda", sep = ""))


# #load SNVdata from github repository
# load(url(paste0("https://github.com/cnieuwoudt/1000-Genomes-Exon-Data/raw/master/Formatted-SNVdata/SNVdata_chrom", chrom_num, ".rda", sep = "")))
