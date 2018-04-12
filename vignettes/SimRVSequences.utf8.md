---
title: "SimRVSequences"
author: "Christina Nieuwoudt and Jinko Graham"
date: "2018-04-11"
output: rmarkdown::pdf_document
setspace: doublespacing
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<style type="text/css">

body{ /* Normal  */
   font-size: 16px;
}
</style>

<style>
    p {line-height: 2em;}
</style>

# Table of contents
1. [Introduction](#introduction)
2. [Create Recombination Map for SLiM 2.0](#SlimMat)
3. [Import SLiM 2.0 Data](#SlimDat)
4. [Select Pool of Causal Variants](#CausalVars)
5. [Import Pedigrees Simulated for Multiple Disease-Affected Relatives](#PedSample)
6. [Simulate Sequence Data](#SimSeq)
7. [References](#Ref)


# 1. Introduction <a name="introduction"></a>

Family-based studies are attractive because they have more power to detect rare variants, require smaller sample sizes, and can more accurately detect sequencing errors than case-control studies.  However, data collection for these studies is both time consuming and expensive.  

`SimRVSequences` provides methods to simualte sequence data for family-based studies.  To simulate sequence data we require: (1) a sample of ascertained pedigrees and (2) single-nucleotide variant (SNV) data from a sample of unrelated individuals.  At present, we streamline use of pedigrees simualted by `SimRVPedigree` and exon-only SNV data produced by SLiM 2.0 [1].  In this vignette we illustrate how to use the methods provided by `SimRVSequences` to accomplish this task.

# 2. Create Recombination Map for Exon-Only Data with SLiM 2.0 <a name="SlimMap"></a>

With SLiM 2.0 users may specify a recombination map to simulate recombination hotspots.  Additionally, the recombination map can be used to simulate mutations over unlinked regions (i.e. in different chromosomes) or in linked but non-contiguous regions (i.e in exon-only data).  The `create_slimMap` function may be used to create a recombination map to simulate exon-only data with SLiM.  We now illustrate this process using the `hg_exons` dataset.


```r
# load the SimRVSequences library
library(SimRVSequences)

# load the hg_exons dataset
data("hg_exons")

# print the first three rows of hg_exons
head(hg_exons, n = 4)
```

```
##   chrom exonStart exonEnd                 geneName
## 2     1     11873   12227              NR_046018.2
## 3     1     12612   12721              NR_046018.2
## 4     1     13220   14829 NR_046018.2, NR_024540.1
## 5     1     14969   15038              NR_024540.1
```

As seen in the output above, the `hg_exons` dataset catalogues the position of each exon in the 22 human autosomes.  The variable `chrom` is the chomosome on which the exon resides, `exonStart` is the position of the first base pair in the exon, and `exonStop` is the position of the last base pair in the exon.  The variable `geneName` is the NCBI reference sequence identifier for the gene that the exon resides on.  In `hg_exons` overlapping exons are combined into a single exon; when this occurs variable `geneName` will contain a list of NCBI reference sequence identifiers.


The `create_slimMap` function has three arguments:

1. `exon_df`: A dataframe that contians the positions for each exon of interest.  This dataframe must contain the variables `chrom`, `exonStart`, and `exonEnd`.  *We expect that `exon_df` does not contain any overlapping segments.  Prior to supplying exon data to `create_slimMap` users must combine overlapping exons into a single observation.*  Furthermore, `exon_df` **must** contain the variables  `chrom`, `exonStart`, and `exonEnd`.  The variable `geneName`, as seen in `hg_exons`, is not required.
2. `mutation_rate`: the per site per generation mutation rate. By default, `mutation_rate= 1E-8`, as in [2].
3. `recomb_rate`: the per site per generation recombination rate.  By default, `recomb_rate= 1E-8`, as in [2]   
 

```r
# create recombination map for exon-only data using the hg_exons dataset 
s_map <- create_slimMap(exon_df = hg_exons)

# print first four rows of s_map 
head(s_map, n = 4)
```

```
##   chrom  dist no  recRate mutRate   type simDist endPos
## 1     1 11872  1 0.00e+00   0e+00 intron       1      1
## 2     1   355  1 1.00e-08   1e-08   exon     355    356
## 3     1   384  2 3.84e-06   0e+00 intron       1    357
## 4     1   110  2 1.00e-08   1e-08   exon     110    467
```


The `create_slimMap` function returns a dataframe with several variables, but only three are required by SLiM 2.0 to simulate exon-only data:

1. `recRate`: the per site per generation recombination rate.  For computational efficiency, introns between exons on the same chromosome are simulated as a single base pair with `rec_rate` = recombination rate * number of base pairs in the intron.  For each chromosome, a single intron is created between the last exon on the previous chromosome and the first exon of the current chromosome.  This exon will have recombination rate 0.5, so that exons on separate chromosomes are unlinked. (cite Harris?)
2. `mutRate`: the per site per generation mutation rate.  Since we are interested in exon-only data, the mutation rate in introns is set to zero.
3. `endPos`: The position of the last base pair of the segment. 

The other variables seen in the output above are used to remap mutations to their correct postions after simulation. 

SLiM 2.0 is written in a scripting language called Eidos. Unlike an `R` array, the first position in an Eidos array is 0.  Therefore, we must shift the variable `endPos` forward 1 unit before suppying this data to SLiM 2.0.


```r
# restrict output to the variables required by SLiM
slimMap <- s_map[, c(4, 5, 8)]

# shift endPos up by one unit
slimMap$endPos <- slimMap$endPos - 1

# print first four rows of slimMap 
head(slimMap)
```

```
##    recRate mutRate endPos
## 1 0.00e+00   0e+00      0
## 2 1.00e-08   1e-08    355
## 3 3.84e-06   0e+00    356
## 4 1.00e-08   1e-08    466
## 5 4.98e-06   0e+00    467
## 6 1.00e-08   1e-08   2077
```

The creators of SLiM provide excellent resources and documentation to slimulate forwards-in-time evoulutionary data, which can be found at the SLiM website: https://messerlab.org/slim/ .

# 3. Import SLiM 2.0 Data <a name="SlimDat"></a>

To import data simulated by SLiM, we provide the `read_slim` function.  The `read_slim` function is only appropriate for data produced by SLiM's outputFull() method. Presently, we do not support output in MS or VCF format.

The `read_slim` function has three arguements:

1. `file_path`: The file path of the .txt output file created by SLiM 2.0.
2. `keep_maf`: The largest allele frequency for retained SNVs. All variants with allele frequency greater than `keep_maf` will be removed. Please note, removing common variants is recommended for large datasets due to the limitations of data allocation in R.
3. `recomb_map`: A recombination map of the same format as the data frame returned by `create_slimMap`.  Users who followed the instructions in Section 2, may simply supply the default output from `create_slimMap` as `recomb_map`.  

To clarify, the argument `recomb_map` is used to remap mutations to their actual locations and chromosomes.  This is necessary when data has been simulated over non-contiguous regions such as exon-only data.

In addition to reducing the size of the data, the argument `keep_maf` has a practible applicability as well.  In family-based studies, common SNVs are generally filtered out prior to analysis.  Users who intend to study common variants in addition to rare variants may need to run analyses separately for different chromosomes to allow for allocation of large data sets.


```r
# Let's suppose the output is named slimOut.txt and is saved in the 
# current working directory.  We import the data using the read_slim function.
s_out <- read_slim(file_path  = slimOut.txt, 
                   recomb_map = create_slimMap(hg_exons))

summary(s_out)
```

As shown above, `s_out` is a list containing two items: a data frame named `Mutations` and a dgCMatrix named `Genomes`.


```r
# view the first 4 observations of the Mutations dataset
head(s_out$Mutations, n = 4)
```

The `Mutations` dataframe is used to catalogue the SNVs in `Genomes`.  The variable `colID` associates the rows of `Mutations` to the columns of `Genomes`, `chrom` identifies the chromosome number, `position` is the position of the SNV in base pairs, `afreq` is the derived allele frequency of each SNV, and `marker` is a unique identified for each SNV.


```r
# view the first 5 rows and columns of Genomes
s_out$Genomes[1:5, 1:5]
```

As expected, `Genomes` is a sparse matrix of class dgCMatrix.  Entries of `1` represent a mutated allele, while the wild type is represented as `.`.  Recall that we retain SNVs with derived allele frequency less than or equal to `keep_maf`; hence, the majority of entries represent the wild type allele.  

# 4. Select Pool of Causal Variants <a name="CausalVars"></a>

**VERY SIMILAR TO PAPER.. EDIT THIS**

Now that we have simulated a wide array of mutations over the exons, we must decide which mutations will be modeled as causal rare variants. Family-based studies are attractive because of their ability to detect rare causal variation.  However, even with family-based studies researchers may not be able to identify a strong association between the disease of interest and a single gene or variant.  When this occurs, the next step is to consider that variation in a pathway, or a set of related genes, may be sufficient to predispose individuals to the disease of interest.  We will focus on the latter approach, and note that implementation of a single causal variant is equivalent to a pathway containing a single gene with a single mutation.

When cells become damaged it is biologically advantageous for these cells to be eliminated.  Apoptosis, also known as programmed cell death, is a biological process that destroys damaged or unnecessary cells (cite ThomThom).  Since cancer involves the uncontrolled replication of a damaged cell, it is possible that deleterious mutations in the apoptosis pathway may increase the risk of developing cancer.  We will demonstrate an approach to select variants from an apoptosis-like pathway.  The TNFSF10 gene is a member of the human apoptosis pathway (cite KEGG).  We defined a pseudo-apoptosis sub-pathway centered about TNFSF10 using the 25 genes that had the highest interaction with this gene (citeUCSCGeneInteraction).  The data for this pseudo-pathway is contained in the `hg_apopPath` data set.


```r
# load the hg_apopPath data
data("hg_apopPath")

#View the first 4 observations of hg_apopPath
head(hg_apopPath, n = 4)
```

```
##   chrom exonStart   exonEnd
## 2     1 155689090 155689243
## 3     1 155709032 155709125
## 4     1 155709772 155709824
## 5     1 155717005 155717128
##                                                                   geneName
## 2 NM_001199851.1, NM_001199850.1, NM_001199849.1, NM_004632.3, NM_033657.2
## 3                                                           NM_001199849.1
## 4 NM_001199851.1, NM_001199850.1, NM_001199849.1, NM_004632.3, NM_033657.2
## 5                 NM_001199850.1, NM_001199849.1, NM_004632.3, NM_033657.2
##   gene
## 2 DAP3
## 3 DAP3
## 4 DAP3
## 5 DAP3
```

The `hg_apopPath` dataset is similar the `hg_exons` dataset dicussed in section 2.  However, `hg_apopPath` only catalogues the postions of exons contained in our pseudo-pathway.  

The `identify_pathwayRVs` function is used to further refine the `Mutations` data set returned by `read_slim`.  The `identify_pathwayRVs` function has three arguments:

1. `markerDF` a dataframe containing SNV data, which must be of the same format as the `Mutations` dataframe returned by `read_slim`. 2. `pathwayDF` A dataframe that contains the positions for each exon in the pathway of interest.  This dataframe must contain the variables `chrom`, `exonStart`, and `exonEnd`.  *We expect that `pathwayDF` does not contain any overlapping segments.*
3. `carrier_prob`  The carrier probability for all causal variants with relative-risk of disease GRR. By default, `carrier_prob = 0.002`.  

for the mutations in the 

Before we simulate sequence data we create an object of class


```r
# identify variants located in exons contained in our pathway
mutDF <- identify_pathwaySNVs(markerDF  = s_out$Mutations,
                              pathwayDF = hg_apopPath)

#View the first 4 observations of mutDF
head(mutDF, n = 4)
```

Upon supplying the `Mutations` dataset and a pathway dataset, such as `hg_apopPath`, to `identify_pathwayRVs`, the variable `possibleRV` is marked FALSE for all SNVs not located in the pathway of interest.  Additionally, `possibleRV` is marked FALSE for any SNVs within our pathway with derived allele frequency greater than `carrier_prob`.

For `SimRVpedigree`'s ascertainment process to be valid we require the cumulative allele frequency of all causal rare variants to be small, i.e $\le$ 0.002 (cite Nieuwoudt 2017).  Hence, for a pool of approximately 20 causal rare variants, each could have derived allele frequency 0.0001. 
 
 
We randomly sample 20 variants with minor allele frequency 0.0001 from the genes contained in our pathway and define these variants to be our pool of candidate disease variants.  We then sample familial causal rare variants from this pool so that different families segregate different rare variants. Upon identifying the familial rare variant we then sample haplotypes for each founder from the distribution of haplotypes conditioned on the founder's rare variant status at the familial disease locus.  This ensures that the causal rare variant is introduced by the correct founder.  With this simulation tool in hand, researchers can investigate a wide variety of methods to identify causal rare variants in ascertained pedigrees. 


# 5. Import Pedigrees Simulated for Multiple Disease-Affected Relatives <a name="PedSample"></a>

The R package `SimRVPedigree` [3] is used to simulate pedigrees ascertained for multiple disease-affected relatives.  For the purpose of illustration we will use the dataset `EgPeds`, which is included with `SimRVPedigree`.  To learn more about simulating pedigrees with `SimRVPedigree` please refer to the vignette included with the `SimRVPedigree` package. 



```r
# load the SimRVPedigree library
library(SimRVPedigree)

# import the EgPeds dataset
data(EgPeds)

# view first 4 obsetvations of EgPeds
head(EgPeds, n = 4)
```

```
##   FamID ID sex dadID momID affected DA1 DA2 birthYr onsetYr deathYr RR
## 1     1  1   0    NA    NA    FALSE   1   0    1910      NA    1948 15
## 2     1  2   1    NA    NA    FALSE   0   0      NA      NA      NA  1
## 3     1  3   0     1     2    FALSE   1   0    1930      NA    1942 15
## 4     1  4   0     1     2    FALSE   1   0    1938      NA    1987 15
##   available Gen proband
## 1      TRUE   1   FALSE
## 2     FALSE   1   FALSE
## 3      TRUE   2   FALSE
## 4      TRUE   2   FALSE
```

*Maybe Perhaps we should note which variables are required to simulate sequence data*

# 6. Simulate Sequence Data <a name="SimSeq"></a>

# 7. References <a name="Ref"></a>

[1] Benjamin C. Haller and Philipp W. Messer (2017). **Slim 2: Flexible, interactive forward genetic simulations**. Molecular
Biology and Evolution; 34(1), pp. 230-240.

[2] Kelley Harris and Rasmus Nielsen (2016). **The genetic cost of neanderthal introgression**. Genetics, 203(2): pp. 881-891.

[3] Christina Nieuwoudt and Jinko Graham (2018). 
  **SimRVPedigree: Simulate Pedigrees Ascertained for a Rare Disease.** 
  *R package version 0.1.0*
   https://CRAN.R-project.org/package=SimRVPedigree.
  
