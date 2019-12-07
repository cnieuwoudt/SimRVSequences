## NEXT UPDATE
 - add catch for same famID in multiple study peds
 - add catch for fully 0 haplotypes

## Re-submission December 5, 2019
In this resubmission we have fixed all improper uses of class() == *, and instead switched to using inherits() or is().
Many thanks to Kurt Hornik for pointing out this deficiency and for sharing the informative blog by Martin Maechler. 



## Re-submission November 3, 2019

In this resubmission we have added several new features to SimRVSequences:
  * We now allow users to model both dominant and recessive disease models
  * We include a new function to import pre-formatted exon-only SNV data
  * We include a new function to help users process SNV data from vcf files
  * We have expanded the vignette to describe how users can format vcf file data for use with SimRVSequences

- Updated the DESCRIPTION file:
    - Updated the version number to 0.2.0 to reflect:
         - the new features listed above, and 
         - the deprecated arguments SNV_map and haplos from the sim_RVstudy function,
    - Fixed vignette yaml header to remove NOTEs on package check page for flavor:
 r-devel-linux-x86_64-debian-clang, r-devel-linux-x86_64-debian-gcc, nd r-devel-windows-ix86+x86_64.

Test Environments prior to resubmission
* https://win-builder.r-project.org/
* local Windows OS install, R 3.5.3
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.0


## Re-submission April 18, 2019
In this resubmission we have corrected the sampling probability for causal rare variants in sim_RVstudy.  Previously, was sampling with equal probability, now sampling according to allele frequency in the population.
  


## Re-submission March 8, 2019
In this re-submission we have 
  * removed examples from unexported functions, and
  * created a working example for the read_slim function.
  
  Previously, the examples for the read_slim function were "dontrun," and only illustrated how the function could be called under certain assumptions (e.g. given file name and path). In this re-submission I have, instead, included a working example.  To accomplish this, I created a github repository which contains publicly-available data for this purpose.

## Test environments
* local Windows OS install, R 3.5.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.2
* win-builder (release)

## R CMD check results
0 errors | 0 warnings | 1 notes

Possibly mis-spelled words in DESCRIPTION:
  Jinko (13:9)
  Nieuwoudt (12:15)
  
These are not mis-spelled words, they are the names of authors.



## Re-submission March 5, 2019
* In response to the reviewer's question: Is there some reference about the 
method you can add in the Description field in the form Authors (year) 
<doi:.....>?

We posted a preprint of our upcoming manuscript, which details the methodology for this R package, on bioRxiv and added the reference (Authors (Year) <doi: >) to the DESCRIPTION.  We plan to submit the manuscript to BioInformatics in the coming weeks.

* In response to the reviewer's question:  We see many Rd files do not have any examples, can you add more, please?

It is true that most of the Rd files do not have examples. The reason for this is because we export 5 datasets, and 1 function to import and format data simulated by SLiM (Haller and Messer, 2017).  For the method that imports data, we have added a "dontrun" example to illustrate how this function could be called, assuming file path, etc.  

Additionally, we have included an extensive discussion of all methods and datasets in the vigentte.

## Test environments
* local Windows OS install, R 3.5.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.2
* win-builder (release)

## R CMD check results
0 errors | 0 warnings | 1 note

Possibly mis-spelled words in DESCRIPTION:
  Jinko (13:9)
  Nieuwoudt (12:15)

These are not mis-spelled words, these are the names of authors.  

## February 9, 2019 Submission
## Test environments
* local Windows OS install, R 3.5.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.2
* win-builder (release)

## R CMD check results
0 errors | 0 warnings | 1 note

Possibly mis-spelled words in DESCRIPTION:
  nucleotide (10:42)

This is the correct spelling for nucleotide, as in "single-nucleotide variant".


## Downstream Dependencies
Currently, there are no downstream dependencies for this package.
