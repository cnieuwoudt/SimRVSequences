## Re-submission February March 5, 2019
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
  nucleotide (10:42)

This is the correct spelling for nucleotide, as in "single-nucleotide variant".


## Downstream Dependencies
Currently, there are no downstream dependencies for this package.
