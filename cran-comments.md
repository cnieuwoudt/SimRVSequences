## Re-submission March 8, 2019
In this re-submission we have 
  * removed examples from unexported functions, and
  * added a working example to the read_slim function.
  
Note: 
  The read_slim function imports genomic data (simulated by SLiM 2.0+) to R.
  Previously, the examples for the read_slim function were "dontrun," and only illustrated how the function could be called under certain assumptions (e.g. given file name and path). 

  In this re-submission I have, instead, included a working example.  To accomplish this, I created a github repository which contains publicly-available data for this purpose.   Unfortunately, the length of the URL prompts in a new note when running checks, which I've displayed below as "Second Note".

## Test environments
* local Windows OS install, R 3.5.0
* ubuntu 14.04.5 LTS (on travis-ci), R 3.5.2
* win-builder (release)

## R CMD check results
0 errors | 0 warnings | 2 notes

#First Note:
Possibly mis-spelled words in DESCRIPTION:
  Jinko (13:9)
  Nieuwoudt (12:15)


#Second Note:
* checking Rd line widths ... NOTE
Rd file 'read_slim.Rd':
  \examples lines wider than 100 characters:
     file_url <-'https://raw.githubusercontent.com/cnieuwoudt/Example--SLiMSim/master/example_SLIMout.txt'

These lines will be truncated in the PDF manual.



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
