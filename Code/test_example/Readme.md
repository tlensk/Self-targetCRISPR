# Test example

This folder contains files necessary to run the analysis on the test example of two strains of _Acetobacteraceae bacterium_:

1. ***Acetobacteraceae bacterium (a-proteobacteria) 868*** <br>
2 replicons: chromosome (CP039460.1) and plasmid (CP039461.1)
2. ***Acetobacteraceae bacterium (a-proteobacteria) 880*** <br>
1 replicon: chromosome (CP039459.1)


File __list.txt__ contains accession numbers for three replicons in these two organisms.<br>
This file is used to check if all required files (__.fasta__ and __.csv__) are provided for the replicons to conduct the analysis.

_For each replicon, the corresponding .fasta and .csv files should contain the replicon accession number as the file name followed by the file extention._


## Input files:
* __export.txt__ <br>
Organism level statistics and information about the analyzed genome assembly downloaded from CRISPRCasdb 
* __20190618_spacer_34.fasta__ <br>
A complete list of spacers downloaded from CRISPRCasdb
* __summary.txt__ <br> 
Additional information about each replicon (whether a replicon is linear or circular)<br>
_This information is missing from fasta format files and needs to be obtained from NCBI separately._
* __.fasta files__ for the replicons <br>
Each .fasta file contains a genomic sequence of the corresponding replicon 
* __.csv files__ for the replicons <br>
Each .csv file contains information about CRISPR arrays and Cas-gene clusters found in the corresponding replicon

## Output files:

### Organism level statistics:
* __organisms.txt__ <br>
The number of spacers and the number of self-targeting spacrs per organism supplemented by organism level statistics from export.txt
* __organisms_timing.txt__ <br>
The processing time per organism

## Spacer level statistics (the copies found for each spacer in a given organism):
* __1_Acetobacte_spinfo.txt__
* __2_Acetobacte_spinfo.txt__ 

_Naming convention for spinfo.txt files is the following (three parts separated by underscrore symbol):_ <br>
(a) the leading number represent the order of a given organism as it appers in export.txt;<br>
(b) ten first letters from the strain name in export.txt;<br>
(c) "spinfo.txt"._


