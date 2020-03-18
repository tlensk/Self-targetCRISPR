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
Organism level statistics and information about the analyzed genome assemblies for the organisms of interest downloaded from CRISPRCasdb 
* __20190618_spacer_34.fasta__ <br>
A complete list of spacers downloaded from CRISPRCasdb
* __summary.txt__ <br> 
Additional information about each replicon (in the first place whether a replicon is linear or circular)<br>
_This information is missing from fasta format files and needs to be obtained from NCBI separately._
* __.fasta files__ for the replicons  list.txt<br>
Each .fasta file contains a genomic sequence of the corresponding replicon 
* __.csv files__ for the replicons <br>
Each .csv file contains information about CRISPR arrays and Cas-gene clusters found in the corresponding replicon

## Output files:

### Organism level statistics:
* __organisms.txt__ <br>
The number of spacers and the number of self-targeting spacers found by our methods in the organisms of interest supplemented by the organism level statistics from export.txt.

### Spacer level statistics:
* __spacer.txt__ <br>
Detailed information about the analysis of spacers found in a given organism by our methods for all the organisms of interest combined. This information also includes spacer length, whether a spacer is located on a plasmid, and whether a CRISPR array bearing a spacer is  functional (whether it is accompanied by at least one cluster of Cas genes).

### Wildcard spacers
* __wildcard_spacers.txt__ <br>
Information about the wildcard spacers that were discarded prior to the analysis.

_Wildcard spacer is a spacer that contains symbols other then {T,C,A,G}, e.g., N, Y, and W._


