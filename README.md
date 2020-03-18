# Self-targetCRISPR 

Tatiana Lenskaia (lensk010@umn.edu) and Daniel Boley (boley@umn.edu)

Computational analysis of CRISPR-Cas systems and their autoimmunity potential in prokaryotes

## Part I. Data acquisition and preprocessing

### Data files from CRISPRCasdb

The following three files were downloaded from CRISPRCasdb<sup>1</sup> accessed in February 2020:

1. __export.xlsx__ <br> 
a xlsx file from Strain list section (filters: CRISPR level 4, Bacteria and Archaea) <br>
We saved this xlsx file as tab-delimited txt file __export.txt__ for the further analysis <br>
https://crisprcas.i2bc.paris-saclay.fr/MainDb/StrainList


2. __20190618_spacer_34.fasta__ <br>
CRISPRCasdb spacers file from Download section <br>
https://crisprcas.i2bc.paris-saclay.fr/Home/Download


3. __20190617_ccpp.sql__ <br>
CRISPRCasdb (SQL Dump) file from Download section <br>
https://crisprcas.i2bc.paris-saclay.fr/Home/Download
_We rebuilt the Postgres database from the dump sql file. For each sequence, we created a separate (comma-delimited) file with the sequence accession number as its name and ".csv" as the extension to store information about the found CRISPR arrays and Cas gene clusters extracted from entity, region, crisprlocus, and clustercas tables (Postgres). We uploaded the .csv files for the test example to this repository (test_example folder)._


From __export.txt__, we extracted information about the organisms that possess CRISPR arrays level 4 (6865 prokaryotes) and constructed the corresponding list of accession numbers for replicons in their genomes (__list.txt__, 13816 replicons in total). 

From __20190618_spacer_34.fasta__, we obtained information about 221397 spacer entries.  We discarded spacers that contained at least one symbol other than {A,C,G,T}, i.e., "wildcard spacers" (__wildcard_spacers.txt__, 290 spacer entries) and considered only "exact" spacers for further analysis (221107 spacer entries). We found 326187 instances of exact spacers in the analyzed genomes (the same spacer string might appear in several genomes and consequently it yields several spacer instances, __Fig.1__). An instance of a spacer indicates that a spacer string is present in at least one CRISPR array in a given genome (if the spacer string is present in several CRISPR arrays and/or several times in the same CRISPR array in a given genome then it yields spacer copies, i.e., spacer duplicates).

![Fig.1](/images/wildcard_spacerinstance.png)

__Fig.1.__ Example of wildcard spacer, spacer entry, and spacer instance.

### Data files from NCBI

We downloaded information about genomic sequences for replicon accession numbers stored in __list.txt__ from the NCBI database using the Batch Entrez interface https://www.ncbi.nlm.nih.gov/sites/batchentrez (both fasta and summary files). We also extracted the information about replicon name, length, type (linear or circular), and base (DNA or RNA) and recorded this information for every replicon accession number in __summary.txt__. The sequences (~ 30 Gb in total) and the corresponding summary information can be obtained from the NCBI database.

_We uploaded fasta files and the corresponding short version of export.txt, list.txt, and summary.txt for the test example organisms to test_example directory in this repository_. 



## Part II. Analysis of self-targeting events in Bacteria and Archaea

### Overview
Instead of using an alignment-based approach, we use an “exact matching” approach inspired by the CRISPR mechanism itself. To search for exact matches and to accurately determine the location of the found matches, we utilized our dictionary-based methods. Our approach is made efficient by using a dictionary (hash table) data structure. To search for self-targeting spacers in a genome, we took information about found CRISPR arrays in the genome and grouped all the found spacers by length. For each of spacer lengths, we created a dictionary of the strings found in a genome with the unique strings of a given length as the keys and the lists of positions of these strings in the genome as the values. Then, we searched each dictionary for all spacers of that length. To find copies on the forward and reverse strands, we searched the dictionary for the spacer itself and its reverse complement. As a result, for each spacer, we recorded all position(s) of the copies found for this spacer on all replicons in the genome and compared these positions to the location of CRISPR arrays in the genome. This helped us accurately distinguish between duplicate spacers and matched strings located outside CRISPR arrays and precisely identify self-targeting spacers. After we benchmarked our methods on the dataset employed by Stern et al.<sup>2</sup>(see Benchmark.xlsx in the Supplementary Materials), we conducted a similar analysis on the latest dataset available at CRISPRCasdb (June, 2019), a successor of CRISPRdb<sup>3</sup>.

### Implementation
The core methods for the analysis of self-targeting events are implemented in Python 3 and stored in __Find_ST.py__ file availabe in Code directory. 

### Results for the CRISPRCasdb dataset

The results of this analysis conducted on the CRISPRCasdb data are available in the Supplementary Materials: 

#### Organism level statistics:
* __CRISPRCasdb_organisms.txt__  <br>
The number of spacers and the number of self-targeting spacers found in the organisms of interest during the analysis supplemented by the organism level statistics from export.txt.

#### Spacer level statistics:
* __CRISPRCasdb_spacers.zip__ <br>
This zip archive contains __CRISPRCasdb_spacer.txt__ (~39 Mb) that provides the detailed information about the analysis of spacers found in a given organism by our methods for all the organisms of interest combined. This information also includes spacer length, whether a spacer is located on a plasmid, and whether CRISPR array bearing a spacer is functional (whether it is accompanied by at least one cluster of Cas genes).
<br>

The calculation of the number of spacers and the number of self-targeting spacers for a toy set of two organisms is illustrated on __Fig.2__.

![Fig.2](/images/Stat_desc.png)
__Fig.2.__ The calculation of the number of spacers and the number of self-targeting spacers for a toy set of organisms.

<br><br><br>
References

1. Pourcel, C., Touchon, M., Villeriot, N., Vernadet, J. P., Couvin, D., Toffano-Nioche, C., & Vergnaud, G. (2020). CRISPRCasdb a successor of CRISPRdb containing CRISPR arrays and Cas genes from complete genome sequences, and tools to download and query lists of repeats and spacers. _Nucleic acids research_, 48(D1), D535-D544.

2. Stern, A., Keren, L., Wurtzel, O., Amitai, G., & Sorek, R. (2010). Self-targeting by CRISPR: gene regulation or autoimmunity? _Trends in genetics_, 26(8), 335-340.

3. Grissa, I., Vergnaud, G., & Pourcel, C. The CRISPRdb database and tools to display CRISPRs and to generate dictionaries of spacers and repeats. _BMC bioinformatics_, 8(1), 172, 2007.
