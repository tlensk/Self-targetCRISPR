# Leaders results

This folder contains information about the leaders among Archaea and Bacteria with respect to the number of self-targeting spacers. 
One archaeal organism _Methanobrevibacter olleyae (euryarchaeotes) YLM1_ (CP014265.1) has 18 self-targeting spacers,
and two strains of _Megasphaera elsdenii DSM 20460_ (CP027570.1 and HE576794.1) have 35 self-targeting spacers per a genome.
Each file contains information about the self-targeting spacers found in the corresponding genome. 

We analyzed gbk files for these organisms to identify genes that bear targets for these self-targeting spacers. Each line in a file corresponds to the found relation between a self-targeting spacer and an annotated entry where the target for this spacer is located. The number of lines in a file can exceed the number of the self-targeting spacers found in the corresponding genome. It happens if at least one of the following conditions is met:

If __a spacer have more than one targets__ in a genome then the corresponding file will contain several lines with the same values for Spacer_start_pos and Spacer_end_pos but values for Target_start_pos and Target_end_pos will be different. 

If __a spacer has duplicates__ in a genome then the corresponding file will contain lines with the same values for Target_start_pos and Target_end_pos, and the values for Spacer_start_pos and Spacer_end_pos will be different for all the targets found in the genome.

If __a target of a self-targeting spacer overlaps with several annotated entries__ in a genome then then the corresponding file will contain a separate line for each of those entries. The values for Target_start_pos, Target_end_pos, Spacer_start_pos, and Spacer_end_pos will be the same for those lines, but the Entry_info, Entry_ID, Entry_start_pos, and Entry_end_pos will be different. 


### For example:
* __Info_CP014265.1.txt__ <br>
One of the self-targeting spacers has a target that overlaps with two genes.

* __Info_CP027570.1.txt__ and __Info_HE576794.1.txt__ <br>
For each of these genomes, one of the self-targeting spacers has duplicates in the same CRISPR array and only one target that they aim to. Also, another self-targeting spacer have two targets in each of the two genomes.




