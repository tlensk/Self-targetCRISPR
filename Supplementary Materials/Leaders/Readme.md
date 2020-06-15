# Leaders results

This folder contains information about the leaders among Archaea and Bacteria with respect to the number of self-targeting spacers. 
One archaeal organism _Methanobrevibacter olleyae (euryarchaeotes) YLM1_ (CP014265.1) - 18 self-targeting spacers
and two strains of _Megasphaera elsdenii DSM 20460_ (CP027570.1 and HE576794.1) - 35 self-targeting spacers found in each strain
Each file contains information about the matches found for the self-targeting found in the corresponding genome. 

We analyzed gbk files for these organisms to identify genes that bear targets for these self-targeting spacers. Each line in a file corresponds to the found relation between a self-targeting spacer and an annotated entry where the target is located. The number of lines in a file can exceed the number of the self-targeting spacers found in the corresponding genome, if at least one of the following conditions is met:

If __a spacer have more than one targets in a genome__ then the file will contain several lines with the same values for Spacer_start_pos and Spacer_end_pos but values for Target_start_pos and Target_end_pos will be different. 

If __a spacer has duplicates in a genome__ then the file will contain lines with the same values for Target_start_pos and Target_end_pos, and the values for Spacer_start_pos and Spacer_end_pos will be different for all the targets found in the genome.

If __a target of a self-targeting spacer overlaps with several annotated entries in a genome__ then then the file will contain a separate line for each of those entries. The values for Target_start_pos, Target_end_pos, Spacer_start_pos, and Spacer_end_pos will be the same for those lines, but the Entry_info, Entry_ID, Entry_start_pos, and Entry_end_pos will be different. 






