#!/bin/bash
# SetBackground.sh
#Code adapted from code on Taylor lab github site

MUT="ROB_PL3_1-A01_COMPILEDJW_final.txt"
export NA=$1
cd ${NA}

# count mutations per isolate
perl ../count_PRUV_mutations_isolate.pl <${MUT}

# get fasta trinuc context for SNVs
perl ../prUV_process_trinuc_formutations.pl <${MUT}

# extract SNVs
perl ../process_prUV_mutations.pl <${MUT}

# extract DNP mutations 
perl ../prUV_process_DNPmutations.pl <${MUT}

# take care of mutation classes per isolate
#perl ../count_trinuc_UVB_mutations_isolate.pl <${MUT}

# count dnp mutations
#perl ../count_DNV_UVBmutations_isolate.pl <${MUT}
