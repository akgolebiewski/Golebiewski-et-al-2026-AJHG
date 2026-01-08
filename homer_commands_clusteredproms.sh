#!/bin/bash

cd /path/erg_rela_clusteredproms/

######### find motifs in clusters 1 and 2 at Pbasal set
findMotifsGenome.pl clusters_1and2_ergrelarep_pbasal.txt hg38 clusters_1and2_ergrelarep_pbasal_pindbackground -size given -mask -bg clusters_1and2_ergrelarep_pinducible.txt

# and the reverse, at Pinducible
findMotifsGenome.pl clusters_1and2_ergrelarep_pinducible.txt hg38 clusters_1and2_ergrelarep_pinducible_pbasbackground -size given -mask -bg clusters_1and2_ergrelarep_pbasal.txt

######### find motifs in clusters 3 6 and 7 at Pinducible set
findMotifsGenome.pl clusters_3and6and7_ergrelaact_pinducible.txt hg38 clusters_3and6and7_ergrelaact_pinducible_pbasbackground -size given -mask -bg clusters_3and6and7_ergrelaact_pbasal.txt

# and the reverse, at basal
findMotifsGenome.pl clusters_3and6and7_ergrelaact_pbasal.txt hg38 clusters_3and6and7_ergrelaact_pbasal_pindbackground -size given -mask -bg clusters_3and6and7_ergrelaact_pinducible.txt
