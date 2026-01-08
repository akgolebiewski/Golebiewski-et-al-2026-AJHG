#!/bin/bash

# find motifs at inducible promoters
findMotifsGenome.pl /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_only.activating_1500bpProm_inducible.txt hg38 /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_only.activating_inducibleProm_1500bp_motifs -size given -mask -bg /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_only.activating_1500bpProm_basal.txt

findMotifsGenome.pl /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.repressing_1500bpProm_inducible.txt hg38 /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.repressing_inducibleProm_1500bp_motifs -size given -mask -bg /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.repressing_1500bpProm_basal.txt

findMotifsGenome.pl /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.activating_1500bpProm_inducible.txt hg38 /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.activating_inducibleProm_1500bp_motifs -size given -mask -bg /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.activating_1500bpProm_basal.txt

## and at basal promoters
findMotifsGenome.pl /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_only.activating_1500bpProm_basal.txt hg38 /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_only.activating_basalProm_1500bp_motifs -size given -mask -bg /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_only.activating_1500bpProm_inducible.txt

findMotifsGenome.pl /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.repressing_1500bpProm_basal.txt hg38 /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.repressing_basalProm_1500bp_motifs -size given -mask -bg /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.repressing_1500bpProm_inducible.txt

findMotifsGenome.pl /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.activating_1500bpProm_basal.txt hg38 /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.activating_basalProm_1500bp_motifs -size given -mask -bg /data3/anna/splicing_analysis_2024/promoterBindingRatios/p65_and_erg.activating_1500bpProm_inducible.txt
