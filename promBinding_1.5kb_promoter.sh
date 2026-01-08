#!/bin/bash

ntpeak="/data3/anna/splicing_analysis_2024/20241203_1500bp_promoter_AFnotx.txt"
il1bpeak="/data3/anna/splicing_analysis_2024/20241203_1500bp_promoter_AFil1b.txt"

# p65
annotatePeaks.pl $ntpeak hg38 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/p65/*.hg38.tagDir -size given -rpkm > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20241203_1500bpProm_AFnotx_byDonor_p65.txt

annotatePeaks.pl $il1bpeak hg38 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/p65/*.hg38.tagDir -size given -rpkm > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20241203_1500bpProm_AFil1b_byDonor_p65.txt

# ATAC
annotatePeaks.pl $ntpeak hg38 -d /data4/vari-gene-final/hg38/atac/tagdirs_hg38/notx/*.hg38.tagDir /data4/vari-gene-final/hg38/atac/tagdirs_hg38/il1b/*.hg38.tagDir -size given -rpkm > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20241203_1500bpProm_AFnotx_byDonor_ATAC.txt

annotatePeaks.pl $il1bpeak hg38 -d /data4/vari-gene-final/hg38/atac/tagdirs_hg38/notx/*.hg38.tagDir /data4/vari-gene-final/hg38/atac/tagdirs_hg38/il1b/*.hg38.tagDir -size given -rpkm > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20241203_1500bpProm_AFil1b_byDonor_ATAC.txt

# H3K27ac
annotatePeaks.pl $ntpeak hg38 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/h3k27ac/notx/*.hg38.tagDir /data4/vari-gene-final/hg38/chip/tagdirs_hg38/h3k27ac/il1b/*.hg38.tagDir -size given -rpkm > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20241203_1500bpProm_AFnotx_byDonor_H3K27ac.txt

annotatePeaks.pl $il1bpeak hg38 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/h3k27ac/notx/*.hg38.tagDir /data4/vari-gene-final/hg38/chip/tagdirs_hg38/h3k27ac/il1b/*.hg38.tagDir -size given -rpkm > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20241203_1500bpProm_AFil1b_byDonor_H3K27ac.txt


# ERG
annotatePeaks.pl $ntpeak hg38 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/erg/notx/*.hg38.tagDir /data4/vari-gene-final/hg38/chip/tagdirs_hg38/erg/il1b/*.hg38.tagDir -size given -rpkm > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20241203_1500bpProm_AFnotx_byDonor_ERG.txt

annotatePeaks.pl $il1bpeak hg38 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/erg/notx/*.hg38.tagDir /data4/vari-gene-final/hg38/chip/tagdirs_hg38/erg/il1b/*.hg38.tagDir -size given -rpkm > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20241203_1500bpProm_AFil1b_byDonor_ERG.txt


## RNA (proof of concept
annotatePeaks.pl $ntpeak hg38 -d /data4/vari-gene-final/hg38/rna/tagdirs/*.tagDir -size given -rpkm > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20241203_1500bpProm_AFnotx_byDonor_RNA.txt

annotatePeaks.pl $il1bpeak hg38 -d /data4/vari-gene-final/hg38/rna/tagdirs/*.tagDir -size given -rpkm > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20241203_11500bpProm_AFil1b_byDonor_RNA.txt


########## all Donors in one tagdir for histograms
# H3K27ac
annotatePeaks.pl $ntpeak hg38 -size 5000 -hist 10 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/h3k27ac/notx/ChIP-notx-hg38-h3k27ac-MERGEDALLDONORS.tagDir /data4/vari-gene-final/hg38/chip/tagdirs_hg38/h3k27ac/il1b/ChIP-il1b-hg38-h3k27ac-MERGEDALLDONORS.tagDir > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20240227_1500bpProm_AFnotx_H3K27ac_size5000_hist10.txt

annotatePeaks.pl $il1bpeak hg38 -size 5000 -hist 10 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/h3k27ac/notx/ChIP-notx-hg38-h3k27ac-MERGEDALLDONORS.tagDir /data4/vari-gene-final/hg38/chip/tagdirs_hg38/h3k27ac/il1b/ChIP-il1b-hg38-h3k27ac-MERGEDALLDONORS.tagDir > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20240227_1500bpProm_AFil1b_H3K27ac_size5000_hist10.txt

# ATAC
annotatePeaks.pl $ntpeak hg38 -size 5000 -hist 10 -d /data4/vari-gene-final/hg38/atac/tagdirs_hg38/notx/ATAC-notx-hg38-MERGEDALLDONORS.tagDir /data4/vari-gene-final/hg38/atac/tagdirs_hg38/il1b/ATAC-il1b-hg38-MERGEDALLDONORS.tagDir > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20240227_1500bpProm_AFnotx_ATAC_size5000_hist10.txt

annotatePeaks.pl $il1bpeak hg38 -size 5000 -hist 10 -d /data4/vari-gene-final/hg38/atac/tagdirs_hg38/notx/ATAC-notx-hg38-MERGEDALLDONORS.tagDir /data4/vari-gene-final/hg38/atac/tagdirs_hg38/il1b/ATAC-il1b-hg38-MERGEDALLDONORS.tagDir > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20240227_1500bpProm_AFil1b_ATAC_size5000_hist10.txt


### ERG
annotatePeaks.pl $ntpeak hg38 -size 5000 -hist 10 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/erg/notx/ChIP-notx-hg38-ERG-MERGEDALLDONORS.tagDir /data4/vari-gene-final/hg38/chip/tagdirs_hg38/erg/il1b/ChIP-il1b-hg38-erg-MERGEDALLDONORS.tagDir > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20240227_1500bpProm_AFnotx_ERG_size5000_hist10.txt

annotatePeaks.pl $il1bpeak hg38 -size 5000 -hist 10 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/erg/notx/ChIP-notx-hg38-ERG-MERGEDALLDONORS.tagDir /data4/vari-gene-final/hg38/chip/tagdirs_hg38/erg/il1b/ChIP-il1b-hg38-erg-MERGEDALLDONORS.tagDir > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20240227_1500bpProm_AFil1b_ERG_size5000_hist10.txt

### p65
annotatePeaks.pl $ntpeak hg38 -size 5000 -hist 10 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/p65/ChIP-Il1b-hg38-p65-MERGEDALLDONORS.tagDir > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20240227_1500bpProm_AFnotx_p65_size5000_hist10.txt

annotatePeaks.pl $il1bpeak hg38 -size 5000 -hist 10 -d /data4/vari-gene-final/hg38/chip/tagdirs_hg38/p65/ChIP-Il1b-hg38-p65-MERGEDALLDONORS.tagDir > /data3/anna/splicing_analysis_2024/promoterBindingRatios/20240227_1500bpProm_AFil1b_p65_size5000_hist10.txt






