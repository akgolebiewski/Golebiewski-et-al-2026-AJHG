library(tidyverse)
library(dplyr)

########################### NOTX ###########################
ntqtls <- read.delim("/path/notx/splicing_QTLs_notx_rna_cis100kbCis_genelevelfdr0.05sig.txt", as.is = T, stringsAsFactors = F)

#adding in the gene symbols
clust_anno<-read.delim("/path/final_hg38_cluster_significance.txt",
                     sep = "\t", header = T, as.is = T, stringsAsFactors = F)
clust_anno$cluster_base_name<-unlist(lapply(strsplit(clust_anno$cluster, ":"), "[[",2))

#ntqtls$GeneSym<-clust_anno$genes[match(ntqtls$Cluster_name, clust_anno$cluster_base_name)]

#Adding in the eQTL data
notx_eQTL<-read.delim("/path/rna/notx/20_01_27_gene.level.FDR_PROCESSED.w.Betas.txt",
                      sep = "\t", header = T, as.is = T, stringsAsFactors = F)
il1b_eQTL<-read.delim("/path/rna/il1b/20_01_27_gene.level.FDR_PROCESSED.w.Betas.txt",
                      sep = "\t", header = T, as.is = T, stringsAsFactors = F)

## add in molQTL data
notx_histQTL<-read.delim("/path/h3k27ac_chip/notx/results_PROCESSED_18_12_12_h3k27ac_notx_peakOnly_given.txt",
                      sep = "\t", header = T, as.is = T, stringsAsFactors = F)
il1b_histQTL<-read.delim("/path/h3k27ac_chip/il1b/results_PROCESSED_18_12_12_h3k27ac_il1b_peakOnly_given.txt",
                      sep = "\t", header = T, as.is = T, stringsAsFactors = F)
il1b_p65QTL <- read.delim("/path/p65_chip/il1b/results_PROCESSED_18_12_12_p65_il1b_peakOnly_given.txt", 
                          sep = "\t", header = T, as.is = T, stringsAsFactors = F)
notx_ergQTL<-read.delim("/path/erg_chip/notx/results_PROCESSED_18_12_12_erg_notx_peakOnly_given.txt",
                         sep = "\t", header = T, as.is = T, stringsAsFactors = F)
il1b_ergQTL<-read.delim("/path/erg_chip/il1b/results_PROCESSED_18_12_12_erg_il1b_peakOnly_given.txt",
                         sep = "\t", header = T, as.is = T, stringsAsFactors = F)

### add caQTLs (ATAC)
notx_caQTL<-read.delim("/path/atac/notx/results_PROCESSED_18_12_12_atac_notx_peakOnly_given.txt",
                      sep = "\t", header = T, as.is = T, stringsAsFactors = F)
il1b_caQTL<-read.delim("/path/atac/il1b/results_PROCESSED_18_12_12_atac_il1b_peakOnly_given.txt",
                      sep = "\t", header = T, as.is = T, stringsAsFactors = F)

#Add in the differential splicing info
leaf <- read.delim("/path/annotatedFINAL_leafcutter_results_53HAECsIL1Bornotx_final_set_w_depth.txt", as.is = T, stringsAsFactors = F)

####Putting in the differential expression data
diff_exp<-read.delim("/path/22_09_07_edgeR_differentialExpression-Treatment_coef_PCA.ancestry_Sex_UniTags.txt",
                     sep = "\t", header = T, as.is = T, stringsAsFactors = F)


#### add information into the big data table
#cluster = ifelse(!is.na(ntqtls$snps), unlist(lapply(strsplit(ntqtls$snp_gene, ":"), "[[", 4)), unlist(lapply(strsplit(ntqtls$gene_IL1B, ":"), "[[", 4)))

ntqtls %>%
  mutate(intron = gene,
         cluster = unlist(lapply(strsplit(intron, ":"), "[[", 4)),
         # add gene symbol
         geneSym = clust_anno$gene[match(cluster, clust_anno$cluster_base_name)],
         # add eQTLs
         eQTL_Beta_NT = notx_eQTL$Beta[match(snps, notx_eQTL$SNP)],
         eQTL_GeneLevelFDR_NT = notx_eQTL$Gene.level.FDR[match(snps, notx_eQTL$SNP)],
         eQTL_Gene_NT = notx_eQTL$Gene_sym[match(snps, notx_eQTL$SNP)],
         eQTL_Beta_IL1B = il1b_eQTL$Beta[match(snps, il1b_eQTL$SNP)],
         eQTL_GeneLevelFDR_IL1B = il1b_eQTL$Gene.level.FDR[match(snps, il1b_eQTL$SNP)],
         eQTL_Gene_IL1B = il1b_eQTL$Gene_sym[match(snps, il1b_eQTL$SNP)],
         # add caQTLs
         caQTL_FDR_NT = notx_caQTL$Q.VALUE[match(snps, notx_caQTL$rs_ID)],
         caQTL_Feature_NT = notx_caQTL$Feature_ID[match(snps, notx_caQTL$rs_ID)],
         caQTL_FDR_IL1B = il1b_caQTL$Q.VALUE[match(snps, il1b_caQTL$rs_ID)],
         caQTL_Feature_IL1B = il1b_caQTL$Feature_ID[match(snps, il1b_caQTL$rs_ID)],
         # add H3K27ac QTLs
         histQTL_FDR_NT = notx_histQTL$Q.VALUE[match(snps, notx_histQTL$rs_ID)],
         histQTL_Feature_NT = notx_histQTL$Feature_ID[match(snps, notx_histQTL$rs_ID)],
         histQTL_FDR_IL1B = il1b_histQTL$Q.VALUE[match(snps, il1b_histQTL$rs_ID)],
         histQTL_Feature_IL1B = il1b_histQTL$Feature_ID[match(snps, il1b_histQTL$rs_ID)],
         # add ERG QTLs
         ergQTL_FDR_NT = notx_ergQTL$Q.VALUE[match(snps, notx_ergQTL$rs_ID)],
         ergQTL_Feature_NT = notx_ergQTL$Feature_ID[match(snps, notx_ergQTL$rs_ID)],
         ergQTL_FDR_IL1B = il1b_ergQTL$Q.VALUE[match(snps, il1b_ergQTL$rs_ID)],
         ergQTL_Feature_IL1B = il1b_ergQTL$Feature_ID[match(snps, il1b_ergQTL$rs_ID)],
         # add RELA QTLs
         relaQTL_FDR_IL1B = il1b_p65QTL$Q.VALUE[match(snps, il1b_p65QTL$rs_ID)],
         relaQTL_Feature_IL1B = il1b_p65QTL$Feature_ID[match(snps, il1b_p65QTL$rs_ID)],
         # add differential tests
         diff_expn_LogFC = diff_exp$LogFC[match(geneSym, diff_exp$Gene.sym)],
         diff_expn_FDR = diff_exp$FDR[match(geneSym, diff_exp$Gene.sym)],
         diff_spl_dPSI = leaf$deltapsi[match(intron, gsub("chr", "", leaf$intron))],
         diff_spl_p.adjust = leaf$p.adjust[match(intron, gsub("chr", "", leaf$intron))],
         spliceType = leaf$spliceType[match(intron, gsub("chr", "", leaf$intron))]
         ) -> combined

################### COMPLETE SPLICE TYPE ANNOTATION--if they were not tested for Diff Spl there isn't an annotation

### add in splice type annotations for ALL introns tested in the DSG analysis
combined %>%
  filter(is.na(spliceType)) %>%
  select(intron) %>% distinct() %>% 
  mutate(junction = vapply(strsplit(intron, ":"), function(x) paste(x[c(2,3)], collapse = "-"), character(1L)),
         chr = vapply(strsplit(intron, ":"), function(x) paste(x[c(1)]), character(1L))) %>%
  tibble() -> introns

#### annotating the results that weren't tested for diff splicing but were tested for sQTLs:
suppa.AF <- read.delim("/path/variGene_diffSplice_SUPPA_AF_RESULTS.txt", as.is = T, stringsAsFactors=F)
suppa.AL <- read.delim("/path/variGene_diffSplice_SUPPA_AL_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.A3 <- read.delim("/path/variGene_diffSplice_SUPPA_A3_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.A5 <- read.delim("/path/variGene_diffSplice_SUPPA_A5_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.MX <- read.delim("/path/variGene_diffSplice_SUPPA_MX_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.RI <- read.delim("/path/variGene_diffSplice_SUPPA_RI_RESULTS.txt", as.is = T, stringsAsFactors=T)
## the names in this file are not correct, fix before annotating
names(suppa.RI) <- c("gene_ID", "gene_sym", "chr", "s1", "strand", "e2", "e1.s2", "psiPerLocalEvent_dPSI", "psiPerLocalEvent_pValue")
suppa.SE <- read.delim("/path/variGene_diffSplice_SUPPA_SE_RESULTS.txt", as.is = T, stringsAsFactors=T)
## SE splices need to add e1-s3 (the "long" splice)
e1 <- lapply(strsplit(suppa.SE$e1.s2, "-"), "[[", 1)
e3 <- lapply(strsplit(suppa.SE$e2.s3, "-"), "[[", 2)
suppa.SE$e1.s3 <- paste(e1, e3, sep = "-")

## add splice type to the intron table
introns$SUPPA_annotation <- ifelse ( introns$junction %in% suppa.AF$e1.s3 | introns$junction %in% suppa.AF$e2.s3, paste("AF"),
                                  
                                  ifelse ( introns$junction %in% suppa.AL$e1.s2 | introns$junction %in% suppa.AL$e1.s3, paste("AL"),
                                           
                                           ifelse ( introns$junction %in% suppa.A3$e1.s2 | introns$junction %in% suppa.A3$e1.s3, paste("A3"), 
                                                    
                                                    ifelse ( introns$junction %in% suppa.A5$e1.s3 | introns$junction %in% suppa.A5$e2.s3, paste("A5"),
                                                             
                                                             ifelse ( introns$junction %in% suppa.MX$e1.s2 | introns$junction %in% suppa.MX$e2.s4 | introns$junction %in% suppa.MX$e1.s3 | introns$junction %in% suppa.MX$e3.s4, paste("MX"),
                                                                      
                                                                      ifelse ( introns$junction %in% suppa.RI$e1.s2, paste("RI"),
                                                                               
                                                                               ifelse ( introns$junction %in% suppa.SE$e1.s2 | introns$junction %in% suppa.SE$e2.s3 | introns$junction %in% suppa.SE$e1.s3, paste("SE"),
                                                                                        NA)))))))
table(introns$SUPPA_annotation)


### impute the remaining annotations based on their cluster annotation
introns %>%
  mutate(cluster = vapply(strsplit(intron, ":"), function(x) paste(x[c(4)]), character(1L)),
         cluster_ann = combined$spliceType[match(cluster, combined$cluster)],
         final_ann = ifelse(is.na(SUPPA_annotation), cluster_ann, SUPPA_annotation),
         final_ann = ifelse(is.na(final_ann), "cryptic: unannotated", final_ann)) -> introns

#save into the final table so we know what was/wasn't imputed
combined$spliceType = ifelse(is.na(combined$spliceType), introns$final_ann[match(combined$intron, introns$intron)], combined$spliceType)

####### save the results
write.table(combined, "/path/2025_06_26_final_sQTL_results_LDpruned_with_molQTLs_eQTLs_diffSPl_diffExpn_GWAS.txt", sep = "\t", col.names = T, quote = F)



########################### IL1B ###########################
ntqtls <- read.delim("/path/il1b/splicing_QTLs_il1b_rna_cis100kbCis_genelevelfdr0.05sig.txt", as.is = T, stringsAsFactors = F)

#adding in the gene symbols
clust_anno<-read.delim("/path/final_hg38_cluster_significance.txt",
                       sep = "\t", header = T, as.is = T, stringsAsFactors = F)
clust_anno$cluster_base_name<-unlist(lapply(strsplit(clust_anno$cluster, ":"), "[[",2))

#ntqtls$GeneSym<-clust_anno$genes[match(ntqtls$Cluster_name, clust_anno$cluster_base_name)]

#Adding in the eQTL data
notx_eQTL<-read.delim("/path/rna/notx/20_01_27_gene.level.FDR_PROCESSED.w.Betas.txt",
                      sep = "\t", header = T, as.is = T, stringsAsFactors = F)
il1b_eQTL<-read.delim("/path/rna/il1b/20_01_27_gene.level.FDR_PROCESSED.w.Betas.txt",
                      sep = "\t", header = T, as.is = T, stringsAsFactors = F)

## add in molQTL data
notx_histQTL<-read.delim("/path/h3k27ac_chip/notx/results_PROCESSED_18_12_12_h3k27ac_notx_peakOnly_given.txt",
                         sep = "\t", header = T, as.is = T, stringsAsFactors = F)
il1b_histQTL<-read.delim("/path/h3k27ac_chip/il1b/results_PROCESSED_18_12_12_h3k27ac_il1b_peakOnly_given.txt",
                         sep = "\t", header = T, as.is = T, stringsAsFactors = F)
il1b_p65QTL <- read.delim("/path/p65_chip/il1b/results_PROCESSED_18_12_12_p65_il1b_peakOnly_given.txt", 
                          sep = "\t", header = T, as.is = T, stringsAsFactors = F)
notx_ergQTL<-read.delim("/path/erg_chip/notx/results_PROCESSED_18_12_12_erg_notx_peakOnly_given.txt",
                        sep = "\t", header = T, as.is = T, stringsAsFactors = F)
il1b_ergQTL<-read.delim("/path/erg_chip/il1b/results_PROCESSED_18_12_12_erg_il1b_peakOnly_given.txt",
                        sep = "\t", header = T, as.is = T, stringsAsFactors = F)

### add caQTLs (ATAC)
notx_caQTL<-read.delim("/path/atac/notx/results_PROCESSED_18_12_12_atac_notx_peakOnly_given.txt",
                       sep = "\t", header = T, as.is = T, stringsAsFactors = F)
il1b_caQTL<-read.delim("/path/atac/il1b/results_PROCESSED_18_12_12_atac_il1b_peakOnly_given.txt",
                       sep = "\t", header = T, as.is = T, stringsAsFactors = F)

#Add in the differential splicing info
leaf <- read.delim("/path/annotatedFINAL_leafcutter_results_53HAECsIL1Bornotx_final_set_w_depth.txt", as.is = T, stringsAsFactors = F)

####Putting in the differential expression data
diff_exp<-read.delim("/path/22_09_07_edgeR_differentialExpression-Treatment_coef_PCA.ancestry_Sex_UniTags.txt",
                     sep = "\t", header = T, as.is = T, stringsAsFactors = F)


#### add information into the big data table
#cluster = ifelse(!is.na(ntqtls$snps), unlist(lapply(strsplit(ntqtls$snp_gene, ":"), "[[", 4)), unlist(lapply(strsplit(ntqtls$gene_IL1B, ":"), "[[", 4)))

ntqtls %>%
  mutate(intron = gene,
         cluster = unlist(lapply(strsplit(intron, ":"), "[[", 4)),
         # add gene symbol
         geneSym = clust_anno$gene[match(cluster, clust_anno$cluster_base_name)],
         # add eQTLs
         eQTL_Beta_NT = notx_eQTL$Beta[match(snps, notx_eQTL$SNP)],
         eQTL_GeneLevelFDR_NT = notx_eQTL$Gene.level.FDR[match(snps, notx_eQTL$SNP)],
         eQTL_Gene_NT = notx_eQTL$Gene_sym[match(snps, notx_eQTL$SNP)],
         eQTL_Beta_IL1B = il1b_eQTL$Beta[match(snps, il1b_eQTL$SNP)],
         eQTL_GeneLevelFDR_IL1B = il1b_eQTL$Gene.level.FDR[match(snps, il1b_eQTL$SNP)],
         eQTL_Gene_IL1B = il1b_eQTL$Gene_sym[match(snps, il1b_eQTL$SNP)],
         # add caQTLs
         caQTL_FDR_NT = notx_caQTL$Q.VALUE[match(snps, notx_caQTL$rs_ID)],
         caQTL_Feature_NT = notx_caQTL$Feature_ID[match(snps, notx_caQTL$rs_ID)],
         caQTL_FDR_IL1B = il1b_caQTL$Q.VALUE[match(snps, il1b_caQTL$rs_ID)],
         caQTL_Feature_IL1B = il1b_caQTL$Feature_ID[match(snps, il1b_caQTL$rs_ID)],
         # add H3K27ac QTLs
         histQTL_FDR_NT = notx_histQTL$Q.VALUE[match(snps, notx_histQTL$rs_ID)],
         histQTL_Feature_NT = notx_histQTL$Feature_ID[match(snps, notx_histQTL$rs_ID)],
         histQTL_FDR_IL1B = il1b_histQTL$Q.VALUE[match(snps, il1b_histQTL$rs_ID)],
         histQTL_Feature_IL1B = il1b_histQTL$Feature_ID[match(snps, il1b_histQTL$rs_ID)],
         # add ERG QTLs
         ergQTL_FDR_NT = notx_ergQTL$Q.VALUE[match(snps, notx_ergQTL$rs_ID)],
         ergQTL_Feature_NT = notx_ergQTL$Feature_ID[match(snps, notx_ergQTL$rs_ID)],
         ergQTL_FDR_IL1B = il1b_ergQTL$Q.VALUE[match(snps, il1b_ergQTL$rs_ID)],
         ergQTL_Feature_IL1B = il1b_ergQTL$Feature_ID[match(snps, il1b_ergQTL$rs_ID)],
         # add RELA QTLs
         relaQTL_FDR_IL1B = il1b_p65QTL$Q.VALUE[match(snps, il1b_p65QTL$rs_ID)],
         relaQTL_Feature_IL1B = il1b_p65QTL$Feature_ID[match(snps, il1b_p65QTL$rs_ID)],
         # add differential tests
         diff_expn_LogFC = diff_exp$LogFC[match(geneSym, diff_exp$Gene.sym)],
         diff_expn_FDR = diff_exp$FDR[match(geneSym, diff_exp$Gene.sym)],
         diff_spl_dPSI = leaf$deltapsi[match(intron, gsub("chr", "", leaf$intron))],
         diff_spl_p.adjust = leaf$p.adjust[match(intron, gsub("chr", "", leaf$intron))],
         spliceType = leaf$spliceType[match(intron, gsub("chr", "", leaf$intron))]
  ) -> combined_il1b

################### COMPLETE SPLICE TYPE ANNOTATION--if they were not tested for Diff Spl there isn't an annotation

### add in splice type annotations for ALL introns tested in the DSG analysis
combined_il1b %>%
  filter(is.na(spliceType)) %>%
  select(intron) %>% distinct() %>% 
  mutate(junction = vapply(strsplit(intron, ":"), function(x) paste(x[c(2,3)], collapse = "-"), character(1L)),
         chr = vapply(strsplit(intron, ":"), function(x) paste(x[c(1)]), character(1L))) %>%
  tibble() -> introns

#### annotating the results that weren't tested for diff splicing but were tested for sQTLs:
suppa.AF <- read.delim("/path/variGene_diffSplice_SUPPA_AF_RESULTS.txt", as.is = T, stringsAsFactors=F)
suppa.AL <- read.delim("/path/variGene_diffSplice_SUPPA_AL_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.A3 <- read.delim("/path/variGene_diffSplice_SUPPA_A3_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.A5 <- read.delim("/path/variGene_diffSplice_SUPPA_A5_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.MX <- read.delim("/path/variGene_diffSplice_SUPPA_MX_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.RI <- read.delim("/path/variGene_diffSplice_SUPPA_RI_RESULTS.txt", as.is = T, stringsAsFactors=T)
## the names in this file are not correct, fix before annotating
names(suppa.RI) <- c("gene_ID", "gene_sym", "chr", "s1", "strand", "e2", "e1.s2", "psiPerLocalEvent_dPSI", "psiPerLocalEvent_pValue")
suppa.SE <- read.delim("/path/variGene_diffSplice_SUPPA_SE_RESULTS.txt", as.is = T, stringsAsFactors=T)
## SE splices need to add e1-s3 (the "long" splice)
e1 <- lapply(strsplit(suppa.SE$e1.s2, "-"), "[[", 1)
e3 <- lapply(strsplit(suppa.SE$e2.s3, "-"), "[[", 2)
suppa.SE$e1.s3 <- paste(e1, e3, sep = "-")

## add splice type to the intron table
introns$SUPPA_annotation <- ifelse ( introns$junction %in% suppa.AF$e1.s3 | introns$junction %in% suppa.AF$e2.s3, paste("AF"),
                                     
                                     ifelse ( introns$junction %in% suppa.AL$e1.s2 | introns$junction %in% suppa.AL$e1.s3, paste("AL"),
                                              
                                              ifelse ( introns$junction %in% suppa.A3$e1.s2 | introns$junction %in% suppa.A3$e1.s3, paste("A3"), 
                                                       
                                                       ifelse ( introns$junction %in% suppa.A5$e1.s3 | introns$junction %in% suppa.A5$e2.s3, paste("A5"),
                                                                
                                                                ifelse ( introns$junction %in% suppa.MX$e1.s2 | introns$junction %in% suppa.MX$e2.s4 | introns$junction %in% suppa.MX$e1.s3 | introns$junction %in% suppa.MX$e3.s4, paste("MX"),
                                                                         
                                                                         ifelse ( introns$junction %in% suppa.RI$e1.s2, paste("RI"),
                                                                                  
                                                                                  ifelse ( introns$junction %in% suppa.SE$e1.s2 | introns$junction %in% suppa.SE$e2.s3 | introns$junction %in% suppa.SE$e1.s3, paste("SE"),
                                                                                           NA)))))))
table(introns$SUPPA_annotation)


### impute the remaining annotations based on their cluster annotation
introns %>%
  mutate(cluster = vapply(strsplit(intron, ":"), function(x) paste(x[c(4)]), character(1L)),
         cluster_ann = combined_il1b$spliceType[match(cluster, combined_il1b$cluster)],
         final_ann = ifelse(is.na(SUPPA_annotation), cluster_ann, SUPPA_annotation),
         final_ann = ifelse(is.na(final_ann), "cryptic: unannotated", final_ann)) -> introns

#save into the final table so we know what was/wasn't imputed
combined_il1b$spliceType = ifelse(is.na(combined_il1b$spliceType), introns$final_ann[match(combined_il1b$intron, introns$intron)], combined_il1b$spliceType)

####### save the results
write.table(combined_il1b, "/path/2025_06_26_final_IL1B_sQTL_results_LDpruned_with_molQTLs_eQTLs_diffSPl_diffExpn_GWAS.txt", sep = "\t", col.names = T, quote = F)

#################################COMBINE THE TWO TREATMENTS ############################
combined_il1b$sqtl_treatment = "IL1B"
combined %>% mutate(sqtl_treatment = "notx") %>%
 rbind(combined_il1b[,names(combined_il1b)[names(combined_il1b) %in% c(names(combined), "sqtl_treatment")]]) %>%
  tibble() -> combined

### save
write.table(combined, "/path/2025_06_30_final_NTandIL1B_sQTL_results_FDR0.05_with_molQTLs_eQTLs_diffSPl_diffExpn.txt", sep = "\t", col.names = T, quote = F, row.names = F)

##################### VISUALIZE RESULTS ##########################
##### add distance to junction and TSS
geno.ann <- read.delim("notx/Genotype_annotations_notx.txt", as.is=T, stringsAsFactors = F)

combined %>% ungroup() %>% 
  mutate(chr = geno.ann$chr[match(snps, geno.ann$id)], 
         pos = geno.ann$pos[match(snps, geno.ann$id)]) -> combined


hg38 <- readGFF("/path/anno_files/hg38_anno/gencode.v26.annotation.gtf")

firstExons <- hg38[hg38$exon_number == "1" & !is.na(hg38$exon_number) & hg38$type == 'exon',]

##### now match up the first exons to the introns

firstExons <- firstExons[firstExons$gene_name %in% combined$geneSym,]

combined %>%
  mutate(start = as.numeric(unlist(lapply(strsplit(gene, split = ":"), "[[", 2))),
         end = as.numeric(unlist(lapply(strsplit(gene, split = ":"), "[[", 3))),
         strand = hg38$strand[match(geneSym, hg38$gene_name)],
         gene = geneSym) %>%
  select(chr, start, end, strand, gene) -> leaf.pos


FE.pos <- data.frame(chr = firstExons$seqid,
                     start = firstExons$start,
                     end = firstExons$end,
                     strand = firstExons$strand,
                     gene = firstExons$gene_name
)

#make granges objects
leaf.grange <- makeGRangesFromDataFrame(leaf.pos)
FE.grange <- makeGRangesFromDataFrame(FE.pos)
# match the nearest FE
nearest <- nearest(leaf.grange, FE.grange, ignore.strand = F)

#add to dataframe
combined$FE_start <- FE.pos$start[match(nearest, row.names(FE.pos))]
combined$FE_end <- FE.pos$end[match(nearest, row.names(FE.pos))]

combined %>% 
  mutate(strand = hg38$strand[match(geneSym, hg38$gene_name)], 
         intron_start = as.numeric(unlist(lapply(strsplit(gene, split = ":"), "[[", 2))),
         intron_end = as.numeric(unlist(lapply(strsplit(gene, split = ":"), "[[", 3))), # add strand and positions
         ### calculate distance to TSS
         distanceToTSS = ifelse(strand == "+" & !is.na(strand), pos - FE_start, pos - FE_end),
         ### calculate a distance to the 5' and 3' end of the intron 
         distanceto5prime = intron_start - pos,
         
         distanceto3prime = intron_end - pos
  ) -> combined_ann

### save annotated result
write.table(combined_ann, "2025_06_26_final_IL1BandNotx_sQTL_results_LDpruned_with_molQTLs_eQTLs_diffSPl_diffExpn_GWAS.txt", sep = "\t", quote = F)
#combined_ann <- read.delim("2025_06_26_final_IL1BandNotx_sQTL_results_LDpruned_with_molQTLs_eQTLs_diffSPl_diffExpn_GWAS.txt", as.is = T, stringsAsFactors = F)

## make a version of this that is with all three in one plot--it'll be easier to see and there really isn't a functional difference between the 5' and 3' junction (at least that I see)
bins = c(seq(0, 1000, by = 100), seq(1000, 5000, by = 1000), seq(5000, 20000, by = 5000), seq(20000, 100000, by = 10000))
bins = bins[!duplicated(bins)]
combined_ann %>% select(snps, distanceToTSS, distanceto5prime, distanceto3prime, sqtl_treatment, spliceType) %>% 
  pivot_longer(names_to = "site", cols = c("distanceToTSS", "distanceto5prime", "distanceto3prime"), values_to = "distance") %>% 
  mutate(site = gsub("to", " to ", site), site = gsub("distanceToTSS", "distance to TSS", site), category = ifelse(site == "distance to TSS", "TSS", "junction")) %>%
  mutate(abs_distance = abs(distance),
         binned_distance = cut(abs_distance, breaks = bins, labels = bins[-1], include.lowest = TRUE)) %>%
  group_by(binned_distance, site,sqtl_treatment, spliceType) %>%
  summarise(count = n()) %>% group_by(site, sqtl_treatment) %>% 
  mutate(percent = count / sum(count) * 100) -> bincounts

# Plot the line chart
bincounts$binned_distance <-as.numeric(as.character(bincounts$binned_distance))

ggplot(bincounts, aes(x = log10(binned_distance), y = count, color = site)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "sQTL distance from TSS and splice junctions",
       x = "basepairs",
       y = "Proportion") +
  scale_x_continuous(breaks = unlist(lapply(c(100, 500, 1000, 2500,  5000, 10000, 20000, 50000, 100000), log10)), labels = c("100", "500", "1000", "2500", "5000", "10000", "20000", "50000", "100000"), transform = "log10") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + facet_grid(cols = vars(sqtl_treatment), rows = vars(spliceType))


## sQTLs are either 1) close to the TSS, 2) close to the junction
combined_ann %>% select(snps, distanceToTSS, distanceto5prime, distanceto3prime, sqtl_treatment) %>%
  pivot_longer(names_to = "site", cols = c("distanceToTSS", "distanceto5prime", "distanceto3prime"), values_to = "distance")%>%
  mutate(abs_distance = abs(distance)) %>% 
  group_by(snps) %>% slice_max(order_by = abs_distance, n = 1) %>%
  mutate(binned_distance = cut(abs_distance, breaks = bins, labels = bins[-1], include.lowest = TRUE)) %>%
  group_by(binned_distance, site, sqtl_treatment) %>%
  summarise(count = n()) %>% group_by(site, sqtl_treatment) %>%
  mutate(percent = count / sum(count) * 100,
         binned_distance = as.numeric(as.character(binned_distance)),
         category = gsub("distanceto", "closest to ", gsub("distanceTo", "distanceto", site))) -> bincounts_grouped


ggplot(bincounts_grouped, aes(x = log10(binned_distance), y = percent, color = category)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "sQTL distance from TSS and splice junctions", subtitle = "by closest site for each SNP",
       x = "basepairs",
       y = "Proportion") +
  scale_x_continuous(breaks = unlist(lapply(c(100, 500, 1000, 2500,  5000, 10000, 20000, 50000, 100000), log10)), labels = c("100", "500", "1000", "2500", "5000", "10000", "20000", "50000", "100000"), transform = "log10") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90),, text = element_text(size =12)) + ylim(0,12.5) + facet_wrap(~sqtl_treatment)

combined_ann %>% select(snps, intron, distanceToTSS, distanceto5prime, distanceto3prime, spliceType, sqtl_treatment) %>%
  pivot_longer(names_to = "site", cols = c("distanceToTSS", "distanceto5prime", "distanceto3prime"), values_to = "distance") %>% mutate(abs_distance = abs(distance)) %>% 
  group_by(snps, intron, sqtl_treatment) %>% slice_max(order_by = abs_distance, n = 1) -> nearest_ann

nearest_ann %>% ungroup() %>% group_by(site, sqtl_treatment) %>% summarise(count = n()) %>%
  mutate(category = gsub("distanceto", "closest to ", gsub("distanceTo", "distanceto", site)), closestCategory = "") %>% ungroup()-> toplotbargraph

ggplot(toplotbargraph, aes(x = category, y = count, fill = category)) + geom_bar(position="dodge", stat="identity") + theme_bw() + ggtitle("sQTL category") + xlab("") + ylab("Number of sQTLs") + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size =12)) + facet_wrap(~sqtl_treatment)

### how many are in the 1.5 kb promoter?
head(combined_ann)
combined_ann %>%
  mutate(unstranded_distToTSS = ifelse(strand == "+" & !is.na(strand), distanceToTSS, -1*distanceToTSS),
         direction = ifelse(unstranded_distToTSS > 0, "upstream", "downstream"),
         SNPInProm = ifelse(direction == "upstream" & abs(unstranded_distToTSS) <=1000 | direction == "downstream" & abs(unstranded_distToTSS) <=500, "yes", "no")) %>%
  group_by(SNPInProm) %>% summarise(count = sum(SNPInProm == "yes"))


######## sQTL and other QTL overlaps
library(UpSetR)

upsetdf <- data.frame(
  sQTL = ifelse(combined$gene_level_FDR[combined$sqtl_treatment == "notx"] < 0.05, 1, 0),
  eQTL = ifelse(combined$eQTL_GeneLevelFDR_NT[combined$sqtl_treatment == "notx"] < 0.05 | combined$eQTL_GeneLevelFDR_IL1B[combined$sqtl_treatment == "notx"] < 0.05, 1,0),
  caQTL = ifelse(combined$caQTL_FDR_NT [combined$sqtl_treatment == "notx"] < 0.05 | combined$caQTL_FDR_IL1B[combined$sqtl_treatment == "notx"] < 0.05, 1,0),
  H3K27acQTL = ifelse(combined$histQTL_FDR_NT[combined$sqtl_treatment == "notx"] < 0.05 | combined$histQTL_FDR_IL1B[combined$sqtl_treatment == "notx"] < 0.05, 1,0),
  ERGbQTL = ifelse(combined$ergQTL_FDR_NT[combined$sqtl_treatment == "notx"] < 0.05 | combined$ergQTL_FDR_IL1B[combined$sqtl_treatment == "notx"] < 0.05, 1,0),
  RELAbQTl = ifelse(combined$relaQTL_FDR_IL1B[combined$sqtl_treatment == "notx"] < 0.05, 1,0)
)

upsetdf <- replace(upsetdf,is.na(upsetdf),0)

upset(upsetdf, order.by = "freq", nsets = 6)

####### venn diagram of overlap between control and il1b
library(ggvenn)
vennlist = list(
  Control = paste(combined$snps[combined$sqtl_treatment == "notx"], combined$gene[combined$sqtl_treatment == "notx"], sep = "_"),
  IL1B = paste(combined$snps[combined$sqtl_treatment == "IL1B"], combined$gene[combined$sqtl_treatment == "IL1B"], sep = "_")
)
ggvenn::ggvenn(vennlist)

# or just by intron-- this makes more sense to compare since we took one signal per intron/locus
vennlist = list(
  Control = combined$gene[combined$sqtl_treatment == "notx"],
  IL1B = combined$gene[combined$sqtl_treatment == "IL1B"]
)

ggvenn::ggvenn(vennlist, fill_color = c("steelblue", "darkred"))


### bargraph of splice types in sQTLs
combined_ann$spliceType_toplot <- factor(ifelse(as.character(combined_ann$spliceType) %in% c("cryptic_threeprime", "cryptic_fiveprime", "cryptic_unanchored", "novel annotated pair", "unknown_strand"), "cryptic", as.character(combined_ann$spliceType)), levels =c("A3", "A5", "AF", "AL", "MX", "RI", "SE", "cryptic"))

palette <- c("#76E1B2", "#E56E63", "#BB46E0", "#DCC7D8", "#D0DEAF", "#DDAB69", "#83A1DA", "#96DADD", "#9272D5", "#85E359", "#DE77B6", "#7E7D73", "#D8E068")


ggplot(combined_ann[combined_ann$spliceType %in% combined_ann$spliceType_toplot,], aes(x = factor(spliceType, levels = c("RI", "MX", "A5", "A3", "AL", "AF", "SE")), fill = spliceType)) + 
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat = "count", aes(label = ..count..), hjust = -0.1) +
  theme_bw() + ylim(0,1550)+
  scale_fill_manual(values = palette) + 
  ylab("# of sIntrons") + xlab("") + ggtitle("sIntrons (from sQTL mapping)")  + theme(text = element_text(size = 15), axis.text = element_text(size = 13), legend.position = "none") +
  coord_flip()  -> sGenes

# plot DSTs too so they're the same size and fit next to each other
leaf <- read.delim("annotatedFINAL_leafcutter_results_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig.txt", as.is = T, stringsAsFactors = F)
leaf$spliceType_toplot <- factor(ifelse(as.character(leaf$spliceType) %in% c("cryptic_threeprime", "cryptic_fiveprime", "cryptic_unanchored", "novel annotated pair", "unknown_strand"), "cryptic", as.character(leaf$spliceType)), levels =c("A3", "A5", "AF", "AL", "MX", "RI", "SE", "cryptic"))

palette <- c("#76E1B2", "#E56E63", "#BB46E0", "#DCC7D8", "#D0DEAF", "#DDAB69", "#83A1DA", "#96DADD", "#9272D5", "#85E359", "#DE77B6", "#7E7D73", "#D8E068")


ggplot(leaf[leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05 & !leaf$spliceType_toplot == "cryptic",], aes(x = factor(spliceType, levels = c("RI", "MX", "A5", "A3", "AL", "AF", "SE")), fill = spliceType)) + 
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat = "count", aes(label = ..count..), hjust = -0.1) +
  theme_bw() + ylim(0,415) +
  scale_fill_manual(values = palette) + 
  ylab("# of DSTs") + xlab("") + ggtitle("DSTs (Control vs IL1B)")  + theme(text = element_text(size = 15), axis.text = element_text(size = 13), legend.position = "none") +
  coord_flip()->  splicetype

require(gridExtra)
grid.arrange(splicetype, sGenes, ncol = 2)

