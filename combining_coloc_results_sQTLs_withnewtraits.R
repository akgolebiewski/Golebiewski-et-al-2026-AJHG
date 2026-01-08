###### combining coloc results (nt-sQTLs)

library(dplyr)
library(tidyverse)
library(data.table)


traits <- data.frame(names = c("biomarkers-30760-both_sexes-irnt", "biomarkers-30690-both_sexes-irnt", "continuous-4079-both_sexes-irnt", "continuous-4080-both_sexes-irnt", "icd10-I25-both_sexes", "continuous-LDLC-both_sexes-medadj_irnt", "categorical-6152-both_sexes-5", "categorical-6152-both_sexes-7"),
                     pheno = c("HDLC", "totalCholesterol", "diastolicBP", "systolicBP", "CAD", "LDLC", "DVT", "PE"),
                     n = c(374709, 409385, 405308, 407904, 32353+407078, 407824, 435946, 434215),
                     type = c("quant", "quant", "quant", "quant", "cc", "quant","cc", "cc"))# cc for case-control or quant for quantitatitve


### get the positions for the SNPs in hg19 to match back to gwas
geno.ann <- read.delim("/path/notx/Genotype_annotations_notx.txt", as.is=T, stringsAsFactors = F)

hg19_ann <- read.delim("/data4/vari-gene-final/hg19/QTL_Testing_files/rna/matrixeqtl/Anno_file.txt")
geno.ann %>% mutate(hg19_pos =  paste(hg19_ann$chr[match(id, hg19_ann$id)], hg19_ann$pos[match(id, hg19_ann$id)], sep = "_"),
                    hg38_pos = paste(chr, pos, sep = "_")) -> geno.ann

ntqtls <- read.delim("/path/splicing_QTLs_notx_cis100kbCis_LDpruned_genelevelfdrsig.txt", as.is = T, stringsAsFactors = F)
il1bqtls <- read.delim("/path/splicing_QTLs_il1b_cis100kbCis_LDpruned_genelevelfdrsig.txt", as.is = T, stringsAsFactors = F)

#### for every trait, read in results and make a df of sig colocalizations
for (tx in c("notx", "il1b")) {

for (i in 1:nrow(traits)) {
  filename = traits[i, "names"]
  name = traits[i, "pheno"]
  gwasname = paste0("PanUKBB_meta_", traits[i, "pheno"])
  
  
  # Import GWAS table
  gwas.full <- read_tsv(paste0("/path/gwas_summary_stats_for_coloc/", filename, ".tsv.bgz"))
  gwas.full %>% mutate(snp = paste0("chr", chr, "_", pos)) -> gwas.full
  
  # import results
  summary <- read.delim(paste0("/path/coloc/20250625_coloc_", gwasname, "_sqtl_", tx, "_summary.txt"), as.is = T, stringsAsFactors = F)
  results <- read.delim(paste0("/path/coloc/20250625_coloc_", gwasname, "_sqtl_", tx, "_results.txt"), as.is = T, stringsAsFactors = F)
  # get significant results (PP.H4 is > 0.8, 80% confidence of colocalization)
  summary %>% filter(PP.H4.abf > 0.8) -> sig_introns
  
  if (name %in% c("DVT", "PE")) {
    results %>% 
      filter(gene %in% row.names(sig_introns)) %>% group_by(gene) %>%
      mutate(rank = rank(-SNP.PP.H4)) %>% filter(rank ==1)-> sig_results
    sig_results %>%
      mutate(rsid = geno.ann$id[match(snp, geno.ann$hg38_pos)],
             hg19_pos = geno.ann$hg19_pos[match(snp, geno.ann$hg38_pos)],
             gwas_beta = gwas.full$beta_meta[match(hg19_pos, gwas.full$snp)],
             gwas_pvalue = 10^-gwas.full$neglog10_pval_meta[match(hg19_pos, gwas.full$snp)],
             intron_PP.H4.abf = summary$PP.H4.abf[match(gene, row.names(summary))],
             gwas = name,
             treatment = tx)  -> sig_results
    sig_results
  
    } else {
  
  # restrict to the top SNP for the sig loci
  results %>% 
    filter(gene %in% row.names(sig_introns)) %>% group_by(gene) %>%
    mutate(rank = rank(-SNP.PP.H4)) %>% filter(rank ==1)-> sig_results
  sig_results %>%
    mutate(rsid = geno.ann$id[match(snp, geno.ann$hg38_pos)],
           hg19_pos = geno.ann$hg19_pos[match(snp, geno.ann$hg38_pos)],
           gwas_beta = gwas.full$beta_meta_hq[match(hg19_pos, gwas.full$snp)],
           gwas_pvalue = 10^-gwas.full$neglog10_pval_meta_hq[match(hg19_pos, gwas.full$snp)],
          intron_PP.H4.abf = summary$PP.H4.abf[match(gene, row.names(summary))],
          gwas = name,
          treatment = tx)  -> sig_results
  sig_results
   }       
  
  assign(paste0("coloc_", name, "_", tx), sig_results)
}  
}

coloc_LDLC_notx$gwas = "LDLC"
coloc_HDLC_notx$gwas = "HDLC"
coloc_diastolicBP_notx$gwas = "diastolicBP"
coloc_systolicBP_notx$gwas = "systolicBP"
coloc_totalCholesterol_notx$gwas = "totalCholesterol"

coloc_LDLC_il1b$gwas = "LDLC"
coloc_HDLC_il1b$gwas = "HDLC"
coloc_diastolicBP_il1b$gwas = "diastolicBP"
coloc_systolicBP_il1b$gwas = "systolicBP"
coloc_totalCholesterol_il1b$gwas = "totalCholesterol"

coloc_LDLC_notx$treatment = "notx"
coloc_HDLC_notx$treatment = "notx"
coloc_diastolicBP_notx$treatment = "notx"
coloc_systolicBP_notx$treatment = "notx"
coloc_totalCholesterol_notx$treatment = "notx"

coloc_LDLC_il1b$treatment = "il1b"
coloc_HDLC_il1b$treatment = "il1b"
coloc_diastolicBP_il1b$treatment = "il1b"
coloc_systolicBP_il1b$treatment = "il1b"
coloc_totalCholesterol_il1b$treatment = "il1b"

listofdfs <- c(ls(pattern = "_notx"), ls(pattern = "_il1b"))

coloc_all <- bind_rows(lapply(listofdfs, get))

#write.table(coloc_all, "/path/coloc/20250625_coloc_HDLC_LDCL_totalchol_bloodpressure_panukbb_bothtreatments_summary.txt", sep = "\t", quote = F)
#coloc_all <- read.delim("/path/coloc/20250625_coloc_HDLC_LDCL_totalchol_bloodpressure_panukbb_summary.txt", as.is = T, stringsAsFactors = F)

### add CAD from VanDerHarst
cad <- read.delim(gzfile("/path/vanDerHast_GWAS/29212778-GCST005194-EFO_0000378.h.tsv.gz"), sep = "\t")
resil1b <- read.delim("/path/20250115_coloc_vanderharstCAD_sQTLil1b_results.txt", as.is = T, stringsAsFactors = F)
sumil1b <- read.delim("/path/20250115_coloc_vanderharstCAD_sQTLil1b_summary.txt", as.is = T, stringsAsFactors = F)
# summarize
sig_events <- row.names(sumil1b)[sumil1b$PP.H4.abf > 0.8]
sigil1b <- resil1b[resil1b$gene %in% sig_events,]
sigil1b %>% group_by(gene) %>% top_n(SNP.PP.H4, n = 1) -> causalil1b
cad_il1b <- as.data.frame(causalil1b[!duplicated(causalil1b$gene),])

# notx results CAD
resnotx <- read.delim("/path/20250115_coloc_vanderharstCAD_sQTLnotx_results.txt", as.is = T, stringsAsFactors = F)
sumnotx <- read.delim("/path/20250115_coloc_vanderharstCAD_sQTLnotx_summary.txt", as.is = T, stringsAsFactors = F)
# summarize
sig_events <- row.names(sumnotx)[sumnotx$PP.H4.abf > 0.8]
signotx <- resnotx[resnotx$gene %in% sig_events,]
signotx %>% group_by(gene) %>% top_n(SNP.PP.H4, n = 1) -> causalnotx
cad_notx <- as.data.frame(causalnotx[!duplicated(causalnotx$gene),])

cad_il1b %>% mutate(treatment = "il1b",
                    intron_PP.H4.abf = sumil1b$PP.H4.abf[match(gene, row.names(sumil1b))]) -> cad_il1b

cad_notx %>% mutate(treatment = "notx", intron_PP.H4.abf = sumnotx$PP.H4.abf[match(gene, row.names(sumnotx))]) %>%
  bind_rows(., cad_il1b) %>%
  mutate(gwas = "CAD",
         gwas_beta = cad$hm_beta[match(snp, cad$variant_id)],
         gwas_pvalue = cad$p_value[match(snp, cad$variant_id)],
         rsid = snp,
         hg19_pos = NA,
         rank = NA) -> cad_both

# add to the big df
coloc_all %>% 
  rbind(., cad_both[,names(cad_both)[names(cad_both) %in% names(coloc_all)]])  -> coloc_all

### add FDR and beta for the sQTLs
ntqtls <- read.delim("/path/notx/splicing_QTLs_notx_rna_cis100kbCis_genelevelfdr0.05sig.txt", as.is = T, stringsAsFactors = F)
il1bqtls <- read.delim("/path/il1b/splicing_QTLs_il1b_rna_cis100kbCis_genelevelfdr0.05sig.txt", as.is = T, stringsAsFactors = F)

ntqtls %>% mutate(snpid = unlist(lapply(strsplit(snps, ":"), "[[", 1)),
                  qtl = paste(snpid, gene, sep = '_')) -> ntqtls
il1bqtls %>% mutate(snpid = unlist(lapply(strsplit(snps, ":"), "[[", 1)),
                  qtl = paste(snpid, gene, sep = '_')) -> il1bqtls

coloc_all %>%
  mutate(snpid = unlist(lapply(strsplit(rsid, ":"), "[[", 1)),
         qtl = paste(snpid, gene, sep = '_'),
         FDR_Control = ntqtls$gene_level_FDR[match(qtl, ntqtls$qtl)],
         Beta_Control = ntqtls$beta[match(qtl, ntqtls$qtl)],
         FDR_IL1B = il1bqtls$gene_level_FDR[match(qtl, il1bqtls$qtl)],
         Beta_IL1B = il1bqtls$beta[match(qtl, il1bqtls$qtl)]) -> coloc_all


# save
write.table(coloc_all, "/path/coloc/20250722_coloc_ALLtraits_summaryofsigresults.txt", sep = "\t", quote = F, row.names = F)


