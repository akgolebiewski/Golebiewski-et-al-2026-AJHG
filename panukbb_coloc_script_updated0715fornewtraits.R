#### coloc with GWAS from Pan-UK UKBB study
library(coloc)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(rtracklayer)


setwd("/path/")

#### set up genotypes and MAF
geno <- read.delim("/path/notx/genotypes_NT_ordered.txt", as.is = T, stringsAsFactors = F)
names(geno)
names(geno)[names(geno) == "X"] <- "SNP"

geno.ann <- read.delim("/path/notx/Genotype_annotations_notx.txt", as.is=T, stringsAsFactors = F)
#put IDs in the cad.sig.siggwas formatting
#refalt <- vapply(strsplit(geno.ann$id, split = ":"), function(x) paste(x[c(3,4)], collapse = ":"), character(1L))
chr <- gsub("chr", "", geno.ann$chr)
geno.ann$variant_id <- paste(chr, geno.ann$pos, sep = "_")
head(geno.ann$variant_id)
#add to geno df
geno$variant_id <- geno.ann$variant_id[match(geno$SNP, geno.ann$id)]


### calculate allele frequencies
geno$test <- rowSums(geno[, 2:54])
geno$af <- geno$test/106
geno$maf <- ifelse(geno$af < 0.5, geno$af, 1-geno$af)

#### set up GWAS info

### GWAS from PanUKBB study (https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit?gid=1450719288#gid=1450719288)
#30760 = HDLC
#30690 = total cholesterol
#4079 = diastolic blood pressure
#4080 = systolic blood pressure

traits <- data.frame(names = c("categorical-6152-both_sexes-5", "categorical-6152-both_sexes-7", "phecode-415.21-both_sexes", "continuous-22671-both_sexes-irnt", "icd10-I70-both_sexes", "icd10-I73-both_sexes"),
                     pheno = c("DVT", "PE", "PPH", "meanCarotidThickness", "atherosclerosis", "PAD"),
                     n = c(435946, 434215, 412289, 45984, 422216, 440821),
                     type = c("cc", "cc", "cc", "quant", "cc", "cc"))# cc for case-control or quant for quantitatitve

#### for every trait, run the coloc test for both notx and IL1B sQTLs  


for (i in 2:nrow(traits)) {
  filename = traits[i, "names"]
  
  gwasname = paste0("PanUKBB_meta_", traits[i, "pheno"])
  n = traits[i, "n"]
  gwas.type = traits[i, "type"] 
  
############################## GWAS ################################

# Import GWAS table
gwas.full <- read_tsv(paste0("gwas_summary_stats_for_coloc/", filename, ".tsv.bgz"))
# sum AF
#gwas.full$af_meta <- rowMeans(gwas.full[, c("af_cases_meta", "af_controls_meta")])
  
## this GWAS was performed in hg19/GRCh37 and we need to lift the positions to hg38
  path = "gwas_summary_stats_for_coloc/hg19ToHg38.over.chain"
  ch = import.chain(path)
  ch
  
  gwas_hg19 <- makeGRangesFromDataFrame(gwas.full, start.field = "pos", end.field = "pos", keep.extra.columns = T)
  
  seqlevelsStyle(gwas_hg19) = "UCSC"  # necessary
  gwas_lifted = liftOver(gwas_hg19, ch)
  gwas_hg38 <- as_tibble(gwas_lifted) 


# pull relevant information with standard names
  if (gwas.type == "cc") {
gwas <- data.frame(
  chr = gwas_hg38$seqnames,
  pos = gwas_hg38$start,
  pvalue = 10^(-gwas_hg38$neglog10_pval_meta_hq),
  beta = gwas_hg38$beta_meta_hq,
  se = gwas_hg38$se_meta_hq,
  af = ifelse(gwas_hg38$af_cases_meta_hq > 0.5, 1-gwas_hg38$af_cases_meta_hq , gwas_hg38$af_cases_meta_hq),
  variant_id = paste(gwas_hg38$seqnames, gwas_hg38$start, sep = "_")
)
  } else {
    gwas <- data.frame(
      chr = gwas_hg38$seqnames,
      pos = gwas_hg38$start,
      pvalue = 10^(-gwas_hg38$neglog10_pval_meta_hq),
      beta = gwas_hg38$beta_meta_hq,
      se = gwas_hg38$se_meta_hq,
      af = ifelse(gwas_hg38$af_meta_hq > 0.5, 1-gwas_hg38$af_meta_hq , gwas_hg38$af_meta_hq),
      variant_id = paste(gwas_hg38$seqnames, gwas_hg38$start, sep = "_")
    )
}

############################### read in QTL results ###########################
#treatments <- c("notx", "il1b")
treatments = c("il1b")
for (tx in treatments) {
  # these are the snps with no pvalue restiriction! why: we want LD for all the signals in the region
  # BUT every intron in here has one sig SNP by gene_level_FDR
  qtl <- read.delim(paste0("/path/", tx, "/splicing_QTLs_", tx, "_rna_cis100kbCis_allIntronsWithASigSNP_forSuSiE.txt"), as.is = T, stringsAsFactors = F)
  
  qtl %>%
    # add SE and fix variant IDs to match gwas
    mutate(se_beta = beta/statistic)  %>%
    # add chr and position
    inner_join(., geno[, c("variant_id", "maf", "SNP")], by = c("snps" = "SNP")) %>%
    mutate(chr = geno.ann$chr[match(variant_id, geno.ann$variant_id)],
           pos = geno.ann$pos[match(variant_id, geno.ann$variant_id)],
           variant_id = paste0("chr", variant_id)) %>%
    #filter for only SNPs in both gwas and sQTL
    filter(variant_id %in% gwas$variant_id) %>% tibble() -> qtl_var
  
  ## also restrict gwas to sites that match
  gwas <- gwas[gwas$variant_id %in% qtl_var$variant_id,]
  
############################### Set up ranges ##################################
# all splicing introns tested that have at least one significant QTL
introns <- unique(qtl_var$gene)

# Define range as within 100 kb of the splice event
starts <- as.numeric(unlist(lapply(strsplit(qtl_var$gene, ":"), "[[", 2)))
ends <- as.numeric(unlist(lapply(strsplit(qtl_var$gene, ":"), "[[", 3)))

gwas.range <- data.frame(chr = qtl_var$chr, start = starts - 100000, end = ends + 100000, gene = qtl_var$gene)
gwas.range <- gwas.range[!duplicated(gwas.range$gene),] # only one per gene


############################### Perform test ###################################
# open lists to write results into
resultslist<- list()
summarylist <- list()

for (x in introns) {
  if (!x == "3:98770956:98771121:clu_6636_NA") {
  print(x)
  snps <- qtl_var$variant_id[qtl_var$gene == x]
  range <- gwas.range[gwas.range$gene == x,]
  
  # Apply range to each dataset
  gwas.small <- gwas %>% filter(gwas$chr == range$chr, gwas$pos >= range$start, gwas$pos <= range$end, !is.na(gwas$beta))
  gwas.small <- gwas.small[!duplicated(gwas.small$variant_id),]
  
  qtl.small <- qtl_var %>% filter(qtl_var$chr == range$chr, qtl_var$pos >= range$start, qtl_var$pos <= range$end, qtl_var$gene == x)
  qtl.small <- qtl.small[!duplicated(qtl.small$variant_id),]
  
  # confirm snps match
  
  gwas.dataset <- list(beta = gwas.small$beta, 
                      snp = gwas.small$variant_id, 
                      varbeta = (gwas.small$se)^2, 
                      type = gwas.type, 
                      #position = gwas.small$hg38pos,
                      N = n, 
                      MAF = gwas.small$af,
                      pvalues = gwas.small$pvalue)
  
  qtl.dataset <- list(beta = qtl.small$beta,
                      snp = qtl.small$variant_id,
                      varbeta = (qtl.small$se_beta)^2, 
                      type = "quant",
                      #position = qtl.small$pos,
                      N = 53,
                      MAF = as.numeric(qtl.small$maf),
                      pvalues = qtl.small$pvalue)
  
  print(gwas.dataset$snp[duplicated(gwas.dataset$snp)])
  print(qtl.dataset$snp[duplicated(qtl.dataset$snp)])
  
  if (sum(qtl.small$variant_id %in% gwas.small$variant_id) > 0) {
  
  # run coloc and save
    res <- coloc.abf(gwas.dataset, qtl.dataset)
    res$results$gene <- x
    
    resultslist[[x]] <- res$results
    summarylist[[x]] <- res$summary
    
  }
}
results <- as.data.frame(do.call(rbind, resultslist))
results<- results[order(results$SNP.PP.H4, decreasing = T),]

summary <- as.data.frame(do.call(rbind, summarylist))
summary<- summary[order(summary$PP.H4.abf, decreasing = T),]



## save results
write.table(results, file = paste0("/path//20260625_coloc_", gwasname, "_sqtl_", tx, "_results.txt"), sep = "\t", quote = F)
write.table(summary, file = paste0("/path//20250625_coloc_", gwasname, "_sqtl_", tx, "_summary.txt"), sep = "\t", quote = F)

}

}
}
