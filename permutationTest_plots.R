########### PLOTTING ENRICHMENT TESTS ###########

library(dplyr)
library(tidyverse)
library(ggplot2)
set.seed(12678)

setwd("/Volumes/data3/leafcutter/hg38_splicing_QTLs/eQTL_molQTL_enrichmentTests/")

sqtl_nt_sumstat <- read.table("/Volumes/data3/leafcutter/hg38_splicing_QTLs/notx/splicing_QTLs_notx_rna_cis1000bpCis_with_gene_level_FDR.txt", sep = "\t", as.is = T, stringsAsFactors = F, header = T)
head(sqtl_nt_sumstat)
sqtl_il1b_sumstat <- read.table("/Volumes/data3/leafcutter/hg38_splicing_QTLs/il1b/splicing_QTLs_il1b_rna_cis1000bpCis_with_gene_level_FDR.txt", sep = "\t", as.is = T, stringsAsFactors = F, header = T)
head(sqtl_nt_sumstat)

################ RELA ##################
## read in results and plotted on local to get them perfectly shaped
results_p65 <- read.delim("/Volumes/data3/leafcutter/hg38_splicing_QTLs/eQTL_molQTL_enrichmentTests/20251002_sQTL_nt_p65_enrichvsrandom1000perm.txt")

results_p65 %>%
  mutate(xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = z_score)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "histQTL Enrichment by in sQTLs",
       x = "cumulative p-values from sQTLs",
       y = "Enrichment Z-score") + theme(text = element_text(size = 13)) + 
  theme_bw()

### or as a line
results_p65 %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = gsub("observed_n", "RELA molQTLs (IL1B)", category),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + geom_line(size = 1) + geom_point(size = 3) + 
  theme_bw() + scale_color_manual(values = c("black", "seagreen"), name = NULL) + 
  theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 13)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping molQTLs and sQTLs (log"[10]*")"))

## hypergeometric test for each point
results_p65 %>% group_by(bin) %>%
  mutate(pvalbin =  10^-(bin),
         nonhit = sum(sqtl_il1b_sumstat$gene_level_FDR <= pvalbin & sqtl_il1b_sumstat$gene_level_FDR > pvalbin/10),
         hypergeom = dhyper(observed_n,  observed_n, nonhit, observed_n)) -> stats_p65
write.table(stats_p65, "20251002_sQTL_nt_p65_enrichvsrandom1000perm_hypergeometricTest.txt")

################ ERG notx ##################
## read in results and plotted on local to get them perfectly shaped
results_ergnt <- read.delim("/Volumes/data3/leafcutter/hg38_splicing_QTLs/eQTL_molQTL_enrichmentTests/20251002_sQTL_nt_erg.notx_enrichvsrandom1000perm.txt")

results_ergnt %>%
  mutate(xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = z_score)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "histQTL Enrichment by in sQTLs",
       x = "cumulative p-values from sQTLs",
       y = "Enrichment Z-score") + theme(text = element_text(size = 13)) + 
  theme_bw()

### or as a line
results_ergnt %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = gsub("observed_n", "ERG molQTLs (NT)", category),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + geom_line(size = 1) + geom_point(size = 3) + 
  theme_bw() + scale_color_manual(values = c("black", "seagreen"), name = NULL) + 
  theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 13)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping molQTLs and sQTLs (log"[10]*")"))

## hypergeometric test for each point
results_ergnt %>% group_by(bin) %>%
  mutate(pvalbin =  10^-(bin),
         nonhit = sum(sqtl_il1b_sumstat$gene_level_FDR <= pvalbin & sqtl_il1b_sumstat$gene_level_FDR > pvalbin/10),
         hypergeom = dhyper(observed_n,  observed_n, nonhit, observed_n)) -> stats_ergnt
write.table(stats_ergnt, "20251002_sQTL_nt_erg.notx_enrichvsrandom1000perm_hypergeometricTest.txt")
################ ERG notx ##################
## read in results and plotted on local to get them perfectly shaped
results_ergil1b <- read.delim("/Volumes/data3/leafcutter/hg38_splicing_QTLs/eQTL_molQTL_enrichmentTests/20251002_sQTL_il1b_erg.il1b_enrichvsrandom1000perm.txt")

results_ergil1b %>%
  mutate(xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = z_score)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "histQTL Enrichment by in sQTLs",
       x = "cumulative p-values from sQTLs",
       y = "Enrichment Z-score") + theme(text = element_text(size = 13)) + 
  theme_bw()

### or as a line
results_ergil1b %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = gsub("observed_n", "ERG molQTLs (NT)", category),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + geom_line(size = 1) + geom_point(size = 3) + 
  theme_bw() + scale_color_manual(values = c("black", "seagreen"), name = NULL) + 
  theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 13)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping molQTLs and sQTLs (log"[10]*")"))

## hypergeometric test for each point
results_ergil1b %>% group_by(bin) %>%
  mutate(pvalbin =  10^-(bin),
         nonhit = sum(sqtl_il1b_sumstat$gene_level_FDR <= pvalbin & sqtl_il1b_sumstat$gene_level_FDR > pvalbin/10),
         hypergeom = dhyper(observed_n,  observed_n, nonhit, observed_n)) -> stats_ergil1b
write.table(stats_ergil1b, "20251002_sQTL_nt_erg.il1b_enrichvsrandom1000perm_hypergeometricTest.txt")
################ ATAC notx ##################
## read in results and plotted on local to get them perfectly shaped
results_atacnt <- read.delim("/Volumes/data3/leafcutter/hg38_splicing_QTLs/eQTL_molQTL_enrichmentTests/20251002_sQTL_nt_atac.notx_enrichvsrandom1000perm.txt")

results_atacnt %>%
  mutate(xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = z_score)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "histQTL Enrichment by in sQTLs",
       x = "cumulative p-values from sQTLs",
       y = "Enrichment Z-score") + theme(text = element_text(size = 13)) + 
  theme_bw()

### or as a line
results_atacnt %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = gsub("observed_n", "ERG molQTLs (NT)", category),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + geom_line(size = 1) + geom_point(size = 3) + 
  theme_bw() + scale_color_manual(values = c("black", "seagreen"), name = NULL) + 
  theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 13)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping molQTLs and sQTLs (log"[10]*")"))

## hypergeometric test for each point
results_atacnt %>% group_by(bin) %>%
  mutate(pvalbin =  10^-(bin),
         nonhit = sum(sqtl_il1b_sumstat$gene_level_FDR <= pvalbin & sqtl_il1b_sumstat$gene_level_FDR > pvalbin/10),
         hypergeom = dhyper(observed_n,  observed_n, nonhit, observed_n)) -> stats_atacnt
write.table(stats_atacnt, "20251002_sQTL_nt_atac.notx_enrichvsrandom1000perm_hypergeometricTest.txt")

################ ATAC IL1B ##################
## read in results and plotted on local to get them perfectly shaped
results_atacil1b <- read.delim("/Volumes/data3/leafcutter/hg38_splicing_QTLs/eQTL_molQTL_enrichmentTests/20251002_sQTL_il1b_atac.il1b_enrichvsrandom1000perm.txt")

results_atacil1b %>%
  mutate(xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = z_score)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "histQTL Enrichment by in sQTLs",
       x = "cumulative p-values from sQTLs",
       y = "Enrichment Z-score") + theme(text = element_text(size = 13)) + 
  theme_bw()

### or as a line
results_atacil1b %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = gsub("observed_n", "ERG molQTLs (NT)", category),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + geom_line(size = 1) + geom_point(size = 3) + 
  theme_bw() + scale_color_manual(values = c("black", "seagreen"), name = NULL) + 
  theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 13)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping molQTLs and sQTLs (log"[10]*")"))

## hypergeometric test for each point
results_atacil1b %>% group_by(bin) %>%
  mutate(pvalbin =  10^-(bin),
         nonhit = sum(sqtl_il1b_sumstat$gene_level_FDR <= pvalbin & sqtl_il1b_sumstat$gene_level_FDR > pvalbin/10),
         hypergeom = dhyper(observed_n,  observed_n, nonhit, observed_n)) -> stats_atacil1b
write.table(stats_atacil1b, "20251002_sQTL_nt_atac.il1b_enrichvsrandom1000perm_hypergeometricTest.txt")


################ h3k27ac notx ##################
## read in results and plotted on local to get them perfectly shaped
results_h3k27acnt <- read.delim("/Volumes/data3/leafcutter/hg38_splicing_QTLs/eQTL_molQTL_enrichmentTests/20251002_sQTL_nt_h3k27ac.notx_enrichvsrandom1000perm.txt")

results_h3k27acnt %>%
  mutate(xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = z_score)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "histQTL Enrichment by in sQTLs",
       x = "cumulative p-values from sQTLs",
       y = "Enrichment Z-score") + theme(text = element_text(size = 13)) + 
  theme_bw()

### or as a line
results_h3k27acnt %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = gsub("observed_n", "H3K27ac molQTLs (NT)", category),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + geom_line(size = 1) + geom_point(size = 3) + 
  theme_bw() + scale_color_manual(values = c("black", "seagreen"), name = NULL) + 
  theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 13)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping molQTLs and sQTLs (log"[10]*")"))

## hypergeometric test for each point
results_h3k27acnt %>% group_by(bin) %>%
  mutate(pvalbin =  10^-(bin),
         nonhit = sum(sqtl_il1b_sumstat$gene_level_FDR <= pvalbin & sqtl_il1b_sumstat$gene_level_FDR > pvalbin/10),
         hypergeom = dhyper(observed_n,  observed_n, nonhit, observed_n)) -> stats_h3k27acnt
write.table(stats_h3k27acnt, "20251002_sQTL_nt_h3k27ac.notx_enrichvsrandom1000perm_hypergeometricTest.txt")

################ h3k27ac IL1B ##################
## read in results and plotted on local to get them perfectly shaped
results_h3k27acil1b <- read.delim("/Volumes/data3/leafcutter/hg38_splicing_QTLs/eQTL_molQTL_enrichmentTests/20251002_sQTL_il1b_h3k27ac.il1b_enrichvsrandom1000perm.txt")

results_h3k27acil1b %>%
  mutate(xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = z_score)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "histQTL Enrichment by in sQTLs",
       x = "cumulative p-values from sQTLs",
       y = "Enrichment Z-score") + theme(text = element_text(size = 13)) + 
  theme_bw()

### or as a line
results_h3k27acil1b %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = gsub("observed_n", "ERG molQTLs (NT)", category),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + geom_line(size = 1) + geom_point(size = 3) + 
  theme_bw() + scale_color_manual(values = c("black", "seagreen"), name = NULL) + 
  theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 13)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping molQTLs and sQTLs (log"[10]*")"))

## hypergeometric test for each point
results_h3k27acil1b %>% group_by(bin) %>%
  mutate(pvalbin =  10^-(bin),
         nonhit = sum(sqtl_il1b_sumstat$gene_level_FDR <= pvalbin & sqtl_il1b_sumstat$gene_level_FDR > pvalbin/10),
         hypergeom = dhyper(observed_n,  observed_n, nonhit, observed_n)) -> stats_h3k27acil1b
write.table(stats_h3k27acil1b, "20251002_sQTL_nt_h3k27ac.il1b_enrichvsrandom1000perm_hypergeometricTest.txt")


################ rna notx ##################
## read in results and plotted on local to get them perfectly shaped
results_rnant <- read.delim("/Volumes/data3/leafcutter/hg38_splicing_QTLs/eQTL_molQTL_enrichmentTests/20251002_sQTL_nt_rna.notx.eqtl_enrichvsrandom1000perm.txt")

results_rnant %>%
  mutate(xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = z_score)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "histQTL Enrichment by in sQTLs",
       x = "cumulative p-values from sQTLs",
       y = "Enrichment Z-score") + theme(text = element_text(size = 13)) + 
  theme_bw()

### or as a line
results_rnant %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = gsub("observed_n", "eQTLs molQTLs (NT)", category),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + geom_line(size = 1) + geom_point(size = 3) + 
  theme_bw() + 
  #scale_color_manual(values = c("black", "seagreen"), name = NULL) + 
  theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 13)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping molQTLs and sQTLs (log"[10]*")"))

## hypergeometric test for each point
results_rnant %>% group_by(bin) %>%
  mutate(pvalbin =  10^-(bin),
         nonhit = sum(sqtl_il1b_sumstat$gene_level_FDR <= pvalbin & sqtl_il1b_sumstat$gene_level_FDR > pvalbin/10),
         hypergeom = dhyper(observed_n,  observed_n, nonhit, observed_n)) -> stats_rnant
write.table(stats_rnant, "20251002_sQTL_nt_rna.notx_enrichvsrandom1000perm_hypergeometricTest.txt")

################ rna IL1B ##################
## read in results and plotted on local to get them perfectly shaped
results_rnail1b <- read.delim("/Volumes/data3/leafcutter/hg38_splicing_QTLs/eQTL_molQTL_enrichmentTests/20251002_sQTL_il1b_rna.il1b.eqtl_enrichvsrandom1000perm.txt")

results_rnail1b %>%
  mutate(xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = z_score)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "histQTL Enrichment by in sQTLs",
       x = "cumulative p-values from sQTLs",
       y = "Enrichment Z-score") + theme(text = element_text(size = 13)) + 
  theme_bw()

### or as a line
results_rnail1b %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = gsub("observed_n", "ERG molQTLs (NT)", category),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + geom_line(size = 1) + geom_point(size = 3) + 
  theme_bw() + scale_color_manual(values = c("black", "seagreen"), name = NULL) + 
  theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 13)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping molQTLs and sQTLs (log"[10]*")"))

## hypergeometric test for each point
results_rnail1b %>% group_by(bin) %>%
  mutate(pvalbin =  10^-(bin),
         nonhit = sum(sqtl_il1b_sumstat$gene_level_FDR <= pvalbin & sqtl_il1b_sumstat$gene_level_FDR > pvalbin/10),
         hypergeom = dhyper(observed_n,  observed_n, nonhit, observed_n)) -> stats_rnail1b
write.table(stats_rnail1b, "20251002_sQTL_nt_rna.il1b_enrichvsrandom1000perm_hypergeometricTest.txt")




###### plot the significant ones (eQTL and histQTL) on one plot together nicely
library(gdata)
library(gridExtra)
library(ggsignif)

results_rnail1b %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = gsub("observed_n", "eQTLs (IL1B)", category),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + geom_line(size = 1) + geom_point(size = 3) + 
  theme_bw() + scale_color_manual(values = c("seagreen", "black"), name = NULL) + 
  #theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 11)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping QTLs (log"[10]*")")) -> p_rnail1b

results_rnant %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = gsub("observed_n", "eQTLs (Control)", category),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + geom_line(size = 1) + geom_point(size = 3) + 
  theme_bw() + scale_color_manual(values = c("seagreen", "black"), name = NULL) + 
  #theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 11)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping QTLs (log"[10]*")")) -> p_rnant

grid.arrange(p_rnant, p_rnail1b)

### plot eQTLs on one plot
gdata::combine(results_rnant, results_rnail1b) -> rna

rna %>%
  mutate(source = gsub("results_rnant", "(Control)", source), source = gsub("results_rnail1b", "(IL1B)", source)) %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = paste(gsub("observed_n", "eQTLs", category), source, sep = " "),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + 
  geom_line(size = 1, aes(linetype = source)) + geom_point(size = 2, position = position_dodge(width = 0.3)) + 
  theme_bw() +
  scale_color_manual(values = c("darkblue", "skyblue", "darkgray", "lightgray"), name = NULL) + 
  #theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 11)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping eQTLs (log"[10]*")"))-> rnaPlots


gdata::combine(results_h3k27acnt, results_h3k27acil1b) -> h3k27ac

h3k27ac %>%
  mutate(source = gsub("results_h3k27acnt", "(Control)", source), source = gsub("results_h3k27acil1b", "(IL1B)", source)) %>%
  pivot_longer(cols = c(observed_n, mean), names_to = "category", values_to = "value") %>%
  mutate(category = gsub("mean", "random (1000 permutations)", category),
         category = paste(gsub("observed_n", "histQTLs", category), source, sep = " "),
         xlab = paste0("1e-", as.character(bin)),
         xlab = gsub("1e-8", "<1e-8", xlab)) %>%
  ggplot(aes(x = factor(xlab, levels = c("<1e-8","1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1e-0")), y = log10(value), color = category, group = category)) + 
  geom_line(size = 1, aes(linetype = source)) + geom_point(size = 2, position = position_dodge(width = 0.3)) + 
  theme_bw() +
  scale_color_manual(values = c("darkorchid", "violet", "darkgray", "lightgray"), name = NULL) + 
  #theme(legend.position = c(0.33, 0.9), legend.box.background = element_rect(color = "black", fill = NA, linewidth = 1), text = element_text(size = 11)) +
  xlab("cumulative p-values from sQTLs") + ylab(expression("Counts of overlapping histQTLs (log"[10]*")"))-> h3k27acPlots

grid.arrange(rnaPlots, h3k27acPlots, ncol = 2)
