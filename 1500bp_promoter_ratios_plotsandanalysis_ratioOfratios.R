setwd("/Volumes/data3/anna/splicing_analysis_2024/promoterBindingRatios/")
set.seed(12568)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(dplyr)

######################## H3K27ac #############
## nt
ntprom <- read.delim("20241203_1500bpProm_AFnotx_byDonor_H3K27ac.txt", header = T, row.names = 1)
head(ntprom)[,1:10]
names(ntprom)[19:ncol(ntprom)] <-vapply(strsplit(names(ntprom)[19:ncol(ntprom)], split = "\\."), function(x) paste(x[c(10,16)], collapse = "_"), character(1L))
ntprom[["gene"]] <- unlist(lapply(strsplit(row.names(ntprom), split = ":"), "[[", 5))

## il1b
il1bprom <- read.delim("20241203_1500bpProm_AFil1b_byDonor_H3K27ac.txt", header = T, row.names = 1)
head(il1bprom)[,1:10]
names(il1bprom)[19:ncol(il1bprom)] <-vapply(strsplit(names(il1bprom)[19:ncol(il1bprom)], split = "\\."), function(x) paste(x[c(10,16)], collapse = "_"), character(1L))
il1bprom[["gene"]] <- unlist(lapply(strsplit(row.names(il1bprom), split = ":"), "[[", 5))


## compare the ratio at the promoters in each treatment condition
both <- merge(ntprom[, 19:ncol(ntprom)], il1bprom[, 19:ncol(il1bprom)], by = "gene", suffixes = c("_NTprom", "_IL1Bprom"))

# separate into the two treatment conditions
il1b_both <- both[grep("il1b_", names(both))]
il1b_both$gene <- both$gene
nt_both <- both[grep("notx_", names(both))]
nt_both$gene <- both$gene

# get list of donors
donors_i <- unlist(lapply(strsplit(names(il1b_both[1:ncol(il1b_both)-1]), "_"), "[[", 2))
donors_n <- unlist(lapply(strsplit(names(nt_both[1:ncol(nt_both)-1]), "_"), "[[", 2))
donors <- donors_i[donors_i %in% donors_n]

## calculate ratio in il1b treatment
il1b_ratio <-list()
for (x in donors) {
  il1b = il1b_both[grep(paste0(x, "_IL1Bprom"), names(il1b_both))]
  nt = il1b_both[grep(paste0(x, "_NTprom"), names(il1b_both))]
  # adjust 0 values
  il1b[,1] <- ifelse(il1b[,1] == 0, NA, il1b[,1])
  nt[,1] <- ifelse(nt[,1] == 0, NA, nt[,1])
  ratio = il1b/nt ## high ratio means MORE binding at the IL1B promoter with IL1B tx
  il1b_ratio[[x]] <- ratio[,1]
}

il1b_ratio <- as.data.frame(il1b_ratio)

il1b_ratio$avg <- rowMeans(il1b_ratio)

il1b_ratio$gene <- il1b_both$gene

ggplot(data = il1b_ratio, aes(x = gene, y = log10(avg))) + 
  geom_point() + geom_hline(yintercept = 0, linetype = "dashed") + ggtitle("IL1B promoter ratios", subtitle = "+ value = more H3K27ac binding at IL1B promoters in IL1B condition")

notx_ratio <-list()
for (x in donors) {
  il1b = nt_both[grep(paste0(x, "_IL1Bprom"), names(nt_both))]
  nt = nt_both[grep(paste0(x, "_NTprom"), names(nt_both))]
  # adjust 0 values
  il1b[,1] <- ifelse(il1b[,1] == 0, NA, il1b[,1])
  nt[,1] <- ifelse(nt[,1] == 0, NA, nt[,1])
  ratio = il1b/nt ## high ratio means MORE binding at the IL1B promoter with IL1B tx
  notx_ratio[[x]] <- ratio[,1]
}


notx_ratio <- as.data.frame(notx_ratio)

notx_ratio$avg <- rowMeans(notx_ratio)

notx_ratio$gene <- nt_both$gene

ggplot(data = notx_ratio, aes(x = gene, y = log10(avg))) + 
  geom_point() + geom_hline(yintercept = 0, linetype = "dashed") + ggtitle("NT ratios", subtitle = "+ value = more H3K27ac binding at IL1B promoters in NT condition")

### now merge and compare the two 
pairs_h3 <- merge(notx_ratio, il1b_ratio, by = "gene", suffixes = c("_nt", "_il1b"))
names(pairs_h3)

# look at ratio of ratios
ratio_il1bnt_H3 <- list()
for (x in donors) {
  il1b = pairs_h3[grep(paste0(x, "_il1b"), names(pairs_h3))]
  nt = pairs_h3[grep(paste0(x, "_nt"), names(pairs_h3))]
  d = il1b/nt
  ratio_il1bnt_H3[[x]] <- d[,1]
}

ratio_il1bnt_H3 <- as.data.frame(ratio_il1bnt_H3)
sum(is.na(ratio_il1bnt_H3)) ## some NA values


ratio_il1bnt_H3$gene <- pairs_h3$gene
ratio_il1bnt_H3 <- ratio_il1bnt_H3[!duplicated(ratio_il1bnt_H3$gene),]
row.names(ratio_il1bnt_H3) <- ratio_il1bnt_H3$gene

ratio_il1bnt_H3$category <- ifelse(rowMeans(ratio_il1bnt_H3[,1:37], na.rm = T) < 1, "Decreased H3K27ac", "Increased H3K27ac")
ratio_il1bnt_H3 <- ratio_il1bnt_H3[order(rowMeans(ratio_il1bnt_H3[,1:37], na.rm = T), decreasing = F),]

ggplot(ratio_il1bnt_H3, aes(x = reorder(gene, rowMeans(ratio_il1bnt_H3[,1:37], na.rm = T), decreasing = T) , y = rowMeans(ratio_il1bnt_H3[,1:37], na.rm = T))) + geom_col() + facet_grid(~category)

## student's T-Test to compare ratios in Il1B to ratios in NT
pairs_h3$gene <- make.unique(pairs_h3$gene)
genes <- pairs_h3$gene

t_H3K27ac <- list()
for (x in genes) {
  il1b = pairs_h3[pairs_h3$gene == x,2:38]
  nt = pairs_h3[pairs_h3$gene == x, 40:76]
  t_H3K27ac[[x]] <- t.test(as.numeric(il1b[1,1:37]), as.numeric(nt[1,1:37]), paired = T)$p.value
}


pairs_h3$pvalue <- as.numeric(t(unlist(t_H3K27ac)))

df_H3K27ac <- melt(pairs_h3[,!names(pairs_h3) %in% c("avg_nt", "avg_il1b", "pvalue")])
df_H3K27ac$treatment <- unlist(lapply(strsplit(as.character(df_H3K27ac$variable), split = "_"), "[[",2))

df_H3K27ac$pvalue <- pairs_h3$pvalue[match(df_H3K27ac$gene, pairs_h3$gene)]

geneorder <- factor(ratio_il1bnt_H3$gene)
#df_H3K27ac <- df_H3K27ac[order(df_H3K27ac$gene, geneorder),]


df_H3K27ac$pvalue <- ifelse(is.na(df_H3K27ac$pvalue), 1, df_H3K27ac$pvalue)
df_H3K27ac$category<- ratio_il1bnt_H3$category[match(df_H3K27ac$gene, ratio_il1bnt_H3$gene)]
df_H3K27ac$category <- ifelse(df_H3K27ac$pvalue > 0.05, "not significant", df_H3K27ac$category)

df_H3K27ac <- df_H3K27ac[order(df_H3K27ac$category),]

##### final plot!!
ggplot(df_H3K27ac, aes(x = reorder(gene, value), y = value, fill = treatment)) + 
  geom_boxplot(outlier.size = 0.1) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 7), plot.title = element_text(hjust = 0.5, face = "bold", size = 15, vjust = 1), plot.subtitle = element_text(hjust = 0.5))  + 
  ggtitle("H3K27ac at AF promoters", subtitle = "positive value indicates more accessibility at the Il1b-favored promoter") + 
  facet_grid(category ~ ., scales = "free", space = "free") + geom_hline(yintercept = 0, linetype = "dashed") + 
  ylim(-3.5,3.5) + ylab("Ratio of H3K27ac signal (il1b/notx)") + xlab("") 


### heatmap

# heatmap of ratio of ratios (I think this is simpler to look at and more difficult to understand)
ror_H3K27ac <- melt(ratio_il1bnt_H3[,!names(ratio_il1bnt_H3) %in% c("category")])
ror_H3K27ac$pvalue <- pairs_h3$pvalue[match(ror_H3K27ac$gene, pairs_h3$gene)]
#ror_H3K27ac$value <- ifelse(is.na(ror_H3K27ac$value), 0.0000000001, ror_H3K27ac$value)
ggplot(na.omit(ror_H3K27ac), aes(x = reorder(gene, -1* value), y = reorder(variable, value), fill = log2(value))) + 
  geom_tile() + scale_fill_gradient2(low = "blue", mid="white", high="red") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 6, hjust = 0.9), axis.text.y = element_text(size = 6), plot.title = element_text(hjust = 0.5, face = "bold", size = 15), plot.subtitle = element_text(hjust = 0.5)) + 
  ggtitle("H3K27ac at AF promoters", subtitle = "positive value indicates more H3K27ac at the inducible promoter with IL1B treatment") + xlab("AF promoters") + ylab("ratio of H3K27ac by donor")

######################### ATAC ##########################
## nt
ntprom <- read.delim("20241203_1500bpProm_AFnotx_byDonor_ATAC.txt", header = T, row.names = 1)
head(ntprom)[,1:10]
names <- gsub("\\_.*","", gsub("ATAC_", "", vapply(strsplit(names(ntprom)[19:ncol(ntprom)], split = "\\."), function(x) paste(x[c(13,14)], collapse = "_"), character(1L))))
tx <- vapply(strsplit(names(ntprom)[19:ncol(ntprom)], split = "\\."), function(x) paste(x[c(9)], collapse = "_"), character(1L))
names(ntprom)[19:ncol(ntprom)] <- paste(names, tx, sep = "_")

ntprom[["gene"]] <- unlist(lapply(strsplit(row.names(ntprom), split = ":"), "[[", 5))

# il1b
il1bprom <- read.delim("20241203_1500bpProm_AFil1b_byDonor_ATAC.txt", header = T, row.names = 1)
head(il1bprom)[,1:10]
names <- gsub("\\_.*","", gsub("ATAC_", "", vapply(strsplit(names(il1bprom)[19:ncol(il1bprom)], split = "\\."), function(x) paste(x[c(13,14)], collapse = "_"), character(1L))))
tx <- vapply(strsplit(names(il1bprom)[19:ncol(il1bprom)], split = "\\."), function(x) paste(x[c(9)], collapse = "_"), character(1L))
names(il1bprom)[19:ncol(il1bprom)] <- paste(names, tx, sep = "_")

il1bprom[["gene"]] <- unlist(lapply(strsplit(row.names(il1bprom), split = ":"), "[[", 5))

## compare the ratio at the promoters in each treatment condition
both <- merge(ntprom[, 19:ncol(ntprom)], il1bprom[, 19:ncol(il1bprom)], by = "gene", suffixes = c("_NTprom", "_IL1Bprom"))

# separate into the two treatment conditions
il1b_both <- both[grep("_il1b_", names(both))]
il1b_both$gene <- both$gene
nt_both <- both[grep("_notx_", names(both))]
nt_both$gene <- both$gene

# get list of donors
donors_i <- unlist(lapply(strsplit(names(il1b_both[1:ncol(il1b_both)-1]), "_"), "[[", 1))
donors_n <- unlist(lapply(strsplit(names(nt_both[1:ncol(nt_both)-1]), "_"), "[[", 1))
donors <- donors_i[donors_i %in% donors_n]

## calculate ratio in il1b treatment
il1b_ratio <-list()
for (x in donors) {
  il1b = il1b_both[grep(paste0(x, "_il1b_IL1Bprom"), names(il1b_both))]
  nt = il1b_both[grep(paste0(x, "_il1b_NTprom"), names(il1b_both))]
  # adjust 0 values
  il1b[,1] <- ifelse(il1b[,1] == 0, NA, il1b[,1])
  nt[,1] <- ifelse(nt[,1] == 0, NA, nt[,1])
  ratio = il1b/nt ## high ratio means MORE binding at the IL1B promoter with IL1B tx
  il1b_ratio[[x]] <- ratio[,1]
}

il1b_ratio <- as.data.frame(il1b_ratio)

il1b_ratio$avg <- rowMeans(il1b_ratio)

il1b_ratio$gene <- il1b_both$gene

ggplot(data = il1b_ratio, aes(x = gene, y = log10(avg))) + 
  geom_point() + geom_hline(yintercept = 0, linetype = "dashed") + ggtitle("IL1B promoter ratios", subtitle = "+ value = more ATAC binding at IL1B promoters in IL1B condition")

notx_ratio <-list()
for (x in donors) {
  il1b = nt_both[grep(paste0(x, "_notx_IL1Bprom"), names(nt_both))]
  nt = nt_both[grep(paste0(x, "_notx_NTprom"), names(nt_both))]
  # adjust 0 values
  il1b[,1] <- ifelse(il1b[,1] == 0, NA, il1b[,1])
  nt[,1] <- ifelse(nt[,1] == 0, NA, nt[,1])
  ratio = il1b/nt ## high ratio means MORE binding at the IL1B promoter with IL1B tx
  notx_ratio[[x]] <- ratio[,1]
}


notx_ratio <- as.data.frame(notx_ratio)

notx_ratio$avg <- rowMeans(notx_ratio)

notx_ratio$gene <- nt_both$gene

ggplot(data = notx_ratio, aes(x = gene, y = log10(avg))) + 
  geom_point() + geom_hline(yintercept = 0, linetype = "dashed") + ggtitle("NT ratios", subtitle = "+ value = more ATAC binding at IL1B promoters in NT condition")

### now merge and compare the two 
pairs_ATAC <- merge(notx_ratio, il1b_ratio, by = "gene", suffixes = c("_nt", "_il1b"))
names(pairs_ATAC)

# look at ratio of ratios
ratio_il1bnt_ATAC <- list()
for (x in donors) {
  il1b = pairs_ATAC[grep(paste0(x, "_il1b"), names(pairs_ATAC))]
  nt = pairs_ATAC[grep(paste0(x, "_nt"), names(pairs_ATAC))]
  d = il1b/nt
  ratio_il1bnt_ATAC[[x]] <- d[,1]
}

ratio_il1bnt_ATAC <- as.data.frame(ratio_il1bnt_ATAC)
sum(is.na(ratio_il1bnt_ATAC)) ## some NA values


ratio_il1bnt_ATAC$gene <- pairs_ATAC$gene
ratio_il1bnt_ATAC <- ratio_il1bnt_ATAC[!duplicated(ratio_il1bnt_ATAC$gene),]
row.names(ratio_il1bnt_ATAC) <- ratio_il1bnt_ATAC$gene

ratio_il1bnt_ATAC$category <- ifelse(rowMeans(ratio_il1bnt_ATAC[,1:37], na.rm = T) < 1, "Decreased ATAC", "Increased ATAC")
ratio_il1bnt_ATAC <- ratio_il1bnt_ATAC[order(rowMeans(ratio_il1bnt_ATAC[,1:37], na.rm = T), decreasing = F),]

ggplot(ratio_il1bnt_ATAC, aes(x = reorder(gene, rowMeans(ratio_il1bnt_ATAC[,1:37], na.rm = T), decreasing = T) , y = rowMeans(ratio_il1bnt_ATAC[,1:37], na.rm = T))) + geom_col() + facet_grid(~category)

## student's T-Test to compare ratios in Il1B to ratios in NT
pairs_ATAC$gene <- make.unique(pairs_ATAC$gene)
genes <- pairs_ATAC$gene

t_ATAC <- list()
for (x in genes) {
  il1b = pairs_ATAC[pairs_ATAC$gene == x,2:38]
  nt = pairs_ATAC[pairs_ATAC$gene == x, 40:76]
  t_ATAC[[x]] <- t.test(as.numeric(il1b[1,1:37]), as.numeric(nt[1,1:37]), paired = T)$p.value
}


pairs_ATAC$pvalue <- as.numeric(t(unlist(t_ATAC)))

df_ATAC <- melt(pairs_ATAC[,!names(pairs_ATAC) %in% c("avg_nt", "avg_il1b", "pvalue")])
df_ATAC$treatment <- unlist(lapply(strsplit(as.character(df_ATAC$variable), split = "_"), "[[",2))

df_ATAC$pvalue <- pairs_ATAC$pvalue[match(df_ATAC$gene, pairs_ATAC$gene)]

geneorder <- factor(ratio_il1bnt_ATAC$gene)
#df_ATAC <- df_ATAC[order(df_ATAC$gene, geneorder),]


df_ATAC$pvalue <- ifelse(is.na(df_ATAC$pvalue), 1, df_ATAC$pvalue)
df_ATAC$category<- ratio_il1bnt_ATAC$category[match(df_ATAC$gene, ratio_il1bnt_ATAC$gene)]
df_ATAC$category <- ifelse(df_ATAC$pvalue > 0.05, "not significant", df_ATAC$category)

df_ATAC <- df_ATAC[order(df_ATAC$category),]

##### final plot!!
ggplot(df_ATAC, aes(x = reorder(gene, value), y = value, fill = treatment)) + 
  geom_boxplot(outlier.size = 0.1) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 7), plot.title = element_text(hjust = 0.5, face = "bold", size = 15, vjust = 1), plot.subtitle = element_text(hjust = 0.5))  + 
  ggtitle("ATAC at AF promoters", subtitle = "positive value indicates more accessibility at the Il1b-favored promoter") + 
  facet_grid(category ~ ., scales = "free", space = "free") + geom_hline(yintercept = 1, linetype = "dashed") + 
  ylim(-3.5,3.5) + ylab("Ratio of ATAC signal (il1b/notx)") + xlab("") 


### heatmap

# heatmap of ratio of ratios (I think this is simpler to look at and more difficult to understand)
ror_ATAC <- melt(ratio_il1bnt_ATAC[,!names(ratio_il1bnt_ATAC) %in% c("category")])
#ror_ATAC$value <- ifelse(is.na(ror_ATAC$value), 0.0000000001, ror_ATAC$value)
ggplot(na.omit(ror_ATAC), aes(x = reorder(gene, -1* value), y = reorder(variable, value), fill = log2(value))) + 
  geom_tile() + scale_fill_gradient2(low = "blue", mid="white", high="red") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 6, hjust = 0.9), axis.text.y = element_text(size = 6), plot.title = element_text(hjust = 0.5, face = "bold", size = 15), plot.subtitle = element_text(hjust = 0.5)) + 
  ggtitle("ATAC at AF promoters", subtitle = "positive value indicates more ATAC at the inducible promoter with IL1B treatment") + xlab("AF promoters") + ylab("ratio of ATAC by donor")


######################### ERG ##########################
## nt
ntprom <- read.delim("20241203_1500bpProm_AFnotx_byDonor_ERG.txt", header = T, row.names = 1)
head(ntprom)[,1:10]
names <-vapply(strsplit(names(ntprom)[19:ncol(ntprom)], split = "\\."), function(x) paste(x[c(16, 10)], collapse = "_"), character(1L))
names(ntprom)[19:ncol(ntprom)] <- names

ntprom[["gene"]] <- unlist(lapply(strsplit(row.names(ntprom), split = ":"), "[[", 5))

# il1b
il1bprom <- read.delim("20241203_1500bpProm_AFil1b_byDonor_ERG.txt", header = T, row.names = 1)
head(il1bprom)[,1:10]
names <-vapply(strsplit(names(il1bprom)[19:ncol(il1bprom)], split = "\\."), function(x) paste(x[c(16, 10)], collapse = "_"), character(1L))
names(il1bprom)[19:ncol(il1bprom)] <- names

il1bprom[["gene"]] <- unlist(lapply(strsplit(row.names(il1bprom), split = ":"), "[[", 5))

## compare the ratio at the promoters in each treatment condition
both <- merge(ntprom[, 19:ncol(ntprom)], il1bprom[, 19:ncol(il1bprom)], by = "gene", suffixes = c("_NTprom", "_IL1Bprom"))

# separate into the two treatment conditions
il1b_both <- both[grep("_il1b_", names(both))]
il1b_both$gene <- both$gene
nt_both <- both[grep("_notx_", names(both))]
nt_both$gene <- both$gene

# get list of donors
donors_i <- unlist(lapply(strsplit(names(il1b_both[1:ncol(il1b_both)-1]), "_"), "[[", 1))
donors_n <- unlist(lapply(strsplit(names(nt_both[1:ncol(nt_both)-1]), "_"), "[[", 1))
donors <- donors_i[donors_i %in% donors_n]

## calculate ratio in il1b treatment
il1b_ratio <-list()
for (x in donors) {
  il1b = il1b_both[grep(paste0(x, "_il1b_IL1Bprom"), names(il1b_both))]
  nt = il1b_both[grep(paste0(x, "_il1b_NTprom"), names(il1b_both))]
  # adjust 0 values
  il1b[,1] <- ifelse(il1b[,1] == 0, NA, il1b[,1])
  nt[,1] <- ifelse(nt[,1] == 0, NA, nt[,1])
  ratio = il1b/nt ## high ratio means MORE binding at the IL1B promoter with IL1B tx
  il1b_ratio[[x]] <- ratio[,1]
}

il1b_ratio <- as.data.frame(il1b_ratio)

il1b_ratio$avg <- rowMeans(il1b_ratio, na.rm = T)

il1b_ratio$gene <- il1b_both$gene

ggplot(data = il1b_ratio, aes(x = gene, y = log10(avg))) + 
  geom_point() + geom_hline(yintercept = 0, linetype = "dashed") + ggtitle("IL1B promoter ratios", subtitle = "+ value = more ERG binding at IL1B promoters in IL1B condition")

notx_ratio <-list()
for (x in donors) {
  il1b = nt_both[grep(paste0(x, "_notx_IL1Bprom"), names(nt_both))]
  nt = nt_both[grep(paste0(x, "_notx_NTprom"), names(nt_both))]
  # adjust 0 values
  il1b[,1] <- ifelse(il1b[,1] == 0, NA, il1b[,1])
  nt[,1] <- ifelse(nt[,1] == 0, NA, nt[,1])
  ratio = il1b/nt ## high ratio means MORE binding at the IL1B promoter with IL1B tx
  notx_ratio[[x]] <- ratio[,1]
}
notx_ratio <- as.data.frame(notx_ratio)

notx_ratio$avg <- rowMeans(notx_ratio, na.rm = T)

notx_ratio$gene <- nt_both$gene

### get p values for ERG just in NT (why, Il1B data quality is poor)
genes <- nt_both$gene
t_ERG <- list()
for (x in genes) {
  il1b = nt_both[grep("_notx_IL1Bprom", names(nt_both))][ nt_both$gene == x,]
  nt = nt_both[grep("_notx_NTprom", names(nt_both))][ nt_both$gene == x,]
  t_ERG[[x]] <- t.test(as.numeric(il1b[1,]), as.numeric(nt[1,]), paired = F)$p.value
}


notx_ratio$pvalue <- as.numeric(t(unlist(t_ERG)))[match(notx_ratio$gene, names(t_ERG))]

write.table(notx_ratio, "/Volumes/data3/anna/splicing_ms_2025/final_data/ERG_ratio_alt_promoters_byDonor_notxOnly.txt", sep = "\t", quote = F)



ggplot(data = notx_ratio, aes(x = gene, y = log10(avg))) + 
  geom_point() + geom_hline(yintercept = 0, linetype = "dashed") + ggtitle("NT ratios", subtitle = "+ value = more ERG binding at IL1B promoters in NT condition")

### now merge and compare the two 
pairs_ERG <- merge(notx_ratio, il1b_ratio, by = "gene", suffixes = c("_nt", "_il1b"))
names(pairs_ERG)

# look at ratio of ratios
ratio_il1bnt_ERG <- list()
for (x in donors) {
  il1b = pairs_ERG[grep(paste0(x, "_il1b"), names(pairs_ERG))]
  nt = pairs_ERG[grep(paste0(x, "_nt"), names(pairs_ERG))]
  d = il1b/nt
  ratio_il1bnt_ERG[[x]] <- d[,1]
}

ratio_il1bnt_ERG <- as.data.frame(ratio_il1bnt_ERG)
sum(is.na(ratio_il1bnt_ERG)) ## some NA values


ratio_il1bnt_ERG$gene <- pairs_ERG$gene
ratio_il1bnt_ERG <- ratio_il1bnt_ERG[!duplicated(ratio_il1bnt_ERG$gene),]
row.names(ratio_il1bnt_ERG) <- ratio_il1bnt_ERG$gene

ratio_il1bnt_ERG$category <- ifelse(rowMeans(ratio_il1bnt_ERG[,1:15], na.rm = T) < 1, "Decreased ERG", "Increased ERG")
ratio_il1bnt_ERG <- ratio_il1bnt_ERG[order(rowMeans(ratio_il1bnt_ERG[,1:15], na.rm = T), decreasing = F),]

ggplot(ratio_il1bnt_ERG, aes(x = reorder(gene, rowMeans(ratio_il1bnt_ERG[,1:15], na.rm = T), decreasing = T) , y = rowMeans(ratio_il1bnt_ERG[,1:15], na.rm = T))) + geom_col() + facet_grid(~category)

## student's T-Test to compare ratios in Il1B to ratios in NT
pairs_ERG$gene <- make.unique(pairs_ERG$gene)
genes <- pairs_ERG$gene

t_ERG <- list()
for (x in genes) {
  il1b = pairs_ERG[pairs_ERG$gene == x,2:16]
  nt = pairs_ERG[pairs_ERG$gene == x, 18:33]
  t_ERG[[x]] <- t.test(as.numeric(il1b[1,1:15]), as.numeric(nt[1,1:15]), paired = T)$p.value
}


pairs_ERG$pvalue <- as.numeric(t(unlist(t_ERG)))

df_ERG <- melt(pairs_ERG[,!names(pairs_ERG) %in% c("avg_nt", "avg_il1b", "pvalue")])
df_ERG$treatment <- unlist(lapply(strsplit(as.character(df_ERG$variable), split = "_"), "[[",2))

df_ERG$pvalue <- pairs_ERG$pvalue[match(df_ERG$gene, pairs_ERG$gene)]

geneorder <- factor(ratio_il1bnt_ERG$gene)
#df_ERG <- df_ERG[order(df_ERG$gene, geneorder),]


df_ERG$pvalue <- ifelse(is.na(df_ERG$pvalue), 1, df_ERG$pvalue)
df_ERG$category<- ratio_il1bnt_ERG$category[match(df_ERG$gene, ratio_il1bnt_ERG$gene)]
df_ERG$category <- ifelse(df_ERG$pvalue > 0.05, "not significant", df_ERG$category)

df_ERG <- df_ERG[order(df_ERG$category),]
order <- pairs_ERG$gene[order(log2(pairs_ERG$avg_nt), decreasing = T)]
df_ERG$order <- match(df_ERG$gene, order)

##### final plot!!
ggplot(df_ERG, aes(x = reorder(gene, order), y = log2(value), color = value > 1)) + 
  geom_boxplot(outlier.size = 0.1) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 6, hjust = 0.9), plot.title = element_text(hjust = 0.5, face = "bold", size = 15, vjust = 1), plot.subtitle = element_text(hjust = 0.5))  + 
  ggtitle("ERG at AF promoters", subtitle = "positive value indicates more accessibility at the Il1b-favored promoter") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ylim(-3.5,3.5) + ylab("Ratio of ERG signal (il1b/notx)") + xlab("") + facet_grid(cols = vars(treatment)) + coord_flip()


### heatmap

# heatmap of ratio of ratios (I think this is simpler to look at and more difficult to understand)
ror_ERG <- melt(ratio_il1bnt_ERG[,!names(ratio_il1bnt_ERG) %in% c("category")])
#ror_ERG$value <- ifelse(is.na(ror_ERG$value), 0.0000000001, ror_ERG$value)
ggplot(na.omit(ror_ERG), aes(x = reorder(gene, -1*value), y = variable, fill = log2(value))) + 
  geom_tile() + scale_fill_gradient2(low = "blue", mid="white", high="red", name = "Ratio") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 6, hjust = 0.9), axis.text.y = element_text(size = 6), plot.title = element_text(hjust = 0.5, face = "bold", size = 15), plot.subtitle = element_text(hjust = 0.5)) + 
  ggtitle("ERG at AF promoters", subtitle = "positive value indicates more ERG at the inducible promoter with IL1B treatment") + xlab("AF promoters") + ylab("ratio of ERG by donor")




######################### p65 ##########################
## nt
ntprom <- read.delim("20241203_1500bpProm_AFnotx_byDonor_p65.txt", header = T, row.names = 1)
head(ntprom)[,19:35]
names <-vapply(strsplit(names(ntprom)[19:ncol(ntprom)], split = "\\."), function(x) paste(x[c(15, 17)], collapse = "_"), character(1L))
names(ntprom)[19:ncol(ntprom)] <- names

ntprom[["gene"]] <- unlist(lapply(strsplit(row.names(ntprom), split = ":"), "[[", 5))

# il1b
il1bprom <- read.delim("20241203_1500bpProm_AFil1b_byDonor_p65.txt", header = T, row.names = 1)
head(il1bprom)[,1:10]
names <-vapply(strsplit(names(il1bprom)[19:ncol(il1bprom)], split = "\\."), function(x) paste(x[c(15, 17)], collapse = "_"), character(1L))
names(il1bprom)[19:ncol(il1bprom)] <- names

il1bprom[["gene"]] <- unlist(lapply(strsplit(row.names(il1bprom), split = ":"), "[[", 5))

## compare the ratio at the promoters in each treatment condition
both_p65 <- merge(ntprom[, 19:ncol(ntprom)], il1bprom[, 19:ncol(il1bprom)], by = "gene", suffixes = c("_NTprom", "_IL1Bprom"))

# separate into the two treatment conditions
il1b_both<- both_p65[grep("_IL1b_", names(both_p65))]
il1b_both$gene <- both_p65$gene


# get list of donors
donors_i <- unlist(lapply(strsplit(names(il1b_both[1:ncol(il1b_both)-1]), "_"), "[[", 1))
#donors_n <- unlist(lapply(strsplit(names(nt_both[1:ncol(nt_both)-1]), "_"), "[[", 1))
#donors <- donors_i[donors_i %in% donors_n]

## calculate ratio in il1b treatment
il1b_ratio_p65 <-list()
for (x in donors_i) {
  il1b = il1b_both[grep(paste0(x, "_IL1b_IL1Bprom"), names(il1b_both))]
  nt = il1b_both[grep(paste0(x, "_IL1b_NTprom"), names(il1b_both))]
  # adjust 0 values
  il1b[,1] <- ifelse(il1b[,1] == 0, NA, il1b[,1])
  nt[,1] <- ifelse(nt[,1] == 0, NA, nt[,1])
  ratio = il1b/nt ## high ratio means MORE binding at the IL1B promoter with IL1B tx
  il1b_ratio_p65[[x]] <- ratio[,1]
}

il1b_ratio_p65 <- as.data.frame(il1b_ratio_p65)

il1b_ratio_p65$avg <- rowMeans(il1b_ratio_p65, na.rm = T)

il1b_ratio_p65$gene <- il1b_both$gene

ggplot(data = il1b_ratio_p65, aes(x = gene, y = log10(avg))) + 
  geom_point() + geom_hline(yintercept = 0, linetype = "dashed") + ggtitle("IL1B promoter ratios", subtitle = "+ value = more p65 binding at IL1B promoters in IL1B condition")

## student's T-Test to compare ratios in Il1B to ratios in NT
il1b_ratio_p65$gene <- make.unique(il1b_ratio_p65$gene)
genes <- both_p65$gene

t_p65 <- list()
for (x in genes) {
  il1b = il1b_both[grep("_IL1b_IL1Bprom", names(il1b_both))][ il1b_both$gene == x,]
  nt = il1b_both[grep("_IL1b_NTprom", names(il1b_both))][ il1b_both$gene == x,]
  t_p65[[x]] <- t.test(as.numeric(il1b[1,]), as.numeric(nt[1,]), paired = T)$p.value
}


il1b_ratio_p65$pvalue <- as.numeric(t(unlist(t_p65)))[match(il1b_ratio_p65$gene, names(t_p65))]

df_p65 <- melt(il1b_ratio_p65[,!names(il1b_ratio_p65) %in% c("avg_nt", "avg_il1b", "pvalue")])
df_p65$treatment <- "IL1B"

df_p65$pvalue <- il1b_ratio_p65$pvalue[match(df_p65$gene, il1b_ratio_p65$gene)]

df_p65$category <- ifelse(df_p65$pvalue > 0.05, "not significant", "significant")

df_p65 <- df_p65[order(df_p65$category),]

##### final plot!!
ggplot(df_p65, aes(x = reorder(gene, value), y = log2(value), fill = pvalue < 0.05)) + 
  geom_boxplot(outlier.size = 0.1) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 6, hjust = 0.9), plot.title = element_text(hjust = 0.5, face = "bold", size = 15, vjust = 1), plot.subtitle = element_text(hjust = 0.5))  + 
  ggtitle("p65 at AF promoters", subtitle = "positive value indicates more accessibility at the Il1b-favored promoter") + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ylim(-3.5,3.5) + ylab("Ratio of p65 signal (il1b/notx)") + xlab("") 



### heatmap
# heatmap of ratio for p65 is just inducible/basal promoters
ggplot(na.omit(df_p65), aes(x = reorder(gene, -1*value), y = reorder(variable, value), fill = log2(value))) + 
  geom_tile() + scale_fill_gradient2(low = "blue", mid="white", high="green") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 6, hjust = 0.9), axis.text.y = element_text(size = 6), plot.title = element_text(hjust = 0.5, face = "bold", size = 15), plot.subtitle = element_text(hjust = 0.5)) + 
  ggtitle("p65 at AF promoters", subtitle = "positive value indicates more p65 at the inducible promoter with IL1B treatment") + xlab("AF promoters") + ylab("ratio of p65 by donor")




######################## FINAL ANALYSIS #########################

# add ID to each df
df_p65$mark <- "P65" ### this dataframe is ONLY inducible:basal in IL1B, make sure that goes in the legend, the rest are ratio of ratios
ror_H3K27ac$mark <- "H3K27ac"
ror_ATAC$mark <- "ATAC"
ror_ERG$mark <- "ERG"


#### saving results and categorizing promoters
##promoters with p65 signal increased/decreased with Il1B at inducible prom
p65up <- na.omit(il1b_ratio_p65$gene[il1b_ratio_p65$pvalue < 0.05 & il1b_ratio_p65$avg > 1])
p65down <- na.omit(il1b_ratio_p65$gene[il1b_ratio_p65$pvalue < 0.05 & il1b_ratio_p65$avg < 1])

##promoters with H3K27ac signal increased/decreased with Il1B at inducible prom
H3K27acup <- na.omit(pairs_h3$gene[pairs_h3$pvalue < 0.05 & pairs_h3$avg_il1b/pairs_h3$avg_nt > 1])
H3K27acdown <- na.omit(pairs_h3$gene[pairs_h3$pvalue < 0.05 & pairs_h3$avg_il1b/pairs_h3$avg_nt < 1])

##promoters with ATAC signal increased/decreased with Il1B at inducible prom
ATACup <- na.omit(pairs_ATAC$gene[pairs_ATAC$pvalue < 0.05 & pairs_ATAC$avg_il1b/pairs_ATAC$avg_nt > 1])
ATACdown <- na.omit(pairs_ATAC$gene[pairs_ATAC$pvalue < 0.05 & pairs_ATAC$avg_il1b/pairs_ATAC$avg_nt < 1])

##promoters with ERG signal increased/decreased with Il1B at inducible prom
ERGup <- na.omit(notx_ratio$gene[notx_ratio$pvalue < 0.05 & notx_ratio$avg > 1])
ERGdown <- na.omit(notx_ratio$gene[notx_ratio$pvalue < 0.05 & notx_ratio$avg < 1])

# genes with various combinations of promoter marks
ind_prom_all4 <- p65up[p65up %in% H3K27acup & p65up %in% ATACup & p65up %in% ERGup]
p65_only <- p65up[!p65up %in% ERGup]
erg_only <- ERGup[!ERGup %in% p65up]
p65_erg_act <- ERGup[ERGup %in% p65up]
p65_erg_rep <- p65down[p65down %in% ERGdown]
erg_rep <- ERGdown[!ERGdown %in% p65_erg_rep]
p65_rep <- p65down[!p65down %in% p65_erg_rep]

##### save results of tests
write.table(pairs_ATAC, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/ATAC_ratio_alt_promoters_byDonor.txt", row.names = F, quote = F, sep = "\t")
write.table(pairs_h3, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/H3K27ac_ratio_alt_promoters_byDonor.txt", row.names = F, quote = F, sep = "\t")
write.table(pairs_ERG, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/ERG_ratio_alt_promoters_byDonor.txt", row.names = F, quote = F, sep = "\t")
write.table(il1b_ratio_p65, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/p65_ratio_alt_promoters_byDonor.txt", row.names = F, quote = F, sep = "\t")

##### save results of ratios
write.table(ror_ATAC, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/ATAC_ratio_alt_promoters.txt", row.names = F, quote = F, sep = "\t")
write.table(ror_H3K27ac, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/H3K27ac_ratio_alt_promoters.txt", row.names = F, quote = F, sep = "\t")
write.table(ror_ERG, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/ERG_ratio_alt_promoters.txt", row.names = F, quote = F, sep = "\t")
write.table(df_ERG, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/ERG_ratio_alt_promoters_sepByTreatment.txt", row.names = F, quote = F, sep = "\t")
write.table(df_p65, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/p65_ratio_alt_promoters.txt", row.names = F, quote = F, sep = "\t")

#ror_ATAC <- read.delim("/Volumes/data3/anna/splicing_ms_2025/final_data/ATAC_ratio_alt_promoters.txt", as.is = T, stringsAsFactors = F)
#ror_H3K27ac <- read.delim("/Volumes/data3/anna/splicing_ms_2025/final_data/H3K27ac_ratio_alt_promoters.txt", as.is = T, stringsAsFactors = F)
#ror_ERG <- read.delim("/Volumes/data3/anna/splicing_ms_2025/final_data/ERG_ratio_alt_promoters.txt", as.is = T, stringsAsFactors = F)
#df_p65<- read.delim("/Volumes/data3/anna/splicing_ms_2025/final_data/p65_ratio_alt_promoters.txt", as.is = T, stringsAsFactors = F)

pairs_ATAC <- read.delim( file = "/Volumes/data3/anna/splicing_ms_2025/final_data/ATAC_ratio_alt_promoters_byDonor.txt")
pairs_h3 <- read.delim(file = "/Volumes/data3/anna/splicing_ms_2025/final_data/H3K27ac_ratio_alt_promoters_byDonor.txt")
pairs_ERG <- read.delim( file = "/Volumes/data3/anna/splicing_ms_2025/final_data/ERG_ratio_alt_promoters_byDonor.txt")
il1b_ratio_p65 <- read.delim( file = "/Volumes/data3/anna/splicing_ms_2025/final_data/p65_ratio_alt_promoters_byDonor.txt")

###### making lists of these genes to pull from our promoter lists for motif analysis
write.table(p65_erg_act, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/p65_and_erg_activating_alt_prom.txt", row.names = F, quote = F)
write.table(p65_erg_rep, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/p65_and_erg_repressing_alt_prom.txt", row.names = F, quote = F)
write.table(p65_only, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/p65_only_activating_alt_prom.txt", row.names = F, quote = F)
write.table(p65_rep, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/p65_only_repressing_alt_prom.txt", row.names = F, quote = F)
write.table(erg_only, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/erg_only_activating_alt_prom.txt", row.names = F, quote = F)
write.table(erg_rep, file = "/Volumes/data3/anna/splicing_ms_2025/final_data/erg_only_repressing_alt_prom.txt", row.names = F, quote = F)

p65_erg_act <- readLines("/Volumes/data3/anna/splicing_ms_2025/final_data/p65_and_erg_activating_alt_prom.txt")
p65_erg_rep <- readLines("/Volumes/data3/anna/splicing_ms_2025/final_data/p65_and_erg_repressing_alt_prom.txt")
p65_only <- readLines("/Volumes/data3/anna/splicing_ms_2025/final_data/p65_only_activating_alt_prom.txt")
p65_rep <- readLines("/Volumes/data3/anna/splicing_ms_2025/final_data/p65_only_repressing_alt_prom.txt")
erg_only <- readLines("/Volumes/data3/anna/splicing_ms_2025/final_data/erg_only_activating_alt_prom.txt")
erg_rep <- readLines("/Volumes/data3/anna/splicing_ms_2025/final_data/erg_only_repressing_alt_prom.txt")

##### make one DF with all of these
df_ERG$mark = "ERG"
toplot <- rbind(ror_ATAC, df_ERG[df_ERG$treatment == "nt", c("gene", "variable", "value", "mark")], ror_H3K27ac[, c("gene", "variable", "value", "mark")], df_p65[, c("gene", "variable", "value", "mark")])
## ERG is only the "nt" ratio of inducible:basal because the ERG IL1B ChIP isn't as high of quality as the ERG nt ChIP so the data is skewed when we make ratios between the two

toplot$category <- factor(ifelse(toplot$gene %in% p65_erg_act, "p65 and ERG activated",
                          ifelse(toplot$gene %in% p65_only, "p65 activated",
                                ifelse(toplot$gene %in% p65_rep, "p65 repressed", 
                                       ifelse(toplot$gene %in% erg_only, "ERG activated",
                                              ifelse(toplot$gene %in% erg_rep, "ERG repressed",
                                                     ifelse(toplot$gene %in% p65_erg_rep, "p65 and ERG repressed",
                                       "not significant")))))), levels = c("p65 and ERG activated", "p65 activated", "ERG activated","p65 and ERG repressed", "p65 repressed", "ERG repressed", "not significant"))
table(toplot$category)

  
ggplot(toplot[!toplot$category == "not significant",], aes(x = reorder(gene, value), y = variable, fill = log2(value))) + 
  geom_tile() + scale_fill_gradient2(low = "blue", mid="white", high="red", name = "Ratio", limits = c(-5,5)) +
  theme_minimal() + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 7.5), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle = element_text(hjust = 0.5), strip.text = element_text(size = 8, face = "bold"), strip.text.y = element_text(angle = 0), strip.background.x = element_rect(fill = "white"), strip.background.y = element_rect(fill = "white"), panel.border = element_rect(linetype = "solid", fill = NA)) + 
  facet_grid(cols = vars(mark), rows = vars(category), scales = "free", space = "free") + coord_flip() + xlab("") + ylab("")



  ##### add deltaPSI to the heatmap:
dpsi <- read.delim(file = "/Volumes/data3/anna/splicing_analysis_2024/promoterBindingRatios/leafcutter_53HAECsIL1Bornotx_deltapsi_byDonor.txt", as.is=T, stringsAsFactors = F)

head(dpsi)
basalintrons <- dpsi$intron[rowMeans(dpsi[1:53]) < 0]
inducintrons <- dpsi$intron[rowMeans(dpsi[1:53]) > 0]

dpsi %>% pivot_longer(cols = names(dpsi)[!names(dpsi) %in% c("intron", "leafcutter_dpsi", "average")]) -> to.plot.dpsi
head(to.plot.dpsi)
names(to.plot.dpsi) <- c("intron", "donor", "value")
to.plot.dpsi$category <- ifelse(to.plot.dpsi$intron %in% inducintrons, "Inducible", "Basal")
to.plot.dpsi$type <- "RNA"

# get gene symbols
key <- data.frame(
  gene = unlist(lapply(strsplit(row.names(ntprom), split = ":"), "[[", 5)),
  intron = unlist(vapply(strsplit(row.names(ntprom), split = ":"), function(x) paste(x[c(1,2,3,4)], collapse = ":"), character(1L))))

to.plot.dpsi$gene <- key$gene[match(unlist(to.plot.dpsi$intron), key$intron)]
to.plot.dpsi$intron <- key$intron[match(unlist(to.plot.dpsi$intron), key$intron)]



to.plot.dpsi %>% group_by(gene) %>% filter(gene %in% toplot$gene) %>% mutate(variable = donor, treatment = "IL1B", pvalue = NA, intron = NULL, mark = "deltaPSI (inducible exon)", donor = NULL, type = NULL, promoter = category, category = NA, treatment = NULL)->to.plot.dpsi.matched

to.plot.dpsi.matched$category <- toplot$category[match(to.plot.dpsi.matched$gene, toplot$gene)]
toplot$promoter = NA

to.plot <- rbind(toplot, to.plot.dpsi.matched[,!names(to.plot.dpsi.matched) == 'pvalue'])
# only plot with the inducible dPSI

## and change the values such that the plot has better colors and the legend can be changed later
to.plot$value_toplot[!to.plot$mark == "deltaPSI (inducible exon)"] <- ifelse(log2(to.plot$value[!to.plot$mark == "deltaPSI (inducible exon)"]) > 5.5, 5.5, ifelse(log2(to.plot$value[!to.plot$mark == "deltaPSI (inducible exon)"]) < -5.5, -5.5, log2(to.plot$value[!to.plot$mark == "deltaPSI (inducible exon)"])))

to.plot$value_toplot[to.plot$mark == "deltaPSI (inducible exon)"] <- ifelse(to.plot$value[to.plot$mark == "deltaPSI (inducible exon)"] > 0.5, 0.5, ifelse(to.plot$value[to.plot$mark == "deltaPSI (inducible exon)"] < -0.5, -0.5, to.plot$value[to.plot$mark == "deltaPSI (inducible exon)"]))

### order the heatmap by the "real" deltapsi (which includes covariates)
leafres <- read.delim("/Volumes/data3/anna/splicing_analysis_2024/20240321_annotatedFINAL_leafcutter_results_53HAECsIL1Bornotx_final_set_w_depth.txt")
leafres$genesym <- lapply(strsplit(leafres$genes, ","), "[[", 1)
to.plot$deltaPSI <- leafres$deltapsi[match(to.plot$gene, leafres$genesym)]
to.plot$gene[is.na(to.plot$deltaPSI)]

to.plot$mark <- gsub(" exon", "", to.plot$mark)
to.plot$mark <- factor(to.plot$mark, levels = c("deltaPSI (inducible)", "P65", 'ERG', "H3K27ac", "ATAC"))

## for p65 and ERG
ggplot(to.plot[!to.plot$category == "not significant", ], aes(x = reorder(gene, deltaPSI * -1), y = variable, fill = value_toplot)) + 
  geom_tile() + scale_fill_gradient2(low = "blue", mid="white", high="red", name = "log2(inducible/basal)") +
  theme_minimal() + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 7.5), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle = element_text(hjust = 0.5), strip.text = element_text(size = 8, face = "bold"), strip.text.y = element_text(angle = 0), strip.background.x = element_rect(fill = "white"), strip.background.y = element_rect(fill = "white"), panel.border = element_rect(linetype = "solid", fill = NA)) + 
  facet_grid(cols = vars(mark), rows = vars(category), scales = "free", space = "free") + coord_flip() + xlab("") + ylab("")

## for ATAC,H3K27ac
ggplot(to.plot[!to.plot$category == "not significant", ], aes(x = reorder(gene, deltaPSI * -1), y = variable, fill = value_toplot)) + 
  geom_tile() + scale_fill_gradient2(low = "blue", mid="white", high="red", name = "log2(inducible ratio / \nbasal ratio)", limits = c(-3.5,3.5)) +
  theme_minimal() + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 7.5), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle = element_text(hjust = 0.5), strip.text = element_text(size = 8, face = "bold"), strip.text.y = element_text(angle = 0), strip.background.x = element_rect(fill = "white"), strip.background.y = element_rect(fill = "white"), panel.border = element_rect(linetype = "solid", fill = NA)) + 
  facet_grid(cols = vars(mark), rows = vars(category), scales = "free", space = "free") + coord_flip() + xlab("") + ylab("")

## for deltaPSI

ggplot(to.plot[!to.plot$category == "not significant", ], aes(x = reorder(gene, deltaPSI * -1), y = variable, fill = value_toplot)) + 
  geom_tile() + scale_fill_gradient2(low = "blue", mid="white", high="red", name = "deltaPSI \n(IL1B vs notx)", limits = c(-0.5,0.5)) +
  theme_minimal() + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 7.5), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), plot.subtitle = element_text(hjust = 0.5), strip.text = element_text(size = 8, face = "bold"), strip.text.y = element_text(angle = 0), strip.background.x = element_rect(fill = "white"), strip.background.y = element_rect(fill = "white"), panel.border = element_rect(linetype = "solid", fill = NA)) + 
  facet_grid(cols = vars(mark), rows = vars(category), scales = "free", space = "free") + coord_flip() + xlab("") + ylab("")


#### upset plot of all these ratios
library(UpSetR)

upset_df <- list(
  ATAC_activated = ATACup,
  ATAC_repressed = ATACdown,
  H3K27ac_activated = H3K27acup,
  H3K27ac_repressed = H3K27acdown,
  ERG_activated = ERGup,
  ERG_repressed = ERGdown,
  p65_activated = p65up,
  p65_repressed = p65down
)

upset(fromList(upset_df), nsets = 8, mainbar.y.label = "# of promoters", sets.x.label = "# of significant results",
      sets = c("p65_activated", "ERG_activated", "H3K27ac_activated", "ATAC_activated", "p65_repressed", "ERG_repressed", "H3K27ac_repressed", "ATAC_repressed"), keep.order = T, order.by = "freq")

 