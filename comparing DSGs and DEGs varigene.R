###### comparing DEG and DSGs with Il1B or notx

setwd("/Volumes/data3/anna/splicing_analysis_2024/")

DSGs <- read.delim("/Volumes/data3/anna/splicing_analysis_2024/20240321_annotatedFINAL_leafcutter_results_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig.txt", as.is = T, stringsAsFactors = F)
DEGs <- read.delim("/Volumes/data3/anna/splicing_analysis_2024/24_04_17_variGeneFinal_edgeR_differentialExpression-Treatment_condenseGenes.txt",
                   as.is = T, stringsAsFactors = T, row.names = 1) #run with condenseGenes!

# positive logFC = UP in Il1B
DEGs$direction <- ifelse(DEGs$logCPM > 0, paste("IL1B"), paste("notx"))

# negative logFC = favored site with IL1B
DSGs$direction <- ifelse(DSGs$deltapsi < 0, paste("IL1B"), paste("notx"))
DSGs$gene.sym <- lapply(strsplit(DSGs$genes, split = ","), "[[", 1)


# plot -logFDR of splicing vs -log FDR of expression
DSGs$gex_logFC <- DEGs$logFC[match(DSGs$gene.sym, DEGs$Gene.sym)]
DSGs$gex_FDR <- DEGs$FDR[match(DSGs$gene.sym, DEGs$Gene.sym)]

##DEGs$deltapsi <- DSGs$deltapsi[match(DEGs$Gene.sym, DSGs$gene.sym)]
##DEGs$DSG_FDR <- DSGs$p.adjust[match(DEGs$Gene.sym, DSGs$gene.sym)]
##DEGs$spliceType <- DSGs$spliceType[match(DEGs$Gene.sym, DSGs$gene.sym)]

# remove "NA" and change to 1 for DEG FDR if there is no significant DEG
sum(is.na(DSGs$gex_FDR))
DSGs$gex_FDR <- ifelse(!is.na(DSGs$gex_FDR), DSGs$gex_FDR, 1)
DSGs$gex_FDR <- as.numeric(DSGs$gex_FDR)

# remove duplicates for plotting
to.plot <- DSGs[, c("gex_FDR", "p.adjust", "gene.sym")]
to.plot<- to.plot[!duplicated(to.plot),]

library(ggplot2)
library(ggrepel)
palette <- c("#76E1B2", "#E56E63", "#BB46E0", "#DCC7D8", "#D0DEAF", "#DDAB69", "#83A1DA", "#96DADD", "#9272D5", "#85E359", "#DE77B6", "#7E7D73", "#D8E068")


ggplot(to.plot, aes(x = -log10(gex_FDR), y = -log10(p.adjust), label = unlist(gene.sym))) +
  geom_jitter() + 
  xlab("DEG -log10(FDR)") +
  ylab("DSG -log10(FDR)") +
  ggtitle("DEGs and DSGs between notx and IL1b") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    panel.background = element_rect(fill='transparent'),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent'),
    axis.ticks = element_line (color = "black"),
    axis.line = element_line (color = "black"),
    text = element_text(size = 12)) +
  #geom_label_repel(data = to.plot[to.plot$p.adjust < 0.05 & to.plot$gex_FDR < 0.05,], force = 10, max.overlaps = 15)+
  scale_color_manual(values = palette)


## stats
sum(DSGs$p.adjust < 0.05 & abs(DSGs$deltapsi) > 0.05) #1224 DSGs
sum(DSGs$gex_FDR < 0.05) #3520 DEGs
sum(DSGs$p.adjust < 0.05 & abs(DSGs$deltapsi) > 0.05 & DSGs$gex_FDR < 0.05) #872 DEG and DSG
sum(DSGs$p.adjust < 0.05 & abs(DSGs$deltapsi) > 0.05 & DSGs$gex_FDR < 0.05)/sum(DSGs$p.adjust < 0.05 & abs(DSGs$deltapsi) > 0.05) #71%

sum(DSGs$p.adjust < 0.05 & abs(DSGs$deltapsi) > 0.05 & DSGs$spliceType == "AF")/sum(DSGs$p.adjust < 0.05 & abs(DSGs$deltapsi) > 0.05) #30% AF
sum(DSGs$p.adjust < 0.05 & abs(DSGs$deltapsi) > 0.05 & DSGs$spliceType == "SE")/sum(DSGs$p.adjust < 0.05 & abs(DSGs$deltapsi) > 0.05) # 30% SE

sum(DSGs$p.adjust < 0.05 & abs(DSGs$deltapsi) > 0.05 & DSGs$spliceType == "AF" & DSGs$gex_FDR < 0.05)/sum(DSGs$p.adjust < 0.05 & abs(DSGs$deltapsi) > 0.05& DSGs$gex_FDR < 0.05) #33% of DSG/DEGs are AF


#### barplot by splice type
DSGs$spliceType <- factor(DSGs$spliceType, levels =c("A3", "A5", "AF", "AL", "MX", "RI", "SE", "cryptic_threeprime", "cryptic_fiveprime", "cryptic_unanchored", "novel annotated pair", "unknown_strand"))


DSGs$DEG <- ifelse(DSGs$gex_FDR < 0.05, "DEG and DSG", "DSG only")

ggplot(DSGs[abs(DSGs$deltapsi) > 0.05,], aes(x=spliceType, fill = DEG)) +
  geom_bar(stat="count", position = "stack") +
  #geom_text(stat = ..count../sum(..count..), position = "stack", vjust = -0.5) +
  labs(y = "Number of siginificant events") +
  theme(
    panel.background = element_rect(fill='transparent'),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent'),
    axis.ticks = element_line (color = "black"),
    axis.line = element_line (color = "black"))+
  scale_fill_manual(values = palette)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("DSG and DEG comparison")


## volcano plot of DEGs colored by DSG status
DSGs.sig <- DSGs[abs(DSGs$deltapsi) > 0.05 & DSGs$p.adjust < 0.05,]

DEGs$spliceType <- DSGs.sig$spliceType[match(unlist(DEGs$Gene.sym), DSGs.sig$gene.sym)]
table(DEGs$spliceType)
DEGs$AF <- ifelse(DEGs$spliceType == "AF", "AF", NA)
table(DEGs$AF)
DEGs$AF <- factor(DEGs$AF, levels = c(NA,"AF"), ordered = T)

DEGs$SE <- ifelse(DEGs$spliceType == "SE", "SE", NA)
table(DEGs$SE)
DEGs$SE <- factor(DEGs$SE, levels = c(NA,"SE"), ordered = T)

pdf("/Volumes/data3/anna/splicing_ms_2025/figures/AF_SE_DSGs_in_DEGs_volcanoplots.pdf", height = 4, width = 5)
ggplot(DEGs) +
  geom_point(aes(x = logFC, y = -log10(FDR)), color = "black", size = 0.5) +
  geom_point(data = DEGs[DEGs$AF == "AF",], aes(x = logFC, y = -log10(FDR), color = "orchid"), size = 1)+
  #geom_label_repel(data = DEGs[DEGs$AF == "AF",], aes(x = logFC, y = -log10(FDR), label = Gene.sym), color = "red", size = 2, force = 15) +
  ggtitle("DEGs with IL1B treatment", subtitle = "AF-DSGs in purple") +
  theme_classic() + scale_color_manual(values = c("orchid", "black")) +
 xlab("Gene Expression Ratio \nlog2(IL1B vs notx)") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgray") + geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "darkgray") + 
 annotate(label = "Up IL1B", x = 3, y = 50, geom="text") + 
  annotate(label = "Down IL1B", x = -2.2, y = 50, geom="text") 




ggplot(DEGs) +
  geom_point(aes(x = logFC, y = -log10(FDR)), color = "black", size = 0.5) +
  geom_point(data = DEGs[DEGs$SE == "SE",], aes(x = logFC, y = -log10(FDR), color = "steelblue"), size = 1)+
  #geom_label_repel(data = DEGs[DEGs$AF == "AF",], aes(x = logFC, y = -log10(FDR), label = Gene.sym), color = "red", size = 2, force = 15) +
  ggtitle("Differentially Expressed Genes with IL1B treatment", subtitle = "SE-DSGs in blue") +
  theme_classic() + scale_color_manual(values = c("steelblue", "black")) +
  xlab("log2(FoldChange)")
dev.off()

SE.genes <- DEGs$Gene.sym[!is.na(DEGs$SE) & DEGs$FDR < 0.05]

### testing for splice type enrichment in all DSTs:
# are DEGs enriched for splice types?
table(DEGs$spliceType)
fisher_list <- list()
for (type in unique(DEGs$spliceType)) {
  mat <- matrix(c(
    sum(DEGs$FDR < 0.05),
    sum(DEGs$FDR < 0.05 & DEGs$spliceType %in% type),
    sum(!DEGs$FDR < 0.05),
    sum(!DEGs$FDR < 0.05 & DEGs$spliceType %in% type)),
    nrow = 2
  )
  fisher_list[[type]] <- fisher.test(mat)$p.value
}

fisher_df <- data.frame(t(as.data.frame(fisher_list)))
fisher_df$spliceType <- unlist(row.names(fisher_df))
fisher_df$spliceType <- ifelse(fisher_df$spliceType %in% c("cryptic_unanchored", "cryptic_threeprime", "cryptic_fiveprime", "novel.annotated.pair", "NA."), "cryptic", fisher_df$spliceType)
names(fisher_df) <- c("p.value", "spliceType")

fisher_df$spliceType <- factor(fisher_df$spliceType, levels = c("AF", "A3", "A5", "SE", "MX", "RI", "AL", "cryptic"))

pdf("/Volumes/data3/anna/splicing_ms_2025/figures/spliceType_enrich_in_DEGs_dotplot.pdf", height = 5, width = 5)
ggplot(fisher_df[!fisher_df$spliceType == "cryptic",], aes(-log10(p.value), fct_rev(spliceType), fill = p.value < 0.05)) + geom_dotplot(dotsize = 3.5) + theme_bw() + theme(axis.line = element_line(color = "black"), axis.text.y = element_text(vjust = -0.3, size = 12), axis.text.x = element_text(size = 12)) + scale_fill_manual(values = c("gray", "red"), name = "Fisher's Exact", labels = c("p.value > 0.05", "p.value < 0.05")) + ylab("") + ggtitle("Splice Type Enrichment in DEGs") + xlab("-log10(p.value)")
dev.off()
