### calculating deltaPSI per donor for the 53 HAECs
library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

setwd("/path/promoterBindingRatios/")


### RNA expression by donor
rna_nt <- read.delim("/path/notx/chrall_NT_splicing_phenotype.ratio.expression.txt", as.is = T, stringsAsFactors = F, row.names = 1)
names(rna_nt) <- lapply(strsplit(names(rna_nt), "_"), "[[", 1)
rna_il1b <- read.delim("/path/il1b/chrall_IL1B_splicing_phenotype.ratio.expression.txt", as.is = T, stringsAsFactors = F, row.names = 1)
names(rna_il1b) <- lapply(strsplit(names(rna_il1b), "_"), "[[", 1)
# put in the same order
rna_nt <- rna_nt[order(names(rna_nt), decreasing = F)]
rna_il1b <- rna_il1b[order(names(rna_il1b), decreasing = F)]
match(names(rna_nt), names(rna_il1b))

# get list of introns
introns <- row.names(rna_nt)[row.names(rna_nt) %in% row.names(rna_il1b)]

## for every intron calculate the dpsi from the psi
dpsi <- list()
for (intron in introns) {
  psi_nt = t(rna_nt[row.names(rna_nt) == intron,])
  psi_il1b = t(rna_il1b[row.names(rna_il1b) == intron,])
  
  dpsi[[intron]] <- psi_nt - psi_il1b
}

dpsi_df <- data.frame(t(data.frame(dpsi)))
dpsi_df$intron <- as.character(gsub("\\.", ":", gsub("X", "chr", unlist(row.names(dpsi_df))))) ## put in the same format as the leafcutter outputs so things can be matched up

# save
write.table(dpsi_df, file = "/path/promoterBindingRatios/leafcutter_53HAECsIL1Bornotx_deltapsi_byDonor.txt", sep = "\t", quote = F)

#dpsi_df <- read.delim("/path/promoterBindingRatios/leafcutter_53HAECsIL1Bornotx_deltapsi_byDonor.txt", as.is = T, stringsAsFactors = F)



## match up with the leafcutter res
leaf <- read.delim("/path/leafcutter_results_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig.txt", as.is = T, stringsAsFactors = F)
leaf %>% filter(p.adjust < 0.05, abs(deltapsi) > 0.05) -> leaf.sig
leaf.sig %>% filter(spliceType == "AF") %>% group_by(genes) %>% filter(n() == 2) -> leaf.af

leaf.af$intron %in% dpsi_df$intron

dpsi.af <- dpsi_df[dpsi_df$intron %in% leaf.af$intron,]
dpsi.af$leafcutter_dpsi <- leaf.af$deltapsi[match(dpsi.af$intron, leaf.af$intron)]
dpsi.af$average <- rowMeans(dpsi.af[,1:53])

plot(dpsi.af$average, dpsi.af$leafcutter_dpsi) + title("dpsi by donor vs dpsi from leafcutter (w/ covariates)") ## strong positive correlation between the averages by donor and the leafcutter values (this means only slight differences based on regressing covariates)
cor(dpsi.af$average, dpsi.af$leafcutter_dpsi) # 0.93

dpsi.af %>% pivot_longer(cols = names(dpsi.af)[!names(dpsi.af) %in% c("intron", "leafcutter_dpsi", "average")]) -> to.plot.dpsi

ggplot(to.plot.dpsi, aes(name, reorder(intron, value), fill = value)) + geom_tile() + scale_fill_gradient2(limits = c(-0.5, 0.5), high = "red", low = "blue", mid = "white", na.value = "white", name = "deltaPSI") +  theme(axis.text = element_blank()) + ggtitle("Splicing at AF exons by donor") + ylab("Intron") + xlab("HAEC donor")

#### boxplots for deltapsi for genes in the oxidative stress pathway
# read in gene ontology results:
gse_res <- read.delim("/path/leafcutter_DSGs_fgsea_geneSetEnrichment_results.txt", as.is = T, stringsAsFactors = F)
oxstressgenes <- c("NCOA7", "RCAN1", "GPX4", "ABL1", "SESN1", "KDM6B")

# add gene symbols to dpsi_df
dpsi_df$geneSymbol <- leaf$genes[match(dpsi_df$intron, leaf$intron)]
dpsi_df$deltapsi <- leaf$deltapsi[match(dpsi_df$intron, leaf$intron)]
dpsi_df$p.adjust <- leaf$p.adjust[match(dpsi_df$intron, leaf$intron)]

dpsi_df %>% 
  filter(geneSymbol %in% oxstressgenes & abs(deltapsi) > 0.05 & p.adjust < 0.05) %>%
  mutate(category = ifelse(deltapsi > 0, "Inducible", "Basal"))-> dpsi_oxstress

oxstress_toplot <- melt(dpsi_oxstress, id.vars = c("geneSymbol", "deltapsi", "p.adjust", "intron", "category"))
head(oxstress_toplot)

# create stat df for plotting
leaf <- as_tibble(leaf)
leaf %>% filter(genes %in% oxstressgenes & deltapsi > 0.05 & p.adjust < 0.05) %>% dplyr::select(genes, p.adjust) %>% mutate(group1 = "Basal", group2 = "Inducible", y.position = 0.95, geneSymbol = genes, category = "Inducible") -> oxstress_pvals


ggplot(oxstress_toplot, aes(category, value, fill = category)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = c("steelblue", "darkred"), name = "Alternate First Exon") + facet_grid(cols = vars(geneSymbol), scales = "free") + ylab("deltaPSI") + xlab("") + ggtitle("Oxidative Stress DSGs") + 
  stat_pvalue_manual(oxstress_pvals, label = "p.adjust", size = 3) + ylim(-1,1)


leaf %>% filter(genes %in% toxicgenes & deltapsi > 0.05 & p.adjust < 0.05) %>% dplyr::select(genes, p.adjust) %>% mutate(group1 = "Basal", group2 = "Inducible", y.position = 0.95, geneSymbol = genes, category = "Inducible") -> toxic_pvals


### for genes in the response to orgacidbiosynth substance pathway
gse_res <- read.delim("leafcutter_DSGs_fgsea_geneSetEnrichment_results.txt", as.is = T, stringsAsFactors = F)
orgacidbiosynthgenes <- as.character(vapply(strsplit(gse_res$core_enrichment[6], split = "/"), function(x) paste(x), character(7L))  )

dpsi_df %>% 
  filter(geneSymbol %in% orgacidbiosynthgenes & abs(deltapsi) > 0.05 & p.adjust < 0.05) %>%
  mutate(category = ifelse(deltapsi > 0, "Inducible", "Basal"))-> dpsi_orgacidbiosynth

orgacidbiosynth_toplot <- melt(dpsi_orgacidbiosynth, id.vars = c("geneSymbol", "deltapsi", "p.adjust", "intron", "category"))
head(orgacidbiosynth_toplot)


leaf %>% filter(genes %in% orgacidbiosynthgenes & abs(deltapsi) > 0.05 & p.adjust < 0.05) %>%  mutate(group1 = intron, group2 = intron, y.position = 0.95, geneSymbol = genes, category = "Inducible") %>% dplyr::select(genes, p.adjust, group1, group2, y.position, geneSymbol, category) -> orgacidbiosynth_pvals


ggplot(orgacidbiosynth_toplot, aes(intron, value, fill = category)) + geom_boxplot() + theme_classic() + scale_fill_manual(values = c("steelblue", "darkred"), name = "Alternative Exon") + facet_grid(cols = vars(geneSymbol), scales = "free") + ylab("deltaPSI") + xlab("") + ggtitle("Organic Acid Biosynthesis DSGs") + geom_text(data = orgacidbiosynth_pvals, aes(label = paste0("FDR=", p.adjust), x = 1.5, y = 0.9)) + 
  theme(axis.text.x = element_text(angle = 90))


#### plot all results separated by pathways
dpsi_sig <- dpsi_df[dpsi_df$intron %in% leaf$intron[abs(leaf$deltapsi) > 0.05 & leaf$p.adjust < 0.05],]
head(dpsi_sig)
dpsi_sig$gene <- leaf$genes[match(dpsi_sig$intron, leaf$intron)]
dpsi_sig$direction <- ifelse(leaf$deltapsi[match(dpsi_sig$intron, leaf$intron)] > 0, "Up IL1B", 'Down IL1B')
dpsi_sig$cluster <- unlist(lapply(strsplit(dpsi_sig$intron, ":"), "[[",4))

pathways <- tibble(
  term = gse_res$Description,
  NES = gse_res$NES,
  p.adjust = gse_res$p.adjust,
  genes = strsplit(gse_res$core_enrichment, "/")
)

dpsi_sig$term <- NULL
for (x in 1:nrow(dpsi_sig)) {
dpsi_sig$term[x] <- paste(pathways$term[grep(dpsi_sig$gene[x], pathways$genes)], collapse = ",")
dpsi_sig$topterm[x] <- pathways$term[grep(dpsi_sig$gene[x], pathways$genes)][1]
}
table(dpsi_sig$firstterm)

dpsi_melt <- melt(dpsi_sig[!is.na(dpsi_sig$term),], id.vars = c("gene", "term", "cluster", "direction"))
dpsi_melt$value <- as.numeric(dpsi_melt$value)

#toppaths <- c("response to toxic substance", "response to oxidative stress", "positive regulation of transcription by RNA polymerase II", "organic acid biosynthetic process", "positive regulation of programmed cell death", "positive regulation of cytokine production", "phosphoric ester hydrolase activity", "RNA splicing, via transesterification reactions") these are redundant
toppaths <- c("response to oxidative stress", "positive regulation of transcription by RNA polymerase II", "organic acid metabolic process", "cellular response to chemical stimulus", "positive regulation of programmed cell death", "RNA splicing", "organophosphate metabolic process", "cytokine production", "lipid biosynthetic process", "cellular response to organonitrogen compound")

ggplot(gse_res[gse_res$Description %in% toppaths & gse_res$p.adjust < 0.05,], aes(x = NES, y = reorder(Description, NES), color = p.adjust, size = setSize)) + geom_point() + scale_color_gradient(low = "red", high = "steelblue") + theme_minimal() + ylab("") + xlab("Normalized Enrichment Score") + ggtitle("Pathway Enrichment in IL1B-DSGs") + theme(axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.75, vjust = 1.5,face = "bold"), axis.text = element_text(size = 13))

### repeat this for every pathway
ggplot(dpsi_melt[grep("cellular response to organonitrogen compound", dpsi_melt$term),] , aes(variable, reorder(gene, value), fill = value)) + geom_tile() + facet_grid(cols = vars(direction), scales = "free", space = "free") + scale_fill_gradient(low = "yellow", high = "blue", limits = c(-0.65,0.65)) + theme_classic() + theme(axis.text.y = element_text(size = 7), axis.text.x = element_blank(), strip.text.y = element_text(angle = 360), axis.ticks.x = element_blank()) + xlab("HAEC donor") + ylab("cellular response to organonitrogen compound") 


# how to make this for everything at once?
top <- pathways[pathways$term %in% toppaths,]

key <- list(
  Response_To_Oxidative_Stress = c("NCOA7", "RCAN1", "GPX4",  "ABL1",  "NET1",  "SESN1", "KDM6B", "HTRA2"),
  positive_regulation_of_transcription_by_RNA_polymerase_II = c( "TCF4", "SBNO2",  "ASH2L",    "AKAP8L", "HAND2",  "RREB1",  "SMAD1",  "ZBTB38", "EPC1",   "TNIP1",  "BPTF"),
  cellular_response_to_chemical_stimulus = c( "INSIG1" ,  "GUK1"  ,   "ABCD4",    "BCAR1" ,    "SP100"),
  positive_regulation_of_programmed_cell_death = c("TNFAIP8",  "TSC22D1",  "CDKN1A" , "CASP7"  ,    "DFFA"  ,   "CTSC"  ,   "TP53INP1", "GRN"  ,    "CTNNA1"  , "DAB2IP" ),
  organic_acid_metabolic_process = c("MTHFD2L", "GPX4",    "DBI",    "MTHFR",   "GSTO1" ,   "BPGM"  ,  "ABCD4" ),
  RNA_splicing  = c("WTAP" ,  "THRAP3" ,"UBL5" ,  "DDX46", "CLK4"  , "PTBP1" , "PRPF39"),
  organophosphate_metabolic_process  = c( "PFKFB3",   "CNP"   ,  "GUK1"    ,"BPGM",    "PNPLA8" , "GLYCTK"),
  cytokine_production = c( "NFATC1" , "TSC22D1", "TANK"  ,  "CD34"  )  ,
  lipid_biosynthetic_process   = c(  "PLD1"  ,  "PIK3C2B", "ABHD2"  ,  "LIPG"  ,  "PNPLA8") ,
  cellular_response_to_organonitrogen_compound  = c("TBC1D4", "CASP7",    "KANK1" , "BAIAP2")
)



dpsi_sig$toppath <- ifelse(dpsi_sig$gene %in% key$Response_To_Oxidative_Stress, "response to oxidative stress",
                           ifelse(dpsi_sig$gene %in% key$positive_regulation_of_transcription_by_RNA_polymerase_II, "positive regulation of transcription by RNA polymerase II",
                                  ifelse(dpsi_sig$gene %in% key$cellular_response_to_chemical_stimulus, "cellular response to chemical stimulus",
                                         ifelse(dpsi_sig$gene %in% key$positive_regulation_of_programmed_cell_death, "positive regulation of programmed cell death",
                                                ifelse(dpsi_sig$gene %in% key$organic_acid_metabolic_process, "organic acid metabolic process",
                                                       ifelse(dpsi_sig$gene %in% key$RNA_splicing, "RNA splicing",
                                                              ifelse(dpsi_sig$gene %in% key$organophosphate_metabolic_process, "organophosphate metabolic process", 
                                                                     ifelse(dpsi_sig$gene %in% key$cytokine_production, "cytokine production",
                                                                            ifelse(dpsi_sig$gene %in% key$lipid_biosynthetic_process, "lipid biosynthetic process",
                                                                                   ifelse(dpsi_sig$gene %in% key$cellular_response_to_organonitrogen_compound, "cellular response to organonitrogen compound", "not in a pathway"))))))))))
table(dpsi_sig$toppath)

top_melt <- melt(dpsi_sig[!dpsi_sig$toppath == "not in a pathway",], id.vars = c("gene", "toppath", "direction"))
top_melt$deltaPSI <- ifelse(as.numeric(top_melt$value) > 0.5, 0.5,
                            ifelse(as.numeric(top_melt$value) < -0.5, -0.5, as.numeric(top_melt$value)))

top_melt$path_fdr <- format(signif(gse_res$p.adjust[match(top_melt$toppath, gse_res$Description)], digits = 3),scientific = T)
top_melt$path_nes <- round(gse_res$NES[match(top_melt$toppath, gse_res$Description)], digits = 2)

top_melt$label <- paste0(top_melt$toppath, ", ", "\nNES=", top_melt$path_nes, ", p.adjust=", top_melt$path_fdr)

#top_melt$toppath <- factor(top_melt$toppath, levels = c("Response to oxidative stress","Positive regulation of transcription by RNA polymerase II","Cellular response to chemical stimulus","Positive regulation of programmed cell death","Organic acid metabolic process","RNA splicing","Organophosphate metabolic process", "Cytokine production","Lipid biosynthetic process","Cellular response to organonitrogen compound"))

pdf("/path/figures/DSGs_in_GO_terms_heatmap_byDPSI.pdf")
ggplot(na.omit(top_melt), aes(variable, reorder(gene, abs(deltaPSI)), fill = deltaPSI)) + geom_tile() + facet_grid(cols = vars(direction), rows = vars(reorder(label, -1*path_nes)), scales = "free", space = "free", switch = "y", labeller = label_wrap_gen(width = 40)) + scale_fill_gradient(low = "yellow", high = "blue") + theme_classic() + theme(axis.text.y = element_text(size = 7), axis.text.x = element_blank(), strip.text.y.left = element_text(angle = 360), axis.ticks.x = element_blank(), axis.title = element_text(size= 7)) + xlab("n = 53             n = 53") + ylab("") + ggtitle("Pathway enrichment in DSGs") + scale_y_discrete(position = "right")
dev.off()

