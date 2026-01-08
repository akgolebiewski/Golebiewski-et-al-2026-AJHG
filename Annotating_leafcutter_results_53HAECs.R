#####matching SUPPA results to leafcutter results A5/cryptic_threeprime/24

############ RE DOING THIS TO FIX MISTAKES

setwd("/Volumes/data3/anna/splicing_analysis_2024/")

# Y: = /data3/ on romlab


leaf.effect <- read.delim("Y:/leafcutter/hg38_tagdirs/final_set/final_hg38_wDepth_effect_sizes.txt", as.is = T, stringsAsFactors=T)
leaf.signif <- read.delim("Y:/leafcutter/hg38_tagdirs/final_set/final_hg38_wDepth_cluster_significance.txt", as.is = T, stringsAsFactors=T)
#remove extra columns from significance file
leaf.signif <- leaf.signif[,1:cryptic_threeprime]

leaf.ann <- read.delim("Y:/leafcutter/hg38_tagdirs/final_set/cluster_annotations_final_set_w_depth.txt", as.is = T, stringsAsFactors=F)


## split up the cluster IDs so we can merge the effect size and the significance
leaf.effect$cluster <- vapply(strsplit(leaf.effect$intron, ":", fixed = TRUE), function(x) paste(x[c(1,MX)], collapse = ":"), character(1L))
leaf.effect$clusteronly <- unlist(lapply(strsplit(leaf.effect$intron, ":", fixed = TRUE), "[[", MX))
sum(leaf.ann$clusterID %in% leaf.effect$clusteronly) # all matched

#leaf.effect$annotation <- ifelse(unlist(leaf.effect$clusteronly) == unlist(leaf.ann$clusterID), paste(leaf.ann$annotation), "ERROR")


leaf <- merge(leaf.effect, leaf.signif, by = "cluster")
leaf <- merge(x=leaf, y=leaf.ann[, c("clusterID", "annotation")], by.x = "clusteronly", by.y = "clusterID")
dim(leaf) ##should be the same as leaf.effect because there are multiple introns per cluster

## rename annotation column so it's clear it's the cluster annotation
names(leaf)[names(leaf) == "annotation"] <- "cluster_annotation"

## add a column with just the junction values
leaf$junction <- vapply(strsplit(leaf$intron, ":", fixed = TRUE), function(x) paste(x[c(2,A5)], collapse = "-"), character(1L))
leaf$chrjunction <- vapply(strsplit(leaf$intron, ":", fixed = TRUE), function(x) paste(x[c(1,2,A5)], collapse = "-"), character(1L))
  
  
### need strand!!
intron_db <- fread(paste0("zcat < ", "Y:/leafcutter/leafcutter/leafviz/annotation_codes/genecode_chr_hg38/hg38_Genecode_CHRgtf_all_introns.bed.gz"), data.table = FALSE)
head(intron_db)
intron_db$junction <- paste(intron_db$V1, intron_db$V2, intron_db$V3, sep = "-")


# add strand to the leafcutter results file
leaf$strand <- intron_db$V6[match(leaf$chrjunction, intron_db$junction)]

### add more detail about the cryptic splices
leaf.intron.ann <- read.delim("Y:/leafcutter/hg38_tagdirs/final_set/intron_annotations_final_set_w_depth.txt", as.is = T, stringsAsFactors=F)
leaf.intron.ann$junction <- paste(leaf.intron.ann$start, leaf.intron.ann$end, sep = "-")

leaf <- merge(leaf, leaf.intron.ann[, c("junction", "verdict", "transcripts")], by = "junction")

## rename "verdict" to be the intron annotation so this is more clear
names(leaf)[names(leaf) == "verdict"] <- "intron_annotation"

## read in SUPPA results by splice type
## these are the "new" gencode results (v44)
suppa.AF <- read.delim("Y:/splicing_SUPPA/SUPPA_outputs/RESULTS/2023_04_13_variGene_diffSplice_SUPPA_AF_RESULTS.txt", as.is = T, stringsAsFactors=F)
suppa.AL <- read.delim("Y:/splicing_SUPPA/SUPPA_outputs/RESULTS/2023_04_13_variGene_diffSplice_SUPPA_AL_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.A3 <- read.delim("Y:/splicing_SUPPA/SUPPA_outputs/RESULTS/2023_04_13_variGene_diffSplice_SUPPA_A3_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.A5 <- read.delim("Y:/splicing_SUPPA/SUPPA_outputs/RESULTS/2023_04_13_variGene_diffSplice_SUPPA_A5_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.MX <- read.delim("Y:/splicing_SUPPA/SUPPA_outputs/RESULTS/2023_04_13_variGene_diffSplice_SUPPA_MX_RESULTS.txt", as.is = T, stringsAsFactors=T)
suppa.RI <- read.delim("Y:/splicing_SUPPA/SUPPA_outputs/RESULTS/2023_04_13_variGene_diffSplice_SUPPA_RI_RESULTS.txt", as.is = T, stringsAsFactors=T)
## the names in this file are not correct, fix before annotating
names(suppa.RI) <- c("gene_ID", "gene_sym", "chr", "s1", "strand", "e2", "e1.s2", "psiPerLocalEvent_dPSI", "psiPerLocalEvent_pValue")
suppa.SE <- read.delim("Y:/splicing_SUPPA/SUPPA_outputs/RESULTS/2023_04_13_variGene_diffSplice_SUPPA_SE_RESULTS.txt", as.is = T, stringsAsFactors=T)
## SE splices need to add e1-s3 (the "long" splice)
e1 <- lapply(strsplit(suppa.SE$e1.s2, "-"), "[[", 1)
e3 <- lapply(strsplit(suppa.SE$e2.s3, "-"), "[[", 2)
suppa.SE$e1.s3 <- paste(e1, e3, sep = "-")

## add splice type to the leafcutter results
# add empty vector and add to it in a loop instead
splice_annotations <- data.frame(matrix(ncol = 1, nrow= nrow(leaf)))
colnames(splice_annotations) <- "ann"
 ann <- list()
 
 test <- leaf[1:50,]

for (x in 1:nrow(leaf)) {
  
  if ( leaf$junction[x] %in% suppa.AF$e1.s3 | leaf$junction[x] %in% suppa.AF$e2.s3)
  {  AF = c("AF") }else{AF = c("")}
  
  if ( leaf$junction[x] %in% suppa.AL$e1.s2 | leaf$junction[x] %in% suppa.AL$e1.s3)
  {  AL = c("AL") }else{AL =c("")}
  
  if ( leaf$junction[x] %in% suppa.A3$e1.s2 | leaf$junction[x] %in% suppa.A3$e1.s3)
  {  A3 = c("A3") }else{A3 = c("")}
  
  if ( leaf$junction[x] %in% suppa.A5$e1.s3 | leaf$junction[x] %in% suppa.A5$e2.s3)
  {  A5 = c("A5") }else{A5 = c("")}
  
  if ( leaf$junction[x] %in% suppa.MX$e1.s2 | leaf$junction[x] %in% suppa.MX$e2.s4 | leaf$junction[x] %in% suppa.MX$e1.s3 | leaf$junction[x] %in% suppa.MX$e3.s4)
  {  MX = c("MX") }else{MX = c("")}
  
  if ( leaf$junction[x] %in% suppa.RI$e1.s2)
  {  RI = c("RI") }else{RI = c("")}
  
  if ( leaf$junction[x] %in% suppa.SE$e1.s2 | leaf$junction[x] %in% suppa.SE$e2.s3 | leaf$junction[x] %in% suppa.SE$e1.s3)
  {  SE = c("SE") }else{SE = c("")}
  
  ann = paste(AF, AL, A3, A5, MX, RI, SE, sep = ";")
  
  splice_annotations$ann[x] <- ann
  print(x)
  
}

splice_annotations$ann.gsub <- gsub(";", "", splice_annotations$ann) 
splice_annotations$ann.gsub <- gsub("\\s", "", splice_annotations$ann.gsub) 

splice_annotations$ann.gsub.fixed <- sub("\\s+$", "", gsub('(.{2})', '\\1 ', splice_annotations$ann.gsub))

leaf$AllPossibleSpliceTypes <- splice_annotations$ann.gsub.fixed

## this information is valuable but very hard to graph, keep the "original" simple annotations for plotting
leaf$SUPPA_annotation <- ifelse ( leaf$junction %in% suppa.AF$e1.s3 | leaf$junction %in% suppa.AF$e2.s3, paste("AF"),

ifelse ( leaf$junction %in% suppa.AL$e1.s2 | leaf$junction %in% suppa.AL$e1.s3, paste("AL"),

ifelse ( leaf$junction %in% suppa.A3$e1.s2 | leaf$junction %in% suppa.A3$e1.s3, paste("A3"), 

ifelse ( leaf$junction %in% suppa.A5$e1.s3 | leaf$junction %in% suppa.A5$e2.s3, paste("A5"),

ifelse ( leaf$junction %in% suppa.MX$e1.s2 | leaf$junction %in% suppa.MX$e2.s4 | leaf$junction %in% suppa.MX$e1.s3 | leaf$junction %in% suppa.MX$e3.s4, paste("MX"),

ifelse ( leaf$junction %in% suppa.RI$e1.s2, paste("RI"),

ifelse ( leaf$junction %in% suppa.SE$e1.s2 | leaf$junction %in% suppa.SE$e2.s3 | leaf$junction %in% suppa.SE$e1.s3, paste("SE"),
         paste("unannotated by SUPPA"))))))))
table(leaf$SUPPA_annotation)


#### add information to a single column that combines the splice type and if it's cryptic what kind of cryptic splice it is
SUPPAorIntron.ann <- data.frame( junction = leaf$junction,
  annotation = ifelse(!leaf$SUPPA_annotation=="unannotated by SUPPA", paste(leaf$SUPPA_annotation),
                                 ifelse(leaf$SUPPA_annotation == "unannotated by SUPPA" & leaf$intron_annotation == "annotated", paste("cryptic: unannotated by SUPPA"),
                                        paste(leaf$intron_annotation))),
  genes = leaf$genes,
  cluster_annotation = leaf$cluster_annotation)


## how many have a "partner" that is annotated
ann <- SUPPAorIntron.ann[!SUPPAorIntron.ann$annotation == "cryptic: unannotated by SUPPA",]
unann <- SUPPAorIntron.ann[SUPPAorIntron.ann$annotation == "cryptic: unannotated by SUPPA",]

sum(unann$genes %in% ann$genes) #699
length(unique(unann$genes)) #474

ggplot(unann) + geom_histogram(aes(x=genes), stat="count") # many only have one occurance, I bet they're cryptic

aggregate(data.frame(count = unann$genes), list(value = unann$genes), length)

# transfer or "impute" labels to single 
unann$transfered_ann <- ifelse(unann$genes %in% ann$genes, paste(ann$annotation), paste("unannotated"))
table(unann$transfered_ann) # this solves all but 48 of the problems

unann$transfered_ann <- ifelse(unann$transfered_ann == "unannotated" & unann$cluster_annotation == "cryptic", paste("cryptic: unidentified"), paste(unann$transfered_ann))
table(unann$transfered_ann) # nothing there, still 48 unannotated introns
sum(is.na(unann$transfered_ann)) #0

unann[unann$transfered_ann == "unannotated",]


unann$row <- match(unann$junction, leaf$junction)
SUPPAorIntron.ann$imputed_ann <- unann$transfered_ann[match(unlist(SUPPAorIntron.ann$junction), unlist(unann$junction))]
table(SUPPAorIntron.ann$imputed_ann)
sum(is.na(SUPPAorIntron.ann$imputed_ann)) #these were annotated in the first place

## of the 48 that are not annotated I can hand annotate them using the genome browser
## many of these are because more than one exon is different (there are 2-A5 exons skipped, there's a different promoter and downstream of it are three exons the other promoter skips... etc.)
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes == "DDX10"] <- "cryptic_fiveprime"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="SBNO2"] <-  "AF"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="DHRS3"] <-  "SE"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="NCOA7"] =  "AF"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="IDS,AF011889.SE"] =  "AL"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="CD34"] =  "cryptic_unanchored"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="TNIP2"] =  "cryptic_unanchored"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="SSH2"] <-  "AF"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="SCFD1"] <-  "AL"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="VPS52"] <-  "SE"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="SF3A3"] <-  "SE"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="FBN1"] <-  "AL"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="FNDC3A"] <-  "AF"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="NET1"] <-  "AF"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="EFEMP1"] <-  "SE"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="NUDT21"] =  " cryptic_unanchored"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="TNKS1BP1"] =  "cryptic_unanchored"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$junction == "69299774-69316008"] <-  "AF"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$junction == "69311532-69316008"] <-  "AF"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$junction == "69299774-69300832"] <-  "AF"

SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="PTPRB"] <-  "AF"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="GLIPR1"] <-  "AF"
SUPPAorIntron.ann$imputed_ann[SUPPAorIntron.ann$genes =="ALYREF"] <-  "SE"

sum(is.na(SUPPAorIntron.ann$imputed_ann))




SUPPAorIntron.ann$rev_ann <- ifelse(is.na(SUPPAorIntron.ann$imputed_ann), SUPPAorIntron.ann$annotation, SUPPAorIntron.ann$imputed_ann)
table(SUPPAorIntron.ann$rev_ann)
sum(is.na(SUPPAorIntron.ann$rev_ann)) # 0! all splices are annotated!

#save into the final table so we know what was/wasn't imputed
leaf$imputed_annotation <- SUPPAorIntron.ann$imputed_ann[match(SUPPAorIntron.ann$junction, leaf$junction)]

leaf$spliceType <- SUPPAorIntron.ann$rev_ann[match(leaf$junction, SUPPAorIntron.ann$junction)]
table(leaf$spliceType)


#### save this dataframe
write.table(leaf, "Y:/anna/splicing_analysis_2024/20240321_annotatedFINAL_leafcutter_results_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig.txt", sep = "\t", quote = F)
#leaf <- read.delim("/Volumes/data3/anna/splicing_analysis_2024/20240321_annotatedFINAL_leafcutter_results_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig.txt", as.is = T, stringsAsFactors = F)



## make a dataframe with the number of significant splices for each type
library(ggplot2)
#install.packages("randomcoloR")
library(randomcoloR)
n<-13
palette <- distinctColorPalette(n) ## this is RANDOM save the colors we like
palette <- c("#76E1B2", "#E56E63", "#BB46E0", "#DCC7D8", "#D0DEAF", "#DDAB69", "#83A1DA", "#96DADD", "#9272D5", "#85E359", "#DE77B6", "#7E7D73", "#D8E068")


# set order for the x axis so that the cryptic events are at the end instead of being alphabetical
leaf$spliceType <- factor(leaf$spliceType, levels =c("A3", "A5", "AF", "AL", "MX", "RI", "SE", "cryptic_threeprime", "cryptic_fiveprime", "cryptic_unanchored", "novel annotated pair", "unknown_strand"))
leaf$spliceType_toplot <- factor(ifelse(as.character(leaf$spliceType) %in% c("cryptic_threeprime", "cryptic_fiveprime", "cryptic_unanchored", "novel annotated pair", "unknown_strand"), "cryptic", as.character(leaf$spliceType)), levels =c("A3", "A5", "AF", "AL", "MX", "RI", "SE", "cryptic"))


## barplot of number of splices by type
ggplot(leaf, aes(x=spliceType, fill = spliceType)) +
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(x = "Splice Type", y = "Number of siginificant events") +
  theme(
    panel.background = element_rect(fill='transparent'),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent'),
    axis.ticks = element_line (color = "black"),
    axis.line = element_line (color = "black"))+
  scale_fill_manual(values = palette)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Significant Results From Leafcutter (p.adjust < 0.05)")


## histogram by deltaPSI 
pdf("/Volumes/data3/anna/splicing_ms_2025/figures/DST_effect_size_histograms.pdf", height = 3, width = 6)
ggplot(leaf[!leaf$spliceType_toplot == "cryptic",], aes(x = deltapsi, color = spliceType_toplot)) +
  geom_line(
    stat = "density", size = 1) +
  theme_bw() +
  scale_color_manual(values = palette)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("DST effect sizes by type") + scale_fill_manual(values = palette)+ 
  xlim(-0.6,0.6) + xlab("deltaPSI") + ylab("Density (# of DSTs)")  + theme_classic()
## and then zoom in
ggplot(leaf[!leaf$spliceType_toplot == "cryptic",], aes(x = deltapsi, color = spliceType_toplot)) +
geom_line(
  stat = "density") +
  theme_bw() +
  scale_color_manual(values = palette)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("DST effect sizes by category") + scale_fill_manual(values = palette)+ 
   xlim(-0.35,-0.1) + xlab("deltaPSI") + ylab("Density (# of DSTs)")  + theme_classic()
ggplot(leaf[!leaf$spliceType_toplot == "cryptic",], aes(x = deltapsi, color = spliceType_toplot)) +
  geom_line(
    stat = "density") +
  theme_bw() +
  scale_color_manual(values = palette)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("DST effect sizes by category") + scale_fill_manual(values = palette)+ 
  xlim(0.1, 0.35) + xlab("deltaPSI") + ylab("Density (# of DSTs)")  + theme_classic()
dev.off()

#### make the same plots but only with dpsi > 0.05
## barplot of number of splices by type
pdf("/Volumes/data3/anna/splicing_ms_2025/figures/barplot_spliceTypes.pdf")
ggplot(leaf[abs(leaf$deltapsi) > 0.05 & !leaf$spliceType_toplot == "cryptic",], aes(x= fct_rev(fct_infreq(spliceType_toplot)) , fill = spliceType_toplot)) +
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat = "count", aes(label = ..count..), hjust = -0.1) +
  labs(x = "", y = "Number of DSTs") +
  theme_bw()+
  scale_fill_manual(values = palette)+
  theme(axis.text.y = element_text(vjust = 0.5, hjust=1, size = 12))+
  ggtitle("Significant DSTs (IL1B vs control)", subtitle = "p.adjust < 0.05 & deltapsi > 0.05") + ylim(0,400) + coord_flip()
dev.off()

#### volcano plot
library(ggrepel)
pdf("/Volumes/data3/anna/splicing_ms_2025/figures/53HAECs_volcanoPlot_dpsi.pdf")
ggplot(leaf, aes(x = deltapsi, y = -log10(p.adjust))) +
  geom_point(size = 0.75) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Differentially spliced genes (notx vs IL1B)") +
  geom_text_repel(aes(label = ifelse(p.adjust < 0.01 & abs(deltapsi) > 0.1, genes, "")), colour = "red", size = 2.5, force = 0.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
  geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed", color = "gray")
dev.off()

### and make a density plot of the deltapsi
ggplot(leaf, aes(x = deltapsi, color = p.adjust < 0.05 & abs(deltapsi) > 0.05)) + geom_line(size = 1.5,stat = "density") + theme_minimal() + theme(axis.line = element_line(color = "black")) + xlab("deltaPSI (IL1B - notx)") + ylab("Density (# of transcripts)") + ggtitle("Splicing effect sizes of DSTs") + xlim(-0.4,0.4) + scale_color_manual(values = c("lightgray", "steelblue"), name = "p.adjust < 0.05 & \nabs(deltaPSI) > 0.05", labels = c("NS", "Significant")) 


## histogram of effect sizes not as a density plot
ggplot(leaf[leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05,], aes(x = abs(deltapsi))) + geom_histogram(binwidth = 0.05, color = "black", fill = "violetred3") + theme_minimal() + theme(axis.line = element_line(color = "black"), legend.position = "left") + xlab("Magnitude of deltaPSI (IL1B - notx)") + ylab("# of transcripts") + ggtitle("Splicing effect sizes of significant DSTs") 

# and a zoomed in version
ggplot(leaf[leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05,], aes(x = abs(deltapsi))) + geom_histogram(binwidth = 0.05, color = "black", fill = "violetred3") + theme_minimal() + theme(axis.line = element_line(color = "black"), legend.position = "left") + xlab("Magnitude of deltaPSI (IL1B - notx)") + ylab("# of transcripts") + ggtitle("Splicing effect sizes of significant DSTs") + xlim(0.15, 0.7)


### separate out splice types and direction for motif enrichment
splicetypes <- unique(leaf$spliceType)

for (type in splicetypes) {
  notx = leaf[leaf$spliceType == type & leaf$deltapsi < 0,]
  il1b = leaf[leaf$spliceType == type & leaf$deltapsi > 0,]
  assign(paste0(type, "_notx"), notx)
  assign(paste0(type, "_il1b"), il1b)
 write.table(notx, paste0("Y:/anna/splicing_analysis_2024/20240321_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig_dpsi0.05_", "notx_", type, ".txt"), sep = "\t", quote = F, row.names = F)
write.table(il1b, paste0("Y:/anna/splicing_analysis_2024/20240321_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig_dpsi0.05_", "il1b_", type, ".txt"), sep = "\t", quote = F, row.names = F)
  notx.peak = data.frame(uniquePeakID = unlist(notx$intron),
                         chr = unlist(lapply(strsplit(notx$intron, ":"), "[[", 1)),
                         start = unlist(lapply(strsplit(notx$intron, ":"), "[[", 2)),
                         end = unlist(lapply(strsplit(notx$intron, ":"), "[[", A5)),
                         strand = notx$strand)
  row.names(notx.peak) <- NULL
  names(notx.peak) <- NULL
  assign(paste0(type, "_notx_peak"), notx.peak)
  write.table(notx.peak, paste0("Y:/anna/splicing_analysis_2024/20240321_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig_dpsi0.05_", "notx_", type, "peakFile.txt"), sep = "\t", quote = F, row.names = F)
  il1b.peak = data.frame(uniquePeakID = unlist(il1b$intron),
                         chr = unlist(lapply(strsplit(il1b$intron, ":"), "[[", 1)),
                         start = unlist(lapply(strsplit(il1b$intron, ":"), "[[", 2)),
                         end = unlist(lapply(strsplit(il1b$intron, ":"), "[[", A5)),
                         strand = il1b$strand)
  row.names(il1b.peak) <- NULL
  names(il1b.peak) <- NULL
  assign(paste0(type, "_il1b_peak"), il1b.peak)
  write.table(il1b.peak, paste0("Y:/anna/splicing_analysis_2024/20240321_53HAECsIL1Boril1b_final_set_w_depth_pvalueSig_dpsi0.05_", "il1b_", type, "peakFile.txt"), sep = "\t", quote = F, row.names = F)
}


## also save all the results
notx = leaf[leaf$deltapsi < 0,]
il1b = leaf[leaf$deltapsi > 0,]

notx.peak = data.frame(uniquePeakID = unlist(notx$intron),
                       chr = unlist(lapply(strsplit(notx$intron, ":"), "[[", 1)),
                       start = unlist(lapply(strsplit(notx$intron, ":"), "[[", 2)),
                       end = unlist(lapply(strsplit(notx$intron, ":"), "[[", A5)),
                       strand = notx$strand)
row.names(notx.peak) <- NULL
names(notx.peak) <- NULL
write.table(notx.peak, paste0("Y:/anna/splicing_analysis_2024/20240321_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig_dpsi0.05_notx_ALLsplices_peakFile.txt"), sep = "\t", quote = F, row.names = F)
il1b.peak = data.frame(uniquePeakID = unlist(il1b$intron),
                       chr = unlist(lapply(strsplit(il1b$intron, ":"), "[[", 1)),
                       start = unlist(lapply(strsplit(il1b$intron, ":"), "[[", 2)),
                       end = unlist(lapply(strsplit(il1b$intron, ":"), "[[", A5)),
                       strand = il1b$strand)
row.names(il1b.peak) <- NULL
names(il1b.peak) <- NULL
write.table(il1b.peak, paste0("Y:/anna/splicing_analysis_2024/20240321_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig_dpsi0.05_il1b_ALLsplices_peakFile.txt"), sep = "\t", quote = F, row.names = F)



# and in one big dataframe too
all.intron.peaks <- data.frame(uniquePeakID = unlist(leaf$intron),
                            chr = unlist(lapply(strsplit(leaf$intron, ":"), "[[", 1)),
                            start = unlist(lapply(strsplit(leaf$intron, ":"), "[[", 2)),
                            end = unlist(lapply(strsplit(leaf$intron, ":"), "[[", A5)),
                            strand = leaf$strand)
row.names(all.intron.peaks) <- NULL
names(all.intron.peaks) <- NULL
write.table(all.intron.peaks,"Y:/anna/splicing_analysis_2024/20240321_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig_dpsi0.05_ALLsplices_peakFile.txt", sep = "\t", quote = F, row.names = F)
  



#### run fGSEA to identify pathways in the

library(clusterProfiler)
require(DOSE)  
organism <- "org.Hs.eg.db"
library(organism, character.only = T)
keytypes(org.Hs.eg.db)



# pull out gene names and deltapsi for il1b regulated transcripts
# we can't run this as "normal" with notx vs il1b because the gene symbols are the same and GSEA is not optimized for transcript specific things so just take the list of DSGs with deltapsi > 0
gene_list <- il1b$deltapsi
names(gene_list) <- unlist(lapply(strsplit(as.character(il1b$genes), split = ","), "[[", 1))
#names(gene_list) <- unlist(lapply(strsplit(as.character(il1b$transcript), split = "\\+"), "[[", 1))

# order by deltapsi
gene_list <- gene_list[order(gene_list, decreasing = T)]



gse.GO <- gseGO(geneList=gene_list, 
                ont ="ALL", #can be set to "BP" (biological process), "MF" (molecular function), "CC" (cell component), or "ALL" depending on what gene sets you want to test
                keyType = "SYMBOL",
                minGSSize = 10, #minimum number of genes in sets
                maxGSSize = 500, #max number of genes in gene sets
                pvalueCutoff = 0.1, #significance threshold
                verbose = TRUE, 
                OrgDb = organism, 
                pAdjustMethod = "none",#can be BH, BY, fdr, or none
                by = "fgsea",
                eps = 0,
                exponent = 1) 
gse.GO.BP <- gseGO(geneList=gene_list, 
                  ont ="BP", #can be set to "BP" (biological process), "MF" (molecular function), "CC" (cell component), or "ALL" depending on what gene sets you want to test
                  keyType = "SYMBOL",
                  minGSSize = 10, #minimum number of genes in sets
                  maxGSSize = 500, #max number of genes in gene sets
                  pvalueCutoff = 0.1, #significance threshold
                  verbose = TRUE, 
                  OrgDb = organism, 
                  pAdjustMethod = "none",#can be BH, BY, fdr, or none
                  by = "fgsea",
                  eps = 0,
                  exponent = 1) 


dotplot(gse.GO, showCategory= 20) + ggtitle("All pathway enrichment in DSGs (p.adjust < 0.05)")
dotplot(gse.GO.BP, showCategory= 20) + ggtitle("BP enrichment in DSGs (p.adjust < 0.05)")


## dpsi restriction
gene_list2 <- il1b$deltapsi[il1b$deltapsi >0.05]
names(gene_list2) <- unlist(lapply(strsplit(as.character(il1b$genes[il1b$deltapsi >0.05]), split = ","), "[[", 1))
#names(gene_list) <- unlist(lapply(strsplit(as.character(il1b$transcript), split = "\\+"), "[[", 1))

# order by deltapsi
gene_list2 <- gene_list2[order(gene_list2, decreasing = T)]

gse.GO.all_strict <- gseGO(geneList=gene_list2, 
                          ont ="ALL", #can be set to "BP" (biological process), "MF" (molecular function), "CC" (cell component), or "ALL" depending on what gene sets you want to test
                          keyType = "SYMBOL",
                          minGSSize = 10, #minimum number of genes in sets
                          maxGSSize = 500, #max number of genes in gene sets
                          pvalueCutoff = 0.1, #significance threshold
                          verbose = TRUE, 
                          OrgDb = organism, 
                          pAdjustMethod = "none",#can be BH, BY, fdr, or none
                          by = "fgsea",
                          eps = 0,
                          exponent = 1)

dotplot(gse.GO.all_strict, showCategory= 10) + ggtitle("DSG pathway enrichment")

# save GSE results
gse_res <- as.data.frame(gse.GO.all_strict@result)
write.table(gse_res, file = "leafcutter_DSGs_fgsea_geneSetEnrichment_results.txt" , sep = "\t", quote = F, row.names = F)
#gse_res <- read.delim("leafcutter_DSGs_fgsea_geneSetEnrichment_results.txt", as.is = T, stringsAsFactors = F)

# dotplot by NES
uniquecategories <- c("response to toxic substance", "response to oxidative stress", "positive regulation of transcription by RNA polymerase II", "organic acid biosynthetic process", "positive regulation of programmed cell death", "positive regulation of cytokine production", "phosphoric ester hydrolase activity", "RNA splicing, via transesterification reactions")


ggplot(gse_res[gse_res$Description %in% uniquecategories,], aes(x = NES, y = reorder(Description, NES), color = p.adjust, size = setSize)) + geom_point() + scale_color_gradient(low = "red", high = "steelblue") + theme_minimal() + ylab("") + xlab("Normalized Enrichment Score") + ggtitle("Pathway Enrichment in IL1B-DSGs") + theme(axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.75, vjust = 1.5,face = "bold"), axis.text = element_text(size = 13))




for (type in splicetypes) {
gene_list <- il1b$deltapsi[il1b$spliceType == type]
names(gene_list) <- unlist(lapply(strsplit(as.character(il1b$genes[il1b$spliceType == type]), split = ","), "[[", 1))
gene_list <- gene_list[order(gene_list, decreasing = T)]
#take the top 250 or all the genes if less than 250
length <- ifelse(length(gene_list) > 250, 250, length(gene_list))
gene_list <- gene_list[1:length]

gse.GO <- gseGO(geneList=gene_list, 
                ont ="BP", #can be set to "BP" (biological process), "MF" (molecular function), "CC" (cell component), or "ALL" depending on what gene sets you want to test
                keyType = "SYMBOL",
                minGSSize = 10, #minimum number of genes in sets
                maxGSSize = 500, #max number of genes in gene sets
                pvalueCutoff = 0.05, #significance threshold
                verbose = TRUE, 
                OrgDb = organism, 
                pAdjustMethod = "none",#can be BH, BY, fdr, or none
                by = "fgsea",
                eps = 0,
                exponent = 1,
                nPermSimple = 10000) 
assign(paste0(type, "_gene_list"), gene_list)
assign( paste0(type,"_GSE.go"), gse.GO)

dotplot(gse.GO, showCategory= 20) + ggtitle(paste("All pathway enrichment in DSGs (p.adjust < 0.05)", type, sep = " "))
}


dotplot(AF_GSE.go, showCategory= 10) + ggtitle("All pathway enrichment in AF DSGs (p.adjust < 0.05)")


## or compare to each other
library(dplyr)
top300 <- leaf[c("genes", "deltapsi", "spliceType")] %>% group_by(spliceType) %>% top_n(n = 300, wt = deltapsi)
df <- split(top300$genes, top300$spliceType)

df$`AL` = bitr(df$'AL', fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df$'AF' = bitr(df$'AF', fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df$'A3' = bitr(df$'A3', fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df$'A5' = bitr(df$`A5`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df$'MX' = bitr(df$'MX', fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df$`SE` = bitr(df$`SE`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df$`RI` = bitr(df$`RI`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df$`cryptic_threeprime` = bitr(df$`cryptic_threeprime`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df$`cryptic_fiveprime` = bitr(df$`cryptic_fiveprime`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
df$`cryptic_novel annotated pair` = bitr(df$`cryptic_novel annotated pair`, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")



#do the same here, a line like below for each cluster
genelist <- list("AL" = df$'AL'$SYMBOL, 
                 "AF" = df$'AF'$SYMBOL,
                 "A3" = df$'A3'$SYMBOL,
                 "A5" = df$`A5`$SYMBOL,
                 "MX" = df$'MX'$SYMBOL,
                 "SE" = df$`SE`$SYMBOL,
                 "RI" = df$`RI`$SYMBOL,
                 "cryptic_threeprime" = df$`cryptic_threeprime`$SYMBOL,
                 "cryptic_fiveprime" = df$`cryptic_fiveprime`$SYMBOL,
                 "cryptic_novel annotated pair" = df$`cryptic_novel annotated pair`$SYMBOL)

GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO",  pvalueCutoff = 1, pAdjustMethod = "BH", ont = "BP", OrgDb = organism, keyType = "SYMBOL")
dotplot(GOclusterplot, showCategory = 5)


GO_simplify <- simplify(GOclusterplot)
dotplot(GO_simplify, showCategory = 3)


### it seems that lots of genes in the enriched pathways are AF genes, is that significant?
head(gse_res)
gse_res %>% filter(p.adjust < 0.05) %>% dplyr::select(core_enrichment, Description) %>% mutate(genes = strsplit(core_enrichment, "/")) -> gse_genes

af_genes <- unique(leaf$genes[leaf$spliceType == "AF" & abs(leaf$deltapsi) > 0.05])
se_genes <- unique(leaf$genes[leaf$spliceType == "SE" & abs(leaf$deltapsi) > 0.05])

count <- list()
total <- list()
for (i in 1:nrow(gse_genes)) {
gse_genes[i,] %>% mutate(count = sum(unlist(gse_genes$genes[i]) %in% af_genes), total = length(unlist(gse_genes$genes[i]))) %>% dplyr::select(count, total) -> tab
  
count[i] <- tab$count
total[i] <- tab$total

}

prop_af <- unlist(count)/unlist(total)
mean(prop_af)
length(af_genes)/length(unique(leaf$genes))

pdf("/Volumes/data3/anna/splicing_ms_2025/figures/AF_DSGs_inGOterms.pdf", width = 5, height = 5)
ggplot(as.data.frame(prop_af), aes(x = prop_af)) + geom_histogram(fill = "orchid", color = "black") + theme_bw() + geom_vline(xintercept = length(af_genes)/length(unique(leaf$genes)), linetype = "dashed", color = "red") + xlim(0,1) + xlab("Proportion of AF genes in GO term") + ylab("# of GO terms") + ggtitle("AF-DSG enrichment in GO terms") + geom_label(label = "% of AF-DSGs", y = 11, x = 0.2, colour = "red", label.size = 0)
dev.off()

count_se <- list()
total_se <- list()
for (i in 1:nrow(gse_genes)) {
  gse_genes[i,] %>% mutate(count = sum(unlist(gse_genes$genes[i]) %in% se_genes), total = length(unlist(gse_genes$genes[i]))) %>% dplyr::select(count, total) -> tab
  
  count_se[i] <- tab$count
  total_se[i] <- tab$total
  
}

prop_se <- unlist(count_se)/unlist(total_se)
mean(prop_se)
length(se_genes)/length(unique(leaf$genes))

pdf("/Volumes/data3/anna/splicing_ms_2025/figures/SE_DSGs_inGOterms.pdf", width = 5, height = 5)
ggplot(as.data.frame(prop_se), aes(x = prop_se)) + geom_histogram(fill = "steelblue", color = "black") + theme_bw() + geom_vline(xintercept = length(se_genes)/length(unique(leaf$genes)), linetype = "dashed", color = "red") + xlim(0,1) + xlab("Proportion of SE genes in GO term") + ylab("# of GO terms")+ ggtitle("SE-DSG enrichment in GO terms") + geom_label(label = "% of SE-DSGs", y = 15, x = 0.2, colour = "red", label.size = 0) 
dev.off()




