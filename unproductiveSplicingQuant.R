##### looking for unproductive vs productive splicing events 

### unproductive is determined by stop codon within 50 bp of the exon start, or the hg38 from gencode actually has an NMD label

library(rtracklayer)
library(dplyr)
library(purrr)
library(GenomicRanges)

hg38 <- readGFF("/path/hg38_anno_files/gencode.v41.annotation.gtf")
names(hg38)
table(hg38$type)
stops <- hg38[hg38$type == "stop_codon",]

## match up with leafcutter results

leaf <- read.delim("/path/leafcutter_results_53HAECsIL1Bornotx_final_set_w_depth.txt", as.is = T, stringsAsFactors = F)

# add start and stops
leaf$start_chr <- vapply(strsplit(leaf$chrjunction, split = "-"), function(x) paste(x[c(1,2)], collapse = ":"), character(1L))
leaf$end_chr <- vapply(strsplit(leaf$chrjunction, split = "-"), function(x) paste(x[c(1,3)], collapse = ":"), character(1L))
leaf$start <- vapply(strsplit(leaf$chrjunction, split = "-"), function(x) paste(x[c(2)], collapse = ":"), character(1L))
leaf$end <- vapply(strsplit(leaf$chrjunction, split = "-"), function(x) paste(x[c(3)], collapse = ":"), character(1L))
leaf$chr <- vapply(strsplit(leaf$chrjunction, split = "-"), function(x) paste(x[c(1)], collapse = ":"), character(1L))
row.names(leaf) <- leaf$intron

# add chr to start and end in the gtf
stops$start_chr <- paste(stops$seqid, stops$start, sep = ":")
stops$end_chr <- paste(stops$seqid, stops$end, sep = ":")

# need to make stops non-redundant
stops.info <- stops[,c("seqid",  "start", "end", "gene_name", "type", "tag", "transcript_type")]
stops.info <- stops.info[!duplicated(stops.info),]
stops.info[grep("PFKFB3", stops.info$gene_name),] ## removed duplicate stop codons by position successfully

distances <- list()
nearestStop <- list()
for (x in 1:length(leaf$intron)){
  print(leaf$intron[x])
  gene <- leaf$genes[x]
  # pull all stop codons for the gene
  stop.df <- stops.info[stops.info$gene_name==gene,]
  
  #find all distances
  dist <- data.frame(
    from3 = as.numeric(leaf$start[x]) - stop.df$start,
    from5 = as.numeric(leaf$end[x]) - stop.df$start,
    stopCodon = paste(stop.df$seqid, stop.df$start, stop.df$end, sep = ':'),
    stopInfo = paste(stop.df$tag, stop.df$transcript_type, sep = ":")
  )
  # find minimum distance (the closest stop codon)
  nearest <- list(
    from3 = unique(dist$from3[abs(dist$from3) == min(abs(dist$from3))]),
    from3_stopCodon = unique(dist$stopCodon[abs(dist$from3) == min(abs(dist$from3))]),
    from3_stopInfo = unique(dist$stopInfo[abs(dist$from3) == min(abs(dist$from3))]),
    from5 = unique(dist$from5[abs(dist$from5) == min(abs(dist$from5))]),
    from5_stopCodon = unique(dist$stopCodon[abs(dist$from5) == min(abs(dist$from5))]),
    from5_stopInfo = unique(dist$stopInfo[abs(dist$from5) == min(abs(dist$from5))])
  )
  intron <- leaf$intron[x]
  # save into a list
  distances[[intron]] <- dist
  nearestStop[[intron]] <- nearest
}

## turn into an easier format
nearestStop <- as.data.frame(t(as.data.frame(do.call(cbind, nearestStop))))
nearestStop$intron <- row.names(nearestStop)

##add into the dataframe
# this is complicated because some values are missing, not NA, missing.
leaf$stop_3prime <- as.character(gsub("character(0)", NA,as.character(nearestStop$from3_stopInfo[match(leaf$intron, nearestStop$intron)])))
leaf$stop_3prime_dist <- as.numeric(gsub("numeric(0)", NA, as.character(nearestStop$from3[match(leaf$intron, nearestStop$intron)])))
leaf$stop_5prime <- as.character(gsub("character(0)", NA,as.character(nearestStop$from5_stopInfo[match(leaf$intron, nearestStop$intron)])))
leaf$stop_5prime_dist <- as.numeric(gsub("numeric(0)", NA, as.character(nearestStop$from5[match(leaf$intron, nearestStop$intron)])))

## add category as a column
nmd_introns3 <- leaf$intron[grep("nonsense_mediated_decay", leaf$stop_3prime)]
leaf$stopCodonCategory_3prime<- ifelse(leaf$intron %in% nmd_introns3, "nonsense_mediated_decay", "protein coding")

nmd_introns5 <- leaf$intron[grep("nonsense_mediated_decay", leaf$stop_5prime)]
leaf$stopCodonCategory_5prime<- ifelse(leaf$intron %in% nmd_introns5, "nonsense_mediated_decay", "protein coding")
table(leaf$stopCodonCategory_3prime[leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05])

sig.nmd <- leaf[leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05 & leaf$stopCodonCategory_3prime == "nonsense_mediated_decay" | leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05 & leaf$stopCodonCategory_5prime == "nonsense_mediated_decay",]

### a negative number indicates that the stop codon is upstream of the intron position, for a + sense gene this means the stop codon is in the exon
# make a stranded copy to plot
leaf$stop_3prime_dist_str <- ifelse(leaf$strand == "-", -1*leaf$stop_3prime_dist, leaf$stop_3prime_dist)
leaf$stop_5prime_dist_str <- ifelse(leaf$strand == "-", -1*leaf$stop_5prime_dist, leaf$stop_5prime_dist)

# group cryptic splices together
leaf$spliceType_toplot <- ifelse(leaf$spliceType %in% c("cryptic_unanchored"," cryptic_unanchored",  "cryptic_threeprime", "cryptic_fiveprime", "unknown_strand", "novel annotated pair"), "cryptic", leaf$spliceType)
table(leaf$spliceType_toplot)
leaf$spliceType_toplot <- factor(leaf$spliceType_toplot, levels = c("AF", "AL", "A3", "A5", "SE","MX", "RI", "cryptic"))

library(ggplot2)
ggplot(leaf, aes(x = stop_3prime_dist_str)) +
         geom_histogram(stat = "density") + xlim(-200,200) + facet_grid(rows = vars(p.adjust < 0.05 & abs(deltapsi) > 0.05)) +
  ggtitle("significant vs nonsignificant DSTs distribution of NMD to protein coding")## similar distribution among sig and nonsig

ggplot(leaf[leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05,], aes(x = stop_3prime_dist_str, fill = stopCodonCategory_3prime)) +
  geom_histogram(stat = "density") + xlim(-200,200) + facet_grid(rows = vars(spliceType_toplot), cols = vars(stopCodonCategory_3prime)) + theme_minimal() + geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("VariGene HAECs NMD splicing summary by 3' splice junction", subtitle = "negative values indicate stop codons in the exon")
## this makes it look super dramatic

### make a barplot of how many events this is
leaf$spliceType_toplot <- factor(leaf$spliceType_toplot, levels = c("SE", "AF", "AL", "A3", "A5", "MX", "RI"))
leaf$direction = ifelse(leaf$deltapsi > 0, "Increased with IL1B", "Decreased with IL1B")
leaf$NMD <- ifelse(leaf$stopCodonCategory_3prime == "nonsense_mediated_decay" & leaf$stop_3prime_dist_str < 50, "NMD", "protein coding")

table(leaf$NMD == 'NMD', leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05)
305/1224 # NMD
912/1224 # protein-coding

# split by up and down in Il1B
ggplot(na.omit(leaf[leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05 & !leaf$spliceType_toplot == "cryptic",]), aes(x = spliceType_toplot, fill = NMD)) + theme_bw() +
  geom_bar(position = "fill") + scale_fill_manual(values = c("skyblue1", "coral2"), name = "Transcript type") + ylab("Proportion of DSTs") + xlab("Splice type") + ggtitle("NMD and protein coding transcripts regulated by IL1B") + 
  facet_grid(cols = vars(direction)) + geom_hline(yintercept = 0.75, linetype = "dashed")

ggplot(na.omit(leaf[leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05 & !leaf$spliceType_toplot == "cryptic",]), aes(x = spliceType_toplot, fill = NMD)) + theme_bw() +
  geom_bar(position = "fill") + scale_fill_manual(values = c("skyblue1", "coral2"), name = "Transcript type") + ylab("Proportion of DSTs") + xlab("Splice type") + ggtitle("NMD and protein coding DSTs")  + geom_hline(yintercept = 0.75, linetype = "dashed") + theme(text = element_text(size = 12))

### save
write.table(leaf, "/path/leafcutter_53HAECsIL1Bornotx_annotated_NMDorProteinCoding.txt", sep = "\t", quote = F)

df <- data.frame(table(sign(leaf$deltapsi[leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05 & leaf$stopCodonCategory_3prime == "nonsense_mediated_decay" & leaf$stop_3prime_dist_str < 50]), leaf$spliceType_toplot[leaf$p.adjust < 0.05 & abs(leaf$deltapsi) > 0.05 & leaf$stopCodonCategory_3prime == "nonsense_mediated_decay" & leaf$stop_3prime_dist_str < 50]))

ggplot(df, aes(x = Var2, y = ifelse(Var1 == 1,Freq, -1*Freq), fill = Var1)) + geom_bar(stat = "identity") + coord_flip() + theme_minimal() + ggtitle("NMD events", subtitle = "Decreased                   Increased") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.3)) + scale_fill_manual(values = c("steelblue", "darkred")) + ylab("number of DSTs") + xlab("")


#### annotating protein-coding DSTs as UTR vs protein-coding exons
hg38 %>% filter(gene_name %in% leaf$genes, type %in% c("exon", "UTR")) %>% tibble() %>%
  mutate(chrstart = paste(seqid, start, sep = "-"), chrend = paste(seqid, end, sep = "-"))-> exons

exons %>% filter(gene_name == "PFKFB3",  start == 6213623, type == "UTR") -> test

leaf %>% 
  filter(spliceType %in% c("AF", "AL"), abs(deltapsi) > 0.05, p.adjust < 0.05) %>%
  # set up start and end with chr and pos
  mutate(start = vapply(strsplit(chrjunction, "-"), function(x) paste(x[c(1,2)], collapse = "-"), character(1L)),
         end = vapply(strsplit(chrjunction, "-"), function(x) paste(x[c(1,3)], collapse = "-"), character(1L)),
         
         ## match to exons on the 5' and 3' ends
        exons5prime = map(start, ~ unique(exons$type[exons$chrend == .x])),
        exons3prime = map(end, ~ unique(exons$type[exons$chrstart == .x])),
        
        exons5prime_transcript = map(start, ~ unique(exons$transcript_id[exons$chrend == .x])),
        exons3prime_transcript = map(end, ~ unique(exons$transcript_id[exons$chrstart == .x]))) -> dsgs_ann

unique(dsgs_ann$exons5prime)



#### and for sQTLs:

sqtl.notx <- read.delim("/path/notx/splicing_QTLs_notx_rna_cis1000bpCis.txt",as.is = T, stringsAsFactors = F)
head(sqtl.notx)

sqtl.notx$start <- vapply(strsplit(sqtl.notx$gene, split = ":"), function(x) paste(x[c(2)], collapse = ":"), character(1L))
sqtl.notx$end <- vapply(strsplit(sqtl.notx$gene, split = ":"), function(x) paste(x[c(3)], collapse = ":"), character(1L))
sqtl.notx$chr <- paste0("chr", vapply(strsplit(sqtl.notx$gene, split = ":"), function(x) paste(x[c(1)], collapse = ":"), character(1L)))

nearestStop <- list()
for (x in 1:length(sqtl.notx$gene)){
  print(sqtl.notx$gene[x])
  chr <- sqtl.notx$chr[x]
  # pull all stop codons for the gene
  stop.df <- stops.info[stops.info$seqid==chr,]
  
  #find all distances
  dist <- data.frame(
    from3 = as.numeric(sqtl.notx$start[x]) - stop.df$start,
    from5 = as.numeric(sqtl.notx$end[x]) - stop.df$start,
    stopCodon = paste(stop.df$seqid, stop.df$start, stop.df$end, sep = ':'),
    stopInfo = paste(stop.df$tag, stop.df$transcript_type, sep = ":")
  )
  # find minimum distance (the closest stop codon)
  nearest <- list(
    from3 = unique(dist$from3[abs(dist$from3) == min(abs(dist$from3))]),
    from3_stopCodon = unique(dist$stopCodon[abs(dist$from3) == min(abs(dist$from3))]),
    from3_stopInfo = unique(dist$stopInfo[abs(dist$from3) == min(abs(dist$from3))]),
    from5 = unique(dist$from5[abs(dist$from5) == min(abs(dist$from5))]),
    from5_stopCodon = unique(dist$stopCodon[abs(dist$from5) == min(abs(dist$from5))]),
    from5_stopInfo = unique(dist$stopInfo[abs(dist$from5) == min(abs(dist$from5))])
  )
  intron <- sqtl.notx$gene[x]
  # save into a list
  distances[[intron]] <- dist
  nearestStop[[intron]] <- nearest
}

## turn into an easier format
nearestStop <- as.data.frame(t(as.data.frame(do.call(cbind, nearestStop))))
nearestStop$intron <- row.names(nearestStop)
head(nearestStop)

##add into the dataframe
# this is complicated because some values are missing, not NA, missing.
sqtl.notx$stop_3prime <- as.character(gsub("character(0)", NA,as.character(nearestStop$from3_stopInfo[match(sqtl.notx$gene, nearestStop$intron)])))
sqtl.notx$stop_3prime_dist <- as.numeric(gsub("numeric(0)", NA, as.character(nearestStop$from3[match(sqtl.notx$gene, nearestStop$intron)])))
sqtl.notx$stop_5prime <- as.character(gsub("character(0)", NA,as.character(nearestStop$from5_stopInfo[match(sqtl.notx$gene, nearestStop$intron)])))
sqtl.notx$stop_5prime_dist <- as.numeric(gsub("numeric(0)", NA, as.character(nearestStop$from5[match(sqtl.notx$gene, nearestStop$intron)])))

## add category as a column
nmd_introns3 <- sqtl.notx$gene[grep("nonsense_mediated_decay", sqtl.notx$stop_3prime)]
sqtl.notx$stopCodonCategory_3prime<- ifelse(sqtl.notx$gene %in% nmd_introns3, "nonsense_mediated_decay", "protein coding")

nmd_introns5 <- sqtl.notx$gene[grep("nonsense_mediated_decay", sqtl.notx$stop_5prime)]
sqtl.notx$stopCodonCategory_5prime<- ifelse(sqtl.notx$gene %in% nmd_introns5, "nonsense_mediated_decay", "protein coding")

## results: looks like the same proportion are NMD as are protein coding
table(sqtl.notx$stopCodonCategory_3prime[sqtl.notx$FDR < 0.05 & abs(sqtl.notx$beta) > 0.1])
table(sqtl.notx$stopCodonCategory_5prime[sqtl.notx$FDR < 0.05 & abs(sqtl.notx$beta) > 0.1])
table(sqtl.notx$stopCodonCategory_3prime[sqtl.notx$FDR > 0.05 & abs(sqtl.notx$beta) > 0.1])
table(sqtl.notx$stopCodonCategory_5prime[sqtl.notx$FDR > 0.05 & abs(sqtl.notx$beta) > 0.1])

### what about the results between ERG and NFkB inducible and repressed?


