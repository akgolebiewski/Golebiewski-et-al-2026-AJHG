######### pull first exons
library(rtracklayer)

# load GTF file used for leafcutter analysis (this is huge, clear quickly after taking what is needed out)
hg38 <- readGFF("/path/hg38_anno/gencode.v26.annotation.gtf")
names(hg38)

firstExons <- hg38[hg38$exon_number == "1" & !is.na(hg38$exon_number) & hg38$type == 'exon',]

#rm(hg38) #make space again

dim(firstExons)
firstExons[firstExons$gene_name == "PFKFB3",]
table(firstExons$transcript_type)


##### now match up the first exons to the introns
# read splices
leaf.all <- read.delim("/path/20240321_annotatedFINAL_leafcutter_results_53HAECsIL1Bornotx_final_set_w_depth_pvalueSig.txt", as.is = T, stringsAsFactors = F)

# restrict to only results with a 5% deltapsi and with known strand--otherwise the bedfile won't work for granges
leaf <- leaf.all[leaf.all$p.adjust < 0.05 & abs(leaf.all$deltapsi) > 0.05 & !is.na(leaf.all$strand),]

leaf$intron_start <- as.numeric(unlist(lapply(strsplit(leaf$junction, split = "-"), "[[", 1)))
leaf$intron_end <- as.numeric(unlist(lapply(strsplit(leaf$junction, split = "-"), "[[", 2)))

# get just one gene symbol
leaf$gene.sym <- unlist(lapply(strsplit(leaf$genes, ","), "[[", 1))
#SEPTIN10 got subbed for "10-Sep" at some point-- fix 
leaf$gene.sym <- sub("10-Sep", "SEPTIN10", leaf$gene.sym)

#restrict to just GOIs
firstExons <- firstExons[firstExons$gene_name %in% leaf$gene.sym,]
table(firstExons$gene_name)
sum(leaf$gene.sym %in% firstExons$gene_name)

library(GenomicRanges)

leaf.pos <- data.frame(chr = unlist(lapply(strsplit(leaf$intron, ":"), "[[", 1)),
                       start = ifelse(leaf$strand == "+",leaf$intron_start,leaf$intron_end),
                       end = ifelse(leaf$strand == "+",leaf$intron_start + 10,leaf$intron_end + 10), ## use this instead of the end of the exon because it's getting confused when the intron is huge
                       strand = leaf$strand
)

FE.pos <- data.frame(chr = firstExons$seqid,
                        start = firstExons$start,
                        end = firstExons$end,
                        strand = firstExons$strand,
                     gene = firstExons$gene_name
                        )

#make granges objects
leaf.grange <- makeGRangesFromDataFrame(leaf.pos)
FE.grange <- makeGRangesFromDataFrame(FE.pos)

nearest <- nearest(leaf.grange, FE.grange, ignore.strand = F)

leaf$FE_start <- FE.pos$start[match(nearest, row.names(FE.pos))]
leaf$FE_end <- FE.pos$end[match(nearest, row.names(FE.pos))]
leaf$FE_gene <- FE.pos$gene[match(nearest, row.names(FE.pos))]

# check work
## for non-AF splices these look great, but there are some AF splices that are obviously wrong

AF <- leaf[leaf$spliceType == "AF",]
summary(AF$intron_start - AF$FE_end) # max?
sum(AF$intron_start - AF$FE_end == 0 & AF$strand == "+") + sum(AF$intron_end - AF$FE_start == 0 & AF$strand == "-")
sum(AF$intron_end - AF$FE_start == 0 & AF$strand == "-")

AF$match <- ifelse(AF$intron_start - AF$FE_end == 0 & AF$strand == "+" | AF$intron_end - AF$FE_start == 0 & AF$strand == "-", "yes", "no")
table(AF$match)
AF$distance <- ifelse(AF$strand == "+", AF$intron_start - AF$FE_end, AF$intron_end - AF$FE_start)
hist(AF$distance, breaks = 50)
table(AF$distance)

# separate the 91 wrong ones and fix them-- the issue is they're not annotated as a first exon in the genome so they didn't get pulled right
wrong <- AF[AF$match == "no",]
FE_2 <- hg38[hg38$gene_name %in% AF$gene.sym & hg38$type == "exon",]

wrong.pos <- data.frame(chr = unlist(lapply(strsplit(wrong$intron, ":"), "[[", 1)),
                       start = ifelse(wrong$strand == "+",wrong$intron_start,wrong$intron_end),
                       end = ifelse(wrong$strand == "+",wrong$intron_start + 10,wrong$intron_end + 10), ## use this instead of the end of the exon because it's getting confused when the intron is huge
                       strand = wrong$strand
)

FE2.pos <- data.frame(chr = FE_2$seqid,
                     start = FE_2$start,
                     end = FE_2$end,
                     strand = FE_2$strand,
                     gene = FE_2$gene_name
)
#make granges objects
wrong.grange <- makeGRangesFromDataFrame(wrong.pos)
FE.grange <- makeGRangesFromDataFrame(FE2.pos)

nearest <- nearest(wrong.grange, FE.grange, ignore.strand = F)

wrong$FE_start <- FE2.pos$start[match(nearest, row.names(FE2.pos))]
wrong$FE_end <- FE2.pos$end[match(nearest, row.names(FE2.pos))]
wrong$FE_gene <- FE2.pos$gene[match(nearest, row.names(FE2.pos))]


wrong$distance <- ifelse(wrong$strand == "+", wrong$intron_start - wrong$FE_end, wrong$intron_end - wrong$FE_start)
sum(wrong$distance == 0) # this fixed 43 of 91
hist(wrong$distance, breaks = 100)

## fix manually the 48 that won't run right because they have the opposite end of the junction matching
wrong$chr <- unlist(lapply(strsplit(wrong$chrjunction, split = "-"), "[[", 1))
wrong$e1end <- ifelse(wrong$strand == "+", paste(wrong$chr, wrong$intron_end, sep = ":"), paste(wrong$chr, wrong$intron_start, sep = ":"))
hg38$position <- ifelse(hg38$strand == "+", paste(hg38$seqid, hg38$start, sep = ":"), paste(hg38$seqid, hg38$end, sep = ":"))

# add back in
wrong$FE_start <- ifelse(wrong$distance == 0, wrong$FE_start, hg38$start[match(wrong$e1end, hg38$position)])
wrong$FE_end <- ifelse(wrong$distance == 0, wrong$FE_end, hg38$end[match(wrong$e1end, hg38$position)])

### add these back into the big dataframe
AF$FE_start <- ifelse(AF$distance==0, AF$FE_start, wrong$FE_start[match(AF$intron, wrong$intron)])
AF$FE_end <- ifelse(AF$distance==0, AF$FE_start, wrong$FE_end[match(AF$intron, wrong$intron)])



### make peak files with good promoter hits
AF.notx <- AF[AF$deltapsi > 0,]
AF.notx.peak <- data.frame(
  uniqueID = paste(AF.notx$intron, AF.notx$gene.sym, AF.notx$spliceType, sep = ":"),
  chromosome = unlist(lapply(strsplit(AF.notx$cluster, split = ":"), "[[", 1)),
  start = ifelse(AF.notx$strand == "+", AF.notx$FE_start - 200, AF.notx$FE_end), # 200 bp UPstream of the 5'UTR for +, the 5' end of the 5' UTR for -
  end = ifelse(AF.notx$strand == "+", AF.notx$FE_start, AF.notx$FE_end + 200), # the opposite
  strand = AF.notx$strand)
#write.table(AF.notx.peak, "/path/20240523_200bp_promoter_AFnotx.txt", sep = "\t", row.names = F, quote = F)

AF.il1b <- AF[AF$deltapsi < 0,]
AF.il1b.peak <- data.frame(
  uniqueID = paste(AF.il1b$intron, AF.il1b$gene.sym, AF.il1b$spliceType, sep = ":"),
  chromosome = unlist(lapply(strsplit(AF.il1b$cluster, split = ":"), "[[", 1)),
  start = ifelse(AF.il1b$strand == "+", AF.il1b$FE_start - 200, AF.il1b$FE_end), # 200 bp UPstream of the 5'UTR for +, the 5' end of the 5' UTR for -
  end = ifelse(AF.il1b$strand == "+", AF.il1b$FE_start, AF.il1b$FE_end + 200), # the opposite
  strand = AF.il1b$strand)
#write.table(AF.il1b.peak, "/path/20240523_200bp_promoter_AFil1b.txt", sep = "\t", row.names = F, quote = F)



### for non AF splices the promoters/first exons will be the same so we don't need to separate by notx/IL1B but we do need to save the promoters
other <- leaf[!leaf$spliceType == 'AF',]
other.peak <- data.frame(
  uniqueID = paste(other$cluster, other$gene.sym, other$spliceType, sep = ":"),
  chromosome = unlist(lapply(strsplit(other$cluster, split = ":"), "[[", 1)),
  start = ifelse(other$strand == "+", other$FE_start - 200, other$FE_end), # 200 bp UPstream of the 5'UTR for +, the 5' end of the 5' UTR for -
  end = ifelse(other$strand == "+", other$FE_start, other$FE_end + 200), # the opposite
  strand = other$strand)

# remove duplicates
other.peak.nodups <- other.peak[!duplicated(other.peak),]

#write.table(other.peak.nodups, "/path/20240523_200bp_promoter_nonAFsplices.txt", sep = "\t", row.names = F, quote = F)


####### making larger promoter regions
## 500 bp downstream
## 1 kb upstream
## why? to capture more signal, when you look in the browser you can see we're missing peaks and only getting some of the promoter signal especially for the ATAC and H3K27ac


### make peak files with good promoter hits
AF.notx <- AF[AF$deltapsi > 0,]
AF.notx.peak <- data.frame(
  uniqueID = paste(AF.notx$intron, AF.notx$gene.sym, AF.notx$spliceType, sep = ":"),
  chromosome = unlist(lapply(strsplit(AF.notx$cluster, split = ":"), "[[", 1)),
  start = ifelse(AF.notx$strand == "+", AF.notx$FE_start - 1000, AF.notx$FE_end - 500), 
  end = ifelse(AF.notx$strand == "+", AF.notx$FE_start + 500, AF.notx$FE_end + 1000),
  strand = AF.notx$strand)
write.table(AF.notx.peak, "/path/20241203_1500bp_promoter_AFnotx.txt", sep = "\t", row.names = F, quote = F)

AF.notx.peak$start - AF.notx.peak$end

AF.il1b <- AF[AF$deltapsi < 0,]
AF.il1b.peak <- data.frame(
  uniqueID = paste(AF.il1b$intron, AF.il1b$gene.sym, AF.il1b$spliceType, sep = ":"),
  chromosome = unlist(lapply(strsplit(AF.il1b$cluster, split = ":"), "[[", 1)),
  start = ifelse(AF.il1b$strand == "+", AF.il1b$FE_start - 1000, AF.il1b$FE_end - 500), # 200 bp UPstream of the 5'UTR for +, the 5' end of the 5' UTR for -
  end = ifelse(AF.il1b$strand == "+", AF.il1b$FE_start + 500, AF.il1b$FE_end + 1000), # the opposite
  strand = AF.il1b$strand)
write.table(AF.il1b.peak, "/path/20241203_1500bp_promoter_AFil1b.txt", sep = "\t", row.names = F, quote = F)

AF.il1b.peak$start - AF.il1b.peak$end
