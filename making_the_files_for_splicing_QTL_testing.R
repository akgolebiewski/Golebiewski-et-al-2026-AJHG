tog<-read.delim("X:/leafcutter/hg38_bamfiles/final_set/Final_set_hg38_perind.counts.gz",
                sep = " ", header = T, as.is = T, stringsAsFactors = F)
chromosomes<-unlist(lapply(strsplit(tog$chrom, split = ":"), "[[",1))
table(chromosomes)

#seperate out the bad chromosomes
tog_filt<-tog[-which(chromosomes%in%c("GL000205.2", "KI270733.1")),]

#seperate by treat
nt<-tog_filt[,grep("NT", names(tog_filt))]
il1b<-tog_filt[,grep("IL1", names(tog_filt))]
nt<-cbind(tog_filt$chrom, nt)
il1b<-cbind(tog_filt$chrom, il1b)


names(nt)[1]<-"chrom"
names(il1b)[1]<-"chrom"

write.table(nt, "X:/leafcutter/hg38_splicing_QTLs/notx/NT_perind.counts",
            sep = " ", col.names = T, row.names = F, quote = F)
write.table(il1b, "X:/leafcutter/hg38_splicing_QTLs/il1b/IL1B_perind.counts",
            sep = " ", col.names = T, row.names = F, quote = F)

####Putting_the_chromosomes_together_for_matrixeQTL
##NT##
nt_list<-c("chr10/NT_perind.counts.gz.phen_chr10", "chr14/NT_perind.counts.gz.phen_chr14",  
           "chr18/NT_perind.counts.gz.phen_chr18", "chr21/NT_perind.counts.gz.phen_chr21",
           "chr4/NT_perind.counts.gz.phen_chr4", "chr8/NT_perind.counts.gz.phen_chr8",
           "chr11/NT_perind.counts.gz.phen_chr11", "chr15/NT_perind.counts.gz.phen_chr15",
           "chr19/NT_perind.counts.gz.phen_chr19", "chr22/NT_perind.counts.gz.phen_chr22",
           "chr5/NT_perind.counts.gz.phen_chr5", "chr9/NT_perind.counts.gz.phen_chr9",
           "chr12/NT_perind.counts.gz.phen_chr12", "chr16/NT_perind.counts.gz.phen_chr16",
           "chr1/NT_perind.counts.gz.phen_chr1", "chr2/NT_perind.counts.gz.phen_chr2",
           "chr6/NT_perind.counts.gz.phen_chr6", "chr13/NT_perind.counts.gz.phen_chr13",
           "chr17/NT_perind.counts.gz.phen_chr17", "chr20/NT_perind.counts.gz.phen_chr20",
           "chr3/NT_perind.counts.gz.phen_chr3", "chr7/NT_perind.counts.gz.phen_chr7")

anno_nt<-as.data.frame(matrix(nrow = 0, ncol = 4))
names(anno_nt)<-c("ID","chr","start","end")
tog_nt<-as.data.frame(matrix(nrow = 0, ncol = 53))
names(tog_nt)<-names(nt)[2:54]
for(the.file in nt_list){
  to.add<-read.delim(paste("X:/leafcutter/hg38_splicing_QTLs/notx/",the.file, sep = ""), 
                     sep = "\t", header = T, as.is = T, stringsAsFactors = F)
  
  anno.to.add<-to.add[,c(4,1,2,3)]
  names(anno.to.add)<-c("ID","chr","start","end")
  
  row.names(to.add)<-to.add$ID
  to.add<-to.add[,-c(1:4)]
  
  tog_nt<-rbind(tog_nt, to.add)
  anno_nt<-rbind(anno_nt, anno.to.add)
  
  print(the.file)
}

##IL1B##
il1b_list<-c("chr10/IL1B_perind.counts.gz.phen_chr10", "chr14/IL1B_perind.counts.gz.phen_chr14",  
             "chr18/IL1B_perind.counts.gz.phen_chr18", "chr21/IL1B_perind.counts.gz.phen_chr21",  
             "chr4/IL1B_perind.counts.gz.phen_chr4", "chr8/IL1B_perind.counts.gz.phen_chr8",
             "chr11/IL1B_perind.counts.gz.phen_chr11", "chr15/IL1B_perind.counts.gz.phen_chr15",  
             "chr19/IL1B_perind.counts.gz.phen_chr19", "chr22/IL1B_perind.counts.gz.phen_chr22",  
             "chr5/IL1B_perind.counts.gz.phen_chr5", "chr9/IL1B_perind.counts.gz.phen_chr9",
             "chr12/IL1B_perind.counts.gz.phen_chr12", "chr16/IL1B_perind.counts.gz.phen_chr16",  
             "chr1/IL1B_perind.counts.gz.phen_chr1", "chr2/IL1B_perind.counts.gz.phen_chr2",  
             "chr6/IL1B_perind.counts.gz.phen_chr6", "chr13/IL1B_perind.counts.gz.phen_chr13", 
             "chr17/IL1B_perind.counts.gz.phen_chr17", "chr20/IL1B_perind.counts.gz.phen_chr20",  
             "chr3/IL1B_perind.counts.gz.phen_chr3", "chr7/IL1B_perind.counts.gz.phen_chr7")

anno_il1b<-as.data.frame(matrix(nrow = 0, ncol = 4))
names(anno_il1b)<-c("ID","chr","start","end")
tog_il1b<-as.data.frame(matrix(nrow = 0, ncol = 53))
names(tog_il1b)<-names(il1b)[2:54]
for(the.file in il1b_list){
  to.add<-read.delim(paste("X:/leafcutter/hg38_splicing_QTLs/il1b/",the.file, sep = ""), 
                     sep = "\t", header = T, as.is = T, stringsAsFactors = F)
  
  anno.to.add<-to.add[,c(4,1,2,3)]
  names(anno.to.add)<-c("ID","chr","start","end")
  
  row.names(to.add)<-to.add$ID
  to.add<-to.add[,-c(1:4)]
  
  tog_il1b<-rbind(tog_il1b, to.add)
  anno_il1b<-rbind(anno_il1b, anno.to.add)
  
  print(the.file)
}

write.table(tog_nt,"X:/leafcutter/hg38_splicing_QTLs/notx/chrall_NT_splicing_phenotype.ratio.expression.txt", 
            sep = "\t", col.names = NA, row.names = T, quote = F)
write.table(tog_il1b,"X:/leafcutter/hg38_splicing_QTLs/il1b/chrall_IL1B_splicing_phenotype.ratio.expression.txt", 
            sep = "\t", col.names = NA, row.names = T, quote = F)
write.table(anno_nt,"X:/leafcutter/hg38_splicing_QTLs/notx/ANNOTATION_NT_splicing_phenotype.ratio.expression.txt", 
            sep = "\t", col.names = T, row.names = F, quote = F)
write.table(anno_il1b,"X:/leafcutter/hg38_splicing_QTLs/il1b/ANNOTATION_IL1B_splicing_phenotype.ratio.expression.txt", 
            sep = "\t", col.names = T, row.names = F, quote = F)



####Genotype files####
geno<-read.delim("W:/vari-gene-final/hg19/QTL_Testing_files/rna/matrixeqtl/notx/GTEx_gene_filter/Genotype_file_correct_donors_ordered_notx_filtered.txt.gz",
                 sep = "\t", header = T, as.is = T, stringsAsFactors = F)
geno_anno<-read.delim("W:/vari-gene-final/hg19/QTL_Testing_files/rna/matrixeqtl/notx/GTEx_gene_filter/Anno_file.txt.gz",
                      sep = "\t", header = T, as.is = T, stringsAsFactors = F)

nt_donor_names<-unlist(lapply(strsplit(names(tog_nt), split = "_"),"[[",1))
il1b_donor_names<-unlist(lapply(strsplit(names(tog_il1b), split = "_"),"[[",1))

names(geno)

geno_nt<-geno[,match(nt_donor_names, names(geno))]
geno_il1b<-geno[,match(il1b_donor_names, names(geno))]

#Make a bedfile to liftover
geno_bedfile<-as.data.frame(matrix(nrow = nrow(geno), ncol = 4))
names(geno_bedfile)<-c("Chrom","Start","End","ID")
geno_bedfile$Chrom<-geno_anno$chr
geno_bedfile$Start<-as.integer(geno_anno$pos)
geno_bedfile$End<-as.integer(geno_anno$pos+1)
geno_bedfile$ID<-geno_anno$id
write.table(geno_bedfile, "X:/leafcutter/geno_bedfile_for_liftover.bed", sep = "\t", col.names = F, row.names = F, quote = F)


lifted<-read.delim("X:/leafcutter/geno_bedfile_hg38-uplifted.bed",
                   sep = "\t", header = F, as.is = T, stringsAsFactors = F)

geno_nt_filt<-geno_nt[which(row.names(geno_nt)%in%lifted$V4),]
geno_il1b_filt<-geno_il1b[which(row.names(geno_il1b)%in%lifted$V4),]


write.table(geno_nt_filt, "X:/leafcutter/hg38_splicing_QTLs/notx/genotypes_NT_ordered.txt",
            sep = "\t", col.names = NA, row.names = T, quote = F)
write.table(geno_il1b_filt, "X:/leafcutter/hg38_splicing_QTLs/il1b/genotypes_IL1B_ordered.txt",
            sep = "\t", col.names = NA, row.names = T, quote = F)

geno_anno_toPrint<-geno_anno[which(geno_anno$id%in%lifted$V4),]
geno_anno_toPrint$pos<-lifted$V2[match(geno_anno_toPrint$id, lifted$V4)]
write.table(geno_anno_toPrint,"X:/leafcutter/hg38_splicing_QTLs/notx/Genotype_annotations_notx.txt",
            sep = "\t", col.names = T, row.names = F, quote = F)

