covar_nt<-read.delim("/path/notx/NT_perind.counts.gz.PCs",
                  sep = "\t", header = T, as.is = T, stringsAsFactors = F)
covar2_nt<-read.delim("/path/bamfiles_and_confounders_final_set.txt",
                      sep = "\t", header = F, as.is = T, stringsAsFactors = F)

notx_cov_file<-as.data.frame(matrix(nrow = 18, ncol = 53))
names(notx_cov_file)<-names(covar_nt)[2:54]

donor.names<-unlist(lapply(strsplit(names(notx_cov_file), split = "_"), "[[", 1))

names(notx_cov_file)<-donor.names

for(i in 1:ncol(notx_cov_file)){
  notx_cov_file[1:8,i]<-t(covar2_nt[grep(names(notx_cov_file)[i], covar2_nt$V1)[1],3:10])
  notx_cov_file[9:18,i]<-covar_nt[,grep(names(notx_cov_file)[i], names(covar_nt))]
  
  print(i)
}

row.names(notx_cov_file)<-paste0("factor_",row.names(notx_cov_file))

write.table(notx_cov_file, "/path/notx/covariates_notx_for_qtl_analysis.txt",
            sep = "\t", col.names = NA, row.names = T, quote = F)

#il1b
covar_il1b<-read.delim("/path/il1b/IL1B_perind.counts.gz.PCs",
                     sep = "\t", header = T, as.is = T, stringsAsFactors = F)
covar2_il1b<-read.delim("/path/bamfiles_and_confounders_final_set.txt",
                      sep = "\t", header = F, as.is = T, stringsAsFactors = F)

il1b_cov_file<-as.data.frame(matrix(nrow = 18, ncol = 53))
names(il1b_cov_file)<-names(covar_il1b)[2:54]

donor.names<-unlist(lapply(strsplit(names(il1b_cov_file), split = "_"), "[[", 1))

names(il1b_cov_file)<-donor.names

for(i in 1:ncol(il1b_cov_file)){
  il1b_cov_file[1:8,i]<-t(covar2_il1b[grep(names(il1b_cov_file)[i], covar2_il1b$V1)[1],3:10])
  il1b_cov_file[9:18,i]<-covar_il1b[,grep(names(il1b_cov_file)[i], names(covar_il1b))]
  
  print(i)
}

row.names(il1b_cov_file)<-paste0("factor_",row.names(il1b_cov_file))

write.table(il1b_cov_file, "/path/il1b/covariates_il1b_for_qtl_analysis.txt",
            sep = "\t", col.names = NA, row.names = T, quote = F)

