#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==6)

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
# base.dir = find.package('MatrixEQTL');
base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
#SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
#snps_location_file_name = paste(base.dir, "/data/snpsloc.txt", sep="");
SNP_file_name = args[1]
snps_location_file_name = args[2]

# Gene expression file name
#expression_file_name = paste(base.dir, "/data/GE.txt", sep="");
#gene_location_file_name = paste(base.dir, "/data/geneloc.txt", sep="");

expression_file_name =  args[3]
gene_location_file_name = args[4]
  
# Covariates file name
# Set to character() for no covariates
#covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");

covariates_file_name = args[5]

#austin's output file name (because he doesn't understand tempfile())
outputFile = args[6]

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 1e-6;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 100000;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
thing = (me$cis$eqtls)
otherBiggerThing = (me$trans$eqtls)
#cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)

write.table(thing,paste(outputFile,'Cis.txt',sep = ''),sep = '\t',quote = F)
write.table(otherBiggerThing,paste(outputFile,'Trans.txt',sep=''),sep='\t',quote=F)
## Plot the Q-Q plot of local and distant p-values

pdf(paste(outputFile,'.pdf',sep = ''))
plot(me)
dev.off()

cat('All done!')
