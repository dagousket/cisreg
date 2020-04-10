#!/usr/bin/env Rscript
library("optparse")
library("tools")
library("reshape2")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="list of count files to use in a txt file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="binom_test", 
              help="file prefix to store the output as a tab-delimited file (2 files per sample) [default= %default]", metavar="character"),
  make_option(c("-x", "--chrx"), type="character", default = NULL,
              help="list of features in input files on chrX", metavar="character"),
  make_option(c("--chrx_treshold"), type="integer", default = 10, 
              help="Minimum total coverage required for feature on chrx", metavar="character"),
  make_option(c("--autosome_treshold"), type="integer", default = 10, 
              help="Minimum total coverage required for feature on autosomes", metavar="character"),
  make_option(c("-n", "--null_hypothesis_mean"), type="character",
              help="Table of mean value of ASE under the H0 hypothesis", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list, description = "An R script to perform a binomial test on count file and extract the imbalanced features.");
opt = parse_args(opt_parser);

##------------------------- Part 1 : Run binomial test

# Step 1 : set up environment and general variable you might want to modify
files = read.delim(opt$input,head=F, colClasses = "character")[,1]
h0_means = read.table(file = opt$null_hypothesis_mean, head = TRUE)
chrXFeat = read.delim(opt$chrx ,head=F, colClasses = "character")[,1]

# Step 2 : construct dataframe and apply binom.test (Warning : different test for chrX and autosomes SNPs)
for (i in seq(length(files),1)){
	print(paste('Processing',files[i]))
	df<-data.frame(read.table(files[i],head=TRUE, fill = TRUE, colClasses = c('character','numeric','numeric'))[,1:3])
	# Compute ASE and keep coverage information
	df$coverage = df[,2] + df[,3]
	df$mASE = df[,2]/(df[,3] + df[,2]) # Assume virginizer is the second column
	df = df[!is.na(df$mASE),]
	# Load H0 hyposthesis means
	mean_ase_autosomes = h0_means$mean_autosome[h0_means$sample == files[i]]
	mean_ase_chrx =  h0_means$mean_chrx[h0_means$sample == files[i]]
	df$chr[df$feature %in% chrXFeat] = as.numeric(mean_ase_chrx)
	df$chr[!(df$feature %in% chrXFeat)] = as.numeric(mean_ase_autosomes)
	print("Header of the file :")
	print(head(df))
	raw_pvals = apply(as.matrix(df[,-1]),1,function(c) binom.test(x = as.numeric(c[1]), n = as.numeric(c[1])+as.numeric(c[2]), p = as.numeric(c[5]), alternative = 'two.sided')$p.value)
	df$pval = raw_pvals
	df$adj_pval = p.adjust(as.numeric(raw_pvals), method = "fdr")
	write.table(df, file = paste(opt$out,file_path_sans_ext(basename(files[i])),'_binom_test.tab',sep = ''), quote = FALSE, sep = "\t", row.names = FALSE)

##------------------------- Part 2 : Save the list of imbalanced SNPs

# Step 1 : find the SNP files and extract imbalanced and balanced SNP names
	dfi = df[(df$coverage >= opt$chrx_treshold & df[,1] %in% chrXFeat) | (df$coverage >= opt$autosome_treshold & !(df[,1] %in% chrXFeat)),]
	dfi = dfi[dfi$adj_pval <= 0.05,]
	the_snps = as.character(dfi$feature[!is.na(dfi$adj_pval)])
	write.table(the_snps, file = paste(opt$out,file_path_sans_ext(basename(files[i])),'_imbalanced_SNP_list.tab',sep = ''), quote = FALSE, row.names = FALSE, col.names = FALSE)


	dfb = df[(df$coverage >= opt$chrx_treshold & df[,1] %in% chrXFeat) | (df$coverage >= opt$autosome_treshold & !(df[,1] %in% chrXFeat)),]
	dfb = dfb[dfb$adj_pval > 0.05,]
	the_snps = as.character(dfb$feature[!is.na(dfb$adj_pval)])
	write.table(the_snps, file = paste(opt$out,file_path_sans_ext(basename(files[i])),'_balanced_SNP_list.tab',sep = ''), quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Tip : the python script 'remove_snp_from_list' wonderfully takes care of the next step, i.e. using this file to create a filtered SNP file for Peaks/Genes/... without the imbalanced SNPs
q('no')