#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="list of count files to use in a txt file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="file name to store the output as a tab-delimited file [default= %default]", metavar="character"),
  make_option(c("-x", "--chrx"), type="character", 
              help="list of features in input files on chrX", metavar="character"),
  make_option(c("-n", "--number_of_best"), type="integer",
              help="Number of best coverage features to take to compute mean ASE", metavar="character"),
  make_option(c("--null_hypothesis"), action = 'store_true',
              help="If set, will always write 0.5 as the null hypothesis for the autosomes", metavar="logical")
); 

opt_parser = OptionParser(option_list=option_list, description = "A simple R script to sort the count files by decreasing total coverage and computing the mean ASE value for chrX and autosomes independantly on the N first best covered features.");
opt = parse_args(opt_parser);

# Step 1 : set up environment and general variable you might want to modify
files = read.delim(opt$input,head=F, colClasses = "character")[,1]

# Step 2 : compute mean ASE for the highest coverage features
ase_df = data.frame(feature = files)
for (j in 1:length(files)){
  print(paste('Processing',files[j]))
  df<-data.frame(read.table(files[j],head=TRUE, fill = TRUE, colClasses = c('character','numeric','numeric'))[,1:3])
  #order df and catch the ASE with the best coverage features
  df = df[order(df[,2] + df[,3], decreasing = TRUE),]
  df$mASE = df[,2]/(df[,3] + df[,2]) # Assume virginizer is the second column
  chrXFeat = read.delim(opt$chrx ,head=F, colClasses = "character")[,1]
  meanChrX = mean(df[df[,1] %in% chrXFeat,'mASE'][1:opt$number_of_best], na.rm = TRUE)
  meanAut = mean(df[!(df[,1] %in% chrXFeat),'mASE'][1:opt$number_of_best], na.rm = TRUE)
  ase_df$meanChrX[ase_df$feature == files[j]] = meanChrX
  if (opt$null_hypothesis == TRUE) {
    ase_df$meanAut[ase_df$feature == files[j]] = 0.5
  }  else {
    ase_df$meanAut[ase_df$feature == files[j]] = meanAut
  }

}

remove(df)

colnames(ase_df) = c('sample','mean_chrx','mean_autosome')
write.table(ase_df, file = opt$out, quote = FALSE, sep = '\t', dec = '.', row.names = FALSE)

q("no")