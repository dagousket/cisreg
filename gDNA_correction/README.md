# gDNA correction

This repository contains scripts to help detect genotyping error based on allele-specific count on genomic DNA data.

## get_empirircal_H0_ratio

This script compute the empirical average allelic ratio for autosomes and chrX

```
Usage: get_empirircal_H0_ratio.R [options]
A simple R script to sort the count files by decreasing total coverage and computing the mean ASE value for chrX and autosomes independantly on the N first best covered features.

example run : Rscript --vanilla  get_empirircal_H0_ratio.R -i snp_list.txt -o output.tab -x snp_in_chrX.tab -n 5 --null_hypothesis FALSE

Options:
	-i CHARACTER, --input=CHARACTER
		list of count files to use in a txt file

	-o CHARACTER, --out=CHARACTER
		file name to store the output as a tab-delimited file [default= out.txt]

	-x CHARACTER, --chrx=CHARACTER
		list of features in input files on chrX

	-n CHARACTER, --number_of_best=CHARACTER
		Number of best coverage features to take to compute mean ASE

	--null_hypothesis
		If set, will always write 0.5 as the null hypothesis for the autosomes

	-h, --help
		Show this help message and exit
```



## test_imbalance.R

This script test for imbalance with a binomial test for all SNPs.

```
Usage: test_imbalance.R [options]
An R script to perform a binomial test on count file and extract the imbalanced features.

example run : Rscript --vanilla test_imbalance.R -i snp_list.txt -o test_output -x snp_in_chrX.tab --autosome_treshold 20 -n output.tab 

Options:
	-i CHARACTER, --input=CHARACTER
		list of count files to use in a txt file

	-o CHARACTER, --out=CHARACTER
		file prefix to store the output as a tab-delimited file (2 files per sample) [default= binom_test]

	-x CHARACTER, --chrx=CHARACTER
		list of features in input files on chrX

	--chrx_treshold=CHARACTER
		Minimum total coverage required for feature on chrx

	--autosome_treshold=CHARACTER
		Minimum total coverage required for feature on autosomes

	-n CHARACTER, --null_hypothesis_mean=CHARACTER
		Table of mean value of ASE under the H0 hypothesis

	-h, --help
		Show this help message and exit
```

## remove_snp_from_list.py

This script efficiently remove SNPs from a feature file listing all the SNPs. This is more efficient than grep -f in this case, as we need to remove a large number of items from the file.

```
usage: remove_snp_from_list.py [-h] --snp_file SNP_FILE --snp_list SNP_LIST
                               --out_file OUT_FILE

example run : ../remove_snp_from_list.py -f snp_feature_file.tab -s output_imbalance_testsnp_count_imbalanced_SNP_list.tab -o filtered_feature_snp_file.tab

Remove snps from a snp feature list 

optional arguments:
  -h, --help            show this help message and exit
  --snp_file SNP_FILE, -f SNP_FILE
                        SNP file to filter
  --snp_list SNP_LIST, -s SNP_LIST
                        SNP selected as bad
  --out_file OUT_FILE, -o OUT_FILE
                        Output file name.
```

