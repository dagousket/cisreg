# CisReg

This repository contains scripts used to perform analysis on F1 allelic ratio data. Scritps are accompanied by toy input file and description of their aim.

The repository is divided into four analysis types :

## gDNA_correction

This repository contains scripts to help detect genotyping error based on allele-specific count on genomic DNA data.

## genome_ranges_toolkit

This repository contains scripts used to create the region overlap assignement for four non-coding layers (ATAC, RNA, ChIP).

## mappability_filter

This repository contains scripts to generate mappability filters (synthetic and genomic) and to apply them to a BAM file in order to trim the reads.

## partial_correlation 

This repository contains scripts to perform a partial correlation analysis, using the package GeneNet and also "from scratch" by computing the residuals from a linear regression.