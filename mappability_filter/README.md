# mappability_filter

This repository contains scripts to generate mappability filters (synthetic and genomic) and to apply them to a BAM file in order to trim the reads. Contrary to Bedtools intersect, the resulting reads are not simply discarded when overlapping the filter. Instead, only the portion of the read overlapping the filter is trimmed.  Deletion are created in the read in cases where the filtered region is present in the midlle of the read.

## Generate filters

* Snakefile_genomic_filter : Create genomic filter, based on genomic DNA sequencing of parental lines.

* Snakefile_synthetic_dna_filter : Create synthetic filter for ATAC and ChIP data, based on simulated genomic reads generated with the script `fragment_genome.py`.

* Snakefile_synthetic_rna_filter : Create synthetic filter for RNA data, based on simulated transcriptomic reads generated with the script `simulate_trancript_reads.sh`.

* Snakefile_apply_filter : Apply any of the filter (BED file) to a BAM file by triming or deleting the portions of the read overlapping the reads. This workflow uses the script `intersect_bam.py` to apply the filter.



