#!/bin/bash
# This script will (for a given genome) :
# create bed12 from a given fasta
# create transcript fasta from bed12
# create chromosome-splitted fasta from transcript
# create reads for each chromosome file on cluster

if [ $# -eq 0 ] || [ "$1" == "-h" ] ; then
    echo 'usage: simulate_transcript_reads.sh [-h] genome gffFile fastaFile directory

Creates simulated fastq files of transcript reads from a given GFF and FASTA file. Detail of the command used in log file.
Should run on cluster.

Arguments :
   -h			show this message and exit
   $1		Name of the genome (used to name files)
   $2		Name of the GFF file
   $3		Name of the FASTA file
   $4		Path to directory where the files will be created'
    exit 0
fi


#Arguments
# genome name
genomeN=$1
# Path gff file
gffFile=$2
# Path to fasta file
fastaFile=$3
# Directory to store results (fastq.gz)
directory=$4


python_bin=/g/furlong/project/36_cisRegVar/src/python/tools
sh_bin=/g/furlong/project/36_cisRegVar/src/sh
bedtools_bin=/g/software/bin/bedtools


echo -e "Log file for the commandline :\nsimulate_transcript_reads.sh $1 $2 $3 $4\n" > $directory/$genomeN.log

echo 'Running gff_to_bed...'
echo -e '--------------\ngff_to_bed.py\n--------------' >> $directory/$genomeN.log

$python_bin/gff_to_bed.py -f $gffFile -b12 -nb6 -id -p $directory -n $genomeN -v -d >> $directory/$genomeN.log

echo 'Running bedtools getfasta...'
echo -e "-----------------\nbedtools getfasta\n-----------------
bedtools getfasta -fi $fastaFile -bed $directory/$genomeN.bed12 -fo $directory$genomeN.tr.fa -split -s\n" >> $directory/$genomeN.log

$bedtools_bin getfasta -fi $fastaFile -bed $directory/$genomeN.bed12 -fo $directory/$genomeN.tr.fa -split -s >> $directory/$genomeN.log

echo 'Running sort_chromosome_fasta...'
echo -e "------------------------\nsort_chromosome_fasta.sh\n------------------------
$sh_bin/sort_chromosome_fasta.sh $directory/$genomeN.tr.fa $directory\n" >> $directory/$genomeN.log

$sh_bin/sort_chromosome_fasta.sh $genomeN $directory/$genomeN.tr.fa $directory

echo "Running fragment_genome..."
echo -e "------------------\nfragment_genome.py\n------------------" >> $directory/$genomeN.log
for f in `cat $directory/$genomeN.chromosome.txt`; do
	echo -e "\n## $f ##" >> $directory/$genomeN.log
	echo -e "$python_bin/fragment_genome.py -i $directory/$f.fa -o $directory/$f.fastq.gz --fastq -f 100 -r" >> $directory/$genomeN.log
	$python_bin/fragment_genome.py -i $directory/$f.$genomeN.fa -o $directory/$f.$genomeN.fastq.gz --fastq -f 100 -r >> $directory/$genomeN.log
done

echo -e "Once the jobs on the cluster are over, run the following command to remove temporary files :\nsrc/sh/simulate_transcript_reads_cleaner.sh $directory">>$directory/$genomeN.log



