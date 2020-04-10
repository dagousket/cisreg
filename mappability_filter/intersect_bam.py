#!/usr/bin/env python
"""
Equivalent of bedtools intersect but with a correct handling of split argument for bam files. Will trim the parts of the reads that fall inside a 'bad' region.
"""

import pysam
import bam_tools
import psl_tools
import argparse
import sys

parser = argparse.ArgumentParser(description='A script for intersecting a BAM file with a bed file while considering the splicing or blocked information.')
parser.add_argument('--in_bam', '-i', required = True, help = 'The path to the BAM file to filter.')
parser.add_argument('--bed_filter', '-b', required = True, help = 'The path to the BED file used as a filter, should contain the regions you want to discard.')
parser.add_argument('--out_bam', '-o', required = True, help = 'The path to the resulting filtered BAM file.')
parser.add_argument('--out_bam_error', '-e', required = False, help = 'The path to the BAM file storing problematic reads (raise KeyError or ValueError, for example the unmapped ones).')
parser.add_argument('--out_bam_fail', '-f', required = False, help = 'The path to the BAM file storing the reads which are entirely discarded by the filter.')

def bed_as_set(bed_filter):
	#Creates a dictionnary for all region to filter out, key = chromosome name, value = set of all (yup, all) the bases corresponding the bad region
	bedDict = {}
	with open(bed_filter) as bed:
		for line in bed :
			line_args = line.rstrip('\r\n').split('\t')
			key = line_args[0]
			coor = set(range(int(line_args[1]),int(line_args[2])))
			if key in bedDict :
				bedDict[key].update(coor)
			else :
				bedDict[key]=coor
	bed.close()
	return bedDict


def main(args):
	#Open all the required files using pysam
	inSam = pysam.Samfile(args.in_bam, 'rb')
	samHeader = inSam.header
	outSam = pysam.Samfile(args.out_bam, 'wb', header = samHeader)
	if args.out_bam_error:
		outSamError = pysam.Samfile(args.out_bam_error, 'wb', header = samHeader)
	if args.out_bam_fail:
		outSamFail = pysam.Samfile(args.out_bam_fail, 'wb', header = samHeader)
	inBed = bed_as_set(args.bed_filter)

	for myRead in inSam :
		try :
			#Get chromosome name for the read
			readchr = inSam.getrname(myRead.tid)
			#Get the base which intersect between the read and the filter
			if readchr in inBed :
				BaseSet = inBed[readchr]
			else:
				BaseSet = set()
			Pairs = myRead.aligned_pairs
			BasePair = set(k[1] for k in Pairs)
			Bad = BaseSet.intersection(BasePair)
			#Based on the previous set, modify the aligned_pairs class accordingly (but first only convert it to trimming)
			newPairs=[(Pairs[k][0],None) if Pairs[k][1] in Bad and Pairs[k][0] != None else Pairs[k] for k in range(0,len(Pairs))]
			#Do not take into account this pair if it doesn't have anything mapping
			#Match and insertion (M=0,I=1,D=2,N=3,S=4)
			preCigar = bam_tools.cigar_tuples_from_aligned_pairs(newPairs)
			if sum([k[1] for k in preCigar if k[0] == 0]) != 0:
				#Correct the aligned_pairs to add a deletion event together with the trimming previously added
				CorPairs = psl_tools.merge_aligned_pairs(Pairs,newPairs)
				#Write the read to a new file once its class values has beed updated
				valid_cig = bam_tools.cigar_tuples_from_aligned_pairs(CorPairs, remove_non_aligned_ref_start = True)
				valid_pos = min([k[1] for k in CorPairs if k[0] is not None and k[1] is not None])
				myRead.cigar = valid_cig
				myRead.pos = valid_pos
				outSam.write(myRead)
			else:
				#Corresponds to reads that do not have any matching bases in the good region
				if args.out_bam_fail:
					outSamFail.write(myRead)
				else:
					pass
		except (ValueError,KeyError) :
			#Corresponds to the reads that raise an error (strange chromosome or unmappable read = empty cigar, usually)
			if args.out_bam_error:
				outSamError.write(myRead)
			else:
				pass

	inSam.close()
	outSam.close()
	if args.out_bam_error:
		outSamError.close()
	if args.out_bam_fail:
		outSamFail.close()



if len(sys.argv) == 1:
    parser.parse_args(['--h'])
else:
    args = parser.parse_args()
    main(args)