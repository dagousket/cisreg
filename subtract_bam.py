#!/usr/bin/env python3
"""
Equivalent of bedtools subtract but with a correct handling of spliced reads for bam files. Will trim the parts of the reads that fall inside a 'bad' region.
"""

import pysam
import argparse
import operator
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

def cigar_tuples_from_aligned_pairs(aligned_pairs, remove_non_aligned_ref_start = False):
    """Generates a simple set of cigar tuples from a set of aligned_pairs.
    No distinction is made between deletions and splices, so you arejust going to end up with Ns here."""
    if remove_non_aligned_ref_start is True:
        minRef = [k for k in aligned_pairs if k[0] is not None and k[1] is not None][0][1]
        aligned_pairs = [k for k in aligned_pairs if k[1] is None or k[1] >= minRef]
    if len(aligned_pairs) == 0:
        return []
    if aligned_pairs[0][0] is None or aligned_pairs[-1][0] is None:
        print(aligned_pairs)
        raise ValueError('Improperly formatted set of aligned_pairs!')
    cigartuples = []
    if aligned_pairs[0][1] is None:
        alignmentType = 1 #an insertion...we will edit first and last insertions into soft-clips at the end
    else:
        alignmentType = 0 #alignment
    alignmentLen = 0
    for myPair in aligned_pairs:
        if myPair[0] is not None and myPair[1] is not None:
            currentType = 0
        elif myPair[0] is None and myPair[1] is not None:
            currentType = 3
        elif myPair[0] is not None and myPair[1] is None:
            currentType = 1
        else:
            print(myPair)
            raise ValueError('You should not have this type of alignment pair')
        if currentType != alignmentType:
            cigartuples.append((alignmentType, alignmentLen))
            alignmentType = currentType
            alignmentLen = 0
        alignmentLen = alignmentLen + 1
    cigartuples.append((alignmentType, alignmentLen))
    if cigartuples[0][0] == 1:
        cigartuples[0] = (4, cigartuples[0][1])
    if cigartuples[-1][0] == 1:
        cigartuples[-1] = (4, cigartuples[-1][1])
    return cigartuples


def collect(l, index):
   return map(operator.itemgetter(index), l)


def insert_pair(myPair, aligned_pairs):
    mySet = collect(aligned_pairs, 0)
    myKey = min([k for k in mySet if k], key = lambda x:abs(x-myPair[0]))
    if myPair[0] < myKey:
        aligned_pairs.insert(mySet.index(myKey), myPair)
    else:
        aligned_pairs.insert(mySet.index(myKey) + 1, myPair)


def merge_aligned_pairs(pairs1, pairs2):
    """This is a simple module designed to return a set of aligned pairs corresponding to read bases that are aligned in BOTH pairs1 and pairs2.
    Importantly, the same base from query1 and query2 must align to the same target base.
    """
    t_range1 = [k[1] for k in pairs1 if k[1] is not None and k[0] is not None]
    t_range2 = [k[1] for k in pairs2 if k[1] is not None and k[0] is not None]
    t_min = min(min(t_range1), min(t_range2))
    t_max = max(max(t_range1), max(t_range2))
    #ideally, these next two values are equal, though hard trimming could screw something up
    q_min = min(pairs1[0][0], pairs2[0][0])
    q_max = max(pairs1[-1][0], pairs2[-1][0])
    if q_min is None or q_max is None:
        raise ValueError('These aligned pairs end or begin with None!')
    pairs1_dict = dict(pairs1)
    pairs2_dict = dict(pairs2)
    aligned_pairs = []
    aligned_q = set()
    aligned_t = set()
    for i in range(q_min, q_max + 1):
        myKey = i
        if myKey in pairs1_dict and myKey in pairs2_dict:
            if pairs1_dict[myKey] is not None and pairs1_dict[myKey] == pairs2_dict[myKey]:
                aligned_pairs.append((i, pairs1_dict[myKey]))
                aligned_q.add(i)
                aligned_t.add(pairs1_dict[myKey])
    for i in range(t_min,t_max):
        if i not in aligned_t: aligned_pairs.append((None, i))
    aligned_pairs.sort(key=operator.itemgetter(1))
    if len(aligned_q) == 0:
        return []
    start_pairs, end_pairs = [], []
    qStart = min(aligned_q)
    qEnd = max(aligned_q)
    for i in range(q_min,q_max + 1):
        if i not in aligned_q:
            if i < qStart: start_pairs.append((i, None))
            elif i > qEnd : end_pairs.append((i,None))
            else: insert_pair((i, None), aligned_pairs)
    return start_pairs + aligned_pairs + end_pairs


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
			preCigar = cigar_tuples_from_aligned_pairs(newPairs)
			if sum([k[1] for k in preCigar if k[0] == 0]) != 0:
				#Correct the aligned_pairs to add a deletion event together with the trimming previously added
				CorPairs = merge_aligned_pairs(Pairs,newPairs)
				#Write the read to a new file once its class values has beed updated
				valid_cig = cigar_tuples_from_aligned_pairs(CorPairs, remove_non_aligned_ref_start = True)
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