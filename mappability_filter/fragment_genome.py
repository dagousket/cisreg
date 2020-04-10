#!/g/software/bin/python2.7
"""Will transform a fasta file into all possible fragments
Remember -- the end(s) of your sequences will be missing fragment starts on each end.
"""

import sys
from Bio import SeqIO
import argparse
import gzip

asFastq = True #if True, will write as fastq rather than fasta


# f_in = sys.stdin
# kmer_len = int(sys.argv[1])
# f_out = sys.stdout


def get_fragments(kmer_len, f_out, asFastq, chromosome, sequence, randName):
    for i in range(len(sequence) - kmer_len):
        if asFastq is True:
            #qualString = '~' * kmer_len
            qualString = 'F' * kmer_len
            outLine = '@' + chromosome + randName + ':' + str(i+1) + '-' + str(i+kmer_len) +  '\n' + sequence[i:i+kmer_len].tostring() + '\n'
            outLine = outLine + '+' + chromosome + randName + ':' + str(i+1) + '-' + str(i+kmer_len) +  '\n' + qualString+ '\n'
            f_out.write(outLine)
        else:
            outLine = '>' + chromosome + randName + ':' + str(i+1) + '-' + str(i+kmer_len) + '\n' + sequence[i:i+kmer_len].tostring() + '\n'
            f_out.write(outLine)

def main(args):
    try:
        args.fragment_size = int(args.fragment_size)
    except:
        raise ValueError('It does not appear that fragment size is an integer')
    if args.inFasta == 'stdin':
        f_in = sys.stdin
    else:
        f_in = open(args.inFasta, 'rU')
    if args.outFile == 'stdout':
        f_out = sys.stdout
    else:
        if args.outFile.split('.')[-1] == 'gz':
            f_out = gzip.open(args.outFile,'w')
        else:
            f_out = open(args.outFile, 'w')
    if args.randomName is True:
        import string
        import random
    else:
        randName = ''
    for record in SeqIO.parse(f_in,"fasta"):
        chromosome =  record.id
        if args.randomName is True:
            randName = '_' + ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(8))
        #sequence = record.seq.tostring()
        sequence = record.seq
        get_fragments(args.fragment_size, f_out, args.fastq, chromosome, sequence, randName)
        if args.both_strands is True:
            rev_seq = sequence.reverse_complement()
            get_fragments(args.fragment_size, f_out, args.fastq, chromosome + '_rev', rev_seq, randName)



parser = argparse.ArgumentParser(description='generates all possible reads from a given fasta file')
parser.add_argument('--inFasta', '-i', required = False, default = 'stdin' ,metavar = 'infile', help = 'Fasta formated input file. Default is stdin')
parser.add_argument('--outFile', '-o', required = False, default = 'stdout', metavar = 'outfile', help = 'Default is stdout')
parser.add_argument('--fastq', required = False, default = False, action = 'store_true', help = 'If set, program will output reads in Fastq format rather than Fasta')
parser.add_argument('--both_strands', '-s', required = False, default = False, action = 'store_true', help = 'If set, program will generate reads from both strands')
parser.add_argument('--fragment_size', '-f', required = True, metavar = 'fragmentsize', help = 'The size of the fragments')
parser.add_argument('--randomName', '-r', required = False, default = False, action = 'store_true' , help = 'If specified, a random name will be added to the sequence')


if len(sys.argv) == 1:
    parser.parse_args(['--h'])
else:
    args = parser.parse_args()
    main(args)

