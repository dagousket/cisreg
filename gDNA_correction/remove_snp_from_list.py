#!/usr/bin/env python2.7
"""
Remove from a SNP file the SNP selected as bad
"""

import sys
import argparse
import gzip

parser = argparse.ArgumentParser(description='Add two count files')
parser.add_argument('--snp_file', '-f', required = True, help = "SNP file to filter")
parser.add_argument('--snp_list', '-s', required = True, help = "SNP selected as bad")
parser.add_argument('--out_file', '-o', required = True, help = 'Output file name.')


def open_and_test_gzip(myFile, readType):
    if myFile.split('.')[-1] == 'gz':
        return gzip.open(myFile, readType)
    else:
        return open(myFile, readType)


def main(args):
    with open_and_test_gzip(args.out_file,'wb') as outfile :
        bad_snp = set(line.rstrip('\r\n').strip() for line in open_and_test_gzip(args.snp_list,'rb'))
        with open_and_test_gzip(args.snp_file,'rb') as infile :
            header = next(infile)
            outfile.write(header)
            for snp in infile :
                atoms = snp.rstrip('\r\n').split('\t')
                if str(atoms[0])+':'+str(atoms[1])+'-'+str(atoms[2]) not in bad_snp :
                    if int(atoms[2])-int(atoms[1]) == 1:
                        outfile.write(snp)
                    else:
                        print(snp)
                        print("Above SNP seems to not be one basepair long")


if len(sys.argv) == 1:
    parser.parse_args(['--h'])
else:
    args = parser.parse_args()
    main(args)
