#!/usr/bin/env python
"""
Convert a BED file to a minimal GTF file
"""

import random, sys
import argparse

class Feature(object):
    """
    Class to represent a feature
    """
    def __init__(self, chromosome, start, end, name):
        self.name = name
        self.chrom = chromosome
        self.start = int(start)
        self.end = int(end)
        self.strand = "+"
    
def makeSeq(name, len):
    """
    Make a fake sequence
    """
    seq = ""
    for i in range(0, len):
        seq += random.choice("ACGT")
    return f">{name}\n{seq}"
if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Convert a BED file to a minimal GTF file")
    args.add_argument("-i", "--bed", help="BED file", required=True)
    
    args.add_argument("-o", "--fasta", help="FASTA file") 
    args.add_argument("--min-len", help="Minimum length of feature [default: %(default)s]", default=0, type=int)
    args = args.parse_args()


    if args.fasta is None:
        out = sys.stdout
    else:
        out = open(args.fasta, 'w')

   
    chrLen = {}
    with open(args.bed, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            chromosome, start, end, name = line.split()

            if not chromosome in chrLen:
                chrLen[chromosome] = max(int(end), int(start), args.min_len)
            else:
                chrLen[chromosome] = max(int(end), int(start), chrLen[chromosome])

    for chrname in chrLen:
        print(makeSeq(chrname, chrLen[chrname]), file=out)