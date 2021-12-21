#!/usr/bin/env python
"""
Convert a BED file to a minimal GTF file
"""

import os, sys
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
    
   
if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Convert a BED file to a minimal GTF file")
    args.add_argument("-i", "--bed", help="BED file", required=True)
    args.add_argument("-o", "--gtf", help="GTF file")
    args.add_argument("--feat", help="Feature type [default: %(default)s]", default="exon")
    args.add_argument("--id", help="ID attribute [default: %(default)s]", default="gene_id")
    args.add_argument("--min-len", help="Minimum length of feature [default: %(default)s]", default=0, type=int)
    args.add_argument("--add-len", help="Add length of feature [default: %(default)s]", default=0, type=int)
    args = args.parse_args()


    if args.gtf is None:
        out = sys.stdout
    else:
        out = open(args.gtf, 'w')

    header = "##gtf-version 3"
    features = []
    chrLen = {}
    print(header, file=out)

    
    with open(args.bed, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            chromosome, start, end, name = line.split()
            features.append(Feature(chromosome, start, end, name))
            if not chromosome in chrLen:
                chrLen[chromosome] = max(int(end) + args.add_len, int(start) + args.add_len, args.min_len)
            else:
                chrLen[chromosome] = max(int(end) + args.add_len, int(start) + args.add_len, chrLen[chromosome])

    for chrname in chrLen:
        print(f"##sequence-region {chrname} 1 {chrLen[chrname]}", file=out)

    for feat in features:
        #NC_001422.1     Prodigal:002006 gene    51      221     .       +       0       gene_id "nbis-gene-1";
        print(f"{feat.chrom}\tbamtocov\t{args.feat}\t{feat.start}\t{feat.end}\t.\t{feat.strand}\t.\t{args.id} \"{feat.name}\"", file=out)