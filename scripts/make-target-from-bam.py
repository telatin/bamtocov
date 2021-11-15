#!/usr/bin/env python
"""
Generate a target covering the full chromosomes from a 
BAM file
"""

import os, sys
import pysam
import argparse

def readHeaderFromBAM(File):
    """
    Read the header from a BAM file
    """
    header = {}
    with pysam.Samfile(File, "rb") as bam:
        header["n_ref"] = bam.nreferences
        header["ref_names"] = bam.references
        header["ref_lengths"] = bam.lengths
    return header

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bam', help='BAM file')
    parser.add_argument('-o', '--out', help='Output file')
    opts = parser.parse_args()

    header = readHeaderFromBAM(opts.bam)

    if opts.out:
        out = open(opts.out, "w")
    else:
        out = sys.stdout

    for chromosome, length in zip(header["ref_names"], header["ref_lengths"]):
        print(chromosome, 0, length, chromosome, sep="\t", file=out)

