#!/usr/bin/env python3
"""
A tool to generate simulated BAM files
starting from a list of coordinates
"""
import os, sys
import argparse
import pysam


def loadSeqsFromFasta(fastaFile):
    """
    Load sequences from a fasta file
    """
    seqs = {}
    with open(fastaFile, 'r') as f:
        for line in f:
            if line[0] == '>':
                seqName = line[1:].strip()
                seqs[seqName] = ''
            else:
                seqs[seqName] += line.strip()
    return seqs

def makeHeader(seqs):
    """
    Make a header for the BAM file
    """
    header = { 'HD': {'VN': '1.0'},
            'SQ': [] }

    for seqName in seqs:
        # {'LN': 1575, 'SN': 'chr1'},
        header['SQ'].append({'LN': len(seqs[seqName]), 'SN': seqName})
    
    return header
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="A tool to generate simulated BAM files")
    #parser.add_argument("-i", "--input", help="Input file with coordinates", required=True)
    parser.add_argument("-f", "--fasta", help="Input fasta file", required=True)
    parser.add_argument("-l", "--length", help="Read length", default=100)
    parser.add_argument("-o", "--output", help="Output BAM file", required=True)
    parser.add_argument("-s", "--seed", help="Random seed", default=2)
    parser.add_argument("--verbose", help="Verbose mode", action="store_true")
    opts = parser.parse_args()

    seqs = loadSeqsFromFasta(opts.fasta)
    header = makeHeader(seqs)

    with pysam.AlignmentFile(opts.output, "wb", header=header) as outf:
        for i, seq in enumerate(seqs):
            n = 0
            if opts.verbose:
                print("Chromosome: {}".format(seq))
            for pos in range(len(seqs[seq]) - int(opts.length)):
                n += 1
                a = pysam.AlignedSegment()
                a.query_name = "read_" + str(pos) + "_" + str(n)
    
                a.query_sequence = seqs[seq][pos:pos+int(opts.length)]
                a.flag = 0
                a.reference_id = i
                a.reference_start = pos
                a.mapping_quality = 20
                a.cigarstring = str(opts.length) + 'M'
                a.next_reference_id = 0
                a.next_reference_start = 0
                a.template_length = 0
                a.query_qualities = [20] * len(seqs[seq][pos:pos+int(opts.length)])
                outf.write(a)
             