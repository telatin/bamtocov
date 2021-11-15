#!/usr/bin/env python
"""
Print a fixed width Wiggle file from a coverage Bed graph stream.
"""

import os, sys
import argparse
import subprocess


def eprint(*args, **kwargs):
    """
    Print to stderr.
    """
    print(*args, file=sys.stderr, **kwargs)



def avgCoverage(bam, bin="bamtocov", out=sys.stdout):
    """
    get the STDOUT from bamtocov and return only the intervals
    between min and max 
    """
    cmd = [bin, bam]
    
    # Process the output line by line with a stream reader
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    sizes = {}
    totcov = {}
    for line in p.stdout:
        
        fields = line.decode().split()
        chrom = fields[0]
        span = int(fields[2]) - int(fields[1])
        coverage = int(fields[3])
        if chrom not in sizes:
            sizes[chrom] = span
            totcov[chrom] = coverage * span
        else:
            sizes[chrom] += span
            totcov[chrom] += coverage * span
    
    totalsize = sum(sizes.values())
    totalcov  = sum(totcov.values())
    for refseq in sizes:
        print(refseq, sizes[refseq], round(totcov[refseq]/sizes[refseq], 2), sep="\t", file=out)

    print("#Total", totalsize, round(totalcov/totalsize, 2), sep="\t", file=out)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bam', help='BAM file')
    parser.add_argument('-o', '--output', help='Output file')
    parser.add_argument('-b', '--bin', default="bamtocov", help='bamtocov binary [default: bamtocov]')
    parser.add_argument('-s', '--step', default=10, help='Step size [default: 10]', type=int)
    opts = parser.parse_args()

    if opts.output is None:
        out = sys.stdout
    else:
        out = open(opts.output, 'w')
    
    avgCoverage(opts.bam, opts.bin, out)