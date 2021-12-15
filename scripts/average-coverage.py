#!/usr/bin/env python
"""
Return statistics of coverage of a BAM file using "bamtocov".
Will print a three columns report: chromosome, length, average coverage
"""

import os, sys
import argparse
import subprocess


def eprint(*args, **kwargs):
    """
    Print to stderr.
    """
    print(*args, file=sys.stderr, **kwargs)



def AvgCov(bam, bin="bamtocov", out=sys.stdout):
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
    avgCov = {}
    for refseq in sizes:
        #print(refseq, sizes[refseq], round(totcov[refseq]/sizes[refseq], 2), sep="\t", file=out)
        avgCov[refseq] = round(totcov[refseq]/sizes[refseq], 2)
    
    #print("#Total", totalsize, round(totalcov/totalsize, 2), sep="\t", file=out)
    avgCov["#Total"] = round(totalcov/totalsize, 2)
    return avgCov

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bam',  help='BAM file', nargs='+')
    parser.add_argument('-o', '--output', help='Output file')
    parser.add_argument('-t', '--total', help='Print total coverage', action='store_true')
    #parser.add_argument('-n', '--normalize', help='Normalize coverage', action='store_true')
    parser.add_argument('-b', '--bin', default="bamtocov", help='bamtocov binary [default: bamtocov]')

    opts = parser.parse_args()

    if opts.output is None:
        out = sys.stdout
    else:
        out = open(opts.output, 'w')
    
    matrix = {}
    files = []
    for bamFile in opts.bam:
        files.append(os.path.basename(bamFile).replace(".bam", ""))
        table = AvgCov(bamFile, opts.bin, out)
        for refseq in table:
            if refseq not in matrix:
                matrix[refseq] = []
            matrix[refseq].append(table[refseq])

    # Print header: files
    print("Chromosome\t", "\t".join(files), file=out)
    for refseq in matrix:
        if len(matrix[refseq]) != len(opts.bam):
            eprint("Warning: missing coverage for", refseq)
            exit(1)
        if not opts.total and refseq == "#Total":
            continue
        print(refseq, *matrix[refseq], sep="\t", file=out)

    if opts.output is not None:
        out.close()
