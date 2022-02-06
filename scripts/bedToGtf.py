#!/usr/bin/env python

import argparse, sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert bed file to gtf file.")
    parser.add_argument("bed", help="Input bed file.")
    parser.add_argument("gtf", help="Output gtf file.")
    args = parser.parse_args()

    if args.gtf:
        outfile = open(args.gtf, "w")
    else:
        outfile = sys.stdout
    with open(args.bed) as f:
        for line in f:
            line = line.strip()
            if line:
       
                fields = line.split()
       
                chrom, start, end, name = fields[0], fields[1], fields[2], fields[3]
                strand = "+"
                #chr2    215593349       215593782       REG1
                #seq1    bamtocov        exon    5       111     .       +       .       gene_id "include_5"
                print(f"{chrom}\tbamtocov\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"{name}\"", file=outfile)
                