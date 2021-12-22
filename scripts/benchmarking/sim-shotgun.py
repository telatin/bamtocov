#!/usr/bin/env python3
"""
Given a list of chromosome names and their lengths, generate a sorted SAM file
of a shotgun.
"""

import os, sys, random
import argparse

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parseSize(size):
    # Size in the format of NUMBER[K|M|G|T]
    size = size.upper()
    try:
        if size.endswith('K'):
            return int(size[:-1]) * 1000
        elif size.endswith('M'):
            return int(size[:-1]) * 1000 * 1000
        elif size.endswith('G'):
            return int(size[:-1]) * 1000 * 1000 * 1000
        elif size.endswith('T'):
            return int(size[:-1]) * 1000 * 1000 * 1000 * 1000
        else:
            return int(size)
    except ValueError:
        eprint("Invalid size:", size)
        sys.exit(1)

def aln(readname, chrname, pos):
    qual = 100
    flag = 0
    cigar = "100M"
    insertLen = 0
    posN = 0
    chrN = "*"
    if args.randomstrand and random.random() < 0.5:
        flag = 16
    return f"{readname}\t{flag}\t{chrname}\t{pos}\t{qual}\t{cigar}\t{chrN}\t{posN}\t{insertLen}\t*\t*"

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("chrlen",   help="Chromosome names and lengths in the NAME:SIZE format or simply SIZE", nargs="+")
    args.add_argument("-s", "--step", help="Step size for the coverage [default: 100]", type=int, default=100)
    args.add_argument("--stack",  help="Depth in region stacks [default: 10]", default=10, type=int)
    args.add_argument("--noise",  help="Random error in position [default: 10]", default=10, type=int)
    args.add_argument("randomstrand",  help="Random strand [default: no]", action="store_true")
    args.add_argument("--out",   help="Output file in SAM format [default: stdout]")
    args = args.parse_args()

    if args.out:
        eprint("Writing to", args.out)
        out = open(args.out, "w")
    else:
        eprint("Writing to stdout")
        out = sys.stdout

    # Read chromosome names and lengths
    chroms = {}
    c = 0
    for chromSize in args.chrlen:
        c += 1
        try:
            chrom, size = chromSize.split(":")
            chroms[chrom] = parseSize(size)
            eprint(f" - {chrom}")
        except ValueError:
            chroms[f"chr{c}"] = parseSize(chromSize)
            
    
    # Sort chroms by size (descending)
    chroms = sorted(chroms.items(), key=lambda x: x[1], reverse=True)

    # Write header
    print("@HD\tVN:1.4\tSO:coordinate", file=out)
    for chrom, length in chroms:
        print(f"@SQ\tSN:{chrom}\tLN:{length}", file=out)
    
    # Write reads
    for chrom, length in chroms:
        for pos in range(1, length, args.step):
            if pos + args.step > length:
                continue
            readname = f"{chrom}-{pos}"
            for i in range(0, args.stack):
                # Random number between -10 and 10
                pos = pos + random.randint(-1 * args.noise, args.noise)
                if pos <= 0:
                    pos = 1
                elif pos >= length:
                    pos = length - args.step
                print( aln(readname, chrom, pos), file=out)