#!/usr/bin/env python

import os, sys

GTFTYPE = "exon"
SUBSET  = 100

if (len(sys.argv) != 3):
    print("Usage: %s <input_dir> <output_file>" % sys.argv[0])
    sys.exit(1)

inputfile = sys.argv[1]
outputfile = sys.argv[2]

outstream = open(outputfile, "w")

with open(inputfile) as f:
    # Read input file line by line
    lines = f.readlines()
    c = 0
    for line in lines:
        if line.startswith("#"):
            print(line, file=outstream)
        else:
            fields = line.strip().split("\t")
            if fields[2] != GTFTYPE:
                continue
            c += 1
            if c % SUBSET == 0:
                print(line, file=outstream)

