#!/usr/bin/env python3
"""
Update the annotation chromosome name based on the original contig
names. Supplying the original contig FASTA files, the renamed FASTA
file and an annotation in GFF or BED format, this will be converted
using the original names.
"""

import os, sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def loadFastaNames(fn):
    fasta_names = []
    with open(fn) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line.strip().split()[0][1:]
                fasta_names.append(name)
    return fasta_names

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ANNOTATION', help='Input annotation in GFF3 or BED format')
    parser.add_argument('-r', '--reference-contigs', help='Reference genome in FASTA format', required=True)
    parser.add_argument('-p', '--prokka-contigs', help='Prokka contigs renamed', required=True)
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')
    args = parser.parse_args()

    originalChrNames = loadFastaNames(args.reference_contigs)

    if args.verbose:
        eprint("Original references:", len(originalChrNames), ":", originalChrNames[:3])
    
    prokkaChrNames = loadFastaNames(args.prokka_contigs)
    if args.verbose:
        eprint("Prokka references:", len(prokkaChrNames), ":", prokkaChrNames[:3])

    
    if len(originalChrNames) != len(prokkaChrNames):
        eprint("ERROR: Number of contigs do not match")
        sys.exit(1)

    with open(args.ANNOTATION) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    # ##sequence-version CHRNAME START END
                    chrName = line.strip().split()[1]
                    # Replace chrName with the original name
                    if chrName in prokkaChrNames:
                        line = line.replace(chrName, originalChrNames[prokkaChrNames.index(chrName)])
                         
                        print(line.strip())
                    else:
                        eprint("ERROR:", chrName, "not found in prokka contigs")
                        sys.exit(1)
                else:
                    print(line.strip())
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")

            # Bed or GFF?
            if len(fields) > 2:
                if fields[0] in prokkaChrNames:
                    fields[0] = originalChrNames[prokkaChrNames.index(fields[0])]
                    print("\t".join(fields))
                else:
                    eprint("ERROR: Contig not found:", fields[0])
                    sys.exit(1)