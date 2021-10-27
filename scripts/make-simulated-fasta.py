#!/usr/bin/env python3
"""
A tool to generate simulated FASTA files.
"""
import os, sys
import argparse
import random

def formatSequence(seq, width=80):
    """
    Format the sequence for FASTA output.
    """
    if width < 1:
        return seq
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))
def makeSequence(name, length, width):
    """
    Generate a random sequence of the specified length.
    """
    seq = ""
    for i in range(length):
        seq += random.choice("ACGT")
    return ">{}\n{}".format(name, formatSequence(seq, width))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="A tool to generate simulated BAM files")
    parser.add_argument("SEQUENCE_LENGTH", type=int, help="The length of the sequence to generate", nargs="+")
    parser.add_argument("-o", "--output", type=str, help="Output file [default: stdout]", default="-")
    parser.add_argument("-n", "--name", type=str, help="Sequence prefix [defaut: seq]", default="seq")
    
    parser.add_argument("-l", "--add-len", action="store_true", help="Add length as comment")
    parser.add_argument("-s", "--seed", type=int, help="Random seed [default: 1]", default=1)
    parser.add_argument("-w", "--line-width", type=int, help="Fasta sequence line width [default: 100]", default=100)
    parser.add_argument("--sep", type=str, help="Sequence prefix and number separator [default: .]", default=".")
    opts = parser.parse_args()

    # Set seed
    random.seed(opts.seed)

    # Open output file
    if opts.output == "-":
        out = sys.stdout
    else:
        out = open(opts.output, "w")

    comment = ""
    for i, length in enumerate(opts.SEQUENCE_LENGTH):
        if opts.add_len:
            comment = " length={}".format(length)
        print(makeSequence(opts.name + opts.sep + str(i) + comment, length, opts.line_width), file=out)
