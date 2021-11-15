#!/usr/bin/env python3
"""
A script to remove the sequeces from a bam file.
"""
import os, sys
import argparse
import pysam
import random

def eprint(*args, **kwargs):
    """
    Print to stderr.
    """
    print(*args, file=sys.stderr, **kwargs)

def shuffle_string(s):
    """
    Shuffle a string.
    """
    if s is None:
        return None
    else:
        return ''.join(random.sample(s, len(s)))
        
def erase_seq(s):
    """
    Erase the sequence from a string.
    """
    if s is None:
        return None
    else:
        # Repeat "N" * len(s) times
        return 'N' * len(s)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Strip sequence and qualities from BAM file')
    parser.add_argument('bam', help='BAM file')
    parser.add_argument("-o", "--output", help="Output BAM file", default=None)
    parser.add_argument("-r", "--randomize-seq", help="Shuffle seq instead removing (quality unchanged)", action="store_true")
    parser.add_argument("-m", "--mask-seq", help="Replace sequence with Ns (quality unchanged)", action="store_true")
    parser.add_argument("-n", "--seq-prefix", help="Rename sequences with a prefix and progressive number", default=None)
    parser.add_argument("-u", "--remove-unmapped", help="Skip unaligned reads", default=None)
    opts = parser.parse_args()

    if opts.output is None:
        opts.output = os.path.splitext(opts.bam)[0] + '.no-seq.bam'
        eprint(f"Output: {opts.output}")
    
    if opts.randomize_seq and opts.mask_seq:
        raise ValueError("Cannot randomize and wipe sequence at the same time.")

    seqnames = {}
    c = 0

    bam = pysam.AlignmentFile(opts.bam, 'rb')
    bam_out = pysam.AlignmentFile(opts.output, 'wb', template=bam)
    
    for read in bam:
        if opts.remove_unmapped is not None and read.is_unmapped:
            continue
        if opts.randomize_seq:
            # Shuffle sequence
            new_qual = read.qual
            read.seq = shuffle_string(read.seq)
            read.qual = new_qual
        elif opts.mask_seq:
            # Replace sequence with Ns
            new_qual = read.qual
            read.seq = erase_seq(read.seq)
            read.qual = new_qual
        else:
            read.seq = None
            read.qual = None

        if opts.seq_prefix is not None:
            if read.qname in seqnames:
                read.qname = seqnames[read.qname]
            else:
                c += 1
                newname = opts.seq_prefix + str(c)
                seqnames[read.qname] = newname
                read.qname = newname
        bam_out.write(read)
    bam.close()
