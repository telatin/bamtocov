#!/usr/bin/env python
"""
Simulate a BAM file with long reads mapped against a hypothetical genome
"""

import os, sys
import argparse
import pysam
import random
import re

class Alignmnent(object):
    """
    Class to represent an alignment
    """
    def __init__(self, ref, pos, cigar, flag):
        self.ref = ref
        self.pos = int(pos)
        self.cigar = cigar
        self.flag = int(flag)

def eprint(*args, **kwargs):
    """
    Print to stderr
    """
    print(*args, file=sys.stderr, **kwargs)
 

def makeHeader(seqs):
    """
    Make a header for the BAM file given
    a dictionary of sequences and their lengths
    """
    header = { 'HD': {'VN': '1.0'},
            'SQ': [] }

    for seqName in seqs:
        # {'LN': 1575, 'SN': 'chr1'},
        header['SQ'].append({'LN': seqs[seqName], 'SN': seqName })
    
    return header
 
def loadBed(file):
    """
    Load alignment coordinates from a BED file: chr, start, end, strand
    """
    alignments = []
    references = {}
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            chromosome, start, end, strand = line.split()
            Cigar = str(int(end) - int(start)) + 'M'
            Flag = 0
            End = int(end) + 100
            if strand == '-':
                Flag = 16
            if chromosome not in references:
                references[chromosome] = End
            else:
                references[chromosome] = max(references[chromosome], End)
            
            alignments.append(Alignmnent(chromosome, int(start) - 1, Cigar, Flag))
    return alignments, references

def cigarToLen(cigar):
    """
    Calculate sequence length from CIGAR string
    """
    # Split "cigar" on capital letters
    span = re.split('[A-Z]', cigar)
    ops   = re.split('[0-9]+', cigar)
    len = 0
    del(span[-1])
    del(ops[0])
    for i, span in enumerate(span):
      if ops[i] in ["M", "I", "S"]:
        len += int(span)
    return len

def getRefId(refs, refname):
    """
    Get the reference ID for a reference name
    """
    for i, ref in enumerate(refs):
        
        if ref == refname:
            return i
    return -1
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate a BAM file with long reads mapped against a hypothetical genome")

    parser.add_argument("-i", "--input",  help="Bed file with alignments", required=True)
    parser.add_argument("-o", "--output", help="Output BAM file", required=True)
    parser.add_argument("-l", "--length", help="Length of the genome [default: %(default)s]", default="100M")
    
    parser.add_argument("-n", "--num-reads",       help="Number of reads [default: %(default)s]", default="1M")
    parser.add_argument("-t", "--target-size",     help="Bases in target [default: %(default)s]", default="10M")
    parser.add_argument("-f", "--target-features", help="Number of features [default: %(default)s]", type=int, default=10)
    
    parser.add_argument("-m", "--min-len", help="Minimum read length [default: 1000]", type=int, default=50)
    parser.add_argument("-M", "--max-len", help="Minimum read length [default: 10000]", type=int,default=300)
    parser.add_argument("-s", "--seed",    help="Random seed [default: 42]", type=int,default=42)
    parser.add_argument("--multiply",      help="Multiply the number of reads by this number", type=int, default=1)
    parser.add_argument("--randomcigar",   help="Use random CIGAR strings", action="store_true")
    parser.add_argument("--progress",      help="Print progress every INT reads [default: 10000]", type=int, default=10000)

    opts = parser.parse_args()

    alignments, references = loadBed(opts.input)
    bamHeader = makeHeader(references)
    print(bamHeader)
    with pysam.AlignmentFile(opts.output, "wb", header=bamHeader) as outf:
        n = 0
        for aln in alignments:
                        n += 1
                        a = pysam.AlignedSegment(header=outf.header)
                        a.reference_id = getRefId(references,aln.ref)
                        
                        a.reference_name = aln.ref
                        a.flag = aln.flag
                        
                        a.reference_start = aln.pos
                        a.mapping_quality = 60
                        a.cigarstring = aln.cigar
                        
                        a.query_name = f"read-{aln.ref}:{getRefId(references,aln.ref)}_" + str(n)
                        #a.next_reference_id = 0
                        #a.next_reference_start = 0
                        #a.template_length = 0
                        a.query_sequence = "T" * cigarToLen(a.cigarstring)
                        #a.query_qualities = [10] * len(a.query_sequence)
                        print(f"{a.query_name} {a.reference_id} {a.reference_start} {a.cigarstring}")
                        outf.write(a)

    exit(0)
    # Reference genome
    try:
        genomeSize = stringSizeToInt(opts.length)
    except:
        parser.error(f"ERROR: Invalid genome size {opts.length}")

    genome = { "chromosome": genomeSize}
    bamHeader = makeHeader(genome)
    target = makeTarget(genomeSize, stringSizeToInt(opts.target_size), opts.target_features)

    bedOutput = opts.output.replace("bam", "bed")
    savetarget(target, bedOutput)
    totalReads = stringSizeToInt(opts.num_reads)
    generatedBases = 0
    n = 0
    seqPosLength = {}
     
    for n in track(range(totalReads), total=totalReads, description="Generating positions..."):
    #while n < totalReads:
        #n += 1
        readLength = int(opts.min_len + (opts.max_len - opts.min_len) * random.random())
        
        # Generate a random position being included in a random interval of the target
        intervalStart = random.choice(list(target.keys()))
        intervalEnd = target[intervalStart]
        # Pos should be >= intervalStart and < intervalEnd - readLength
        pos = random.randint(intervalStart, intervalEnd - readLength)
        if pos in seqPosLength:
            seqPosLength[pos].append(readLength)
        else:
            seqPosLength[pos] = [readLength]


    
    
    # Sort seqPosLength keys ascending
    seqPosLength = {k: v for k, v in sorted(seqPosLength.items())}
    # Generate reads
    with pysam.AlignmentFile(opts.output, "wb", header=bamHeader) as outf:
        n = 0
        for index, pos in  track(enumerate(seqPosLength), total=len(seqPosLength), description="Writing BAM...         "):
            
            
            for length in seqPosLength[pos]:
                
                
                for clone in range(opts.multiply):
                    n += 1
                    # Generate a list of length elements equal to 20
                    qualityList = [9] * length
                    seqString   = "A" * length
                    try:
                        a = pysam.AlignedSegment()
                        #a.query_sequence = "chr1"
                        a.flag = 0
                        a.reference_id = 0
                        a.reference_start = pos
                        a.mapping_quality = 60
                        if not opts.randomcigar:
                            a.cigarstring = str(length) + 'M'
                        else:
                            a.cigarstring = randomCIGAR(length)
                            if a.infer_read_length() != length or pos + length >= genomeSize:
                                print(f"Invalid CIGAR: {a.cigarstring} {a.infer_read_length()}/{length}")
                                exit(1)
                        
                        a.query_name = "pos_" + str(pos) + "_len" + str(cigarToLen(a.cigarstring)) + "_end" + str(pos + cigarToLen(a.cigarstring)) + "_" + str(clone)
                        a.next_reference_id = 0
                        a.next_reference_start = 0
                        a.template_length = 0
                        a.query_sequence = seqString
                        if pos + length > genomeSize:
                            raise Exception("Read extends past end of genome")
                        outf.write(a)


                    except Exception as e:
                        eprint("ERROR:", sys.exc_info(), a)
                        eprint(f"{n} {pos} {seqPosLength[pos]}")
                        eprint(f"{len(seqString)},{seqString}")
                        eprint(f"{len(qualityList)},{qualityList}")
                        sys.exit(1)