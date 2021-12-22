#!/usr/bin/env python3
"""
Simulate a BAM file with long reads mapped against a hypothetical genome
"""

import os, sys
import argparse
import pysam
import random
import re
from rich.progress import track

def eprint(*args, **kwargs):
    """
    Print to stderr
    """
    print(*args, file=sys.stderr, **kwargs)

def stringSizeToInt(stringSize):
    """
    Convert suffixes like Mbp Kpb to the int
    """
    suffixes = {'G': 1000000000, 'M': 1000000, 'K': 1000}
    if stringSize[-1] in suffixes:
        return int(stringSize[:-1]) * suffixes[stringSize[-1]]
    else:
        return int(stringSize)

def randomCIGAR(length):
    """
    Generate a random CIGAR string
    """
    total = 0
    cigar = ""
    ops = ['M', 'I', 'D']
    while total < length:
         
        oplen = random.randint(1, length)
        if total == 0:
            op = "M"
        else:
            op = random.choice(['M', 'I', 'D'])
            if op == cigar[-1]:
                while op == cigar[-1]:
                    op = random.choice(['M', 'I', 'D'])
                     
            if op == cigar[-1]:
                raise "Invalid CIGAR"
        if oplen + total > length:
            oplen = length - total
            
        if op in ['M', 'I', 'S']:
            total += oplen     
    
        cigar += str(oplen) + op
     

    if cigar[-1] != "M":
        cigar = cigar[:-1] + "S"
    
    if cigarToLen(cigar) != length:
        print(f"Invalid CIGAR: {cigar} {cigarToLen(cigar)}/{length}")
        exit(1)
    return cigar

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


def makeHeader(seqs):
    """
    Make a header for the BAM file
    """
    header = { 'HD': {'VN': '1.0'},
            'SQ': [] }

    for seqName in seqs:
        # {'LN': 1575, 'SN': 'chr1'},
        header['SQ'].append({'LN': seqs[seqName], 'SN': seqName })
    
    return header

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate a BAM file with long reads mapped against a hypothetical genome")
    parser.add_argument("-o", "--output", help="Output BAM file", required=True)
    parser.add_argument("-l", "--length", help="Length of the genome [default: 10M]", default="10M")
    parser.add_argument("-c", "--coverage", help="Approximate coverage [default: 50]", type=float, default=50)
    parser.add_argument("-m", "--min-len", help="Minimum read length [default: 1000]", type=int, default=1000)
    parser.add_argument("-M", "--max-len", help="Minimum read length [default: 10000]", type=int,default=10000)
    parser.add_argument("-x", "--multiply", help="Clone each read N times", type=int,default=1)
    parser.add_argument("-s", "--seed", help="Random seed [default: 42]", type=int,default=42)
    parser.add_argument("-r", "--randomcigar", help="Generate a random CIGAR", action="store_true")
    parser.add_argument("--progress", help="Print progress every INT reads [default: 10000]", type=int, default=10000)

    
    
    opts = parser.parse_args()

    # Reference genome
    try:
        genomeSize = stringSizeToInt(opts.length)
    except:
        parser.error(f"ERROR: Invalid genome size {opts.length}")

    genome = { "chromosome": genomeSize}
    bamHeader = makeHeader(genome)
    totalBases = genomeSize * opts.coverage
    totalReads = 0
    generatedBases = 0
    
    seqPosLength = {}
    
    eprint(f"Generating positions to cover {totalBases} bp")
    while generatedBases < totalBases:
        readLength = int(opts.min_len + (opts.max_len - opts.min_len) * random.random())
        generatedBases += readLength
        totalReads += 1
        pos = int(random.random() * (genomeSize - readLength - opts.max_len))
        if pos in seqPosLength:
            seqPosLength[pos].append(readLength)
        else:
            seqPosLength[pos] = [readLength]

    
    eprint(f"Generated positions for {totalReads} reads (from {len(seqPosLength)} positions)")

    # Sort seqPosLength keys ascending
    seqPosLength = {k: v for k, v in sorted(seqPosLength.items())}
    eprint(f"Generating {totalReads} reads, {generatedBases} bp")
    # Generate reads
    with pysam.AlignmentFile(opts.output, "wb", header=bamHeader) as outf:
        n = 0
        for index, pos in  track(enumerate(seqPosLength), total=len(seqPosLength), description="Writing BAM..."):
            
            
            for length in seqPosLength[pos]:
                
                
                for clone in range(opts.multiply):
                    n += 1
                    # Generate a list of length elements equal to 20
                    qualityList = [9] * length
                    seqString   = "A" * length
                    try:
                        a = pysam.AlignedSegment()
                        
            
                        a.query_sequence = "chr1"
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