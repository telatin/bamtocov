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
    Make a header for the BAM file given
    a dictionary of sequences and their lengths
    """
    header = { 'HD': {'VN': '1.0'},
            'SQ': [] }

    for seqName in seqs:
        # {'LN': 1575, 'SN': 'chr1'},
        header['SQ'].append({'LN': seqs[seqName], 'SN': seqName })
    
    return header

def makeTarget(genomeSize, targetSize, numFeatures=1):
    """
    Make a target genome with total lengt of targetSize
    being a subset of genomeSize. Return a dictionary of
    start and end positions for each feature
    """
    target = {}
    featureSize = int (targetSize / numFeatures)
    while len(target) < numFeatures:
        start =  round(random.randint(0, genomeSize - featureSize))
        end = start + featureSize
        # Check if this feature overlaps with any other
        overlap = False
        for f in target:
            if start < f and end > target[f]:
                overlap = True
                break
        if not overlap:
            target[start] = end
    return target

def savetarget(target, file):
    """
    Save target in bed format
    """
    bed = open(file, "w")
    chrname = "chromosome"
    
    for t in target:
        name = f"chromosome:{t}-{target[t]}"
        print(f"{chrname}\t{t}\t{target[t]}\t{name}", file=bed)
    

    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate a BAM file with long reads mapped against a hypothetical genome")
    parser.add_argument("-o", "--output", help="Output BAM file", required=True)
    parser.add_argument("-l", "--length", help="Length of the genome [default: %(default)s]", default="100M")
    
    parser.add_argument("-n", "--num-reads", help="Number of reads [default: %(default)s]", default="1M")
    parser.add_argument("-t", "--target-size", help="Bases in target [default: %(default)s]", default="10M")
    parser.add_argument("-f", "--target-features", help="Number of features [default: %(default)s]", type=int, default=10)
    
    parser.add_argument("-m", "--min-len", help="Minimum read length [default: 1000]", type=int, default=50)
    parser.add_argument("-M", "--max-len", help="Minimum read length [default: 10000]", type=int,default=300)
    parser.add_argument("-s", "--seed", help="Random seed [default: 42]", type=int,default=42)
    parser.add_argument("--multiply", help="Multiply the number of reads by this number", type=int, default=1)
    parser.add_argument("--randomcigar", help="Use random CIGAR strings", action="store_true")
    parser.add_argument("--progress", help="Print progress every INT reads [default: 10000]", type=int, default=10000)

    
    
    opts = parser.parse_args()

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