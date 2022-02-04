#!/usr/bin/env python3
"""
Print some statistics (total seqs, max len, n50) for both the reference and the
query sequences of a BAM file
"""
import pysam
import argparse
import sys
import os
import math
def dict_stats(sizes):
    stats = {
        'tot': 0,
        'sum': 0,
        'n50': 0,
        'min': 0,
        'max': 0,
        'avg': 0,
        'half': 0,
    }
    stats["sum"] = sum([x*sizes[x] for x in sizes])
    stats["half"] = stats["sum"]/2
    
    # sort dictionary by keys descending
    sorted_sizes = sorted(sizes.items(), key=lambda x: x[0], reverse=True)

    partial = 0
    stats["max"] = sorted_sizes[0][0]
    stats["min"] = sorted_sizes[-1][0]
    for key, value in sorted_sizes:
        stats["tot"] += value
        partial += key*value
        if partial > stats["half"] and stats["n50"] == 0:
            stats["n50"] = key
    stats["avg"] = stats["sum"]/stats["tot"]
    return f"{stats['tot']}\t{stats['max']}\t{stats['n50']}"  
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam", 
                            help="BAM file to calculate stats for.")
    
    opts = parser.parse_args()

    bam = pysam.AlignmentFile(opts.bam, "rb")

    # Total number of reads
    total_mapped_reads = 0
    read_lengths = {}
    ref_lengths = {}
    
    for i in bam.references:
        if bam.header.get_reference_length(i) in ref_lengths:
            ref_lengths[bam.header.get_reference_length(i)] += 1
        else:
            ref_lengths[bam.header.get_reference_length(i)] = 1
    
    min = 1000000000
    max = 0
    for read in bam:
        if not read.is_unmapped:
            total_mapped_reads += 1 
            readLen = read.infer_query_length()
            if not readLen is None :
                if readLen in read_lengths:
                    read_lengths[readLen] += 1
                else:
                    read_lengths[readLen] = 1

    print("Reference\t", dict_stats(ref_lengths))
    print("Reads\t", dict_stats(read_lengths))
# ['__class__', 
# '__copy__', 
# '__deepcopy__', 
# '__delattr__', 
# '__dir__', 
# '__doc__', 
# '__eq__', 
# '__format__', 
# '__ge__', 
# '__getattribute__', 
# '__gt__', 
# '__hash__', 
# '__init__', 
# '__init_subclass__', 
# '__le__', 
# '__lt__', 
# '__ne__', 
# '__new__', 
# '__pyx_vtable__', 
# '__reduce__', 
# '__reduce_ex__', 
# '__repr__', 
# '__setattr__', 
# '__setstate__', 
# '__sizeof__', 
# '__str__', 
# '__subclasshook__', 
# 'aend', 
# 'alen', 
# 'aligned_pairs', 
# 'bin', 
# 'blocks', 
# 'cigar', 
# 'cigarstring', 
# 'cigartuples', 
# 'compare', 
# 'flag', 
# 'from_dict', 
# 'fromstring', 
# 'get_aligned_pairs', 
# 'get_blocks', 
# 'get_cigar_stats', 
# 'get_forward_qualities', 
# 'get_forward_sequence', 
# 'get_overlap', 
# 'get_reference_positions', 
# 'get_reference_sequence', 
# 'get_tag', 
# 'get_tags', 
# 'has_tag', 
# 'header', 
# 'infer_query_length', 
# 'infer_read_length', 
# 'inferred_length', 
# 'is_duplicate', 
# 'is_paired', 
# 'is_proper_pair', 
# 'is_qcfail', 
# 'is_read1', 
# 'is_read2', 
# 'is_reverse', 
# 'is_secondary', 
# 'is_supplementary', 
# 'is_unmapped', 
# 'isize', 
# 'mapping_quality', 
# 'mapq', 
# 'mate_is_reverse', 
# 'mate_is_unmapped', 
# 'mpos', 
# 'mrnm', 
# 'next_reference_id', 
# 'next_reference_name', 
# 'next_reference_start', 
# 'opt', 
# 'overlap', 
# 'pnext', 
# 'pos', 
# 'positions', 
# 'qend', 
# 'qlen', 
# 'qname', 
# 'qqual', 
# 'qstart', 
# 'qual', 
# 'query', 
# 'query_alignment_end', 
# 'query_alignment_length', 
# 'query_alignment_qualities', 
# 'query_alignment_sequence', 
# 'query_alignment_start', 
# 'query_length', 
# 'query_name', 
# 'query_qualities', 
# 'query_sequence', 
# 'reference_end', 
# 'reference_id', 
# 'reference_length', 
# 'reference_name', 
# 'reference_start', 
# 'rlen', 
# 'rname', 
# 'rnext', 
# 'seq', 
# 'setTag', 
# 'set_tag', 
# 'set_tags', 
# 'tags', 
# 'template_length', 
# 'tid', 
# 'tlen', 
# 'to_dict', 
# 'to_string', 
# 'tostring']