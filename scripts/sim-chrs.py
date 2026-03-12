#!/usr/bin/env python3
"""
sim-chrs.py — Generate synthetic BAM files with highly variable per-contig read counts.

No annotation file required: contigs are generated synthetically.
Writes N BAM files and a ground-truth {outdir}/contigs.tsv.

Usage example:
  python sim-chrs.py -c 12 -n 3 -o sim_out/ \\
      --mean-reads 500 --read-length 100 \\
      --seed 42 --suppl-ratio 0.05

Then check coverage with bamtocov:
  bamtocov sim_out/sample*.bam
"""

import argparse
import os
import random
import csv
import sys
import math
import pysam


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Simulate BAM files with variable per-contig coverage."
    )
    p.add_argument("-c", "--contig-num", type=int, default=10,
                   help="Number of synthetic contigs/chromosomes (default: 10)")
    p.add_argument("-n", "--num-samples", type=int, default=3,
                   help="Number of BAM files to generate (default: 3)")
    p.add_argument("-o", "--outdir", required=True,
                   help="Output directory")

    # Contig size
    p.add_argument("--min-contig-len", type=int, default=10_000,
                   help="Minimum contig length in bp (default: 10000)")
    p.add_argument("--max-contig-len", type=int, default=500_000,
                   help="Maximum contig length in bp (default: 500000)")

    # Read simulation
    p.add_argument("--mean-reads", type=int, default=500,
                   help="Geometric mean of reads per contig per sample (default: 500)")
    p.add_argument("--read-length", type=int, default=100,
                   help="Read length in bp (default: 100)")
    p.add_argument("--variability", type=float, default=2.5,
                   help="Log-normal sigma for count variability — higher = more uneven (default: 2.5)")
    p.add_argument("--paired", action="store_true",
                   help="Simulate paired-end reads (default: single-end)")

    # Noise flags (same conventions as sim-counts.py)
    p.add_argument("--suppl-ratio", type=float, default=0.0,
                   help="Ratio of supplementary alignments per primary read (default: 0)")
    p.add_argument("--secondary-ratio", type=float, default=0.0,
                   help="Ratio of secondary alignments per primary read (default: 0)")
    p.add_argument("--duplicate-ratio", type=float, default=0.0,
                   help="Ratio of PCR duplicate reads (default: 0)")
    p.add_argument("--unmapped-ratio", type=float, default=0.0,
                   help="Ratio of unmapped reads sprinkled in (default: 0)")
    p.add_argument("--low-mapq-ratio", type=float, default=0.0,
                   help="Ratio of reads given MAPQ=0 (default: 0)")

    # Reproducibility
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducibility")

    return p.parse_args()


# ---------------------------------------------------------------------------
# Contig generation
# ---------------------------------------------------------------------------

def generate_contigs(n, min_len, max_len, rng):
    """Return list of (name, length) for n synthetic contigs."""
    contigs = []
    for i in range(1, n + 1):
        name = f"chr{i}"
        length = rng.randint(min_len, max_len)
        contigs.append((name, length))
    return contigs


def lognormal_count(mean_reads, sigma, rng):
    """Draw a read count from a log-normal distribution with given geometric mean and sigma."""
    mu = math.log(mean_reads)
    z = rng.gauss(0.0, 1.0)
    return max(1, int(round(math.exp(mu + sigma * z))))


# ---------------------------------------------------------------------------
# BAM header
# ---------------------------------------------------------------------------

def build_header(contigs):
    sq_list = [{"SN": name, "LN": length} for name, length in contigs]
    return pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": sq_list,
        "PG": [{"ID": "sim_chrs", "PN": "sim-chrs", "VN": "1.0"}],
    })


# ---------------------------------------------------------------------------
# Read builders (mirrors sim-counts.py)
# ---------------------------------------------------------------------------

def _base_read(name, chrom, pos0, read_length, mapq, flag, header):
    a = pysam.AlignedSegment(header)
    a.query_name      = name
    a.flag            = flag
    a.reference_id    = header.get_tid(chrom)
    a.reference_start = pos0
    a.mapping_quality = mapq
    a.cigar           = [(0, read_length)]
    a.query_sequence  = "A" * read_length
    a.query_qualities = pysam.qualitystring_to_array("I" * read_length)
    a.set_tag("NM", 0)
    return a


def make_primary(name, chrom, pos0, read_length, header, mapq=60):
    return _base_read(name, chrom, pos0, read_length, mapq, 0, header)


def make_supplementary(name, chrom, pos0, read_length, header):
    return _base_read(name, chrom, pos0, read_length, 0, 0x800, header)


def make_secondary(name, chrom, pos0, read_length, header):
    return _base_read(name, chrom, pos0, read_length, 0, 0x100, header)


def make_duplicate(name, chrom, pos0, read_length, header):
    return _base_read(name, chrom, pos0, read_length, 60, 0x400, header)


def make_low_mapq(name, chrom, pos0, read_length, header):
    return _base_read(name, chrom, pos0, read_length, 0, 0, header)


def make_unmapped(name, header):
    a = pysam.AlignedSegment(header)
    a.query_name      = name
    a.flag            = 0x4
    a.reference_id    = -1
    a.reference_start = 0
    a.mapping_quality = 0
    a.cigar           = []
    a.query_sequence  = "A" * 50
    a.query_qualities = pysam.qualitystring_to_array("I" * 50)
    return a


def make_paired_primary(name, chrom, pos0, read_length, header, frag_size=300, mapq=60):
    """Return (read1, read2) for a paired-end fragment."""
    r2_pos = max(0, pos0 + frag_size - read_length)

    f1 = pysam.AlignedSegment(header)
    f1.query_name      = name
    f1.flag            = 0x1 | 0x2 | 0x40 | 0x20   # paired, proper, read1, mate-reverse
    f1.reference_id    = header.get_tid(chrom)
    f1.reference_start = pos0
    f1.mapping_quality = mapq
    f1.cigar           = [(0, read_length)]
    f1.query_sequence  = "A" * read_length
    f1.query_qualities = pysam.qualitystring_to_array("I" * read_length)
    f1.set_tag("NM", 0)
    f1.next_reference_id    = header.get_tid(chrom)
    f1.next_reference_start = r2_pos
    f1.template_length = frag_size

    f2 = pysam.AlignedSegment(header)
    f2.query_name      = name
    f2.flag            = 0x1 | 0x2 | 0x80 | 0x10   # paired, proper, read2, reverse
    f2.reference_id    = header.get_tid(chrom)
    f2.reference_start = r2_pos
    f2.mapping_quality = mapq
    f2.cigar           = [(0, read_length)]
    f2.query_sequence  = "T" * read_length
    f2.query_qualities = pysam.qualitystring_to_array("I" * read_length)
    f2.set_tag("NM", 0)
    f2.next_reference_id    = header.get_tid(chrom)
    f2.next_reference_start = pos0
    f2.template_length = -frag_size

    return f1, f2


# ---------------------------------------------------------------------------
# Core simulation per sample
# ---------------------------------------------------------------------------

def simulate_sample(sample_id, contigs, args, header, rng):
    """
    Simulate one BAM file.
    Returns dict {contig_name: true_primary_count}.
    """
    bam_path = os.path.join(args.outdir, f"sample{sample_id}.bam")
    counts   = {}
    read_idx = 0

    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:

        for chrom, chrom_len in contigs:
            # Very variable count: log-normal draw
            read_count = lognormal_count(args.mean_reads, args.variability, rng)
            true_count = 0

            max_pos = max(0, chrom_len - args.read_length)

            for _ in range(read_count):
                read_idx += 1
                rname = f"read_{sample_id}_{read_idx}"
                pos0  = rng.randint(0, max_pos) if max_pos > 0 else 0

                if args.paired:
                    r1, r2 = make_paired_primary(
                        rname, chrom, pos0, args.read_length, header)
                    bam.write(r1)
                    bam.write(r2)
                    true_count += 1
                else:
                    mapq = 0 if rng.random() < args.low_mapq_ratio else 60
                    r = make_primary(rname, chrom, pos0, args.read_length,
                                     header, mapq=mapq)
                    bam.write(r)
                    if mapq > 0:
                        true_count += 1

                # noise: supplementary
                if rng.random() < args.suppl_ratio:
                    read_idx += 1
                    bam.write(make_supplementary(
                        f"read_{sample_id}_{read_idx}",
                        chrom, min(pos0 + 5, max_pos), args.read_length, header))

                # noise: secondary
                if rng.random() < args.secondary_ratio:
                    read_idx += 1
                    bam.write(make_secondary(
                        f"read_{sample_id}_{read_idx}",
                        chrom, pos0, args.read_length, header))

                # noise: PCR duplicate
                if rng.random() < args.duplicate_ratio:
                    read_idx += 1
                    bam.write(make_duplicate(
                        f"read_{sample_id}_{read_idx}",
                        chrom, pos0, args.read_length, header))

                # noise: extra low-mapq copy
                if rng.random() < args.low_mapq_ratio:
                    read_idx += 1
                    bam.write(make_low_mapq(
                        f"read_{sample_id}_{read_idx}",
                        chrom, pos0, args.read_length, header))

            counts[chrom] = true_count

            # unmapped noise proportional to this contig's read budget
            n_unmapped = int(read_count * args.unmapped_ratio)
            for _ in range(n_unmapped):
                read_idx += 1
                bam.write(make_unmapped(f"read_{sample_id}_{read_idx}", header))

    # Sort and index
    sorted_path = bam_path.replace(".bam", ".sorted.bam")
    pysam.sort("-o", sorted_path, bam_path)
    os.replace(sorted_path, bam_path)
    pysam.index(bam_path)
    print(f"  [✓] {bam_path}")
    return counts


# ---------------------------------------------------------------------------
# Write ground-truth table
# ---------------------------------------------------------------------------

def write_table(outdir, contigs, all_counts, n_samples):
    """Write {outdir}/contigs.tsv with contig x sample matrix."""
    table_path = os.path.join(outdir, "contigs.tsv")
    sample_names = [f"sample{i+1}" for i in range(n_samples)]
    contig_names = [name for name, _ in contigs]

    with open(table_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["contig"] + sample_names)
        for name in contig_names:
            row = [name] + [all_counts[s].get(name, 0) for s in range(n_samples)]
            writer.writerow(row)

    print(f"\n[✓] Ground-truth table → {table_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    rng  = random.Random(args.seed)

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Generating {args.contig_num} synthetic contigs "
          f"({args.min_contig_len}–{args.max_contig_len} bp each)")
    contigs = generate_contigs(
        args.contig_num, args.min_contig_len, args.max_contig_len, rng)
    for name, length in contigs:
        print(f"  {name}: {length:,} bp")

    header = build_header(contigs)

    print(f"\nSimulating {args.num_samples} sample(s) → {args.outdir}/")
    print(f"  read_length={args.read_length}  paired={args.paired}")
    print(f"  mean_reads={args.mean_reads}  variability(sigma)={args.variability}")
    print(f"  noise: suppl={args.suppl_ratio}  secondary={args.secondary_ratio}  "
          f"dup={args.duplicate_ratio}  unmapped={args.unmapped_ratio}  "
          f"low_mapq={args.low_mapq_ratio}\n")

    all_counts = []
    for i in range(args.num_samples):
        counts = simulate_sample(i + 1, contigs, args, header, rng)
        all_counts.append(counts)

    write_table(args.outdir, contigs, all_counts, args.num_samples)

    bam_files = " ".join(
        os.path.join(args.outdir, f"sample{i+1}.bam")
        for i in range(args.num_samples)
    )
    paired_flag = "-p" if args.paired else ""
    print("\n── bamtocov command ───────────────────────────────────────────────")
    print(f"bamtocov {bam_files}")
    print(f"\n── samtools flagstat (first sample) ───────────────────────────────")
    print(f"samtools flagstat {os.path.join(args.outdir, 'sample1.bam')}")
    print("────────────────────────────────────────────────────────────────────\n")


if __name__ == "__main__":
    main()
