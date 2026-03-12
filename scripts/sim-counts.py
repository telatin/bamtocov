#!/usr/bin/env python3
"""
simulate_bam.py — Generate synthetic BAM files with known counts for featureCounts testing.

Reads a BED/GFF/GTF annotation file, simulates reads per feature (with edge cases
and optional noise), writes N BAM files, and a ground-truth table.tsv.

Usage example:
  python simulate_bam.py -t features.gtf -n 3 -o sim_out/ \\
      --reads-per-feature 50 --read-length 75 \\
      --suppl-ratio 0.05 --secondary-ratio 0.03 \\
      --edge-ratio 0.10 --unmapped-ratio 0.02 \\
      --seed 42

Then test with featureCounts:
  featureCounts -a features.gtf -o fc_counts.txt sim_out/sample*.bam
"""

import argparse
import os
import random
import csv
import sys
import pysam


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Simulate BAM files with known counts for featureCounts testing."
    )
    p.add_argument("-t", "--target", required=True,
                   help="Annotation file: BED, GFF, or GTF")
    p.add_argument("-n", "--num-samples", type=int, default=3,
                   help="Number of BAM files to generate (default: 3)")
    p.add_argument("-o", "--outdir", required=True,
                   help="Output directory")

    # Read simulation
    p.add_argument("--reads-per-feature", type=int, default=50,
                   help="Mean reads fully inside each feature per sample (default: 50)")
    p.add_argument("--read-length", type=int, default=100,
                   help="Read length in bp (default: 100)")
    p.add_argument("--paired", action="store_true",
                   help="Simulate paired-end reads (default: single-end)")

    # Noise flags
    p.add_argument("--suppl-ratio", type=float, default=0.0,
                   help="Ratio of supplementary alignments added per feature read (default: 0)")
    p.add_argument("--secondary-ratio", type=float, default=0.0,
                   help="Ratio of secondary alignments added per feature read (default: 0)")
    p.add_argument("--duplicate-ratio", type=float, default=0.0,
                   help="Ratio of PCR duplicate reads added (default: 0)")
    p.add_argument("--unmapped-ratio", type=float, default=0.0,
                   help="Ratio of unmapped reads to sprinkle in (default: 0)")
    p.add_argument("--edge-ratio", type=float, default=0.10,
                   help="Ratio of reads placed at feature edges / partially overlapping (default: 0.10)")
    p.add_argument("--multi-overlap-ratio", type=float, default=0.0,
                   help="Ratio of reads spanning two adjacent features (default: 0)")
    p.add_argument("--low-mapq-ratio", type=float, default=0.0,
                   help="Ratio of reads given MAPQ=0 (multi-mapper simulation, default: 0)")

    # Reproducibility
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducibility")

    return p.parse_args()


# ---------------------------------------------------------------------------
# Annotation parsers  →  list of (chrom, start0, end, name, strand)
# ---------------------------------------------------------------------------

def parse_bed(path):
    features = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            name   = cols[3] if len(cols) > 3 else f"{chrom}:{start}-{end}"
            strand = cols[5] if len(cols) > 5 else "+"
            features.append((chrom, start, end, name, strand))
    return features


def _gff_attr(attr_str, key, fmt="gtf"):
    """Extract a single attribute value from GFF/GTF attribute string."""
    if fmt == "gtf":
        # key "value";
        for part in attr_str.split(";"):
            part = part.strip()
            if part.startswith(key + " "):
                return part.split('"')[1] if '"' in part else part.split()[1]
    else:
        # key=value
        for part in attr_str.split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                if k.strip() == key:
                    return v.strip()
    return None


def parse_gff_gtf(path):
    features = []
    fmt = "gtf" if path.endswith(".gtf") else "gff"
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            feature_type = cols[2]
            # For GTF keep "exon" rows; for GFF keep "gene" or "mRNA" or anything
            if fmt == "gtf" and feature_type not in ("exon", "gene", "transcript"):
                continue
            chrom  = cols[0]
            start0 = int(cols[3]) - 1   # GFF/GTF is 1-based
            end    = int(cols[4])
            strand = cols[6] if cols[6] in ("+", "-") else "+"
            attr   = cols[8]
            name   = (_gff_attr(attr, "gene_id", fmt) or
                      _gff_attr(attr, "ID",      fmt) or
                      f"{chrom}:{start0}-{end}")
            features.append((chrom, start0, end, name, strand))
    # Deduplicate by name — keep first occurrence
    seen, unique = set(), []
    for f in features:
        if f[3] not in seen:
            seen.add(f[3])
            unique.append(f)
    return unique


def load_features(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".bed":
        return parse_bed(path)
    elif ext in (".gff", ".gff3", ".gtf"):
        return parse_gff_gtf(path)
    else:
        # Guess from content
        with open(path) as fh:
            head = fh.read(2000)
        if "\t" in head and head.count("\t") > head.count(" "):
            cols = head.split("\n")[0].split("\t")
            if len(cols) >= 9:
                return parse_gff_gtf(path)
        return parse_bed(path)


# ---------------------------------------------------------------------------
# BAM header / chromosome length helpers
# ---------------------------------------------------------------------------

def build_sq_list(features):
    """Build {chrom: max_end} from features to use as BAM header."""
    lengths = {}
    for chrom, start, end, *_ in features:
        lengths[chrom] = max(lengths.get(chrom, 0), end + 5000)
    return [{"SN": c, "LN": l} for c, l in sorted(lengths.items())]


# ---------------------------------------------------------------------------
# Read / alignment builders
# ---------------------------------------------------------------------------

def _base_read(name, chrom, pos0, read_length, mapq, flag, strand, header):
    """Return a configured AlignedSegment (not yet added to file)."""
    a = pysam.AlignedSegment(header)
    a.query_name       = name
    a.flag             = flag
    a.reference_id     = header.get_tid(chrom)
    a.reference_start  = pos0          # 0-based
    a.mapping_quality  = mapq
    a.cigar            = [(0, read_length)]   # simple match
    a.query_sequence   = "A" * read_length
    a.query_qualities  = pysam.qualitystring_to_array("I" * read_length)
    a.set_tag("NM", 0)
    # Strand: reverse bit (0x10)
    if strand == "-":
        a.flag |= 0x10
    return a


def make_primary(name, chrom, pos0, read_length, strand, header, mapq=60):
    """Standard primary alignment."""
    return _base_read(name, chrom, pos0, read_length, mapq, 0, strand, header)


def make_supplementary(name, chrom, pos0, read_length, strand, header):
    """Supplementary alignment (chimeric / split-read)."""
    return _base_read(name, chrom, pos0, read_length, 0, 0x800, strand, header)


def make_secondary(name, chrom, pos0, read_length, strand, header):
    """Secondary alignment (multi-mapper)."""
    return _base_read(name, chrom, pos0, read_length, 0, 0x100, strand, header)


def make_duplicate(name, chrom, pos0, read_length, strand, header):
    """PCR duplicate (FLAG 0x400)."""
    return _base_read(name, chrom, pos0, read_length, 60, 0x400, strand, header)


def make_low_mapq(name, chrom, pos0, read_length, strand, header):
    """Multi-mapper simulation: MAPQ=0."""
    return _base_read(name, chrom, pos0, read_length, 0, 0, strand, header)


def make_unmapped(name, header):
    """Unmapped read."""
    a = pysam.AlignedSegment(header)
    a.query_name      = name
    a.flag            = 0x4              # unmapped
    a.reference_id    = -1
    a.reference_start = 0
    a.mapping_quality = 0
    a.cigar           = []
    a.query_sequence  = "A" * 50
    a.query_qualities = pysam.qualitystring_to_array("I" * 50)
    return a


# ---------------------------------------------------------------------------
# Paired-end helpers
# ---------------------------------------------------------------------------

def make_paired_primary(name, chrom, pos0, read_length, strand, header,
                        frag_size=300, mapq=60):
    """Return (read1, read2) for a paired-end fragment."""
    # Read1
    f1 = pysam.AlignedSegment(header)
    f1.query_name       = name
    f1.flag             = 0x1 | 0x2 | 0x40      # paired, proper, read1
    f1.reference_id     = header.get_tid(chrom)
    f1.reference_start  = pos0
    f1.mapping_quality  = mapq
    f1.cigar            = [(0, read_length)]
    f1.query_sequence   = "A" * read_length
    f1.query_qualities  = pysam.qualitystring_to_array("I" * read_length)
    f1.set_tag("NM", 0)
    f1.next_reference_id     = header.get_tid(chrom)
    f1.next_reference_start  = pos0 + frag_size - read_length
    f1.template_length  = frag_size
    if strand == "-":
        f1.flag |= 0x10          # R1 on minus
        f1.flag |= 0x20          # R2 on plus (mate)
    else:
        f1.flag |= 0x20          # R2 on minus (mate)

    # Read2 — reverse strand, downstream
    r2_pos = pos0 + frag_size - read_length
    f2 = pysam.AlignedSegment(header)
    f2.query_name       = name
    f2.flag             = 0x1 | 0x2 | 0x80 | 0x10   # paired, proper, read2, reverse
    f2.reference_id     = header.get_tid(chrom)
    f2.reference_start  = max(0, r2_pos)
    f2.mapping_quality  = mapq
    f2.cigar            = [(0, read_length)]
    f2.query_sequence   = "T" * read_length
    f2.query_qualities  = pysam.qualitystring_to_array("I" * read_length)
    f2.set_tag("NM", 0)
    f2.next_reference_id    = header.get_tid(chrom)
    f2.next_reference_start = pos0
    f2.template_length  = -frag_size

    return f1, f2


# ---------------------------------------------------------------------------
# Core simulation per sample
# ---------------------------------------------------------------------------

def simulate_sample(sample_id, features, args, header, rng):
    """
    Simulate one BAM file.
    Returns: dict {feature_name: true_count} (primary non-dup reads inside feature)
    """
    bam_path = os.path.join(args.outdir, f"sample{sample_id}.bam")
    counts   = {}
    read_idx = 0

    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:

        for chrom, feat_start, feat_end, feat_name, strand in features:
            feat_len   = feat_end - feat_start
            read_count = max(1, rng.randint(
                int(args.reads_per_feature * 0.7),
                int(args.reads_per_feature * 1.3)
            ))
            true_count = 0   # only primary, non-dup, non-secondary, non-suppl

            for i in range(read_count):
                read_idx += 1
                rname = f"read_{sample_id}_{read_idx}"

                # ---- placement ----
                # Determine if this read is an "edge" read
                is_edge = rng.random() < args.edge_ratio
                if is_edge:
                    # Overhang: read starts just before or ends just after feature
                    overhang = rng.randint(1, args.read_length - 1)
                    if rng.random() < 0.5:
                        # Hang off the left edge (partially overlaps)
                        pos0 = feat_start - overhang
                    else:
                        # Hang off the right edge
                        pos0 = feat_end - overhang
                else:
                    # Fully inside the feature
                    max_start = feat_end - args.read_length
                    if max_start < feat_start:
                        pos0 = feat_start       # feature shorter than read
                    else:
                        pos0 = rng.randint(feat_start, max_start)

                pos0 = max(0, pos0)   # clamp to chrom start

                # ---- emit primary alignment ----
                if args.paired:
                    r1, r2 = make_paired_primary(
                        rname, chrom, pos0, args.read_length, strand, header)
                    bam.write(r1)
                    bam.write(r2)
                else:
                    mapq = 0 if rng.random() < args.low_mapq_ratio else 60
                    r = make_primary(rname, chrom, pos0, args.read_length,
                                     strand, header, mapq=mapq)
                    bam.write(r)

                # Count only proper primary reads (not edge, MAPQ>0 unless low_mapq)
                if not is_edge:
                    true_count += 1

                # ---- noise: supplementary ----
                if rng.random() < args.suppl_ratio:
                    read_idx += 1
                    s = make_supplementary(
                        f"read_{sample_id}_{read_idx}",
                        chrom, pos0 + 5, args.read_length, strand, header)
                    bam.write(s)

                # ---- noise: secondary ----
                if rng.random() < args.secondary_ratio:
                    read_idx += 1
                    sec = make_secondary(
                        f"read_{sample_id}_{read_idx}",
                        chrom, pos0, args.read_length, strand, header)
                    bam.write(sec)

                # ---- noise: duplicates ----
                if rng.random() < args.duplicate_ratio:
                    read_idx += 1
                    dup = make_duplicate(
                        f"read_{sample_id}_{read_idx}",
                        chrom, pos0, args.read_length, strand, header)
                    bam.write(dup)

                # ---- noise: low-MAPQ extra copy ----
                if rng.random() < args.low_mapq_ratio:
                    read_idx += 1
                    lmq = make_low_mapq(
                        f"read_{sample_id}_{read_idx}",
                        chrom, pos0, args.read_length, strand, header)
                    bam.write(lmq)

            counts[feat_name] = true_count

            # ---- intergenic / unmapped noise (per-feature budget) ----
            n_unmapped = int(read_count * args.unmapped_ratio)
            for _ in range(n_unmapped):
                read_idx += 1
                u = make_unmapped(f"read_{sample_id}_{read_idx}", header)
                bam.write(u)

    # Sort and index
    sorted_path = bam_path.replace(".bam", ".sorted.bam")
    pysam.sort("-o", sorted_path, bam_path)
    os.replace(sorted_path, bam_path)
    pysam.index(bam_path)
    print(f"  [✓] {bam_path}  (reads per feature mean ~{args.reads_per_feature})")
    return counts


# ---------------------------------------------------------------------------
# Write ground-truth table
# ---------------------------------------------------------------------------

def write_table(outdir, all_counts, n_samples):
    """Write {outdir}/table.tsv with feature x sample matrix."""
    table_path = os.path.join(outdir, "table.tsv")
    sample_names = [f"sample{i+1}" for i in range(n_samples)]
    all_features = list(all_counts[0].keys())   # consistent order

    with open(table_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["feature"] + sample_names)
        for feat in all_features:
            row = [feat] + [all_counts[s].get(feat, 0) for s in range(n_samples)]
            writer.writerow(row)

    print(f"\n[✓] Ground-truth table → {table_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    rng  = random.Random(args.seed)

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading features from: {args.target}")
    features = load_features(args.target)
    if not features:
        sys.exit("ERROR: No features parsed — check file format.")
    print(f"  {len(features)} features loaded.")

    # Build BAM header from feature extents
    sq_list = build_sq_list(features)
    header  = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": sq_list,
        "PG": [{"ID": "simulate_bam", "PN": "simulate_bam", "VN": "1.0"}],
    })

    print(f"\nSimulating {args.num_samples} sample(s) → {args.outdir}/")
    print(f"  read_length={args.read_length}  paired={args.paired}")
    print(f"  noise: suppl={args.suppl_ratio}  secondary={args.secondary_ratio}  "
          f"dup={args.duplicate_ratio}  unmapped={args.unmapped_ratio}  "
          f"edge={args.edge_ratio}  low_mapq={args.low_mapq_ratio}\n")

    all_counts = []
    for i in range(args.num_samples):
        counts = simulate_sample(i + 1, features, args, header, rng)
        all_counts.append(counts)

    write_table(args.outdir, all_counts, args.num_samples)

    # Print featureCounts command hints
    bam_files = " ".join(
        os.path.join(args.outdir, f"sample{i+1}.bam")
        for i in range(args.num_samples)
    )
    ext = os.path.splitext(args.target)[1].lower()
    fc_fmt = "-F GTF" if ext == ".gtf" else ("-F GFF" if ext in (".gff", ".gff3") else "")
    paired_flag = "-p" if args.paired else ""

    print("\n── featureCounts command ──────────────────────────────────────────")
    print(f"featureCounts {fc_fmt} {paired_flag} \\")
    print(f"  -a {args.target} \\")
    print(f"  -o {os.path.join(args.outdir, 'fc_counts.txt')} \\")
    print(f"  {bam_files}")
    print("\n  Tips for noise flags:")
    print("  --primary          (ignore secondary alignments)")
    print("  --ignoreDup        (ignore duplicate-flagged reads)")
    print("  -Q 10              (minimum MAPQ to exclude low-mapq reads)")
    print("  --fracOverlap 0.5  (require 50 %% of read to overlap feature)")
    print("  --minOverlap 10    (min bp overlap — useful for edge reads)")
    print("────────────────────────────────────────────────────────────────────\n")


if __name__ == "__main__":
    main()
