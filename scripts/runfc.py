#!/usr/bin/env python3
"""
fc_get.py — featureCounts wrapper with GTF auto-detection, BAM auto-indexing,
             and paired-end detection.
"""

import argparse
import collections
import os
import re
import subprocess
import sys
import tempfile


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def eprint(*args, **kwargs):
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def run(cmd, dry_run=False, check=True):
    """Print and optionally execute a shell command (list form)."""
    eprint(f"[cmd] {' '.join(str(c) for c in cmd)}")
    if dry_run:
        return None
    result = subprocess.run(cmd, check=check)
    return result


# ---------------------------------------------------------------------------
# GTF auto-detection
# ---------------------------------------------------------------------------

def detect_gtf_feature_and_id(gtf_path, max_lines=100):
    """
    Scan up to *max_lines* non-comment lines of a GTF and return
    (feature_type, id_attribute) as the most common values found.
    """
    feature_counter = collections.Counter()
    attr_counter = collections.Counter()

    lines_read = 0
    opener = open  # plain text; add gzip support if needed
    if gtf_path.endswith(".gz"):
        import gzip
        opener = gzip.open

    with opener(gtf_path, "rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            feature_counter[parts[2]] += 1
            # parse attribute column: key "value"; pairs
            for m in re.finditer(r'(\w+)\s+"[^"]+"', parts[8]):
                attr_counter[m.group(1)] += 1
            lines_read += 1
            if lines_read >= max_lines:
                break

    if not feature_counter:
        eprint("[warn] Could not detect feature type from GTF; falling back to 'exon'")
        feature_type = "exon"
    else:
        feature_type = feature_counter.most_common(1)[0][0]

    # Prefer gene_id / transcript_id / ID in that priority if present,
    # otherwise take the most common attribute key.
    preferred = ["gene_id", "transcript_id", "ID", "gene_name"]
    id_attr = None
    for p in preferred:
        if p in attr_counter:
            id_attr = p
            break
    if id_attr is None:
        if attr_counter:
            id_attr = attr_counter.most_common(1)[0][0]
        else:
            id_attr = "gene_id"

    eprint(f"[info] GTF auto-detected: feature type = '{feature_type}', "
           f"id attribute = '{id_attr}'")
    return feature_type, id_attr


# ---------------------------------------------------------------------------
# BAM helpers
# ---------------------------------------------------------------------------

def ensure_indexed(bam_path, dry_run=False):
    """Index the BAM file if .bai is missing."""
    bai = bam_path + ".bai"
    if not os.path.exists(bai):
        eprint(f"[info] Index not found for {bam_path}, running samtools index …")
        run(["samtools", "index", bam_path], dry_run=dry_run)
    else:
        eprint(f"[info] Index found: {bai}")


def is_paired(bam_path, dry_run=False, n_reads=1000):
    if dry_run:
        eprint("[info] Dry-run: assuming single-end for command preview")
        return False
    try:
        import pysam
    except ImportError:
        eprint("[warn] pysam not available, falling back to samtools flag check")
        return _is_paired_samtools(bam_path, n_reads)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        paired_count = 0
        total = 0
        for read in bam.head(n_reads):
            if read.flag & 0x1:
                paired_count += 1
            total += 1

    paired = paired_count > 0
    eprint(f"[info] {bam_path}: paired={'yes' if paired else 'no'} "
           f"({paired_count}/{total} reads flagged as paired)")
    return paired


def _is_paired_samtools(bam_path, n_reads=1000):
    """Fallback using samtools view -f 1 (portable, no awk needed)."""
    cmd = ["samtools", "view", "-f", "1", "-c", bam_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    paired_count = int(result.stdout.strip() or 0)
    paired = paired_count > 0
    eprint(f"[info] {bam_path}: paired={'yes' if paired else 'no'} "
           f"({paired_count} reads with paired flag)")
    return paired
# ---------------------------------------------------------------------------
# Output filtering
# ---------------------------------------------------------------------------

def filter_counts(raw_path, out_fh):
    skip_cols = {"Chr", "Start", "End", "Strand", "Length"}
    keep_idx = None

    with open(raw_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if keep_idx is None:
                keep_idx = [i for i, h in enumerate(parts) if h.strip() not in skip_cols]
            out_fh.write("\t".join(parts[i] for i in keep_idx) + "\n")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Wrapper around featureCounts with GTF auto-detection, "
                    "BAM indexing, and paired-end handling.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-b", "--bam", metavar="BAMFILE", nargs="+", required=True,
                   help="One or more BAM files")
    p.add_argument("-a", "--annotation", metavar="FILE", required=True,
                   help="GTF annotation file")
    p.add_argument("-f", "--feature", metavar="FEATURE",
                   help="Feature type to count (e.g. exon, CDS). "
                        "Auto-detected from GTF if omitted.")
    p.add_argument("-t", "--type", metavar="TYPE", dest="id_attr",
                   help="Attribute used as feature ID (e.g. gene_id). "
                        "Auto-detected from GTF if omitted.")
    p.add_argument("-o", "--output", metavar="FILE", default=None,
                   help="Output counts file (default: stdout)")
    p.add_argument("--dry-run", action="store_true",
                   help="Print commands without executing them")
    p.add_argument("-T", "--threads", type=int, default=1,
                   help="Number of threads for featureCounts")
    return p.parse_args()


def main():
    args = parse_args()

    # ------------------------------------------------------------------
    # 1. Resolve feature type and id attribute
    # ------------------------------------------------------------------
    if args.feature and args.id_attr:
        feature_type = args.feature
        id_attr = args.id_attr
        eprint(f"[info] Using provided: feature = '{feature_type}', "
               f"id attribute = '{id_attr}'")
    else:
        feature_type, id_attr = detect_gtf_feature_and_id(args.annotation)
        # Allow partial override
        if args.feature:
            feature_type = args.feature
            eprint(f"[info] Feature type overridden by -f: '{feature_type}'")
        if args.id_attr:
            id_attr = args.id_attr
            eprint(f"[info] ID attribute overridden by -t: '{id_attr}'")

    # ------------------------------------------------------------------
    # 2. Index BAMs if needed
    # ------------------------------------------------------------------
    for bam in args.bam:
        ensure_indexed(bam, dry_run=args.dry_run)

    # ------------------------------------------------------------------
    # 3. Detect paired-end (use first BAM as representative)
    # ------------------------------------------------------------------
    paired = is_paired(args.bam[0], dry_run=args.dry_run)

    # ------------------------------------------------------------------
    # 4. Run featureCounts
    # ------------------------------------------------------------------
    with tempfile.NamedTemporaryFile(
        suffix=".featureCounts.txt", delete=False, mode="w"
    ) as tmp:
        tmp_path = tmp.name

    fc_cmd = [
        "featureCounts",
        "-a", args.annotation,
        "-o", tmp_path,
        "-t", feature_type,
        "-g", id_attr,
        "-T", str(args.threads),
    ]
    if paired:
        fc_cmd.append("-p")

    fc_cmd += args.bam

    run(fc_cmd, dry_run=args.dry_run)

    # ------------------------------------------------------------------
    # 5. Filter and emit output
    # ------------------------------------------------------------------
    if args.dry_run:
        eprint("[info] Dry-run complete; no output written.")
        return

    if args.output:
        with open(args.output, "w") as out_fh:
            filter_counts(tmp_path, out_fh)
        eprint(f"[info] Counts written to {args.output}")
    else:
        filter_counts(tmp_path, sys.stdout)

    # Clean up temp file
    try:
        os.unlink(tmp_path)
        summary = tmp_path + ".summary"
        if os.path.exists(summary):
            os.unlink(summary)
    except OSError:
        pass


if __name__ == "__main__":
    main()
