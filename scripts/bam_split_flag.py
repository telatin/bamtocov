#!/usr/bin/env python3
import argparse
import os
import sys
from collections import defaultdict

import pysam


def pair_flag(flag: int) -> int:
    """
    Given a SAM FLAG integer for one read, return the FLAG you would expect
    on its mate by swapping:
      - READ1 (0x40) <-> READ2 (0x80)
      - REVERSE (0x10) <-> MREVERSE (0x20)
      - UNMAP (0x04) <-> MUNMAP (0x08)

    Other bits are preserved.
    """
    if not (flag & 0x1):
        # Unpaired read: no mate to swap with; return unchanged
        return flag

    swaps = [(0x40, 0x80),  # READ1 <-> READ2
             (0x10, 0x20),  # REVERSE <-> MREVERSE
             (0x04, 0x08)]  # UNMAP <-> MUNMAP

    out = flag
    for a, b in swaps:
        a_set = bool(out & a)
        b_set = bool(out & b)
        out &= ~(a | b)
        if a_set:
            out |= b
        if b_set:
            out |= a
    return out


def resolve_bucket_key(flag: int) -> int:
    """
    For paired reads: use the READ1 version of the flag as the bucket key,
    so both mates land in the same file (e.g., 99/147 -> key 99).
    For unpaired reads: use the flag as-is.
    """
    if not (flag & 0x1):
        return flag
    # If this is READ1, keep it; if READ2, map to mate (which will be READ1).
    return flag if (flag & 0x40) else pair_flag(flag)


def split_bam_by_flag(in_bam_path: str, out_dir: str, progress_step: int = 100000) -> None:
    if not os.path.exists(in_bam_path):
        raise FileNotFoundError(f"Input BAM not found: {in_bam_path}")

    os.makedirs(out_dir, exist_ok=True)

    base = os.path.splitext(os.path.basename(in_bam_path))[0]

    writers = {}  # key_flag -> pysam.AlignmentFile
    total = 0
    counts = defaultdict(int)

    with pysam.AlignmentFile(in_bam_path, "rb") as bam_in:
        header = bam_in.header
        # Include all records even if unsorted/unmapped
        for rec in bam_in.fetch(until_eof=True):
            flag = rec.flag
            key = resolve_bucket_key(flag)

            if key not in writers:
                out_path = os.path.join(out_dir, f"{base}_{key}.bam")
                writers[key] = pysam.AlignmentFile(out_path, "wb", header=header)

            writers[key].write(rec)
            counts[key] += 1
            # progress every N reads
            if progress_step and total % progress_step == 0:
                sys.stderr.write(f"Processed {total:,} reads...\r")

    # Close all writers
    for w in writers.values():
        w.close()

    # Print a short summary to stderr
    sys.stderr.write(f"Wrote {len(writers)} BAM file(s) to {out_dir}\n")
    for key in sorted(counts):
        sys.stderr.write(f"  {base}_{key}.bam : {counts[key]} records\n")


def main():
    ap = argparse.ArgumentParser(
        description="Split a BAM into multiple BAMs by FLAG, grouping mates "
                    "by the first-in-pair flag (READ1)."
    )
    ap.add_argument("-i", "--input", required=True, dest="input_bam",
                    help="Input BAM file")
    ap.add_argument("-o", "--outdir", dest="outdir", default=None,
                    help="Output directory (default: same directory as input)")
    ap.add_argument("--progress", type=int, default=100000,
                    help="Print progress every N reads (default: 100000, 0 disables)")
    args = ap.parse_args()

    in_dir = os.path.dirname(os.path.abspath(args.input_bam))
    out_dir = args.outdir if args.outdir is not None else in_dir
    split_bam_by_flag(args.input_bam, out_dir, args.progress)


if __name__ == "__main__":
    main()
