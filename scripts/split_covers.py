#!/usr/bin/env python3
import argparse
import csv
import gzip
import io
import os
import re
import sys
from typing import Dict, List, Tuple

def open_maybe_gzip(path: str, mode: str = "rt", **kwargs):
    """
    Open a file that may be gzipped based on its extension.
    Text mode by default; pass mode='rb' for bytes.
    """
    if path.endswith(".gz"):
        # For text mode, wrap with TextIOWrapper to ensure newline handling
        if "b" in mode:
            return gzip.open(path, mode)
        gz = gzip.open(path, "rb")
        return io.TextIOWrapper(gz, encoding=kwargs.get("encoding", "utf-8"), newline=kwargs.get("newline", ""))
    return open(path, mode, **kwargs)

def normalize_metric_name(metric: str) -> str:
    """
    Convert header metric names to snake_case filenames.
    Special-cases common CoverM headers to the shorter names requested.
    """
    m = metric.strip().lower()
    # Map known CoverM labels to requested forms
    mapping = {
        "read count": "count",
        "mean": "mean",
        "tpm": "tpm",
        "rpkm": "rpkm",
        "covered fraction": "covered_fraction",
        "covered bases": "covered_bases",
    }
    if m in mapping:
        return mapping[m]
    # Generic fallback: snake_case, alnum + underscore only
    m = m.replace("/", " ").replace("-", " ")
    m = re.sub(r"\s+", "_", m)
    m = re.sub(r"[^a-z0-9_]", "", m)
    return m

def parse_header(header: List[str]) -> Tuple[List[str], List[str], Dict[str, List[int]]]:
    """
    Parse the TSV header.

    Returns:
        samples: list of sample names in the order encountered
        metrics: list of unique metrics (normalized), in the order encountered
        metric_to_indices: map metric_name -> list of column indices for that metric across samples (sample order)
    """
    if not header or header[0] != "Contig":
        raise ValueError("First column must be 'Contig'.")

    samples_order: List[str] = []
    metrics_order: List[str] = []
    metric_to_indices: Dict[str, List[int]] = {}

    # We’ll also keep track of which sample order we’ve established as we see them
    seen_samples = set()

    # From col 1 onward: "<sample> <metric label>"
    for idx, col in enumerate(header[1:], start=1):
        col = col.strip()
        if not col:
            continue
        if " " not in col:
            # If the column doesn't have a space, we can't split sample vs metric reliably
            # Treat the whole thing as a metric with an implicit single sample
            sample_name, raw_metric = "sample", col
        else:
            sample_name, raw_metric = col.split(" ", 1)

        if sample_name not in seen_samples:
            samples_order.append(sample_name)
            seen_samples.add(sample_name)

        metric_norm = normalize_metric_name(raw_metric)
        if metric_norm not in metrics_order:
            metrics_order.append(metric_norm)
        metric_to_indices.setdefault(metric_norm, []).append(idx)

    # Sanity check: each metric should have as many indices as samples
    for m in metrics_order:
        if len(metric_to_indices[m]) != len(samples_order):
            # Not fatal; warn to stderr
            print(f"Warning: metric '{m}' has {len(metric_to_indices[m])} columns but {len(samples_order)} samples.",
                  file=sys.stderr)

    return samples_order, metrics_order, metric_to_indices

def main():
    ap = argparse.ArgumentParser(
        description="Split a CoverM TSV by metric into separate TSVs, one per metric."
    )
    ap.add_argument("-i", "--input", dest="input_file", required=True,
                    help="Input CoverM TSV (optionally gzipped: .gz).")
    ap.add_argument("-o", "--output-basename", dest="basename", required=True,
                    help="Output basename/prefix for generated files.")
    args = ap.parse_args()

    # Read header first
    with open_maybe_gzip(args.input_file, mode="rt", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            print("Error: empty input file.", file=sys.stderr)
            sys.exit(1)

        samples, metrics, metric_to_indices = parse_header(header)

        # Prepare writers: one file per metric
        writers: Dict[str, csv.writer] = {}
        files: Dict[str, io.TextIOBase] = {}

        try:
            for metric in metrics:
                out_path = f"{args.basename}_{metric}.tsv"
                f = open(out_path, "w", newline="", encoding="utf-8")
                files[metric] = f
                w = csv.writer(f, delimiter="\t", lineterminator="\n")
                writers[metric] = w
                # Header: Contig + one column per sample
                w.writerow(["Contig"] + samples)

            # Process rows and write one line per metric
            for row in reader:
                if not row:
                    continue
                contig = row[0]
                for metric in metrics:
                    indices = metric_to_indices.get(metric, [])
                    values = []
                    for idx in indices:
                        # If row is shorter than expected, pad with empty
                        values.append(row[idx] if idx < len(row) else "")
                    writers[metric].writerow([contig] + values)

        finally:
            for f in files.values():
                try:
                    f.close()
                except Exception:
                    pass

    # Report what we wrote
    for metric in metrics:
        print(f"Wrote {args.basename}_{metric}.tsv")

if __name__ == "__main__":
    main()
