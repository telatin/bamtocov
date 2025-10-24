#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
simulate_bam_coverage.py
Generate synthetic BAM files to evaluate coverage over target regions.

Features
- Per-sample chromosome allocation via Dirichlet (length-weighted).
- Global on/off-target fraction per sample, with on-target reads distributed across features per chromosome.
- Single-end or paired-end (FR) reads with fixed read length and insert size.
- Perfect-match CIGARs, fixed MAPQ, coordinate-sorted + indexed BAMs.
- Outputs per-feature counts table across all samples.

Inputs
- Genome TSV: chrom\tlength
- BED targets: chr  start  end  [name]
"""

from __future__ import annotations
import argparse
import os
import sys
import math
import json
import random
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional

import numpy as np
from rich.console import Console
from rich.table import Table
from rich.progress import track
from rich import box
import pysam

console = Console()

# ----------------------------- Helpers & Data ----------------------------- #

@dataclass
class Feature:
    chrom: str
    start: int
    end: int
    name: str  # unique name assigned

def add_required_tags(rec: pysam.AlignedSegment, sample_name: str, rg_id: str, read_len: int):
    """
    Attach RG/SM, and coverage-critical tags NM and MD.
    NM=0 (no edits), MD=<read_len> (all matches).
    """
    # Collect any existing tags you want to keep; here we overwrite RG/SM/NM/MD explicitly
    tags = {k: v for k, v in rec.get_tags() if k not in {"RG", "SM", "NM", "MD"}}
    tags.update({
        "RG": rg_id,
        "SM": sample_name,
        "NM": 0,
        "MD": str(read_len),
    })
    rec.set_tags(list(tags.items()))

def parse_genome_tsv(path: str) -> Dict[str, int]:
    genome = {}
    with open(path, "r") as fh:
        for ln, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"[genome] Line {ln}: expected 2 columns, got {len(parts)}")
            chrom, length_s = parts[0], parts[1]
            try:
                length = int(length_s)
            except ValueError:
                raise ValueError(f"[genome] Line {ln}: length is not an integer: {length_s}")
            if length <= 0:
                raise ValueError(f"[genome] Line {ln}: length must be >0")
            if chrom in genome:
                console.print(f"[yellow]Warning:[/yellow] duplicate chromosome '{chrom}' in genome TSV; overriding previous value.")
            genome[chrom] = length
    if not genome:
        raise ValueError("Genome TSV is empty.")
    return genome

def parse_bed(path: str, genome: Dict[str, int]) -> List[Feature]:
    features: List[Feature] = []
    used_names: Dict[str, int] = {}
    with open(path, "r") as fh:
        for ln, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"[BED] Line {ln}: expected at least 3 columns (chr start end)")
            chrom, start_s, end_s = parts[0], parts[1], parts[2]
            if chrom not in genome:
                raise ValueError(f"[BED] Line {ln}: chromosome '{chrom}' not found in genome TSV.")
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                raise ValueError(f"[BED] Line {ln}: start or end not integer.")
            if not (0 <= start <= end <= genome[chrom]):
                raise ValueError(f"[BED] Line {ln}: interval {start}-{end} out of bounds for {chrom} (len {genome[chrom]}).")
            base_name = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"
            # ensure unique names (keep distinct overlapping features)
            cnt = used_names.get(base_name, 0)
            if cnt == 0:
                name = base_name
            else:
                name = f"{base_name}#{cnt+1}"
            used_names[base_name] = cnt + 1
            features.append(Feature(chrom=chrom, start=start, end=end, name=name))
    return features

def build_complements(genome: Dict[str, int], feats_by_chr: Dict[str, List[Feature]]) -> Dict[str, List[Tuple[int, int]]]:
    """
    For each chromosome, compute non-target intervals as a list of (start, end) 0-based half-open.
    Merge overlapping features first per chromosome to get a simple complement.
    """
    complements: Dict[str, List[Tuple[int, int]]] = {}
    for chrom, length in genome.items():
        feats = feats_by_chr.get(chrom, [])
        if not feats:
            complements[chrom] = [(0, length)]
            continue
        # merge intervals
        ivals = sorted([(f.start, f.end) for f in feats])
        merged = []
        cur_s, cur_e = ivals[0]
        for s, e in ivals[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                merged.append((cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((cur_s, cur_e))
        comp = []
        prev = 0
        for s, e in merged:
            if prev < s:
                comp.append((prev, s))
            prev = e
        if prev < length:
            comp.append((prev, length))
        complements[chrom] = comp
    return complements

def parse_tot(arg: str) -> int:
    """
    Parse total reads/fragments: INT or INT with suffix K/M/G (case-insensitive).
    """
    s = arg.strip().lower().replace("_", "")
    mult = 1
    if s.endswith("k"):
        mult = 1_000
        s = s[:-1]
    elif s.endswith("m"):
        mult = 1_000_000
        s = s[:-1]
    elif s.endswith("g"):
        mult = 1_000_000_000
        s = s[:-1]
    try:
        base = int(s)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid --tot value: {arg}")
    if base <= 0:
        raise argparse.ArgumentTypeError("--tot must be > 0")
    return base * mult

def choose_start_in_region(region: Tuple[int, int], needed_len: int, rng: np.random.Generator) -> Optional[int]:
    """
    Return a random start so that [start, start+needed_len) fits within region (half-open).
    If region too short, return None.
    """
    rlen = region[1] - region[0]
    if rlen < needed_len:
        return None
    max_start = region[1] - needed_len
    start = int(rng.integers(region[0], max_start + 1))
    return start

def multinomial_int(total: int, probs: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Integer allocation using multinomial; handles edge cases for tiny totals."""
    if total <= 0:
        return np.zeros_like(probs, dtype=int)
    probs = probs / probs.sum()
    return rng.multinomial(total, probs)

def dirichlet_weights(alpha: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Safely sample Dirichlet weights; fall back to uniform if alpha invalid."""
    alpha = np.array(alpha, dtype=float)
    if not np.all(alpha > 0) or np.isnan(alpha).any() or np.isinf(alpha).any():
        alpha = np.ones_like(alpha)
    return rng.dirichlet(alpha)

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def bam_header_from_genome(genome: Dict[str, int]) -> Dict:
    sq = [{"LN": ln, "SN": chrom} for chrom, ln in genome.items()]
    return {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}

def build_tid_map(header: Dict) -> Dict[str, int]:
    return {sq["SN"]: i for i, sq in enumerate(header["SQ"])}

# ----------------------------- Core Simulation ----------------------------- #

@dataclass
class SamplePlan:
    sample_name: str
    total_fragments: int
    per_chr_total: Dict[str, int]
    per_chr_on: Dict[str, int]
    per_chr_off: Dict[str, int]
    per_chr_feat_alloc: Dict[str, Dict[str, int]]  # chrom -> feature_name -> counts

def allocate_reads(
    genome: Dict[str, int],
    feats_by_chr: Dict[str, List[Feature]],
    total_fragments: int,
    target_frac: float,
    rng: np.random.Generator,
) -> SamplePlan:
    chroms = list(genome.keys())
    lengths = np.array([genome[c] for c in chroms], dtype=float)
    # Dirichlet alpha proportional to length, scaled for mild smoothness
    alpha = (lengths / lengths.sum()) * 10.0 + 1e-6
    w_chr = dirichlet_weights(alpha, rng)
    per_chr_total_counts = multinomial_int(total_fragments, w_chr, rng)
    per_chr_total = {c: int(n) for c, n in zip(chroms, per_chr_total_counts)}

    total_on = int(round(target_frac * total_fragments))
    total_off = total_fragments - total_on

    # Distribute "on" only to chromosomes that actually have features
    chroms_with_feats = [c for c in chroms if len(feats_by_chr.get(c, [])) > 0]
    if chroms_with_feats:
        weights_on = np.array([per_chr_total[c] for c in chroms_with_feats], dtype=float)
        if weights_on.sum() == 0:
            weights_on = np.ones(len(chroms_with_feats))
        per_chr_on_counts = multinomial_int(total_on, weights_on / weights_on.sum(), rng)
        per_chr_on = {c: 0 for c in chroms}
        for c, n in zip(chroms_with_feats, per_chr_on_counts):
            per_chr_on[c] = int(n)
    else:
        # No features anywhere -> all off-target
        per_chr_on = {c: 0 for c in chroms}
        total_off = total_fragments
        total_on = 0

    # Distribute "off" across all chromosomes, proportional to their totals
    weights_off = np.array([per_chr_total[c] for c in chroms], dtype=float)
    if weights_off.sum() == 0:
        weights_off = np.ones(len(chroms))
    per_chr_off_counts = multinomial_int(total_off, weights_off / weights_off.sum(), rng)
    per_chr_off = {c: int(n) for c, n in zip(chroms, per_chr_off_counts)}

    # Within each chromosome, distribute on-target across its features via Dirichlet
    per_chr_feat_alloc: Dict[str, Dict[str, int]] = {}
    for c in chroms:
        n_on = per_chr_on.get(c, 0)
        feats = feats_by_chr.get(c, [])
        if n_on == 0 or not feats:
            per_chr_feat_alloc[c] = {}
            continue
        alpha_f = np.ones(len(feats))  # equal Dirichlet among features
        wf = dirichlet_weights(alpha_f, rng)
        alloc = multinomial_int(n_on, wf, rng)
        per_chr_feat_alloc[c] = {f.name: int(n) for f, n in zip(feats, alloc)}

    return SamplePlan(
        sample_name="",
        total_fragments=total_fragments,
        per_chr_total=per_chr_total,
        per_chr_on=per_chr_on,
        per_chr_off=per_chr_off,
        per_chr_feat_alloc=per_chr_feat_alloc,
    )

def simulate_sample_bam(
    out_bam: str,
    header: Dict,
    tid_map: Dict[str, int],
    genome: Dict[str, int],
    feats_by_chr: Dict[str, List[Feature]],
    complements: Dict[str, List[Tuple[int, int]]],
    plan: SamplePlan,
    read_len: int,
    paired_end: bool,
    insert_size: int,
    sample_name: str,
    mapq: int,
    rng: np.random.Generator,
    rg_id: str,
):
    """
    Generate alignments according to plan and write a coordinate-sorted BAM.
    """
    # Validate insert size
    if paired_end and insert_size < 2 * read_len:
        console.print(f"[yellow]Warning:[/yellow] --ins ({insert_size}) < 2*read_len; adjusted to {2*read_len}.")
        insert_size = 2 * read_len

    alignments = []  # list of pysam.AlignedSegment to be sorted and written

    cigar = [(0, read_len)]  # 0 = BAM_CMATCH

    for chrom, total_chr in plan.per_chr_total.items():
        # On-target
        feat_counts = plan.per_chr_feat_alloc.get(chrom, {})
        for feat_name, n in feat_counts.items():
            if n <= 0:
                continue
            # find feature interval
            feat = next(f for f in feats_by_chr[chrom] if f.name == feat_name)
            needed = insert_size if paired_end else read_len
            for _ in range(n):
                start = choose_start_in_region((feat.start, feat.end), needed, rng)
                if start is None:
                    # if feature too short, skip placing here (rare) -> try a few fallback tries, else drop
                    fallback = False
                    for _try in range(5):
                        start = choose_start_in_region((feat.start, feat.end), needed, rng)
                        if start is not None:
                            fallback = True
                            break
                    if not fallback:
                        continue
                # Build one fragment (1 or 2 reads)
                if paired_end:
                    r1 = pysam.AlignedSegment()
                    r2 = pysam.AlignedSegment()
                    r1.query_name = f"{sample_name}:on:{chrom}:{feat_name}:{rng.integers(1e12)}"
                    r2.query_name = r1.query_name
                    r1.query_sequence = "A" * read_len
                    r2.query_sequence = "A" * read_len
                    r1.flag = 99   # read1, proper pair, mate reverse, FR
                    r2.flag = 147  # read2, proper pair, reverse, mate forward
                    r1.reference_id = tid_map[chrom]
                    r2.reference_id = tid_map[chrom]
                    r1.reference_start = start
                    r2.reference_start = start + insert_size - read_len
                    r1.mapping_quality = mapq
                    r2.mapping_quality = mapq
                    r1.cigar = cigar
                    r2.cigar = cigar
                    # mate info
                    r1.next_reference_id = r2.reference_id
                    r2.next_reference_id = r1.reference_id
                    r1.next_reference_start = r2.reference_start
                    r2.next_reference_start = r1.reference_start
                    tlen = (r2.reference_start + read_len) - r1.reference_start
                    r1.template_length = tlen
                    r2.template_length = -tlen
                    # tags
                    #r1.set_tags([("RG", rg_id), ("SM", sample_name)])
                    #r2.set_tags([("RG", rg_id), ("SM", sample_name)])
                    add_required_tags(r1, sample_name, rg_id, read_len)
                    add_required_tags(r2, sample_name, rg_id, read_len)
                    alignments.extend([r1, r2])
                else:
                    r = pysam.AlignedSegment()
                    r.query_name = f"{sample_name}:on:{chrom}:{feat_name}:{rng.integers(1e12)}"
                    r.query_sequence = "A" * read_len
                    # Random strand 50/50
                    is_rev = bool(rng.integers(0, 2))
                    r.flag = 16 if is_rev else 0
                    r.reference_id = tid_map[chrom]
                    r.reference_start = start
                    r.mapping_quality = mapq
                    r.cigar = cigar
                    #r.set_tags([("RG", rg_id), ("SM", sample_name)])
                    add_required_tags(r, sample_name, rg_id, read_len)
                    alignments.append(r)

        # Off-target
        n_off = plan.per_chr_off.get(chrom, 0)
        if n_off <= 0:
            continue
        regions = complements.get(chrom, [])
        if not regions:
            # Nowhere to place off-target on this chromosome; skip
            continue
        needed = insert_size if paired_end else read_len
        # Precompute region lengths for weighted picking
        reg_lengths = np.array([max(0, e - s) for s, e in regions], dtype=float)
        eligible_idx = np.where(reg_lengths >= needed)[0]
        if len(eligible_idx) == 0:
            continue
        w = reg_lengths[eligible_idx] / reg_lengths[eligible_idx].sum()
        for _ in range(n_off):
            idx = int(rng.choice(eligible_idx, p=w))
            region = regions[idx]
            start = choose_start_in_region(region, needed, rng)
            if start is None:
                # try a couple retries across regions
                placed = False
                for _try in range(5):
                    idx2 = int(rng.choice(eligible_idx, p=w))
                    region2 = regions[idx2]
                    start = choose_start_in_region(region2, needed, rng)
                    if start is not None:
                        placed = True
                        break
                if not placed:
                    continue
            if paired_end:
                r1 = pysam.AlignedSegment()
                r2 = pysam.AlignedSegment()
                r1.query_name = f"{sample_name}:off:{chrom}:{rng.integers(1e12)}"
                r2.query_name = r1.query_name
                r1.query_sequence = "A" * read_len
                r2.query_sequence = "A" * read_len
                r1.flag = 99
                r2.flag = 147
                r1.reference_id = tid_map[chrom]
                r2.reference_id = tid_map[chrom]
                r1.reference_start = start
                r2.reference_start = start + insert_size - read_len
                r1.mapping_quality = mapq
                r2.mapping_quality = mapq
                r1.cigar = cigar
                r2.cigar = cigar
                r1.next_reference_id = r2.reference_id
                r2.next_reference_id = r1.reference_id
                r1.next_reference_start = r2.reference_start
                r2.next_reference_start = r1.reference_start
                tlen = (r2.reference_start + read_len) - r1.reference_start
                r1.template_length = tlen
                r2.template_length = -tlen
                add_required_tags(r1, sample_name, rg_id, read_len)
                add_required_tags(r2, sample_name, rg_id, read_len)
                alignments.extend([r1, r2])
            else:
                r = pysam.AlignedSegment()
                r.query_name = f"{sample_name}:off:{chrom}:{rng.integers(1e12)}"
                r.query_sequence = "A" * read_len
                is_rev = bool(rng.integers(0, 2))
                r.flag = 16 if is_rev else 0
                r.reference_id = tid_map[chrom]
                r.reference_start = start
                r.mapping_quality = mapq
                r.cigar = cigar
                add_required_tags(r, sample_name, rg_id, read_len)
                alignments.append(r)

    # Sort alignments by (tid, pos, flag to keep read1 before read2)
    alignments.sort(key=lambda a: (a.reference_id, a.reference_start, a.flag))

    # Write BAM
    tmp_bam = out_bam + ".unsorted.tmp.bam"
    with pysam.AlignmentFile(tmp_bam, "wb", header=header) as outf:
        for a in alignments:
            outf.write(a)

    # Sort & index (pysam.sort creates BAM)
    pysam.sort("-o", out_bam, tmp_bam)
    os.remove(tmp_bam)
    try:
        pysam.index(out_bam)
    except Exception as e:
        console.print(f"[yellow]Warning:[/yellow] Could not index {out_bam}: {e}")

# ----------------------------- Counts Aggregation ----------------------------- #

def init_counts(features: List[Feature], num_samples: int) -> Tuple[List[str], Dict[str, List[int]]]:
    names = [f.name for f in features]
    names.append("offtarget")
    table = {name: [0] * num_samples for name in names}
    return names, table

def update_counts_for_plan(table: Dict[str, List[int]], sample_idx: int, plan: SamplePlan):
    # On-target features
    for chrom, feat_alloc in plan.per_chr_feat_alloc.items():
        for feat_name, n in feat_alloc.items():
            table[feat_name][sample_idx] += int(n)
    # Off-target (global)
    off_total = sum(plan.per_chr_off.values())
    table["offtarget"][sample_idx] += int(off_total)

def write_counts_table(path: str, feature_names: List[str], num_samples: int, table: Dict[str, List[int]]):
    with open(path, "w") as fh:
        header = ["feature"] + [f"sample_{i+1}" for i in range(num_samples)]
        fh.write("\t".join(header) + "\n")
        for fname in feature_names:
            row = [fname] + [str(table[fname][i]) for i in range(num_samples)]
            fh.write("\t".join(row) + "\n")

# ----------------------------- CLI ----------------------------- #

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate synthetic BAM datasets with on-/off-target reads to evaluate coverage.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-g", "--genome", required=True, help="Genome TSV file: chrom<tab>length")
    p.add_argument("-t", "--target", required=True, help="BED file with target intervals (0-based, half-open).")
    p.add_argument("-o", "--outdir", required=True, help="Output directory")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("-f", "--target-frac", dest="target_frac", type=float, default=0.8,
                   help="Fraction of reads (fragments) on target per sample (0..1)")
    p.add_argument("--tot", type=parse_tot, default=parse_tot("1M"),
                   help="Total number of reads/fragments per sample; INT or INT+suffix (e.g., 250K, 2M)")
    p.add_argument("-n", "--num_samples", type=int, required=True, help="Number of samples to generate")
    p.add_argument("-l", "--len", dest="read_len", type=int, default=100, help="Read length")
    p.add_argument("--pe", action="store_true", help="Generate paired-end reads in FR orientation")
    p.add_argument("--ins", type=int, default=300, help="Insert size (fragment length) for paired-end")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return p

def main():
    args = build_argparser().parse_args()

    if not (0.0 <= args.target_frac <= 1.0):
        console.print("[red]Error:[/red] --target must be in [0,1].", style="bold")
        sys.exit(2)
    if args.read_len <= 0:
        console.print("[red]Error:[/red] Read length must be > 0.", style="bold")
        sys.exit(2)
    if args.num_samples <= 0:
        console.print("[red]Error:[/red] --num_samples must be > 0.", style="bold")
        sys.exit(2)

    ensure_dir(args.outdir)

    if args.verbose:
        console.rule("[bold]Simulation Parameters")
        t = Table(box=box.SIMPLE_HEAVY)
        t.add_column("Param", style="bold")
        t.add_column("Value")
        t.add_row("Genome TSV", args.genome)
        t.add_row("BED targets", args.target)
        t.add_row("Outdir", args.outdir)
        t.add_row("Seed", str(args.seed))
        t.add_row("Target fraction", f"{args.target_frac:.3f}")
        t.add_row("Total per sample", str(args.tot))
        t.add_row("Num samples", str(args.num_samples))
        t.add_row("Read length", str(args.read_len))
        t.add_row("Paired-end", "yes" if args.pe else "no")
        if args.pe:
            t.add_row("Insert size", str(args.ins))
        console.print(t)

    try:
        genome = parse_genome_tsv(args.genome)
    except Exception as e:
        console.print(f"[red]Error parsing genome TSV:[/red] {e}")
        sys.exit(2)

    try:
        features = parse_bed(args.target, genome)
    except Exception as e:
        console.print(f"[red]Error parsing BED:[/red] {e}")
        sys.exit(2)

    feats_by_chr: Dict[str, List[Feature]] = {}
    for f in features:
        feats_by_chr.setdefault(f.chrom, []).append(f)

    complements = build_complements(genome, feats_by_chr)

    header = bam_header_from_genome(genome)
    tid_map = build_tid_map(header)

    # Prepare counts matrix
    feat_names, counts_table = init_counts(features, args.num_samples)

    # RNG
    base_rng = np.random.default_rng(args.seed)

    # Generate samples
    for i in track(range(1, args.num_samples + 1), description="Generating samples"):
        sample_name = f"sample_{i}"
        rng = np.random.default_rng(base_rng.integers(0, 2**63 - 1))
        plan = allocate_reads(
            genome=genome,
            feats_by_chr=feats_by_chr,
            total_fragments=args.tot,
            target_frac=args.target_frac,
            rng=rng,
        )
        plan.sample_name = sample_name

        # Write BAM
        out_bam = os.path.join(args.outdir, f"{sample_name}.bam")
        rg_id = sample_name
        simulate_sample_bam(
            out_bam=out_bam,
            header=header,
            tid_map=tid_map,
            genome=genome,
            feats_by_chr=feats_by_chr,
            complements=complements,
            plan=plan,
            read_len=args.read_len,
            paired_end=args.pe,
            insert_size=args.ins,
            sample_name=sample_name,
            mapq=60,
            rng=rng,
            rg_id=rg_id,
        )

        # Update counts
        update_counts_for_plan(counts_table, i - 1, plan)

        # Write per-sample meta
        meta = {
            "sample": sample_name,
            "seed": int(rng.integers(0, 2**31 - 1)),
            "total_fragments": plan.total_fragments,
            "per_chr_total": plan.per_chr_total,
            "per_chr_on": plan.per_chr_on,
            "per_chr_off": plan.per_chr_off,
        }
        with open(os.path.join(args.outdir, f"meta_{sample_name}.json"), "w") as mf:
            json.dump(meta, mf, indent=2)

    # Write counts table
    out_table = os.path.join(args.outdir, "table.tsv")
    write_counts_table(out_table, feat_names, args.num_samples, counts_table)

    console.rule("[bold green]Done")
    console.print(f"[bold]BAMs:[/bold] {args.outdir}/sample_{{1..{args.num_samples}}}.bam")
    console.print(f"[bold]Counts table:[/bold] {out_table}")

if __name__ == "__main__":
    main()
