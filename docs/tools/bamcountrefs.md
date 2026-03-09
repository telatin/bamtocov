---
title: BamCountRefs
parent: Tools
nav_order: 3
---

# BamCountsRefs

A program to build a count table from multiple BAM
files (having the same reference sequence).

```text
BamCountRefs 2.9.0

  Usage: bamcountrefs [options]  <BAM-or-CRAM>...

Arguments:

  <BAM-or-CRAM>  the alignment file for which to calculate depth

BAM/CRAM processing options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -W, --workers <workers>      Number of parallel file processors [default: auto]
  -r, --fasta <fasta>          FASTA file for use with CRAM files [default: ].
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]

Output options:
  -o, --output <BASENAME>      Output file basename (generates multiple files:  <BASENAME>_counts.tsv, <BASENAME>_rpkm.tsv, etc. depending on the metrics selected)

                               If not specified, outputs counts to stdout in TSV format
  -n                           [DEPRECATED: use --rpkm] Output RPKM values
  --rpkm                       Calculate RPKM (reads per kilobase per million mapped reads)
  --tpm                        Calculate TPM (transcripts per million)
  --mean                       Calculate mean coverage depth (approximate method)
  --covered-bases              Calculate number of bases with coverage > 0 *
  --covered-ratio              Calculate coverage breadth (fraction of reference covered) *
  --length                     Output reference sequence lengths ()
  --all-metrics                Enable all available metrics

Other options:
  --tag STR                    First column name [default: Contig]
  --multiqc                    Print output as MultiQC table (stdout only)
  --debug                      Enable diagnostics
  -h, --help                   Show help
  
[*] requires extra memory compared to base metrics
```

## Examples

### Basic Usage (stdout)

Output counts to stdout:

```bash
bin/bamcountrefs --tag "Chrom" input/mini.bam input/mini2.bam
```

Output:

```text
Chrom   mini    mini2
seq0    0       1
seq1    15      15
seq2    10      10
```

### Multi-file Output

Generate separate files for different metrics:

```bash
bin/bamcountrefs --output results/sample --rpkm --tpm --mean input/mini.bam input/mini2.bam
```

This creates:
- `results/sample_counts.tsv` - Raw read counts
- `results/sample_rpkm.tsv` - RPKM normalized values
- `results/sample_tpm.tsv` - TPM normalized values
- `results/sample_mean.tsv` - Mean coverage depth (approximate)

### All Metrics at Once

Generate all available metrics with a single command:

```bash
bin/bamcountrefs --output results/sample --all-metrics input/*.bam
```

This creates all output files:
- `results/sample_counts.tsv` - Raw read counts
- `results/sample_rpkm.tsv` - RPKM normalized values
- `results/sample_tpm.tsv` - TPM normalized values
- `results/sample_mean.tsv` - Mean coverage depth (approximate)
- `results/sample_covered_bases.tsv` - Number of bases with coverage > 0
- `results/sample_covered_fraction.tsv` - Fraction of reference covered (breadth)

### Coverage Breadth Metrics

Calculate coverage breadth (what fraction of each reference is covered):

```bash
bin/bamcountrefs --output results/sample --covered-bases --covered-ratio input/*.bam
```

Note: Breadth metrics require tracking per-base coverage, which uses additional memory proportional to reference length.