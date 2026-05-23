---
title: BamCountRefs
parent: Tools
nav_order: 3
---


A program to build a count table from multiple BAM
files (having the same reference sequence).

```text
BamCountRefs 

  Usage: bamcountrefs [options]  <BAM-or-CRAM>...

Arguments:

  <BAM-or-CRAM>  the alignment file for which to calculate depth

BAM/CRAM processing options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -W, --workers <workers>      Number of parallel file processors [default: auto]
  -r, --fasta <fasta>          FASTA file for use with CRAM files [default: ].
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 3844]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]
  -P, --proper-pairs           If paired flag is set, then also proper pair must be set

Output options:
  -o, --output <BASENAME>      Output file basename (generates multiple files: <BASENAME>_counts.tsv, etc.)
                               If not specified, outputs counts to stdout in TSV format
  -n                           [DEPRECATED: use --rpkm] Output RPKM values
  --rpkm                       Calculate RPKM (reads per kilobase per million mapped reads)
  --tpm                        Calculate TPM (transcripts per million)
  --mean                       Calculate mean coverage depth (approximate method, no extra memory)
  --trimmed-mean               Calculate trimmed mean coverage (robust against outliers) [requires extra memory]
  --trim-min <FRACTION>        Remove this smallest fraction of positions when calculating trimmed_mean [default: 5]
  --trim-max <FRACTION>        Maximum fraction for trimmed_mean calculations [default: 95]
  --covered-bases              Calculate number of bases with coverage > 0 [requires extra memory]
  --covered-ratio              Calculate coverage breadth (fraction of reference covered) [requires extra memory]
  --variance                   Calculate variance of coverage depth [requires extra memory]
  --reads-per-base             Calculate reads per base (count / length, normalized read density)
  --discov                     Calculate Distribution of Coverage score [requires extra memory]
  --discov-window <INT>        Window length for DisCov spread component [default: 1000]
  --discov-fold-lower <FLOAT>  Lower fold range around median nonzero coverage [default: 0.5]
  --discov-fold-upper <FLOAT>  Upper fold range around median nonzero coverage [default: 2.0]
  --discov-alpha <FLOAT>       Weight on DisCov spread component [default: 0.5]
  --discov-formula <FORMULA>   DisCov formula: linear or geometric [default: linear]
  --length                     Output reference sequence lengths
  -a, --all-metrics            Enable all available metrics

Other options:
  --tag STR                    First column name [default: Contig]
  --multiqc                    Print output as MultiQC table (stdout, or <BASENAME>.mqc.txt with -o)
  --debug                      Enable diagnostics
  -h, --help                   Show help
```

## Memory Requirements

Different metrics have different memory requirements:

**Low memory** (no extra memory per reference):

- counts, rpkm, tpm, mean, reads-per-base, length

**High memory** (requires per-base tracking):

- covered-bases, covered-ratio, variance, trimmed-mean, discov

For large reference sequences or many samples, high-memory metrics will require RAM proportional to reference length. The algorithm implements several optimizations:

- Zero-coverage references are detected early and skip expensive computations
- Depth arrays are shared between variance and trimmed-mean calculations
- Processing is parallelized across multiple BAM files

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
bin/bamcountrefs --output results/sample --rpkm --tpm --mean --variance input/mini.bam input/mini2.bam
```

This creates:
- `results/sample_counts.tsv` - Raw read counts
- `results/sample_rpkm.tsv` - RPKM normalized values
- `results/sample_tpm.tsv` - TPM normalized values
- `results/sample_mean.tsv` - Mean coverage depth (approximate)
- `results/sample_variance.tsv` - Variance of coverage depth

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
- `results/sample_variance.tsv` - Variance of coverage depth
- `results/sample_trimmed_mean.tsv` - Trimmed mean coverage (robust statistic)
- `results/sample_reads_per_base.tsv` - Reads per base (normalized read density)
- `results/sample_discov.tsv` - Distribution of Coverage score
- `results/sample_covered_bases.tsv` - Number of bases with coverage > 0
- `results/sample_covered_fraction.tsv` - Fraction of reference covered (breadth)
- `results/sample_length.tsv` - Reference sequence lengths

### Coverage Breadth Metrics

Calculate coverage breadth (what fraction of each reference is covered):

```bash
bin/bamcountrefs --output results/sample --covered-bases --covered-ratio input/*.bam
```

Note: Breadth metrics require tracking per-base coverage, which uses additional memory proportional to reference length.

### Trimmed Mean Coverage

Calculate trimmed mean coverage for robust statistics that are less sensitive to outliers:

```bash
bin/bamcountrefs --output results/sample --trimmed-mean --trim-min 10 --trim-max 90 input/*.bam
```

The trimmed mean removes extreme values before calculating the mean:
- `--trim-min 10` removes the bottom 10% of coverage positions
- `--trim-max 90` removes the top 10% of coverage positions
- Default values are 5 and 95 (removing 5% from each tail)

This is particularly useful for:
- Metagenomics data with variable coverage
- Detecting regions with consistently high/low coverage
- Robust coverage estimation in the presence of PCR duplicates or mapping artifacts

### Variance and Statistical Metrics

Calculate variance to understand coverage uniformity:

```bash
bin/bamcountrefs --output results/sample --variance --mean input/*.bam
```

The variance metric indicates how evenly reads are distributed:
- Low variance: uniform coverage across the reference
- High variance: uneven coverage with peaks and valleys

The reads-per-base metric provides a length-normalized read density:

```bash
bin/bamcountrefs --output results/sample --reads-per-base input/*.bam
```

This is equivalent to `count / length` and useful for comparing references of different lengths.

### DisCov

Calculate Distribution of Coverage scores to summarize both spatial coverage spread and evenness of nonzero depth:

```bash
bin/bamcountrefs --output results/sample --discov input/*.bam
```

DisCov uses non-overlapping windows for the spread component and the fraction of nonzero-covered bases within a fold range around the median nonzero coverage for the evenness component. The default score is the linear formula with `--discov-alpha 0.5`, `--discov-window 1000`, `--discov-fold-lower 0.5`, and `--discov-fold-upper 2.0`. See [DisCov]({% link _docs/discov.md %}) for the metric details.
