---
sort: 4
---

# CovToTarget

This tool produce a coverage per feature report based on
a *target* (annotation file in BED or GFF3 format) and the output
of *[covtobed 1.0](https://github.com/telatin/covtobed)*.

:warning: Note that while [bamtocov](bamtocov.md) produces an
identical coverage output, it also includes the built-in
feature to restrict the analysis to the target.

## Help screen

```text
covToTarget $version

  Usage: covtotarget [options] <Target> [<covtobed-output>]

Arguments:                                                                                                                                                 

  <Target>           the BED (or GFF) file containing regions in which to count reads
  <covtobed-output>  covtobed output, or STDIN if not provided

Options:

  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]
  -l, --norm-len               Normalize by gene length
  -b, --bed-output             Output format is BED-like (default is feature_name [tab] counts)
  -h, --help                   Show help
```


## Usage

This program extends `covtobed`, so it can be used in a stream:

```bash
covtobed input/mini.bam | covtotarget input/mini.bed > output/counts.tsv
```

## Output format

A table of features and counts per feature:

```text
target1_8X      699
target2_0X      0
target3_1X      50
for_rev_10Xa    690
for_rev_10Xb    0
for_rev_10Xc    0
```