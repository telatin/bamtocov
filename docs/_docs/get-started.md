---
title: Get started
description: "Install Bamtocov and learn the basics"
tags: [seqfu, installation, usage, fasta, fastq]
---

This page covers installation and basic usage of Bamtocov.

## Installation

You can install _bamtocov_ from BioConda, if you have 
[_conda_](https://docs.conda.io/en/latest/miniconda.html) installed:

```bash
conda install -c conda-forge -c bioconda bamtocov
```

### Compiling from source

1. Install [nim](https://nim-lang.org/) and nimble.
1. Install [hts-lib](https://www.htslib.org)
1. Compile with `nimble build`


## Quick start

```bash
bamtocov alignment.bam > coverage.bed
```

will produce a coverage BED file from the alignment file.

## File formats

### BED files

A BED file (.bed) is a tab-delimited text file that defines a feature track. In this context the magnitude
refers to the _nucleotide coverage_ of the interval.

The columns are _chromosome name_, _start position_ (inclusive, zero-based), _end position_ 
(non-inclusive, zero-based) and _coverage_.
An example is:

```text
seq1    0       9       0
seq1    9       109     5
seq1    109     189     0
seq1    189     200     2
```

### Target statistics

:warning: this format is not final.

For each sample, 5 columns are printed:

* `bam_bases`
* `bam_mean`
* `bam_min`
* `bam_max`
* `bam_length`

| interval     | bam_bases | bam_mean | bam_min | bam_max | bam_length |
| ------------ | --------: | -------: | ------: | ------: | ---------: |
| target1_8X   |       699 |    3.495 |       1 |       6 |        200 |
| target2_0X   |         0 |      0.0 |       0 |       0 |         50 |
| target3_1X   |         . |        . |       . |       . |          . |
| for_rev_10Xa |       100 |     10.0 |      10 |      10 |         10 |
| for_rev_10Xb |       100 |     10.0 |      10 |      10 |         10 |
| for_rev_10Xc |         . |        . |       . |       . |          . |
