---
title: First steps
description: "Some bamtocov examples to get you started"
date created: "2025-12-29"
last updated: "2025-06-22"
---

## Produce the coverage BED from a BAM file

```bash
bamtocov alignments.bam > coverage.bed
```

## Use a target file to count the coverage in each region

```bash
bamtocov --regions target.bed --report target-stats.tsv alignments.bam > coverage.bed 
```

## Generate a matrix from multiple outputs

The output of multiple files in BED format makes no sense, so at this moment we require to
explicitly specifying a `--skip-output` flag. The report matrix is saved with `--report` (or `-o`).

```bash
bamtocov --regions target.bed --report target-stats.tsv --skip-output *.bam 
```

## Count the coverage per target using legacy "covtobed" output

```bash
covtobed file.bam | covtotarget ...
```
