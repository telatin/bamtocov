---
sort: 3
permalink: /usage
---

# Usage

## Tools

**[bamtocov]({{ site.baseurl }}/tools/bamtocov.html)** will produce a _coverage BED_ from a single BAM file, or a count matrix from a set of alignments and a target (in BED, GTF or GFF format).

Used without a target, it is a drop-in replacement for [covtobed](https://github.com/telatin/covtobed), but discarding invalid alignments by default.
When providing the target, it can produce coverage statistics for each region in the target, also with multiple BAM files.

**[bamtocounts]({{ site.baseurl }}/tools/bamtocounts.html)** will count the number of reads covering each target region, rather than the nucleotidic coverage

**[bamcountrefs]({{ site.baseurl }}/tools/bamcountrefs.html)** is a shortcut to count the number of reads per chromosome, with filters on the read flags, length and quality

**[covtotarget]({{ site.baseurl }}/tools/covtotarget.html)** is an utility to create a count table from the _output_ of the original [covtobed](https://github.com/telatin/covtobed) program.

 
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
