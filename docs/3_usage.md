---
sort: 3
permalink: /usage
---

# Usage

* **bamtocov** will produce a _coverage BED_ from a single BAM file, or a count matrix from a set of alignments and a target (in BED, GTF or GFF format)
* **bamtocounts** will count the number of reads covering each target region, rather than the nucleotidic coverage
* **bamcountrefs** is a shortcut to count the number of reads per chromosome, with filters on the read flags, length and quality
* **covtotarget** is an utility to create a count table from the _output_ of the original [covtobed](https://github.com/telatin/covtobed) program.

Each tool has a [dedicated page](tools).
