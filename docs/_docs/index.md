---
title: BamToCov 3.0
description: "A suite of utilities for BAM file coverage analysis"
tags: [bamtocov, bam, coverage, bioinformatics, sequence]
permalink: /
---
 
 
😺 **[Repository](https://github.com/telatin/bamtocov)** | 📦 **[Releases](https://github.com/telatin/bamtocov/releases)** | 📃 **[Paper](https://academic.oup.com/bioinformatics/article/38/9/2617/6535233)**

![bamtocov logo]({% link assets/img/bamtocov-banner.png %})

BamToCov is inspired by the [UNIX Phylosophy](https://en.wikipedia.org/wiki/Unix_philosophy) and the tools are designed for efficient computation
of a very specific task. Integration of multiple samples and specific tasks can be achieved [with scripts](https://telatin.github.io/bamtocov/scripts/)
and we provide a set to demonstrate the process.

**[bamtocov]({% link _docs/commands/bamtocov.md %})** will produce a _coverage BED_ from a single BAM file, or a count matrix from a set of alignments and a target (in BED, GTF or GFF format).
Used without a target, it is a drop-in replacement for [covtobed](https://github.com/telatin/covtobed), but discarding invalid alignments by default.
When providing the target, it can produce coverage statistics for each region in the target, also with multiple BAM files.

**[bamtocounts]({% link _docs/commands/bamtocounts.md %})** will count the number of reads covering each target region,
rather than the nucleotidic coverage

**[bamcountrefs]({% link _docs/commands/bamcountrefs.md %})** is a shortcut to count the number of reads per chromosome,
with filters on the read flags, length and quality

**[covtotarget]({% link _docs/commands/covtotarget.md %})** (_legacy_) is an utility to create a count table from the _output_ of
the original [covtobed](https://github.com/telatin/covtobed) program.
 