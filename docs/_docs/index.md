---
title: BamToCov 3.0
description: "A suite of utilities for BAM file coverage analysis"
tags: [bamtocov, bam, coverage, bioinformatics, sequence]
permalink: /
---
 
 
😺 **[Repository](https://github.com/telatin/bamtocov)** | 📦 **[Releases](https://github.com/telatin/bamtocov/releases)** | 📃 **[Paper](https://academic.oup.com/bioinformatics/article/38/9/2617/6535233)**

![bamtocov logo]({{ '/assets/img/bamtocov-banner.png' | relative_url }})

BamToCov is inspired by the [UNIX Phylosophy](https://en.wikipedia.org/wiki/Unix_philosophy) and the tools are designed for efficient computation
of a very specific task. Integration of multiple samples and specific tasks can be achieved [with scripts]({{ '/scripts/' | relative_url }})
and we provide a set to demonstrate the process.

**[bamtocov]({ /_docs/commands/bamtocov.md  | relative_url})** will produce a _coverage BED_ from a single BAM file, or a count matrix from a set of alignments and a target (in BED, GTF or GFF format).
Used without a target, it is a drop-in replacement for [covtobed](https://github.com/telatin/covtobed), but discarding invalid alignments by default.
When providing the target, it can produce coverage statistics for each region in the target, also with multiple BAM files.

**[bamtocounts]({ /_docs/commands/bamtocounts.md  | relative_url})** will count the number of reads covering each target region,
rather than the nucleotidic coverage

**[bamcountrefs]({ /_docs/commands/bamcountrefs.md  | relative_url})** is a shortcut to count the number of reads per chromosome,
with filters on the read flags, length and quality

**[covtotarget]({ /_docs/commands/covtotarget.md  | relative_url})** (_legacy_) is an utility to create a count table from the _output_ of
the original [covtobed](https://github.com/telatin/covtobed) program.
 
