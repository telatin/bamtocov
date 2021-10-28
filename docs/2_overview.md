---
sort: 2
permalink: /overview
---

# Overview

In bioinformatics, a common task is to align several (usually short) DNA sequences 
against a reference sequence (e. g. a complete genome of the organism).

![coverage in a genome browser](https://camo.githubusercontent.com/340658f1bbad9a553dd9b3e8d943f6696cb4a825840ee22b65a5fe754cda2afa/68747470733a2f2f6d69726f2e6d656469756d2e636f6d2f6d61782f323733342f312a68616f374a4d4c516c6f71626d68412d6535793141512e706e67)

<sub>The screenshot of  shows the alignment of several short DNA sequences (*reads*)
against a reference genome. The program graphically displays the *coverage track*</sub>

[IGV](https://software.broadinstitute.org/software/igv/) itself displays a coverage track,
but for computational reasons only if the zoom is below 50
[kbp](https://en.wikipedia.org/wiki/Base_pair).
On the other hand, it is possible to visualize a
[BED track](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
spanning the whole chromosome.

**bamtocov** reads an alignment file in BAM format and prints a BED track with the
nucleotide coverage and can be used to select the regions with a coverage of
interest
(*e. g.* uncovered regions with 0X coverage, low coverage regions <5X, etc.).

## Features

:white_check_mark: can read sorted BAM files even without the index

:white_check_mark: can read sorte BAM files from stream

:white_check_mark: can calculate the **physical coverage**

:white_check_mark: can calculate the **per-strand coverage**

## Further reading

For a more detailed introduction to coverage analysis see: [Coverage Analysis from the Command Line](https://medium.com/ngs-sh/coverage-analysis-from-the-command-line-542ef3545e2c).
