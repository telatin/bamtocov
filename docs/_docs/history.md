---
title: History
description: "The history of BamToCov development"
tags: [bamtocov, history, development]
---

This project extends [covtobed](https://github.com/telatin/covtobed),
reimplementing the core algorithm in Nim.

* 2.11.0
  - `bamtocounts`: new options to match featureCounts behaviour
    - `--exclude-multimapping` drops reads with an `NH` tag > 1
    - `-B/--both-mapped` (with `--paired`) requires the mate to be mapped
    - `--fraction` splits an ambiguous fragment's credit 1/n across the n regions it overlaps (like `featureCounts -O --fraction`)
    - `--exclude-ambiguous` drops fragments overlapping more than one region (featureCounts' default); mutually exclusive with `--fraction`
  - `bamcountrefs`: `--debug` prints the version and build features (e.g. libdeflate) of the htslib loaded at run time
  - `bamcountrefs`: much lower memory use with identical output: tables are streamed row by row, per-sample results are freed after merging, and only the metrics actually written are computed (e.g. 1 M contigs × 12 BAMs: 5.4 GB → 2.6 GB peak)
  - `bamcountrefs`: `-W/--workers` is now honoured (it was parsed but ignored, starting one thread per BAM)
  - Bug fixes
    - `bamcountrefs`: **RPKM** denominator is now the number of reads passing the filters (`-F`, `-Q`, `-P`), not the index "mapped" total; values now match coverm. RPKM changes wherever filters remove reads; other metrics are unchanged
    - `bamcountrefs`: **CRAM** counts were all zero (missing CRAM index statistics were read as "no reads")
    - `bamcountrefs`: `-P/--proper-pairs` and `--mean` were ignored by the index-count shortcut used with `-F 4 -Q 0`
    - `bamcountrefs`: references are queried by id, so a reference named like a region (e.g. `seq1:1-10`) is no longer reported as 0
    - `bamcountrefs`: unreadable, unindexed or corrupted input files stop the run with an error naming the file (instead of a silent column of zeros)
    - `bamcountrefs`: warns when samples have different references; duplicate sample names are made unique (`name_2`, ...)
    - `bamcountrefs`: invalid numeric options are rejected (e.g. `-Q 300`, `--trim-min` >= `--trim-max`)

* 2.10.0
  - `bamcountrefs` now supports multithreading
  - Added `--length` output option to `bamcountrefs`
  - `bamtocov` now handles Ctrl-C/Ctrl-D gracefully in streams
  - Optimized `bamcountrefs` and refactored `bamtocounts`
  - Added **DisCov** (merenlab) integration to `bamcountrefs`
  - Consolidated Nim configuration (`nim.cfg`)
  - Fixed paired-end fragment counting in `bamtocounts`
  - Bug fixes since 2.9 (grouped under 2.10.0)
    - `bamtocounts`: fixed overlap predicate for exact matches and boundary cases
    - `bamtocounts`: fixed RPKM denominator to account for MAPQ/flag/paired-read filtering
    - `bamcountrefs`: fixed RPKM denominator to count only reads passing the filters (it used the index mapped total, including supplementary alignments)
    - `bamtocov`: fails gracefully when target contigs are missing from BAM/CRAM headers
    - `bamtocov`: preserves target-file order in reports
    - `bamtocov`: fixed extra empty column in stranded quantized BED output
    - `bamtocov`: fixed WIG span state leakage across contigs in stranded mode
    - `bamtocov`: replaced deprecated output destructor with explicit flush
    - Replaced `--op` assert with CLI validation

* 2.9.0
  * **bamcountrefs** is now multithreading
  * bug fixes
* 2.8.0
  * **bamcountrefs** bugfix for multiple files; optimised
  * refactoring: unified imports, added tests, fixed Make test to use the environment's shell
* 2.7.0
  * **bamtocounts** now supports `--paired` to count fragments
  * **bamtocounts** fixed a bug that ignored the user flag to exclude
  * [**EXPERIMENTAL**] testing `--extendRead INT` in bamtocov #6
* 2.6.0
  * **bamtocounts** algorithm rewritten to accept unsorted and non-indexed files. the requirement for indexed files will be permanently dropped, for consistency across the suite, while we recommend using sorted files, and, in the future, this requirement might be added again to allow for performance gains.
* 2.5.0
  * **BamToCounts** rewritten with the target engine of BamToCov
* 2.4.0
  * Standed analysis is now supported with Wig-like output 
* 2.3.0:
  * Improved Wig support
  * Improved GTF format detection
  * Expanded test suite
* 2.2.0:
  * Support for wiggle output (`--wig STEP` and `--op FUNCT`)
* 2.1.1:
  * Help screen for all the tools updated to reflect the new tool names
  * Added `--tag STRING` in `bamcountrefs` to customise the reference column name
  * Initial support for **Wig** files
  * Expanded test suite and memory/speed benchmarks
* 2.1.0:
  * Bug Fix: an assertion caused the program to fail if an alignment ended at the end of a chromosome.
  * Updated compilation instructions for Nimble and Bioconda builds
  * Speed improvements
* 2.0.4:
  * Added new internal classes, like _output\_t_ and _coverage\_t_
* 2.0.2:
  * Added `covtocounts`
  * Added target support
  * Performance improvement
* 2.0.0:
  * Initial release: Nim porting of the original C++ code
