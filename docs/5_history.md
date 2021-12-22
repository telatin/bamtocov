---
sort: 5
permalink: /history
---
# History

This project extends [covtobed](https://github.com/telatin/covtobed),
reimplementing the core algorithm in Nim.

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
