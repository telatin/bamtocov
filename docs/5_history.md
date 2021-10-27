---
sort: 5
permalink: /history
---
# History

This project extends [covtobed](https://github.com/telatin/covtobed),
reimplementing the core algorithm in Nim.

* 2.2.0:
  * Help screen for all the tools updated to reflect the new tool names
  * Added `--tag STRING` in `bamcountrefs` to customise the reference column name
* 2.1.0:
  * Fixed a bug: an assertion was too stringed and caused the program to fail if an alignment ended at the end of a chromosome.
  * Updated compilation instructions
* 2.0.4:
  * Added new internal classes, like _output\_t_ and _coverage\_t_
  * This release is not yet optimized for speed.
* 2.0.2:
  * Added target support. Performance improvement.
* 2.0.0:
  * Initial release

