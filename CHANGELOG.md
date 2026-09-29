# Changelog

All notable changes to BamToCov are documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Releases before 2.10.0 are listed in [docs/_docs/history.md](docs/_docs/history.md).

## [Unreleased]

### Added

- `bamtocounts`: `--exclude-multimapping` drops reads with an `NH` tag > 1,
  matching the default multi-mapping filter of featureCounts (no effect if the
  aligner does not write `NH`).

### Fixed

- `bamcountrefs`: memory no longer grows with references × samples × metrics.
  Output tables are streamed row by row instead of being built in memory, and
  per-sample results are freed as soon as they are merged. On a metagenome
  co-assembly (3.59 M contigs × 21 BAMs, `--all-metrics`) the previous build
  was killed at over 60 GB; it now completes with a 30–35 GB peak, and on 4
  BAMs it went from 107 s / 17 GB to 64 s / 10 GB. Output is byte-identical.
- `bamcountrefs`: `-W/--workers` is now honoured. It was parsed but ignored,
  and one thread was started for every input BAM at once.
- `bamcountrefs`: the RPKM denominator is now the number of reads passing the
  filters (`-F`, `-Q`, `-P`). It was the BAM index "mapped" total, which
  includes supplementary, secondary, duplicate and QC-fail alignments, so RPKM
  was too low whenever filters removed reads (~1.3% on a minimap2 co-assembly
  BAM). Values now match coverm.

### Changed

- `bamcountrefs`: RPKM values change for any BAM where the filters remove
  reads (see above). TPM and all other metrics are unchanged.
- `bamcountrefs`: with many BAMs and a low `-W`, fewer files are processed in
  parallel than before, since the old build ignored `-W`. Raise `-W` to trade
  memory for parallelism.

## [2.10.0]

- `bamcountrefs` now supports multithreading
- Added `--length` output option to `bamcountrefs`
- `bamtocov` now handles Ctrl-C/Ctrl-D gracefully in streams
- Optimized `bamcountrefs` and refactored `bamtocounts`
- Added **DisCov** (merenlab) integration to `bamcountrefs`
- Consolidated Nim configuration (`nim.cfg`)
- Fixed paired-end fragment counting in `bamtocounts`
- Bug fixes since 2.9:
  - `bamtocounts`: fixed overlap predicate for exact matches and boundary cases
  - `bamtocounts`: fixed RPKM denominator to account for MAPQ/flag/paired-read filtering
  - `bamtocov`: fails gracefully when target contigs are missing from BAM/CRAM headers
  - `bamtocov`: preserves target-file order in reports
  - `bamtocov`: fixed extra empty column in stranded quantized BED output
  - `bamtocov`: fixed WIG span state leakage across contigs in stranded mode
  - `bamtocov`: replaced deprecated output destructor with explicit flush
  - Replaced `--op` assert with CLI validation
