# bamtocounts: Multi-BAM parallel processing with TPM support and improved paired-end handling

## Summary

Major refactor of `bamtocounts` to support simultaneous processing of multiple BAM/CRAM files
using parallel worker threads, add TPM normalization, improve paired-end semantics, and
significantly optimize the read counting algorithm.

---

## New Features

### Multi-BAM parallel processing

- Accepts multiple BAM/CRAM files as positional arguments (previously only one was supported).
- Introduces a thread-per-file model via `WorkerTask` / `processOneFile` using Nim's `Thread[T]`.
- New `--workers` / `-W` option controls the number of parallel file processors; defaults to
  `min(numBamFiles, countProcessors())` for automatic tuning.
- Each worker independently opens the BAM file, parses the target file, cooks the target table,
  and streams reads — results are written into a `SampleResult` struct and aggregated after all
  threads finish.
- Output columns are ordered one sample per column group (counts, optional RPKM/TPM/norm-len).

### TPM normalization (`--tpm` / `-m`)

- Adds `rpk()` and `tpm()` helper procs implementing the standard two-pass TPM algorithm:
  1. Compute RPK (Reads Per Kilobase) for every feature.
  2. Sum all RPK values to get a scaling factor.
  3. Divide each RPK by the sum and multiply by 1,000,000.
- This guarantees that TPM values across all features in a sample sum to exactly 1,000,000.
- Computed after all threads join, once per sample, using the aggregated result tables.

### Improved paired-end handling

- `--paired` semantics clarified: counts fragments by processing only R1 reads, skipping R2 to
  avoid double-counting. Previously the flag also required `proper_pair`, conflating two distinct
  filters.
- New `--proper-pairs` flag independently requires that reads from paired experiments be properly
  paired (SAM flag 0x2 set when 0x1 is set). Can be used with or without `--paired`.

---

## Algorithmic Improvements

### Optimized counting in `processOneFile`

- **Index-based chromosome skipping**: when a BAM index is available (`bam.idx != nil`), only
  chromosomes that appear in the cooked target are queried. Chromosomes with zero mapped reads
  are additionally skipped using `stats(bam.idx, index).mapped`.
- **Per-chromosome querying**: uses `bam.query(chrName)` instead of iterating the entire file,
  which avoids decompressing irrelevant data for targeted assays.
- **Simplified overlap check**: replaced the complex three-condition overlap expression with the
  canonical half-open interval test `readStop > region.start and readStart < region.stop`.
- **Pre-initialised feature table**: all target features are seeded with zero counts before
  read iteration, eliminating `hasKey` checks in the inner loop and ensuring zero-coverage
  features always appear in the output.
- **Cached read properties**: `readStart`, `readStop`, and `isReverse` are extracted once per
  read, outside the region inner loop.
- A fallback full-BAM iteration path is retained for stdin, `/dev/stdin`, and `/dev/fd/*` inputs
  where index access is unavailable.

---

## Refactored / Removed Code

- `makeCountsTable` and `alignments_count` procs removed; their logic is consolidated into the
  new `processOneFile` thread proc.
- `FeatureMetrics` object removed (was unused in the previous version).
- `counts_t` type alias removed (also unused).
- The single-BAM main loop is replaced by a thread-dispatch + join + aggregation pipeline.
- Header generation updated to emit one column group per sample, supporting stranded output and
  all normalization columns simultaneously.
- Target-error diagnostic that printed the first non-comment line of the target file on failure
  has been simplified; callers should use `--debug` for verbose diagnostics.

---

## New Dependencies

- `std/sequtils` — `toSeq` for converting table keys to sequences.
- `std/cpuinfo` — `countProcessors()` for automatic worker count.
- `std/algorithm` — `sort` for deterministic chromosome ordering in output.

---

## Compatibility Notes

- Output column order changes when multiple BAM files are supplied: columns are now grouped by
  sample rather than by metric. Single-BAM invocations produce identical output to before
  (except for the new optional `--tpm` column).
- `--paired` no longer implies `--proper-pairs`; scripts that relied on the old implicit
  proper-pair filter should add `--proper-pairs` explicitly.
