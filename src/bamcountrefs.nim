# Standard library
import std/[os, strutils, tables, cpuinfo, algorithm, dynlib]

# External dependencies
import docopt, hts
# hts_idx_get_stat is not re-exported by hts; the wrapper's stats() discards
# its return status, which is needed to tell "no reads" from "no statistics"
from hts/private/hts_concat import hts_idx_t, hts_idx_get_stat, hts_version,
    libname

import ./discov

const NimblePkgVersion {.strdefine.} = "prerelease"

# Named constants for clarity
const
  DEFAULT_EXCLUDE_FLAGS = 3844 # Unmapped(4) + Secondary(256) + QC-fail(512) + Duplicate(1024) + supplementary(2048) = 3844
  READS_PER_MILLION = 1_000_000
  BASES_PER_KILOBASE = 1000

let
  version = NimblePkgVersion

type
  EKeyboardInterrupt = object of CatchableError

  ReferenceMetrics = object
    refName: string
    order: int # Preserve BAM header order
    length: int
    # Per-sample raw counts (one per BAM file)
    sampleCounts: seq[int]
    # Per-sample depth tracking for mean coverage (approximate method)
    sampleTotalDepth: seq[int64]
    # Per-sample breadth tracking (bases with coverage > 0)
    sampleCoveredBases: seq[int]
    # Per-sample sum of squared depths for variance calculation
    sampleSumSquaredDepth: seq[int64]
    # Per-sample values computed in the workers (derived metrics such as
    # RPKM or TPM are computed row by row while writing the output)
    sampleTrimmedMean: seq[float]
    sampleDiscov: seq[float]
    # Only the sequences listed in the run's set[StoredValue] are filled;
    # the others stay empty

  StoredValue = enum
    ## Per-sample values kept for the output tables (see storedValues)
    svCounts, svTotalDepth, svCoveredBases, svSumSquaredDepth, svTrimmedMean,
    svDiscov

  Metric = enum
    ## One output table; the order is the order of the output files
    mCounts, mRPKM, mTPM, mMean, mTrimmedMean, mCoveredBases, mCoveredRatio,
    mVariance, mReadsPerBase, mDiscov, mLength

  # Concurrency-safe types for parallel processing
  RefAggWithName = object
    ## Per-reference aggregated metrics with name (thread-local)
    name: string
    length: int
    order: int
    count: int
    totalDepth: int64
    coveredBases: int
    sumSquaredDepth: int64
    trimmedMean: float
    discov: float

  SampleResult = object
    ## Worker return value containing all metrics for one sample
    ## Uses sequences instead of Tables for thread-safety
    perRef: seq[RefAggWithName] # All reference metrics
    countedReads: float # Reads passing filters, for RPKM
    error: string # Non-empty if the worker failed; perRef is then unusable

  WorkerOpts = object
    ## Immutable options passed to each worker thread
    mapq: uint8
    eflag: uint16
    properPairs: bool # If true, require proper pair flag for paired reads
    sumSpans: bool # Sum alignment spans for the mean (when depths are not tracked)
    trackBreadth: bool
    trackDepths: bool # Track per-base depths for variance and/or trimmed mean
    doVariance: bool
    doTrimmedMean: bool
    doDiscov: bool
    trimMin: float    # Minimum percentile for trimmed mean (0-100)
    trimMax: float    # Maximum percentile for trimmed mean (0-100)
    discovParams: DiscovParams
    threads: int
    fasta: string

  WorkerTask = object
    ## Task description for a worker thread
    bamPath: string
    sampleIndex: int
    opts: WorkerOpts
    result: ptr SampleResult # Pointer to pre-allocated result slot

var
  # Main data structure: ordered table to preserve BAM order
  metricsTable = initOrderedTable[string, ReferenceMetrics]()

proc handler() {.noconv.} =
  raise newException(EKeyboardInterrupt, "Keyboard Interrupt")

setControlCHook(handler)


var
  debug = false

type
  OutputContext = object
    ## Per-sample totals needed by the normalized metrics
    totalCountedReads: seq[float] # RPKM denominator (reads passing filters)
    totalRPK: seq[float]         # TPM denominator (sum of reads per kilobase)

const
  metricFileSuffix: array[Metric, string] = ["_counts.tsv", "_rpkm.tsv",
      "_tpm.tsv", "_mean.tsv", "_trimmed_mean.tsv", "_covered_bases.tsv",
      "_covered_fraction.tsv", "_variance.tsv", "_reads_per_base.tsv",
      "_discov.tsv", "_length.tsv"]

proc newOutputContext(totalCountedReads: seq[float], doTPM: bool): OutputContext =
  result.totalCountedReads = totalCountedReads
  result.totalRPK = newSeq[float](totalCountedReads.len)
  if doTPM:
    # TPM = RPK / sum(RPK) * 1e6, summed in reference order for each sample
    for metrics in metricsTable.mvalues:
      let refLengthKb = metrics.length.float / BASES_PER_KILOBASE.float
      for sampleIdx in 0 ..< totalCountedReads.len:
        result.totalRPK[sampleIdx] += metrics.sampleCounts[sampleIdx].float / refLengthKb

proc formatValue(metrics: ReferenceMetrics, metric: Metric, sampleIdx: int,
    ctx: OutputContext): string =
  ## Format one cell of an output table
  case metric
  of mCounts:
    $metrics.sampleCounts[sampleIdx]
  of mCoveredBases:
    $metrics.sampleCoveredBases[sampleIdx]
  of mLength:
    # Length is a property of the reference, so it is repeated for each sample
    $metrics.length
  else:
    let length = metrics.length.float
    let value = case metric
      of mRPKM:
        # RPKM = (reads × 1,000,000 × 1,000) / (total_counted_reads × reference_length)
        let refLengthKb = length / BASES_PER_KILOBASE.float
        let countedReadsMillions = ctx.totalCountedReads[sampleIdx] /
            READS_PER_MILLION.float
        metrics.sampleCounts[sampleIdx].float / (refLengthKb * countedReadsMillions)
      of mTPM:
        let rpk = metrics.sampleCounts[sampleIdx].float /
            (length / BASES_PER_KILOBASE.float)
        if ctx.totalRPK[sampleIdx] > 0: (rpk / ctx.totalRPK[sampleIdx]) *
            READS_PER_MILLION.float else: 0.0
      of mMean:
        # Approximate method: total_depth = sum of alignment lengths
        metrics.sampleTotalDepth[sampleIdx].float / length
      of mTrimmedMean:
        metrics.sampleTrimmedMean[sampleIdx]
      of mCoveredRatio:
        metrics.sampleCoveredBases[sampleIdx].float / length
      of mVariance:
        # E[X^2] - (E[X])^2; clamp floating point errors below zero
        let meanDepth = metrics.sampleTotalDepth[sampleIdx].float / length
        let variance = metrics.sampleSumSquaredDepth[sampleIdx].float / length -
            meanDepth * meanDepth
        if variance >= 0: variance else: 0.0
      of mReadsPerBase:
        metrics.sampleCounts[sampleIdx].float / length
      of mDiscov:
        metrics.sampleDiscov[sampleIdx]
      of mCounts, mCoveredBases, mLength:
        0.0 # handled above
    formatFloat(value, ffDecimal, 6)

proc writeTable(file: File, samples: seq[string], metric: Metric,
    ctx: OutputContext) =
  ## Write header and one row per reference, formatting each row as it is
  ## written so memory does not grow with references × samples
  file.writeLine(samples.join("\t"))
  var line = newStringOfCap(1024)
  # mvalues borrows each row; pairs would copy it, with all its sequences
  for metrics in metricsTable.mvalues:
    line.setLen(0)
    line.add(metrics.refName)
    for sampleIdx in 0 ..< ctx.totalCountedReads.len:
      line.add('\t')
      line.add(formatValue(metrics, metric, sampleIdx, ctx))
    file.writeLine(line)

proc openOutput(filename: string): File =
  if not open(result, filename, fmWrite):
    stderr.writeLine("ERROR: Unable to write to file: ", filename)
    quit(1)

proc multiqcMetric(doRPKM, doTPM, doMean, doTrimmedMean, doCoveredBases,
    doCoveredRatio, doVariance, doReadsPerBase, doDiscov: bool): Metric =
  ## MultiQC output holds a single table: the first requested metric
  if doRPKM: mRPKM
  elif doTPM: mTPM
  elif doMean: mMean
  elif doTrimmedMean: mTrimmedMean
  elif doCoveredRatio: mCoveredRatio
  elif doCoveredBases: mCoveredBases
  elif doVariance: mVariance
  elif doReadsPerBase: mReadsPerBase
  elif doDiscov: mDiscov
  else: mCounts

proc outputToStdout(samples: seq[string], multiqc: bool, metric: Metric,
    ctx: OutputContext) =
  # Legacy stdout output for backward compatibility
  if multiqc:
    stdout.writeLine("# plot_type: 'table'")
    stdout.writeLine("# section_name: 'BamToCov count'")
    stdout.writeLine("# description: 'Feature table: counts of mapped reads against predicted viral sequences'")
  writeTable(stdout, samples, metric, ctx)

proc outputToFiles(basename: string, samples: seq[string], metrics: set[Metric],
    multiqc: bool, multiqcMetric: Metric, ctx: OutputContext) =
  # Multi-file output: one file per requested metric (counts are always written)
  for metric in metrics:
    let filename = basename & metricFileSuffix[metric]
    var file = openOutput(filename)
    writeTable(file, samples, metric, ctx)
    file.close()
    if debug:
      stderr.writeLine("[debug] Wrote output to: ", filename)

  if multiqc:
    let filename = basename & ".mqc.txt"
    var file = openOutput(filename)
    file.writeLine("# plot_type: 'table'")
    file.writeLine("# section_name: 'BamToCov counts'")
    file.writeLine("# description: 'Feature table: counts of mapped reads against contigs'")
    writeTable(file, samples, multiqcMetric, ctx)
    file.close()
    if debug:
      stderr.writeLine("[debug] Wrote output to: ", filename)

proc indexStats(idx: ptr hts_idx_t, tid: int): tuple[ok: bool, mapped,
    unmapped: uint64] =
  ## Read counts for a reference from the index. ok is false when the index
  ## has no statistics for it: always for CRAM indices, and for BAM/CSI
  ## references with no indexed reads or an index without the pseudo-bin.
  ## Callers must then scan the alignments rather than assume zero reads.
  result.ok = hts_idx_get_stat(idx, cint(tid), addr result.mapped,
      addr result.unmapped) >= 0

proc passesFilters(aln: Record, opts: WorkerOpts): bool {.inline.} =
  if aln.mapping_quality < opts.mapq: return false
  if (aln.flag and opts.eflag) != 0: return false
  # If --proper-pairs is set, skip paired reads that aren't properly paired
  # Flag 0x1 (1) = read is paired
  # Flag 0x2 (2) = read is properly paired
  if opts.properPairs and (aln.flag and 1) != 0 and (aln.flag and 2) == 0:
    return false
  true

proc addAlignment(agg: var RefAggWithName, aln: Record, opts: WorkerOpts,
    covered: var seq[bool], depths: var seq[int]) {.inline.} =
  inc(agg.count)

  # Approximate mean coverage: sum alignment lengths. aln.stop walks the
  # CIGAR, so skip it unless the mean needs it
  if opts.sumSpans:
    agg.totalDepth += int64(aln.stop - aln.start)

  # Optimization: Combined tracking of breadth and depths in single pass
  # When both metrics are requested, this avoids iterating twice through all positions
  if opts.trackBreadth or opts.trackDepths:
    # Allocated on the first passing read: references without index
    # statistics are scanned even when empty, and most contigs of a
    # fragmented assembly have no reads
    if opts.trackBreadth and covered.len == 0:
      covered = newSeq[bool](agg.length)
    if opts.trackDepths and depths.len == 0:
      depths = newSeq[int](agg.length)
    for pos in aln.start ..< aln.stop:
      if pos < agg.length:
        if opts.trackBreadth and not covered[pos]:
          covered[pos] = true
          inc(agg.coveredBases)
        if opts.trackDepths:
          inc(depths[pos])

proc processReferences(bam: Bam, opts: WorkerOpts, perRef: var seq[RefAggWithName]) =
  # hdr.targets builds a new seq (with name strings) on every call
  let targets = bam.hdr.targets
  perRef = newSeq[RefAggWithName](targets.len)

  # Index counts can replace the alignment scan only when they describe
  # exactly the reads that would pass the filters (-F 4 -Q 0, no -P), and no
  # metric needs the alignments themselves
  let indexCountsSuffice = (opts.mapq == 0'u8) and (opts.eflag == 4'u16) and
      (not opts.properPairs) and (not opts.sumSpans) and
      (not opts.trackBreadth) and (not opts.trackDepths)
  # Unmapped reads placed on a reference pass a filter without bit 4
  let countUnmapped = (opts.eflag and 4'u16) == 0

  # Process each reference sequence
  for t in targets:
    var agg: RefAggWithName
    agg.name = t.name
    agg.length = int(t.length)
    agg.order = t.tid

    let st = indexStats(bam.idx, t.tid)
    let indexedReads = st.mapped + (if countUnmapped: st.unmapped else: 0'u64)

    if st.ok and indexedReads == 0:
      # Index confirms no reads on this reference: all metrics are zero
      perRef[t.tid] = agg
      continue

    if st.ok and indexCountsSuffice:
      # Fast path: use pre-computed index statistics
      agg.count = int(st.mapped)
    else:
      # Standard path: iterate and apply filters
      var covered: seq[bool] # allocated by addAlignment
      var depths: seq[int]

      if agg.length <= int(high(int32)):
        # Query by tid: no region-string parsing, which can misread names
        # containing ':' or '-'
        for aln in bam.query(t.tid, 0, agg.length):
          if passesFilters(aln, opts):
            addAlignment(agg, aln, opts, covered, depths)
      else:
        # The wrapper's numeric query takes 32-bit coordinates; braces make
        # htslib treat the whole string as a reference name
        for aln in bam.query("{" & t.name & "}"):
          if passesFilters(aln, opts):
            addAlignment(agg, aln, opts, covered, depths)

      # Calculate metrics from per-base depths; without passing reads they
      # were never allocated and all metrics stay zero
      if opts.trackDepths and depths.len > 0:
        # totalDepth from actual per-base depths (more accurate than the
        # alignment length sum, which is not computed when depths are tracked)
        # Optimization: Consolidated loop for depth metrics
        # Always calculate totalDepth (needed for mean)
        # Conditionally calculate sumSquaredDepth (only for variance)
        for depth in depths:
          agg.totalDepth += int64(depth)
          if opts.doVariance:
            agg.sumSquaredDepth += int64(depth * depth)

        # Calculate trimmed mean if requested
        if opts.doTrimmedMean:
          # Sort a copy: DisCov below needs depths in positional order
          var sortedDepths = depths
          algorithm.sort(sortedDepths)

          let length = sortedDepths.len
          if length > 0:
            # Calculate indices for trimming
            let minIdx = int(float(length) * opts.trimMin / 100.0)
            let maxIdx = int(float(length) * opts.trimMax / 100.0)

            # Calculate mean of trimmed values
            var sum = 0'i64
            var count = 0
            for i in minIdx ..< min(maxIdx, length):
              sum += int64(sortedDepths[i])
              count += 1

            if count > 0:
              agg.trimmedMean = float(sum) / float(count)

        if opts.doDiscov:
          agg.discov = computeDiscov(depths, opts.discovParams).score

    perRef[t.tid] = agg

proc processOneFile(task: WorkerTask) {.thread, gcsafe.} =
  ## Worker thread procedure: process a single BAM file and write results
  ## This runs in a separate thread, so it must not touch global state.
  ## Failures are reported through result.error: an exception escaping a
  ## thread would abort the whole process without naming the file.
  try:
    var bam: Bam
    var fastaCs: cstring = nil
    if task.opts.fasta.len > 0:
      fastaCs = cstring(task.opts.fasta)
    if not open(bam, cstring(task.bamPath), threads = task.opts.threads,
        index = true, fai = fastaCs):
      task.result[].error = "unable to open alignment file"
      return

    if bam.idx == nil:
      task.result[].error = "index not found (run: samtools index)"
      return

    processReferences(bam, task.opts, task.result[].perRef)
  except CatchableError as e:
    task.result[].perRef = @[]
    task.result[].error = e.msg
    return

  # RPKM denominator: reads that passed the filters (-F, -Q, -P), not the
  # index "mapped" total, which includes e.g. supplementary alignments
  var countedReads = 0
  for agg in task.result[].perRef:
    countedReads += agg.count
  task.result[].countedReads = float(countedReads)

proc applySample(res: SampleResult, sampleIdx: int, sampleName: string,
    stored: set[StoredValue], totalCountedReads: var seq[float]) =
  ## Merge a worker's results into the global metricsTable
  ## This runs in the main thread only - no concurrency issues
  totalCountedReads[sampleIdx] = res.countedReads
  var
    missingRefs = 0  # in the first sample's header, not in this one
    matchedRefs = 0
    matchedReads = 0

  # Fast path: when the BAM header has the same references in the same order
  # as the first sample, perRef[i] is row i. Otherwise fall back to a name
  # lookup (built once, on the first mismatch).
  var refLookup: Table[string, int]
  var useLookup = false

  # Process each reference in the global table order
  var rowIdx = 0
  for metrics in metricsTable.mvalues:
    template name: untyped = metrics.refName
    var aggIdx = -1
    if not useLookup and rowIdx < res.perRef.len and res.perRef[rowIdx].name == name:
      aggIdx = rowIdx
    else:
      if not useLookup:
        useLookup = true
        refLookup = initTable[string, int](res.perRef.len)
        for i in 0 ..< res.perRef.len:
          refLookup[res.perRef[i].name] = i
      aggIdx = refLookup.getOrDefault(name, -1)
    inc(rowIdx)

    if aggIdx >= 0:
      template agg: untyped = res.perRef[aggIdx]
      inc(matchedRefs)
      matchedReads += agg.count
      if svCounts in stored: metrics.sampleCounts.add(agg.count)
      if svTotalDepth in stored: metrics.sampleTotalDepth.add(agg.totalDepth)
      if svCoveredBases in stored: metrics.sampleCoveredBases.add(agg.coveredBases)
      if svSumSquaredDepth in stored: metrics.sampleSumSquaredDepth.add(agg.sumSquaredDepth)
      if svTrimmedMean in stored: metrics.sampleTrimmedMean.add(agg.trimmedMean)
      if svDiscov in stored: metrics.sampleDiscov.add(agg.discov)
    else:
      # Reference not present in this sample - fill with zeros
      inc(missingRefs)
      if svCounts in stored: metrics.sampleCounts.add(0)
      if svTotalDepth in stored: metrics.sampleTotalDepth.add(0)
      if svCoveredBases in stored: metrics.sampleCoveredBases.add(0)
      if svSumSquaredDepth in stored: metrics.sampleSumSquaredDepth.add(0)
      if svTrimmedMean in stored: metrics.sampleTrimmedMean.add(0.0)
      if svDiscov in stored: metrics.sampleDiscov.add(0.0)

  # Rows come from the first sample's header, so a sample with different
  # references gets zeros for the missing ones and loses its extra ones
  let extraRefs = res.perRef.len - matchedRefs
  if missingRefs > 0:
    stderr.writeLine("WARNING: ", sampleName, ": ", missingRefs,
        " reference(s) of the first sample are not in its header (reported as 0)")
  if extraRefs > 0:
    stderr.writeLine("WARNING: ", sampleName, ": ", extraRefs,
        " reference(s) not in the first sample's header are not reported (",
        int(res.countedReads) - matchedReads, " reads)")

proc storedValues(metrics: set[Metric]): set[StoredValue] =
  ## Per-sample values the output tables read (see formatValue)
  for metric in metrics:
    case metric
    of mCounts, mRPKM, mTPM, mReadsPerBase: result.incl(svCounts)
    of mMean: result.incl(svTotalDepth)
    of mVariance: result.incl({svTotalDepth, svSumSquaredDepth})
    of mTrimmedMean: result.incl(svTrimmedMean)
    of mCoveredBases, mCoveredRatio: result.incl(svCoveredBases)
    of mDiscov: result.incl(svDiscov)
    of mLength: discard

proc htslibFeatures(): string =
  ## Build features of the loaded htslib (e.g. libdeflate=yes), which explain
  ## most speed differences between installations. Looked up at run time, not
  ## linked: hts_feature_string needs htslib >= 1.11, and static builds have
  ## no shared library to query.
  let lib = loadLib(libname)
  if lib == nil:
    return "unavailable"
  defer: unloadLib(lib)
  let sym = lib.symAddr("hts_feature_string")
  if sym == nil:
    return "unavailable (htslib < 1.11)"
  $cast[proc(): cstring {.cdecl, gcsafe.}](sym)()

proc parseIntOption(value, name: string, minValue: int,
    maxValue = high(int)): int =
  ## Parse an integer option, exiting with an error if invalid or out of range
  ## (unchecked conversions such as uint8(mapq) would silently wrap)
  try:
    result = parseInt(value)
  except ValueError:
    stderr.writeLine("ERROR: ", name, " must be an integer, got: ", value)
    quit(1)
  if result < minValue or result > maxValue:
    if maxValue == high(int):
      stderr.writeLine("ERROR: ", name, " must be >= ", minValue, ", got: ", value)
    else:
      stderr.writeLine("ERROR: ", name, " must be between ", minValue, " and ",
          maxValue, ", got: ", value)
    quit(1)

proc parseFloatOption(value, name: string): float =
  try:
    result = parseFloat(value)
  except ValueError:
    stderr.writeLine("ERROR: ", name, " must be a number, got: ", value)
    quit(1)

proc main(argv: var seq[string]): int =
  let env_fasta = getEnv("REF_PATH")
  let doc = format("""
  BamCountRefs $version

  Usage: bamcountrefs [options]  <BAM-or-CRAM>...

Arguments:

  <BAM-or-CRAM>  the alignment file for which to calculate depth

BAM/CRAM processing options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -W, --workers <workers>      Number of parallel file processors [default: auto]
  -r, --fasta <fasta>          FASTA file for use with CRAM files [default: $env_fasta].
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: $default_flags]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]
  -P, --proper-pairs           If paired flag is set, then also proper pair must be set

Output options:
  -o, --output <BASENAME>      Output file basename (generates multiple files: <BASENAME>_counts.tsv, etc.)
                               If not specified, outputs counts to stdout in TSV format
  -n                           [DEPRECATED: use --rpkm] Output RPKM values
  --rpkm                       Calculate RPKM (reads per kilobase per million reads passing filters)
  --tpm                        Calculate TPM (transcripts per million)
  --mean                       Calculate mean coverage depth (approximate method, no extra memory)
  --trimmed-mean               Calculate trimmed mean coverage (robust against outliers) [requires extra memory]
  --trim-min <FRACTION>        Remove this smallest fraction of positions when calculating trimmed_mean [default: 5]
  --trim-max <FRACTION>        Maximum fraction for trimmed_mean calculations [default: 95]
  --covered-bases              Calculate number of bases with coverage > 0 [requires extra memory]
  --covered-ratio              Calculate coverage breadth (fraction of reference covered) [requires extra memory]
  --variance                   Calculate variance of coverage depth [requires extra memory]
  --reads-per-base             Calculate reads per base (count / length, normalized read density)
  --discov                     Calculate Distribution of Coverage score [requires extra memory]
  --discov-window <INT>        Window length for DisCov spread component [default: 1000]
  --discov-fold-lower <FLOAT>  Lower fold range around median nonzero coverage [default: 0.5]
  --discov-fold-upper <FLOAT>  Upper fold range around median nonzero coverage [default: 2.0]
  --discov-alpha <FLOAT>       Weight on DisCov spread component [default: 0.5]
  --discov-formula <FORMULA>   DisCov formula: linear or geometric [default: linear]
  --length                     Output reference sequence lengths
  -a, --all-metrics            Enable all available metrics

Other options:
  --tag STR                    First column name [default: Contig]
  --multiqc                    Print output as MultiQC table (stdout, or <BASENAME>.mqc.txt with -o)
  --debug                      Enable diagnostics
  -h, --help                   Show help
  """ % ["version", version, "env_fasta", env_fasta, "default_flags",
      $DEFAULT_EXCLUDE_FLAGS])

  let args = docopt(doc, version = version, argv = argv)
  let
    mapq = parseIntOption($args["--mapq"], "--mapq", 0, int(high(uint8)))
    columnName = $args["--tag"]

  # Parse output options
  let
    outputBasename = if $args["--output"] != "nil": $args["--output"] else: ""
    useStdout = outputBasename == ""

  # Parse --workers option (default: auto = min(numFiles, cpuCount))
  let numBamFiles = len(@(args["<BAM-or-CRAM>"]))
  var numWorkers: int
  if $args["--workers"] == "nil" or $args["--workers"] == "auto":
    numWorkers = min(numBamFiles, countProcessors())
  else:
    numWorkers = parseIntOption($args["--workers"], "--workers", 1)

  # Handle deprecated -n flag and metric options
  var
    doRPKM = false
    doTPM = false
    doMean = false
    doTrimmedMean = false
    doCoveredBases = false
    doCoveredRatio = false
    doVariance = false
    doReadsPerBase = false
    doDiscov = false
    doLength = false

  # Trimmed mean parameters
  var
    trimMin = 5.0  # Default: trim bottom 5%
    trimMax = 95.0 # Default: keep up to 95th percentile
    discovParams = DefaultDiscovParams

  if args["-n"]:
    stderr.writeLine("WARNING: -n flag is deprecated, use --rpkm instead")
    doRPKM = true

  if args["--rpkm"]:
    doRPKM = true

  if args["--tpm"]:
    doTPM = true

  if args["--mean"]:
    doMean = true

  if args["--trimmed-mean"]:
    doTrimmedMean = true

  if $args["--trim-min"] != "nil":
    trimMin = parseFloatOption($args["--trim-min"], "--trim-min")

  if $args["--trim-max"] != "nil":
    trimMax = parseFloatOption($args["--trim-max"], "--trim-max")

  if not (trimMin >= 0.0 and trimMin < trimMax and trimMax <= 100.0):
    stderr.writeLine("ERROR: --trim-min and --trim-max must satisfy 0 <= min < max <= 100, got: ",
        trimMin, " and ", trimMax)
    quit(1)

  if args["--covered-bases"]:
    doCoveredBases = true

  if args["--covered-ratio"]:
    doCoveredRatio = true

  if args["--variance"]:
    doVariance = true

  if args["--reads-per-base"]:
    doReadsPerBase = true

  if args["--discov"]:
    doDiscov = true

  discovParams = DiscovParams(
    windowLength: parseIntOption($args["--discov-window"], "--discov-window", 1),
    foldLower: parseFloatOption($args["--discov-fold-lower"], "--discov-fold-lower"),
    foldUpper: parseFloatOption($args["--discov-fold-upper"], "--discov-fold-upper"),
    alpha: parseFloatOption($args["--discov-alpha"], "--discov-alpha"),
    formula: parseDiscovFormula($args["--discov-formula"])
  )
  validateDiscovParams(discovParams, 1)

  if args["--length"]:
    doLength = true

  if args["--all-metrics"]:
    doRPKM = true
    doTPM = true
    doMean = true
    doTrimmedMean = true
    doCoveredBases = true
    doCoveredRatio = true
    doVariance = true
    doReadsPerBase = true
    doDiscov = true
    doLength = true

  debug = args["--debug"]

  var fastaPath = ""
  if $args["--fasta"] != "nil":
    fastaPath = $args["--fasta"]

  var
    eflag = uint16(parseIntOption($args["--flag"], "--flag", 0, int(high(uint16))))
    threads = parseIntOption($args["--threads"], "--threads", 0)
    properPairs = args["--proper-pairs"]

  var
    samples = @[columnName]

  if debug:
    stderr.writeLine("[debug] htslib ", $hts_version(), ": ", htslibFeatures())
    stderr.writeLine("[debug] Processing ", numBamFiles, " BAM file(s) with ",
        numWorkers, " worker(s)")

  # Tables that will be written: stdout holds a single table, so metrics
  # requested but not printed there are not computed
  let singleMetric = multiqcMetric(doRPKM, doTPM, doMean, doTrimmedMean,
      doCoveredBases, doCoveredRatio, doVariance, doReadsPerBase, doDiscov)
  var outputMetrics = {singleMetric}
  if not useStdout:
    outputMetrics.incl(mCounts) # counts are always written
    if doRPKM: outputMetrics.incl(mRPKM)
    if doTPM: outputMetrics.incl(mTPM)
    if doMean: outputMetrics.incl(mMean)
    if doTrimmedMean: outputMetrics.incl(mTrimmedMean)
    if doCoveredBases: outputMetrics.incl(mCoveredBases)
    if doCoveredRatio: outputMetrics.incl(mCoveredRatio)
    if doVariance: outputMetrics.incl(mVariance)
    if doReadsPerBase: outputMetrics.incl(mReadsPerBase)
    if doDiscov: outputMetrics.incl(mDiscov)
    if doLength: outputMetrics.incl(mLength)
  let stored = storedValues(outputMetrics)

  # Determine if we need to track breadth (requires per-base coverage tracking)
  let trackBreadth = outputMetrics * {mCoveredBases, mCoveredRatio} != {}
  # The mean comes from per-base depths when variance, trimmed mean or DisCov
  # are requested, otherwise from alignment spans, which differ for reads
  # extending past the reference end. Keep that choice based on the requested
  # metrics so a printed mean does not depend on which tables are written.
  let meanFromDepths = doVariance or doTrimmedMean or doDiscov
  # Determine if we need to track depths (required for variance, trimmed mean, or DisCov)
  let trackDepths = outputMetrics * {mVariance, mTrimmedMean, mDiscov} != {} or
      (mMean in outputMetrics and meanFromDepths)

  # Prepare worker options
  let workerOpts = WorkerOpts(
    mapq: uint8(mapq),
    eflag: eflag,
    properPairs: properPairs,
    sumSpans: mMean in outputMetrics and not trackDepths,
    trackBreadth: trackBreadth,
    trackDepths: trackDepths,
    doVariance: mVariance in outputMetrics,
    doTrimmedMean: mTrimmedMean in outputMetrics,
    doDiscov: mDiscov in outputMetrics,
    trimMin: trimMin,
    trimMax: trimMax,
    discovParams: discovParams,
    threads: threads,
    fasta: fastaPath
  )

  var results = newSeq[SampleResult](numBamFiles)

  # Create worker tasks
  var bamFiles = @(args["<BAM-or-CRAM>"])
  var workerThreads = newSeq[Thread[WorkerTask]](numBamFiles)
  var tasks = newSeq[WorkerTask](numBamFiles)

  for i, bamFile in bamFiles:
    var sampleName = extractFilename(bamFile).split('.')[0]
    if sampleName in samples:
      # Identical column names would make the table ambiguous downstream
      var suffix = 2
      while sampleName & "_" & $suffix in samples:
        inc(suffix)
      let uniqueName = sampleName & "_" & $suffix
      stderr.writeLine("WARNING: duplicate sample name \"", sampleName,
          "\" for ", bamFile, ", renamed to \"", uniqueName, "\"")
      sampleName = uniqueName
    samples.add(sampleName)

    tasks[i] = WorkerTask(
      bamPath: bamFile,
      sampleIndex: i,
      opts: workerOpts,
      result: addr results[i]
    )

  # Allocate total counted reads array (RPKM denominators)
  var totalCountedReads = newSeq[float](numBamFiles)

  # Run at most numWorkers threads at a time (sliding window): results are
  # merged in input order as each file finishes, then freed, so memory for
  # per-sample results is bounded by the number of workers
  if debug:
    stderr.writeLine("[debug] Starting worker threads...")
  var nextToStart = 0
  while nextToStart < min(numWorkers, numBamFiles):
    createThread(workerThreads[nextToStart], processOneFile, tasks[nextToStart])
    inc(nextToStart)

  if debug:
    stderr.writeLine("[debug] Waiting for workers to complete...")
    stderr.writeLine("[debug] Merging results from all workers...")
  for i in 0 ..< numBamFiles:
    if debug:
      stderr.writeLine("[debug]   Opening BAM/CRAM file ", i)
    joinThread(workerThreads[i])
    if results[i].error.len > 0:
      stderr.writeLine("ERROR: ", bamFiles[i], ": ", results[i].error)
      quit(1)
    if nextToStart < numBamFiles:
      createThread(workerThreads[nextToStart], processOneFile, tasks[nextToStart])
      inc(nextToStart)

    if i == 0:
      # Use first result to establish master reference order
      if debug:
        stderr.writeLine("[debug] Establishing reference order from first sample...")
        stderr.writeLine("[debug] First result has ", results[0].perRef.len, " references")

      if results[0].perRef.len == 0:
        stderr.writeLine("ERROR: First worker returned no references")
        quit(1)

      # perRef is indexed by tid, so it is already in BAM header order
      if debug:
        stderr.writeLine("[debug] Initializing global metrics table with ",
            results[0].perRef.len, " references")
      metricsTable = initOrderedTable[string, ReferenceMetrics](results[0].perRef.len)
      for agg in results[0].perRef:
        var row = ReferenceMetrics(
          refName: agg.name,
          order: agg.order,
          length: agg.length
        )
        if svCounts in stored: row.sampleCounts = newSeqOfCap[int](numBamFiles)
        if svTotalDepth in stored: row.sampleTotalDepth = newSeqOfCap[int64](numBamFiles)
        if svCoveredBases in stored: row.sampleCoveredBases = newSeqOfCap[int](numBamFiles)
        if svSumSquaredDepth in stored: row.sampleSumSquaredDepth = newSeqOfCap[int64](numBamFiles)
        if svTrimmedMean in stored: row.sampleTrimmedMean = newSeqOfCap[float](numBamFiles)
        if svDiscov in stored: row.sampleDiscov = newSeqOfCap[float](numBamFiles)
        metricsTable[agg.name] = row

    if debug:
      stderr.writeLine("[debug]    merging results for sample: \"", samples[
          i+1], "\"")
    applySample(results[i], i, samples[i+1], stored, totalCountedReads)
    results[i] = SampleResult() # free this sample's per-reference results

  # Output results
  let ctx = newOutputContext(totalCountedReads, mTPM in outputMetrics)
  if useStdout:
    # Legacy stdout output (backward compatible): a single metric
    outputToStdout(samples, args["--multiqc"], singleMetric, ctx)
  else:
    # Multi-file output
    outputToFiles(outputBasename, samples, outputMetrics, args["--multiqc"],
        singleMetric, ctx)

  return 0


when isMainModule:
  var args = commandLineParams()
  try:
    discard main(args)
  except EKeyboardInterrupt:
    stderr.writeLine("Quitting.")
    stderr.writeLine( getCurrentExceptionMsg() )
    quit(1)   
