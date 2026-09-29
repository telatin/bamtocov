# Standard library
import std/[os, strutils, tables, cpuinfo, algorithm]

# External dependencies
import docopt, hts

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
    mappedTotal: float # Total mapped reads for RPKM

  WorkerOpts = object
    ## Immutable options passed to each worker thread
    mapq: uint8
    eflag: uint16
    properPairs: bool # If true, require proper pair flag for paired reads
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
    totalMappedReads: seq[float] # RPKM denominator
    totalRPK: seq[float]         # TPM denominator (sum of reads per kilobase)

const
  metricFileSuffix: array[Metric, string] = ["_counts.tsv", "_rpkm.tsv",
      "_tpm.tsv", "_mean.tsv", "_trimmed_mean.tsv", "_covered_bases.tsv",
      "_covered_fraction.tsv", "_variance.tsv", "_reads_per_base.tsv",
      "_discov.tsv", "_length.tsv"]

proc newOutputContext(totalMappedReads: seq[float], doTPM: bool): OutputContext =
  result.totalMappedReads = totalMappedReads
  result.totalRPK = newSeq[float](totalMappedReads.len)
  if doTPM:
    # TPM = RPK / sum(RPK) * 1e6, summed in reference order for each sample
    for metrics in metricsTable.values:
      let refLengthKb = metrics.length.float / BASES_PER_KILOBASE.float
      for sampleIdx in 0 ..< metrics.sampleCounts.len:
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
        # RPKM = (reads × 1,000,000 × 1,000) / (total_mapped_reads × reference_length)
        let refLengthKb = length / BASES_PER_KILOBASE.float
        let mappedReadsMillions = ctx.totalMappedReads[sampleIdx] /
            READS_PER_MILLION.float
        metrics.sampleCounts[sampleIdx].float / (refLengthKb * mappedReadsMillions)
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
  for refName, metrics in metricsTable.pairs:
    line.setLen(0)
    line.add(refName)
    for sampleIdx in 0 ..< metrics.sampleCounts.len:
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

proc processOneFile(task: WorkerTask) {.thread, gcsafe.} =
  ## Worker thread procedure: process a single BAM file and write results
  ## This runs in a separate thread, so it must not touch global state
  var bam: Bam
  var fastaCs: cstring = nil
  if task.opts.fasta.len > 0:
    fastaCs = cstring(task.opts.fasta)
  if not open(bam, cstring(task.bamPath), threads = task.opts.threads,
      index = true, fai = fastaCs):
    stderr.writeLine("ERROR: Unable to open BAM file in worker: ", task.bamPath)
    return

  if bam.idx == nil:
    stderr.writeLine("ERROR: BAM file requires index: ", task.bamPath)
    return

  # hdr.targets builds a new seq (with name strings) on every call
  let targets = bam.hdr.targets

  # Calculate total mapped reads for RPKM denominator
  var totalMapped = 0'u64
  for t in targets:
    totalMapped += stats(bam.idx, t.tid).mapped
  task.result[].mappedTotal = float(totalMapped)

  # Initialize per-reference sequence
  task.result[].perRef = newSeq[RefAggWithName](targets.len)

  # Determine if we can use fast index statistics
  let canUseIndexStats = (task.opts.mapq == 0'u8) and (task.opts.eflag ==
      4'u16) and (not task.opts.trackBreadth) and (not task.opts.trackDepths)

  # Process each reference sequence
  for t in targets:
    var agg: RefAggWithName
    agg.name = t.name
    agg.length = int(t.length)
    agg.order = t.tid

    # Optimization: Check for zero coverage using BAM index
    # If no mapped reads, skip all expensive processing
    let indexMappedCount = stats(bam.idx, t.tid).mapped

    if indexMappedCount == 0:
      # Zero coverage for this reference - skip computation entirely
      agg.count = 0
      agg.totalDepth = 0
      agg.coveredBases = 0
      agg.sumSquaredDepth = 0
      agg.trimmedMean = 0.0
      agg.discov = 0.0
      task.result[].perRef[t.tid] = agg
      continue # Skip to next reference

    if canUseIndexStats:
      # Fast path: use pre-computed index statistics
      agg.count = int(indexMappedCount)
      # Note: depth and breadth cannot be calculated from index stats alone
    else:
      # Standard path: iterate and apply filters
      var covered: seq[bool]
      var depths: seq[int]
      if task.opts.trackBreadth:
        covered = newSeq[bool](agg.length)
      if task.opts.trackDepths:
        depths = newSeq[int](agg.length)

      for aln in bam.query(t.name):
        if aln.mapping_quality < task.opts.mapq: continue
        if (aln.flag and task.opts.eflag) != 0: continue

        # If --proper-pairs is set, skip paired reads that aren't properly paired
        # Flag 0x1 (1) = read is paired
        # Flag 0x2 (2) = read is properly paired
        if task.opts.properPairs:
          if (aln.flag and 1) != 0: # Read is paired
            if (aln.flag and 2) == 0: # But NOT properly paired
              continue # Skip this read

        inc(agg.count)

        # Approximate mean coverage: sum alignment lengths
        agg.totalDepth += int64(aln.stop - aln.start)

        # Optimization: Combined tracking of breadth and depths in single pass
        # When both metrics are requested, this avoids iterating twice through all positions
        if task.opts.trackBreadth or task.opts.trackDepths:
          for pos in aln.start ..< aln.stop:
            if pos < agg.length:
              if task.opts.trackBreadth and not covered[pos]:
                covered[pos] = true
                inc(agg.coveredBases)
              if task.opts.trackDepths:
                inc(depths[pos])

      # Calculate metrics from per-base depths
      if task.opts.trackDepths:
        # Recalculate totalDepth from actual per-base depths (more accurate than alignment length sum)
        agg.totalDepth = 0 # Reset to use accurate calculation

        # Optimization: Consolidated loop for depth metrics
        # Always calculate totalDepth (needed for mean)
        # Conditionally calculate sumSquaredDepth (only for variance)
        for depth in depths:
          agg.totalDepth += int64(depth)
          if task.opts.doVariance:
            agg.sumSquaredDepth += int64(depth * depth)

        # Calculate trimmed mean if requested
        if task.opts.doTrimmedMean:
          # Sort depths to find percentiles
          var sortedDepths = depths
          algorithm.sort(sortedDepths)

          let length = sortedDepths.len
          if length > 0:
            # Calculate indices for trimming
            let minIdx = int(float(length) * task.opts.trimMin / 100.0)
            let maxIdx = int(float(length) * task.opts.trimMax / 100.0)

            # Calculate mean of trimmed values
            var sum = 0'i64
            var count = 0
            for i in minIdx ..< min(maxIdx, length):
              sum += int64(sortedDepths[i])
              count += 1

            if count > 0:
              agg.trimmedMean = float(sum) / float(count)
            else:
              agg.trimmedMean = 0.0
          else:
            agg.trimmedMean = 0.0

        if task.opts.doDiscov:
          agg.discov = computeDiscov(depths, task.opts.discovParams).score

    task.result[].perRef[t.tid] = agg

proc applySample(res: SampleResult, sampleIdx: int, totalMappedReads: var seq[float]) =
  ## Merge a worker's results into the global metricsTable
  ## This runs in the main thread only - no concurrency issues
  totalMappedReads[sampleIdx] = res.mappedTotal

  # Fast path: when the BAM header has the same references in the same order
  # as the first sample, perRef[i] is row i. Otherwise fall back to a name
  # lookup (built once, on the first mismatch).
  var refLookup: Table[string, int]
  var useLookup = false

  # Process each reference in the global table order
  var rowIdx = 0
  for name, metrics in metricsTable.mpairs:
    var aggIdx = -1
    if not useLookup and rowIdx < res.perRef.len and res.perRef[rowIdx].name == name:
      aggIdx = rowIdx
    else:
      if not useLookup:
        useLookup = true
        refLookup = initTable[string, int](res.perRef.len)
        for i, agg in res.perRef:
          refLookup[agg.name] = i
      aggIdx = refLookup.getOrDefault(name, -1)
    inc(rowIdx)

    if aggIdx >= 0:
      template agg: untyped = res.perRef[aggIdx]
      metrics.sampleCounts.add(agg.count)
      metrics.sampleTotalDepth.add(agg.totalDepth)
      metrics.sampleCoveredBases.add(agg.coveredBases)
      metrics.sampleSumSquaredDepth.add(agg.sumSquaredDepth)
      metrics.sampleTrimmedMean.add(agg.trimmedMean)
      metrics.sampleDiscov.add(agg.discov)
    else:
      # Reference not present in this sample - fill with zeros
      metrics.sampleCounts.add(0)
      metrics.sampleTotalDepth.add(0)
      metrics.sampleCoveredBases.add(0)
      metrics.sampleSumSquaredDepth.add(0)
      metrics.sampleTrimmedMean.add(0.0)
      metrics.sampleDiscov.add(0.0)

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
  --rpkm                       Calculate RPKM (reads per kilobase per million mapped reads)
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
    mapq = parse_int($args["--mapq"])
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
    numWorkers = parse_int($args["--workers"])
    if numWorkers < 1:
      numWorkers = 1

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
    trimMin = parseFloat($args["--trim-min"])

  if $args["--trim-max"] != "nil":
    trimMax = parseFloat($args["--trim-max"])

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
    windowLength: parse_int($args["--discov-window"]),
    foldLower: parseFloat($args["--discov-fold-lower"]),
    foldUpper: parseFloat($args["--discov-fold-upper"]),
    alpha: parseFloat($args["--discov-alpha"]),
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
    eflag = uint16(parse_int($args["--flag"]))
    threads = parse_int($args["--threads"])
    properPairs = args["--proper-pairs"]

  var
    samples = @[columnName]

  if debug:
    stderr.writeLine("[debug] Processing ", numBamFiles, " BAM file(s) with ",
        numWorkers, " worker(s)")

  # Determine if we need to track breadth (requires per-base coverage tracking)
  let trackBreadth = doCoveredBases or doCoveredRatio
  # Determine if we need to track depths (required for variance, trimmed mean, or DisCov)
  let trackDepths = doVariance or doTrimmedMean or doDiscov

  # Prepare worker options
  let workerOpts = WorkerOpts(
    mapq: uint8(mapq),
    eflag: eflag,
    properPairs: properPairs,
    trackBreadth: trackBreadth,
    trackDepths: trackDepths,
    doVariance: doVariance,
    doTrimmedMean: doTrimmedMean,
    doDiscov: doDiscov,
    trimMin: trimMin,
    trimMax: trimMax,
    discovParams: discovParams,
    threads: threads,
    fasta: fastaPath
  )

  # Pre-allocate results array with proper initialization
  var results = newSeq[SampleResult](numBamFiles)
  for i in 0 ..< numBamFiles:
    results[i] = SampleResult(
      perRef: @[],
      mappedTotal: 0.0
    )

  # Create worker tasks
  var bamFiles = @(args["<BAM-or-CRAM>"])
  var workerThreads = newSeq[Thread[WorkerTask]](numBamFiles)
  var tasks = newSeq[WorkerTask](numBamFiles)

  for i, bamFile in bamFiles:
    var sampleName = extractFilename(bamFile)
    let sampleBaseName = sampleName.split('.')[0]
    samples.add(sampleBaseName)

    tasks[i] = WorkerTask(
      bamPath: bamFile,
      sampleIndex: i,
      opts: workerOpts,
      result: addr results[i]
    )

  # Allocate total mapped reads array
  var totalMappedReads = newSeq[float](numBamFiles)

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
        metricsTable[agg.name] = ReferenceMetrics(
          refName: agg.name,
          order: agg.order,
          length: agg.length,
          sampleCounts: newSeqOfCap[int](numBamFiles),
          sampleTotalDepth: newSeqOfCap[int64](numBamFiles),
          sampleCoveredBases: newSeqOfCap[int](numBamFiles),
          sampleSumSquaredDepth: newSeqOfCap[int64](numBamFiles),
          sampleTrimmedMean: newSeqOfCap[float](numBamFiles),
          sampleDiscov: newSeqOfCap[float](numBamFiles)
        )

    if debug:
      stderr.writeLine("[debug]    merging results for sample: \"", samples[
          i+1], "\"")
    applySample(results[i], i, totalMappedReads)
    results[i] = SampleResult() # free this sample's per-reference results

  # Output results
  let
    ctx = newOutputContext(totalMappedReads, doTPM)
    singleMetric = multiqcMetric(doRPKM, doTPM, doMean, doTrimmedMean,
        doCoveredBases, doCoveredRatio, doVariance, doReadsPerBase, doDiscov)
  if useStdout:
    # Legacy stdout output (backward compatible): a single metric
    outputToStdout(samples, args["--multiqc"], singleMetric, ctx)
  else:
    # Multi-file output
    var metrics = {mCounts}
    if doRPKM: metrics.incl(mRPKM)
    if doTPM: metrics.incl(mTPM)
    if doMean: metrics.incl(mMean)
    if doTrimmedMean: metrics.incl(mTrimmedMean)
    if doCoveredBases: metrics.incl(mCoveredBases)
    if doCoveredRatio: metrics.incl(mCoveredRatio)
    if doVariance: metrics.incl(mVariance)
    if doReadsPerBase: metrics.incl(mReadsPerBase)
    if doDiscov: metrics.incl(mDiscov)
    if doLength: metrics.incl(mLength)
    outputToFiles(outputBasename, samples, metrics, args["--multiqc"],
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
