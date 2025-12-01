# Standard library
import std/[os, strutils, tables, sequtils, cpuinfo, algorithm]

# External dependencies
import docopt, hts

const NimblePkgVersion {.strdefine.} = "prerelease"

# Named constants for clarity
const
  DEFAULT_EXCLUDE_FLAGS = 1796 # Unmapped(4) + Secondary(256) + QC-fail(512) + Duplicate(1024) = 1796
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
    # Per-sample normalized values (one per BAM file)
    sampleRPKM: seq[float]
    sampleTPM: seq[float]
    sampleMean: seq[float]
    sampleCoveredRatio: seq[float]
    sampleVariance: seq[float]
    sampleTrimmedMean: seq[float]
    sampleReadsPerBase: seq[float]

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

  SampleResult = object
    ## Worker return value containing all metrics for one sample
    ## Uses sequences instead of Tables for thread-safety
    perRef: seq[RefAggWithName] # All reference metrics
    mappedTotal: float # Total mapped reads for RPKM

  WorkerOpts = object
    ## Immutable options passed to each worker thread
    mapq: uint8
    eflag: uint16
    trackBreadth: bool
    trackDepths: bool  # Track per-base depths for variance and/or trimmed mean
    doVariance: bool
    doTrimmedMean: bool
    trimMin: float  # Minimum percentile for trimmed mean (0-100)
    trimMax: float  # Maximum percentile for trimmed mean (0-100)
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

proc calculateRPKM(totalMappedReads: seq[float]) =
  # Calculate RPKM for all references and all samples
  # RPKM = (reads × 1,000,000 × 1,000) / (total_mapped_reads × reference_length)
  for refName, metrics in metricsTable.mpairs:
    metrics.sampleRPKM = newSeqOfCap[float](metrics.sampleCounts.len)
    let refLengthKb = metrics.length.float / BASES_PER_KILOBASE.float

    for sampleIdx in 0 ..< metrics.sampleCounts.len:
      let reads = metrics.sampleCounts[sampleIdx].float
      let mappedReadsMillions = totalMappedReads[sampleIdx] /
          READS_PER_MILLION.float
      let rpkm = reads / (refLengthKb * mappedReadsMillions)
      metrics.sampleRPKM.add(rpkm)

proc calculateTPM() =
  # Calculate TPM for all references and all samples (two-pass algorithm)
  # Pass 1: Calculate RPK (reads per kilobase) for each reference
  # Pass 2: Normalize by sum of RPK and scale to million

  let numSamples = if metricsTable.len > 0: metricsTable.values.toSeq[
      0].sampleCounts.len else: 0
  if numSamples == 0:
    return

  # For each sample, calculate TPM independently
  for sampleIdx in 0 ..< numSamples:
    # Pass 1: Calculate RPK for all references in this sample
    var rpkValues = newSeq[float](metricsTable.len)
    var totalRPK = 0.0
    var refIdx = 0

    for refName, metrics in metricsTable.pairs:
      let reads = metrics.sampleCounts[sampleIdx].float
      let refLengthKb = metrics.length.float / BASES_PER_KILOBASE.float
      let rpk = reads / refLengthKb
      rpkValues[refIdx] = rpk
      totalRPK += rpk
      refIdx += 1

    # Pass 2: Normalize and scale to million
    refIdx = 0
    for refName, metrics in metricsTable.mpairs:
      # Initialize TPM sequence on first sample
      if sampleIdx == 0:
        metrics.sampleTPM = newSeqOfCap[float](numSamples)

      let tpm = if totalRPK > 0: (rpkValues[refIdx] / totalRPK) *
          READS_PER_MILLION.float else: 0.0
      metrics.sampleTPM.add(tpm)
      refIdx += 1

proc calculateMean() =
  # Calculate mean coverage for all references and all samples
  # Mean = total_depth / reference_length
  # Using approximate method: total_depth = sum of alignment lengths
  for refName, metrics in metricsTable.mpairs:
    metrics.sampleMean = newSeqOfCap[float](metrics.sampleTotalDepth.len)

    for sampleIdx in 0 ..< metrics.sampleTotalDepth.len:
      let totalDepth = metrics.sampleTotalDepth[sampleIdx].float
      let meanCoverage = totalDepth / metrics.length.float
      metrics.sampleMean.add(meanCoverage)

proc calculateCoveredRatio() =
  # Calculate coverage breadth (covered bases ratio) for all references and all samples
  # Ratio = covered_bases / reference_length
  for refName, metrics in metricsTable.mpairs:
    metrics.sampleCoveredRatio = newSeqOfCap[float](
        metrics.sampleCoveredBases.len)

    for sampleIdx in 0 ..< metrics.sampleCoveredBases.len:
      let coveredBases = metrics.sampleCoveredBases[sampleIdx].float
      let coveredRatio = coveredBases / metrics.length.float
      metrics.sampleCoveredRatio.add(coveredRatio)

proc calculateVariance() =
  # Calculate coverage variance for all references and all samples
  # Variance = E[X^2] - (E[X])^2
  # Where E[X^2] = sum_squared_depth / length
  # And E[X] = total_depth / length (mean coverage)
  for refName, metrics in metricsTable.mpairs:
    metrics.sampleVariance = newSeqOfCap[float](metrics.sampleSumSquaredDepth.len)

    for sampleIdx in 0 ..< metrics.sampleSumSquaredDepth.len:
      let sumSquaredDepth = metrics.sampleSumSquaredDepth[sampleIdx].float
      let totalDepth = metrics.sampleTotalDepth[sampleIdx].float
      let length = metrics.length.float

      # E[X^2] - (E[X])^2
      let meanSquared = (totalDepth / length) * (totalDepth / length)
      let meanOfSquares = sumSquaredDepth / length
      let variance = meanOfSquares - meanSquared

      # Variance should never be negative, but due to floating point errors it might be
      metrics.sampleVariance.add(if variance >= 0: variance else: 0.0)

proc calculateReadsPerBase() =
  # Calculate reads per base for all references and all samples
  # ReadsPerBase = count / length (normalized read density)
  for refName, metrics in metricsTable.mpairs:
    metrics.sampleReadsPerBase = newSeqOfCap[float](metrics.sampleCounts.len)

    for sampleIdx in 0 ..< metrics.sampleCounts.len:
      let count = metrics.sampleCounts[sampleIdx].float
      let readsPerBase = count / metrics.length.float
      metrics.sampleReadsPerBase.add(readsPerBase)

proc writeTableToFile(filename: string, samples: seq[string], values: seq[seq[string]]) =
  # Write a table to file with header
  var file: File
  if not open(file, filename, fmWrite):
    stderr.writeLine("ERROR: Unable to write to file: ", filename)
    quit(1)

  # Write header
  file.writeLine(samples.join("\t"))

  # Write data rows (preserving order from metricsTable which is OrderedTable)
  var rowIdx = 0
  for refName in metricsTable.keys:
    file.writeLine(refName, "\t", values[rowIdx].join("\t"))
    rowIdx += 1

  file.close()
  if debug:
    stderr.writeLine("[debug] Wrote output to: ", filename)

proc outputToStdout(samples: seq[string], multiqc: bool,
                    doRPKM: bool, doTPM: bool, doMean: bool, doTrimmedMean: bool,
                    doCoveredBases: bool, doCoveredRatio: bool, doVariance: bool,
                    doReadsPerBase: bool, doLength: bool,
                    totalMappedReads: seq[float]) =
  # Legacy stdout output for backward compatibility
  if multiqc:
    echo "# plot_type: 'table'"
    echo "# section_name: 'CovTools count'"
    echo "# description: 'Feature table: counts of mapped reads against predicted viral sequences'"

  echo samples.join("\t")

  if doRPKM:
    # Output RPKM to stdout (both --rpkm and deprecated -n flag)
    calculateRPKM(totalMappedReads)
    for refName, metrics in metricsTable.pairs:
      var rpkmStrings = newSeq[string](metrics.sampleRPKM.len)
      for i in 0 ..< metrics.sampleRPKM.len:
        rpkmStrings[i] = formatFloat(metrics.sampleRPKM[i], ffDecimal, 3)
      echo refName, "\t", rpkmStrings.join("\t")
  elif doTPM:
    calculateTPM()
    for refName, metrics in metricsTable.pairs:
      var tpmStrings = newSeq[string](metrics.sampleTPM.len)
      for i in 0 ..< metrics.sampleTPM.len:
        tpmStrings[i] = formatFloat(metrics.sampleTPM[i], ffDecimal, 3)
      echo refName, "\t", tpmStrings.join("\t")
  elif doMean:
    calculateMean()
    for refName, metrics in metricsTable.pairs:
      var meanStrings = newSeq[string](metrics.sampleMean.len)
      for i in 0 ..< metrics.sampleMean.len:
        meanStrings[i] = formatFloat(metrics.sampleMean[i], ffDecimal, 6)
      echo refName, "\t", meanStrings.join("\t")
  elif doTrimmedMean:
    # Trimmed mean values are already calculated in worker threads
    for refName, metrics in metricsTable.pairs:
      var trimmedMeanStrings = newSeq[string](metrics.sampleTrimmedMean.len)
      for i in 0 ..< metrics.sampleTrimmedMean.len:
        trimmedMeanStrings[i] = formatFloat(metrics.sampleTrimmedMean[i], ffDecimal, 6)
      echo refName, "\t", trimmedMeanStrings.join("\t")
  elif doCoveredRatio:
    calculateCoveredRatio()
    for refName, metrics in metricsTable.pairs:
      var ratioStrings = newSeq[string](metrics.sampleCoveredRatio.len)
      for i in 0 ..< metrics.sampleCoveredRatio.len:
        ratioStrings[i] = formatFloat(metrics.sampleCoveredRatio[i], ffDecimal, 6)
      echo refName, "\t", ratioStrings.join("\t")
  elif doCoveredBases:
    for refName, metrics in metricsTable.pairs:
      var basesStrings = newSeq[string](metrics.sampleCoveredBases.len)
      for i in 0 ..< metrics.sampleCoveredBases.len:
        basesStrings[i] = $metrics.sampleCoveredBases[i]
      echo refName, "\t", basesStrings.join("\t")
  elif doVariance:
    calculateVariance()
    for refName, metrics in metricsTable.pairs:
      var varianceStrings = newSeq[string](metrics.sampleVariance.len)
      for i in 0 ..< metrics.sampleVariance.len:
        varianceStrings[i] = formatFloat(metrics.sampleVariance[i], ffDecimal, 6)
      echo refName, "\t", varianceStrings.join("\t")
  elif doReadsPerBase:
    calculateReadsPerBase()
    for refName, metrics in metricsTable.pairs:
      var rpbStrings = newSeq[string](metrics.sampleReadsPerBase.len)
      for i in 0 ..< metrics.sampleReadsPerBase.len:
        rpbStrings[i] = formatFloat(metrics.sampleReadsPerBase[i], ffDecimal, 6)
      echo refName, "\t", rpbStrings.join("\t")
  else:
    # Default: output counts
    for refName, metrics in metricsTable.pairs:
      var countStrings = newSeq[string](metrics.sampleCounts.len)
      for i in 0 ..< metrics.sampleCounts.len:
        countStrings[i] = $metrics.sampleCounts[i]
      echo refName, "\t", countStrings.join("\t")

proc outputToFiles(basename: string, samples: seq[string], doRPKM: bool,
    doTPM: bool, doMean: bool, doTrimmedMean: bool, doCoveredBases: bool, doCoveredRatio: bool,
    doVariance: bool, doReadsPerBase: bool, doLength: bool, totalMappedReads: seq[float]) =
  # Multi-file output: always write counts, optionally write RPKM, TPM, mean, trimmed mean, covered bases, covered ratio, variance, reads-per-base, and length

  # Always write counts
  var countsValues = newSeq[seq[string]](metricsTable.len)
  var idx = 0
  for refName, metrics in metricsTable.pairs:
    countsValues[idx] = newSeq[string](metrics.sampleCounts.len)
    for i in 0 ..< metrics.sampleCounts.len:
      countsValues[idx][i] = $metrics.sampleCounts[i]
    idx += 1

  writeTableToFile(basename & "_counts.tsv", samples, countsValues)

  # Optionally write RPKM
  if doRPKM:
    calculateRPKM(totalMappedReads)
    var rpkmValues = newSeq[seq[string]](metricsTable.len)
    idx = 0
    for refName, metrics in metricsTable.pairs:
      rpkmValues[idx] = newSeq[string](metrics.sampleRPKM.len)
      for i in 0 ..< metrics.sampleRPKM.len:
        rpkmValues[idx][i] = formatFloat(metrics.sampleRPKM[i], ffDecimal, 3)
      idx += 1

    writeTableToFile(basename & "_rpkm.tsv", samples, rpkmValues)

  # Optionally write TPM
  if doTPM:
    calculateTPM()
    var tpmValues = newSeq[seq[string]](metricsTable.len)
    idx = 0
    for refName, metrics in metricsTable.pairs:
      tpmValues[idx] = newSeq[string](metrics.sampleTPM.len)
      for i in 0 ..< metrics.sampleTPM.len:
        tpmValues[idx][i] = formatFloat(metrics.sampleTPM[i], ffDecimal, 3)
      idx += 1

    writeTableToFile(basename & "_tpm.tsv", samples, tpmValues)

  # Optionally write mean coverage
  if doMean:
    calculateMean()
    var meanValues = newSeq[seq[string]](metricsTable.len)
    idx = 0
    for refName, metrics in metricsTable.pairs:
      meanValues[idx] = newSeq[string](metrics.sampleMean.len)
      for i in 0 ..< metrics.sampleMean.len:
        meanValues[idx][i] = formatFloat(metrics.sampleMean[i], ffDecimal, 6)
      idx += 1

    writeTableToFile(basename & "_mean.tsv", samples, meanValues)

  # Optionally write trimmed mean coverage
  if doTrimmedMean:
    # Trimmed mean values are already calculated in worker threads
    var trimmedMeanValues = newSeq[seq[string]](metricsTable.len)
    idx = 0
    for refName, metrics in metricsTable.pairs:
      trimmedMeanValues[idx] = newSeq[string](metrics.sampleTrimmedMean.len)
      for i in 0 ..< metrics.sampleTrimmedMean.len:
        trimmedMeanValues[idx][i] = formatFloat(metrics.sampleTrimmedMean[i], ffDecimal, 6)
      idx += 1

    writeTableToFile(basename & "_trimmed_mean.tsv", samples, trimmedMeanValues)

  # Optionally write covered bases (raw count)
  if doCoveredBases:
    var coveredBasesValues = newSeq[seq[string]](metricsTable.len)
    idx = 0
    for refName, metrics in metricsTable.pairs:
      coveredBasesValues[idx] = newSeq[string](metrics.sampleCoveredBases.len)
      for i in 0 ..< metrics.sampleCoveredBases.len:
        coveredBasesValues[idx][i] = $metrics.sampleCoveredBases[i]
      idx += 1

    writeTableToFile(basename & "_covered_bases.tsv", samples, coveredBasesValues)

  # Optionally write covered ratio (breadth as fraction)
  if doCoveredRatio:
    calculateCoveredRatio()
    var coveredRatioValues = newSeq[seq[string]](metricsTable.len)
    idx = 0
    for refName, metrics in metricsTable.pairs:
      coveredRatioValues[idx] = newSeq[string](metrics.sampleCoveredRatio.len)
      for i in 0 ..< metrics.sampleCoveredRatio.len:
        coveredRatioValues[idx][i] = formatFloat(metrics.sampleCoveredRatio[i],
            ffDecimal, 6)
      idx += 1

    writeTableToFile(basename & "_covered_fraction.tsv", samples, coveredRatioValues)

  # Optionally write variance
  if doVariance:
    calculateVariance()
    var varianceValues = newSeq[seq[string]](metricsTable.len)
    idx = 0
    for refName, metrics in metricsTable.pairs:
      varianceValues[idx] = newSeq[string](metrics.sampleVariance.len)
      for i in 0 ..< metrics.sampleVariance.len:
        varianceValues[idx][i] = formatFloat(metrics.sampleVariance[i], ffDecimal, 6)
      idx += 1

    writeTableToFile(basename & "_variance.tsv", samples, varianceValues)

  # Optionally write reads per base
  if doReadsPerBase:
    calculateReadsPerBase()
    var readsPerBaseValues = newSeq[seq[string]](metricsTable.len)
    idx = 0
    for refName, metrics in metricsTable.pairs:
      readsPerBaseValues[idx] = newSeq[string](metrics.sampleReadsPerBase.len)
      for i in 0 ..< metrics.sampleReadsPerBase.len:
        readsPerBaseValues[idx][i] = formatFloat(metrics.sampleReadsPerBase[i], ffDecimal, 6)
      idx += 1

    writeTableToFile(basename & "_reads_per_base.tsv", samples, readsPerBaseValues)

  # Optionally write reference lengths
  if doLength:
    var lengthValues = newSeq[seq[string]](metricsTable.len)
    idx = 0
    for refName, metrics in metricsTable.pairs:
      # Length is a property of the reference, so repeat it for each sample
      lengthValues[idx] = newSeq[string](metrics.sampleCounts.len)
      for i in 0 ..< metrics.sampleCounts.len:
        lengthValues[idx][i] = $metrics.length
      idx += 1

    writeTableToFile(basename & "_length.tsv", samples, lengthValues)

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

  # Calculate total mapped reads for RPKM denominator
  var totalMapped = 0'u64
  for t in bam.hdr.targets:
    totalMapped += stats(bam.idx, t.tid).mapped
  task.result[].mappedTotal = float(totalMapped)

  # Initialize per-reference sequence
  task.result[].perRef = newSeq[RefAggWithName](bam.hdr.targets.len)

  # Determine if we can use fast index statistics
  let canUseIndexStats = (task.opts.mapq == 0'u8) and (task.opts.eflag ==
      4'u16) and (not task.opts.trackBreadth) and (not task.opts.trackDepths)

  # Process each reference sequence
  for t in bam.hdr.targets:
    var agg: RefAggWithName
    agg.name = t.name
    agg.length = int(t.length)
    agg.order = t.tid

    if canUseIndexStats:
      # Fast path: use pre-computed index statistics
      agg.count = int(stats(bam.idx, t.tid).mapped)
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
        inc(agg.count)

        # Approximate mean coverage: sum alignment lengths
        agg.totalDepth += int64(aln.stop - aln.start)

        # Track covered positions for breadth calculation
        if task.opts.trackBreadth:
          for pos in aln.start ..< aln.stop:
            if pos < agg.length and not covered[pos]:
              covered[pos] = true
              inc(agg.coveredBases)

        # Track per-base depths for variance and/or trimmed mean calculation
        if task.opts.trackDepths:
          for pos in aln.start ..< aln.stop:
            if pos < agg.length:
              inc(depths[pos])

      # Calculate metrics from per-base depths
      if task.opts.trackDepths:
        # Recalculate totalDepth from actual per-base depths (more accurate than alignment length sum)
        agg.totalDepth = 0  # Reset to use accurate calculation

        # Calculate sum of squared depths for variance if requested
        if task.opts.doVariance:
          for depth in depths:
            agg.totalDepth += int64(depth)
            agg.sumSquaredDepth += int64(depth * depth)
        else:
          # Just calculate total depth
          for depth in depths:
            agg.totalDepth += int64(depth)

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

    task.result[].perRef[t.tid] = agg

proc applySample(res: SampleResult, sampleIdx: int, totalMappedReads: var seq[float]) =
  ## Merge a worker's results into the global metricsTable
  ## This runs in the main thread only - no concurrency issues
  totalMappedReads[sampleIdx] = res.mappedTotal

  # Create a lookup table for fast access
  var refLookup = initTable[string, RefAggWithName]()
  for agg in res.perRef:
    refLookup[agg.name] = agg

  # Process each reference in the global table order
  for name in metricsTable.keys:
    if refLookup.hasKey(name):
      let agg = refLookup[name]
      metricsTable[name].sampleCounts.add(agg.count)
      metricsTable[name].sampleTotalDepth.add(agg.totalDepth)
      metricsTable[name].sampleCoveredBases.add(agg.coveredBases)
      metricsTable[name].sampleSumSquaredDepth.add(agg.sumSquaredDepth)
      metricsTable[name].sampleTrimmedMean.add(agg.trimmedMean)
    else:
      # Reference not present in this sample - fill with zeros
      metricsTable[name].sampleCounts.add(0)
      metricsTable[name].sampleTotalDepth.add(0)
      metricsTable[name].sampleCoveredBases.add(0)
      metricsTable[name].sampleSumSquaredDepth.add(0)
      metricsTable[name].sampleTrimmedMean.add(0.0)

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
  --length                     Output reference sequence lengths
  --all-metrics                Enable all available metrics

Other options:
  --tag STR                    First column name [default: ViralSequence]
  --multiqc                    Print output as MultiQC table (stdout only)
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
    doLength = false

  # Trimmed mean parameters
  var
    trimMin = 5.0  # Default: trim bottom 5%
    trimMax = 95.0  # Default: keep up to 95th percentile

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
    doLength = true

  debug = args["--debug"]

  var fastaPath = ""
  if $args["--fasta"] != "nil":
    fastaPath = $args["--fasta"]

  var
    eflag = uint16(parse_int($args["--flag"]))
    threads = parse_int($args["--threads"])

  var
    samples = @[columnName]

  if debug:
    stderr.writeLine("[debug] Processing ", numBamFiles, " BAM file(s) with ",
        numWorkers, " worker(s)")

  # Determine if we need to track breadth (requires per-base coverage tracking)
  let trackBreadth = doCoveredBases or doCoveredRatio
  # Determine if we need to track depths (required for variance and/or trimmed mean)
  let trackDepths = doVariance or doTrimmedMean

  # Prepare worker options
  let workerOpts = WorkerOpts(
    mapq: uint8(mapq),
    eflag: eflag,
    trackBreadth: trackBreadth,
    trackDepths: trackDepths,
    doVariance: doVariance,
    doTrimmedMean: doTrimmedMean,
    trimMin: trimMin,
    trimMax: trimMax,
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

  # Create worker threads
  var bamFiles = @(args["<BAM-or-CRAM>"])
  var workerThreads = newSeq[Thread[WorkerTask]](numBamFiles)
  var tasks = newSeq[WorkerTask](numBamFiles)

  if debug:
    stderr.writeLine("[debug] Starting worker threads...")

  # Initialize tasks and start threads
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
    createThread(workerThreads[i], processOneFile, tasks[i])

  # Wait for all threads to complete
  if debug:
    stderr.writeLine("[debug] Waiting for workers to complete...")
  for i in 0 ..< numBamFiles:
    if debug:
      stderr.writeLine("[debug]   Opening BAM/CRAM file ", i)
    joinThread(workerThreads[i])

  # Use first result to establish master reference order
  if debug:
    stderr.writeLine("[debug] Establishing reference order from first sample...")
    stderr.writeLine("[debug] First result has ", results[0].perRef.len, " references")

  if results[0].perRef.len == 0:
    stderr.writeLine("ERROR: First worker returned no references")
    quit(1)

  let firstRes = results[0]

  # Sort references by BAM order (tid) and initialize global metricsTable
  var sortedRefs = firstRes.perRef
  algorithm.sort(sortedRefs, proc(a, b: RefAggWithName): int = cmp(a.order, b.order))

  if debug:
    stderr.writeLine("[debug] Initializing global metrics table with ",
        sortedRefs.len, " references")
  metricsTable = initOrderedTable[string, ReferenceMetrics]()
  for agg in sortedRefs:
    metricsTable[agg.name] = ReferenceMetrics(
      refName: agg.name,
      order: agg.order,
      length: agg.length,
      sampleCounts: newSeqOfCap[int](numBamFiles),
      sampleTotalDepth: newSeqOfCap[int64](numBamFiles),
      sampleCoveredBases: newSeqOfCap[int](numBamFiles),
      sampleSumSquaredDepth: newSeqOfCap[int64](numBamFiles),
      sampleRPKM: @[],
      sampleTPM: @[],
      sampleMean: @[],
      sampleCoveredRatio: @[],
      sampleVariance: @[],
      sampleTrimmedMean: @[],
      sampleReadsPerBase: @[]
    )

  # Allocate total mapped reads array
  var totalMappedReads = newSeq[float](numBamFiles)

  # Apply all sample results
  if debug:
    stderr.writeLine("[debug] Merging results from all workers...")
  for i in 0 ..< numBamFiles:
    if debug:
      stderr.writeLine("[debug]    merging results for sample: \"", samples[
          i+1], "\"")
    applySample(results[i], i, totalMappedReads)

  # Output results
  if useStdout:
    # Legacy stdout output (backward compatible)
    # Support both --rpkm and deprecated -n flag for stdout RPKM output
    outputToStdout(samples, args["--multiqc"], doRPKM, doTPM, doMean, doTrimmedMean,
        doCoveredBases, doCoveredRatio, doVariance, doReadsPerBase, doLength, totalMappedReads)
  else:
    # Multi-file output
    outputToFiles(outputBasename, samples, doRPKM, doTPM, doMean, doTrimmedMean,
        doCoveredBases, doCoveredRatio, doVariance, doReadsPerBase, doLength, totalMappedReads)

  return 0


when isMainModule:
  var args = commandLineParams()
  try:
    discard main(args)
  except EKeyboardInterrupt:
    stderr.writeLine("Quitting.")
  except:
    stderr.writeLine(getCurrentExceptionMsg())
    quit(1)
