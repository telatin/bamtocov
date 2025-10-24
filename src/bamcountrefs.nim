# Standard library
import std/[os, strutils, tables, sequtils]

# External dependencies
import docopt, hts

const NimblePkgVersion {.strdefine.} = "prerelease"

# Named constants for clarity
const
  DEFAULT_EXCLUDE_FLAGS = 1796  # Unmapped(4) + Secondary(256) + QC-fail(512) + Duplicate(1024) = 1796
  READS_PER_MILLION = 1_000_000
  BASES_PER_KILOBASE = 1000

let
  version = NimblePkgVersion

type
  EKeyboardInterrupt = object of CatchableError

  ReferenceMetrics = object
    refName: string
    order: int              # Preserve BAM header order
    length: int
    # Per-sample raw counts (one per BAM file)
    sampleCounts: seq[int]
    # Per-sample depth tracking for mean coverage (approximate method)
    sampleTotalDepth: seq[int64]
    # Per-sample normalized values (one per BAM file)
    sampleRPKM: seq[float]
    sampleTPM: seq[float]
    sampleMean: seq[float]

var
  # Main data structure: ordered table to preserve BAM order
  metricsTable = initOrderedTable[string, ReferenceMetrics]()
 
proc handler() {.noconv.} =
  raise newException(EKeyboardInterrupt, "Keyboard Interrupt")
 
setControlCHook(handler)


var
  debug = false

proc get_alignments_per_million(bam:Bam): float =
  # Calculate total mapped reads and normalize to "per million"
  for i in bam.hdr.targets:
    result += float(stats(bam.idx,i.tid).mapped)
  result /= READS_PER_MILLION

proc count_alignments_per_ref(bam:Bam, mapq:uint8, eflag:uint16, sampleIndex: int) =
  # Performance optimization: check if we can use BAM index statistics
  # Index stats give us mapped read counts (reads without flag 4 = unmapped)
  # We can only use them when no additional filtering is required:
  #   - mapq == 0: no mapping quality filter
  #   - eflag == 4: only exclude unmapped reads (already excluded by index)
  # For default eflag=DEFAULT_EXCLUDE_FLAGS (excludes secondary, qc-fail, duplicate), we must iterate
  let canUseIndexStats = (mapq == 0) and (eflag == 4)

  if canUseIndexStats and debug:
    stderr.writeLine("[debug] Using fast index statistics (no additional filtering needed)")
  elif debug:
    stderr.writeLine("[debug] Iterating through alignments (mapq=", mapq, ", eflag=", eflag, ")")

  for chromosome in bam.hdr.targets:
    let refName = chromosome.name
    let refLength = int(chromosome.length)
    let refOrder = chromosome.tid
    var rawCounts = 0
    var totalDepth: int64 = 0

    if canUseIndexStats:
      # Fast path: use pre-computed index statistics
      # This is 10-100x faster for large BAM files
      rawCounts = int(stats(bam.idx, chromosome.tid).mapped)
      # Cannot calculate approximate depth with index stats alone
      # Would need to iterate anyway, so depth calculation is skipped for fast path
      # This means --mean will force iteration even if mapq=0 and eflag=4
    else:
      # Standard path: iterate and apply filters
      for aln in bam.query(refName):
        if aln.mapping_quality < mapq: continue
        if (aln.flag and eflag) != 0: continue
        rawCounts += 1
        # Approximate mean coverage: sum alignment lengths
        # This slightly overestimates when reads overlap, but requires no extra memory
        totalDepth += int64(aln.stop - aln.start)

    # Initialize or update metrics table entry
    if refName notin metricsTable:
      # First time seeing this reference - initialize structure
      var metrics = ReferenceMetrics(
        refName: refName,
        order: refOrder,
        length: refLength,
        sampleCounts: @[],
        sampleTotalDepth: @[],
        sampleRPKM: @[],
        sampleTPM: @[],
        sampleMean: @[]
      )
      metricsTable[refName] = metrics

    # Add count and depth for this sample
    metricsTable[refName].sampleCounts.add(rawCounts)
    metricsTable[refName].sampleTotalDepth.add(totalDepth)

proc calculateRPKM(totalMappedReads: seq[float]) =
  # Calculate RPKM for all references and all samples
  # RPKM = (reads × 1,000,000 × 1,000) / (total_mapped_reads × reference_length)
  for refName, metrics in metricsTable.mpairs:
    metrics.sampleRPKM = newSeqOfCap[float](metrics.sampleCounts.len)
    let refLengthKb = metrics.length.float / BASES_PER_KILOBASE.float

    for sampleIdx in 0 ..< metrics.sampleCounts.len:
      let reads = metrics.sampleCounts[sampleIdx].float
      let mappedReadsMillions = totalMappedReads[sampleIdx] / READS_PER_MILLION.float
      let rpkm = reads / (refLengthKb * mappedReadsMillions)
      metrics.sampleRPKM.add(rpkm)

proc calculateTPM() =
  # Calculate TPM for all references and all samples (two-pass algorithm)
  # Pass 1: Calculate RPK (reads per kilobase) for each reference
  # Pass 2: Normalize by sum of RPK and scale to million

  let numSamples = if metricsTable.len > 0: metricsTable.values.toSeq[0].sampleCounts.len else: 0
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

      let tpm = if totalRPK > 0: (rpkValues[refIdx] / totalRPK) * READS_PER_MILLION.float else: 0.0
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

proc outputToStdout(samples: seq[string], multiqc: bool, doRPKM: bool, totalMappedReads: seq[float]) =
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
  else:
    # Default: output counts
    for refName, metrics in metricsTable.pairs:
      var countStrings = newSeq[string](metrics.sampleCounts.len)
      for i in 0 ..< metrics.sampleCounts.len:
        countStrings[i] = $metrics.sampleCounts[i]
      echo refName, "\t", countStrings.join("\t")

proc outputToFiles(basename: string, samples: seq[string], doRPKM: bool, doTPM: bool, doMean: bool, totalMappedReads: seq[float]) =
  # Multi-file output: always write counts, optionally write RPKM, TPM, and mean

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

proc main(argv: var seq[string]): int =
  let env_fasta = getEnv("REF_PATH")
  let doc = format("""
  BamCountRefs $version

  Usage: bamcountrefs [options]  <BAM-or-CRAM>...

Arguments:

  <BAM-or-CRAM>  the alignment file for which to calculate depth

BAM/CRAM processing options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
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
  --all-metrics                Enable all available metrics (RPKM, TPM, and mean)

Other options:
  --tag STR                    First column name [default: ViralSequence]
  --multiqc                    Print output as MultiQC table (stdout only)
  --debug                      Enable diagnostics
  -h, --help                   Show help
  """ % ["version", version, "env_fasta", env_fasta, "default_flags", $DEFAULT_EXCLUDE_FLAGS])

  let args = docopt(doc, version=version, argv=argv)
  let
    mapq = parse_int($args["--mapq"])
    columnName = $args["--tag"]

  # Parse output options
  let
    outputBasename = if $args["--output"] != "nil": $args["--output"] else: ""
    useStdout = outputBasename == ""

  # Handle deprecated -n flag and metric options
  var
    doRPKM = false
    doTPM = false
    doMean = false

  if args["-n"]:
    stderr.writeLine("WARNING: -n flag is deprecated, use --rpkm instead")
    doRPKM = true

  if args["--rpkm"]:
    doRPKM = true

  if args["--tpm"]:
    doTPM = true

  if args["--mean"]:
    doMean = true

  if args["--all-metrics"]:
    doRPKM = true
    doTPM = true
    doMean = true

  debug = args["--debug"]

  var fasta: cstring
  if $args["--fasta"] != "nil":
    fasta = cstring($args["--fasta"])

  var
    eflag = uint16(parse_int($args["--flag"]))
    threads = parse_int($args["--threads"])
    bam:Bam

  var
    samples = @[columnName]

  # Pre-calculate the number of BAM files for memory pre-allocation optimization
  # This allows us to allocate sequences with the correct capacity upfront,
  # avoiding expensive reallocations as we process each sample
  let numBamFiles = len(@(args["<BAM-or-CRAM>"]))

  if debug:
    stderr.writeLine("[debug] Processing ", numBamFiles, " BAM file(s)")

  var sampleIndex = 0
  for bamFile in @(args["<BAM-or-CRAM>"]):
    var sampleName = extractFilename(bamFile)
    # Cache the base name to avoid duplicate string split operations
    let sampleBaseName = sampleName.split('.')[0]
    samples.add(sampleBaseName)
    try:
      open(bam, cstring(bamFile), threads=threads, index=true, fai=fasta)
      if debug:
        stderr.writeLine("Opening BAM/CRAM file: ", bamFile)
    except:
      stderr.writeLine("Unable to open BAM file: ", bamFile )


    if bam.idx == nil:
      stderr.write_line("ERROR: requires BAM/CRAM index")
      quit(1)

    # Count alignments for this sample
    count_alignments_per_ref(bam, uint8(mapq), eflag, sampleIndex)

    sampleIndex += 1

  # Calculate total mapped reads per sample for RPKM calculation
  var totalMappedReads = newSeq[float](numBamFiles)
  var bamIndex = 0
  for bamFile in @(args["<BAM-or-CRAM>"]):
    try:
      open(bam, cstring(bamFile), threads=threads, index=true, fai=fasta)
      totalMappedReads[bamIndex] = bam.get_alignments_per_million() * READS_PER_MILLION.float
      bamIndex += 1
    except:
      stderr.writeLine("Unable to reopen BAM file for totals: ", bamFile)
      quit(1)

  # Output results
  if useStdout:
    # Legacy stdout output (backward compatible)
    # Support both --rpkm and deprecated -n flag for stdout RPKM output
    outputToStdout(samples, args["--multiqc"], doRPKM, totalMappedReads)
  else:
    # Multi-file output
    outputToFiles(outputBasename, samples, doRPKM, doTPM, doMean, totalMappedReads)

  return 0

 
when isMainModule:
  var args = commandLineParams()
  try:
    discard main(args)
  except EKeyboardInterrupt:
    stderr.writeLine( "Quitting.")
  except:
    stderr.writeLine( getCurrentExceptionMsg() )
    quit(1)   