# Standard library
import std/[os, strutils, tables]

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

var
  tableCounts = initTable[string, seq[int]]()
  tableValues = initTable[string, seq[float]]()
  tableLengths = initTable[string, int]()  # Store reference lengths for RPKM calculation

type EKeyboardInterrupt = object of CatchableError
 
proc handler() {.noconv.} =
  raise newException(EKeyboardInterrupt, "Keyboard Interrupt")
 
setControlCHook(handler)


var
  debug = false
  do_rpkm = false

type
  referenceCounts = tuple[refName: string, order: int, length: int, counts: int, value: float]

proc get_alignments_per_million(bam:Bam): float =
  # Calculate total mapped reads and normalize to "per million"
  for i in bam.hdr.targets:
    result += float(stats(bam.idx,i.tid).mapped)
  result /= READS_PER_MILLION

proc count_alignments_per_ref(bam:Bam, mapq:uint8, eflag:uint16, factor: float, expectedSamples: int): seq[referenceCounts] =
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
    var
      chromCounts : referenceCounts
      rawCounts = 0
    chromCounts.refName = chromosome.name
    chromCounts.order   = chromosome.tid
    chromCounts.length   = int(chromosome.length)

    if canUseIndexStats:
      # Fast path: use pre-computed index statistics
      # This is 10-100x faster for large BAM files
      rawCounts = int(stats(bam.idx, chromosome.tid).mapped)
    else:
      # Standard path: iterate and apply filters
      for aln in bam.query(chromosome.name):
        if aln.mapping_quality < mapq: continue
        if (aln.flag and eflag) != 0: continue
        rawCounts += 1

    chromCounts.counts = rawCounts
    # Calculate normalized value once to avoid duplicate float conversion and division
    let normalizedValue = float(rawCounts) / factor
    chromCounts.value = normalizedValue

    # Initialize table entries on first access with pre-allocated capacity
    # This avoids expensive reallocations when processing multiple BAM files
    discard tableCounts.hasKeyOrPut(chromosome.name, newSeqOfCap[int](expectedSamples))
    discard tableValues.hasKeyOrPut(chromosome.name, newSeqOfCap[float](expectedSamples))
    # Store reference length for RPKM calculation (only needs to be set once)
    if chromosome.name notin tableLengths:
      tableLengths[chromosome.name] = chromCounts.length
    tableCounts[chromosome.name].add( rawCounts )
    tableValues[chromosome.name].add( normalizedValue )

    result.add(chromCounts)

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
  -n, --rpkm                   Output RPKM values (reads per kilobase per million)

Other options:
  --tag STR                    First column name [default: ViralSequence]
  --multiqc                    Print output as MultiQC table
  --debug                      Enable diagnostics    
  -h, --help                   Show help
  """ % ["version", version, "env_fasta", env_fasta, "default_flags", $DEFAULT_EXCLUDE_FLAGS])

  let args = docopt(doc, version=version, argv=argv)
  let
    mapq = parse_int($args["--mapq"])
    columnName = $args["--tag"]
  do_rpkm = args["--rpkm"]
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

    # Calculate alignments per million for THIS BAM file
    # IMPORTANT: Each BAM file has its own read count, so normalization must be per-file
    let currentAlignmentsPerMillion = bam.get_alignments_per_million()

    if debug:
      stderr.writeLine("[debug] Sample: ", sampleBaseName,
                       " - alignments per million: ", formatFloat(currentAlignmentsPerMillion, ffDecimal, 3))

    discard count_alignments_per_ref(bam, uint8(mapq), eflag, currentAlignmentsPerMillion, numBamFiles)
    
  
  if args["--multiqc"]:
    echo "# plot_type: 'table'"
    echo "# section_name: 'CovTools count'"
    echo "# description: 'Feature table: counts of mapped reads against predicted viral sequences'"

  echo samples.join("\t")
  for reference in tableCounts.keys:
    if do_rpkm:
      # Calculate proper RPKM: reads per kilobase per million mapped reads
      let refLengthKb = float(tableLengths[reference]) / float(BASES_PER_KILOBASE)
      var rpkmValues: seq[string] = @[]
      for normalizedValue in tableValues[reference]:
        let rpkm = normalizedValue / refLengthKb
        rpkmValues.add(formatFloat(rpkm, ffDecimal, 3))
      echo reference, "\t", rpkmValues.join("\t")
    else:
      echo reference, "\t", tableCounts[reference].join("\t")

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