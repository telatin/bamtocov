# Standard library
import std/[os, strutils, tables]

# External dependencies
import docopt, hts
when defined(threads):
  import taskpools

# Local modules
import ./covutils

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

  # Basic type definitions for internal use
  stranded_counts_t = tuple[fwd, rev: int]
  feature_coords_t  = tuple[chrom, starts, stops, name: string, length: int]
  counts_t*         = TableRef[interval_t[string], int]
  frozen_region_t   = object
    chrom: string
    start: int
    stop: int
    name: string
  frozen_target_t   = seq[tuple[chrom: string, intervals: seq[frozen_region_t]]]

  FeatureMetrics = object
    ## Metrics for a single feature across all samples
    featureName: string
    chrom: string
    starts: string      # Semicolon-separated for multi-exon features
    stops: string       # Semicolon-separated for multi-exon features
    length: int         # Total length (sum of all exons)
    # Per-sample raw counts (one per BAM file)
    sampleCounts: seq[stranded_counts_t]
    # Per-sample normalized values
    sampleRPKM: seq[float]
    sampleNormLen: seq[float]  # counts per base (reads per base)

  ProcessingOptions = object
    ## Command-line options affecting BAM processing
    mapq: uint8
    eflag: uint16
    threads: int
    jobs: int          # Number of BAM files to process concurrently
    fasta: string
    strict: bool        # Require full containment within feature
    paired: bool        # Count pairs instead of individual reads
    stranded: bool      # Track forward/reverse strand counts separately
    do_rpkm: bool       # Calculate RPKM
    do_norm: bool       # Calculate counts/length
    do_coords: bool     # Output feature coordinates
    debug: bool
    verbose: bool

  TargetOptions = object
    parseAsGff: bool
    gffIdentifier: string
    gffSeparator: string
    gffField: string

  SampleResult = object
    samplePath: string
    counts: OrderedTable[string, stranded_counts_t]
    perMillion: float
    error: string

proc handler() {.noconv.} =
  raise newException(EKeyboardInterrupt, "Keyboard Interrupt")

setControlCHook(handler)

# Helper procedures for stranded counts
proc inc(c: var stranded_counts_t, reverse=false) =
  if reverse == false:
    c.fwd += 1
  else:
    c.rev += 1

proc counts(c: stranded_counts_t): int =
  c.fwd + c.rev

proc countsToString(c: stranded_counts_t, stranded: bool): string =
  if stranded:
    $(c.fwd) & "\t" & $(c.rev)
  else:
    $(c.fwd + c.rev)

proc makeCountsTable(table: var OrderedTable[string, stranded_counts_t], bam:Bam, mapq:uint8, eflag:uint16, regions: target_t, strict = false): float =
  var total: float = 0
  for read in bam:
    total += 1
    if read.tid notin regions:
      continue
    if read.mapping_quality < mapq or ( (read.flag and eflag) != 0):
      continue

    for region in regions[read.tid]:
      
      if (read.start < region.start   and read.stop > region.stop) or (read.stop > region.start   and read.stop < region.stop) or (read.start > region.start  and read.start < region.stop):
        if strict and ( read.start < region.start or read.stop > region.stop ):
            continue
        table[region.label].inc(read.flag.reverse)
  return total / READS_PER_MILLION.float


proc alignments_count(table: var OrderedTable[string, stranded_counts_t], bam:Bam, mapq:uint8, eflag:uint16, regions: target_t, strict = false, paired = false): float =
  var total: float = 0
  for read in bam:
    total += 1
    if read.tid notin regions:
      continue
    if read.mapping_quality < mapq or ( (read.flag and eflag) != 0):
      continue
    if paired and (read.flag.proper_pair == false or read.flag.read2 == true):
      continue


    for region in regions[read.tid]:
      
      if (read.start < region.start   and read.stop > region.stop) or (read.stop > region.start   and read.stop < region.stop) or (read.start > region.start  and read.start < region.stop):
        if strict and ( read.start < region.start or read.stop > region.stop ):
            continue
        table[region.label].inc(read.flag.reverse)
  return total / READS_PER_MILLION.float

proc add(s: var feature_coords_t, z: feature_coords_t) =
  # feature_coords_t = tuple[chrom, starts, stops, name: string, length: int]
  s.chrom &= ";" & z.chrom
  s.starts &= ";" & z.starts
  s.stops &= ";" & z.stops
  s.length += z.length

proc rpkm(f: feature_coords_t, alignmentsPerMillion: float, c: stranded_counts_t): float =
  # RPKM = (reads × 1,000,000 × 1,000) / (total_mapped_reads × feature_length)
  let kb: float = f.length.float / BASES_PER_KILOBASE.float
  if alignmentsPerMillion == 0:
    return 0.0
  float(counts(c)) / alignmentsPerMillion / kb

proc sampleName(path: string): string =
  extractFilename(path)

proc loadTargetTable(targetPath: string, targetOpts: TargetOptions): TableRef[string, seq[region_t]] =
  if targetOpts.parseAsGff:
    result = gff_to_table(targetPath, targetOpts.gffField, targetOpts.gffSeparator, targetOpts.gffIdentifier)
  else:
    result = bed_to_table(targetPath)

proc freezeTarget(targetTable: TableRef[string, seq[region_t]]): frozen_target_t =
  ## Convert parser-owned ref regions into plain immutable values once in the
  ## main thread so worker threads only need BAM-header-specific target cooking.
  for chrom, intervals in targetTable:
    var frozenIntervals = newSeq[frozen_region_t](intervals.len)
    for i, interval in intervals:
      frozenIntervals[i] = frozen_region_t(
        chrom: interval.chrom,
        start: interval.start,
        stop: interval.stop,
        name: interval.name
      )
    result.add((chrom: chrom, intervals: frozenIntervals))

proc collectFeatureMetadata(targetPath: string, targetTable: TableRef[string, seq[region_t]], targetOpts: TargetOptions): tuple[order: seq[string], coords: Table[string, feature_coords_t]] =
  var seen = initTable[string, bool]()
  # Preserve input order so multi-sample output remains deterministic.
  for name in target_names_in_order(
      targetPath,
      formatGff = targetOpts.parseAsGff,
      formatGtf = false,
      gffField = targetOpts.gffField,
      gffSeparator = targetOpts.gffSeparator,
      gffIdentifier = targetOpts.gffIdentifier):
    if name notin seen:
      result.order.add(name)
      seen[name] = true

  result.coords = initTable[string, feature_coords_t]()
  for chrom, intervals in targetTable:
    for interval in intervals:
      let featureName = if interval.name.len > 0: interval.name else: interval.chrom & ":" & $interval.start & "-" & $interval.stop
      let coords: feature_coords_t = (
        chrom: chrom,
        starts: $interval.start,
        stops: $interval.stop,
        name: featureName,
        length: interval.stop - interval.start
      )
      if featureName notin result.coords:
        result.coords[featureName] = coords
      else:
        result.coords[featureName].add(coords)

proc cookFrozenTarget(orig: frozen_target_t, bam: Bam): target_t =
  var chromMap = newTable[string, chrom_t]()
  for t in bam.hdr.targets:
    chromMap[t.name] = t.tid

  var cooked = newTable[chrom_t, seq[interval_t[string]]]()
  for entry in orig:
    doAssert(not (":" in entry.chrom), "bad target")
    if entry.chrom notin chromMap:
      raise newException(ValueError, "Target contig not found in BAM/CRAM header: " & entry.chrom)

    let chrom = chromMap[entry.chrom]
    cooked[chrom] = @[]
    var lastStart = 0
    for interval in entry.intervals:
      doAssert(interval.chrom == entry.chrom, "bad target")
      doAssert(interval.start >= lastStart)
      let name =
        if interval.name.len == 0:
          interval.chrom & ":" & $interval.start & "-" & $interval.stop
        else:
          interval.name
      cooked[chrom].add((pos_t(interval.start), pos_t(interval.stop), name))
      lastStart = interval.start

  cooked

proc processSample(bamPath: string, frozenTarget: frozen_target_t, opts: ProcessingOptions): SampleResult =
  result.samplePath = bamPath

  var
    bam: Bam
    fasta: cstring = nil

  if opts.fasta.len > 0:
    fasta = cstring(opts.fasta)

  if not open(bam, cstring(bamPath), threads=opts.threads, fai=fasta):
    result.error = "Unable to open BAM file: " & bamPath
    return

  if bam.hdr.isNil:
    result.error = "Invalid or empty BAM/CRAM file: " & bamPath
    return

  let cookedTarget = try:
      cookFrozenTarget(frozenTarget, bam)
    except CatchableError as e:
      result.error = e.msg
      return

  for _, intervals in cookedTarget:
    for interval in intervals:
      result.counts[interval.label] = (fwd: 0, rev: 0)

  result.perMillion = result.counts.alignments_count(
    bam,
    opts.mapq,
    opts.eflag,
    cookedTarget,
    opts.strict,
    opts.paired
  )
  
proc main(argv: var seq[string]): int =
  let env_fasta = getEnv("REF_PATH")
  let doc = format("""
  BamToCounts $version

  Usage: bamtocounts [options] <Target> <BAM-or-CRAM>...

Arguments:                                                                                                                                                 

  <Target>       the BED (or GFF) file containing regions in which to count reads
  <BAM-or-CRAM>  the alignment file for which to calculate depth

Options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -j, --jobs <jobs>            BAM files to process concurrently; 0 uses available CPUs [default: 0]
  -r, --fasta <fasta>          FASTA file for use with CRAM files [default: $env_fasta]
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: $default_flags]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 1]
  --paired                     Count read pairs rather than single reads
  --strict                     Read must be contained, not just overlap, with feature
  --stranded                   Print strand-specific counts
  --coords                     Also print coordinates of each feature

  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]

  -n, --rpkm                   Add a RPKM column
  -l, --norm-len               Add a counts/length column (after RPKM when both used)
  -p, --precision INT          Digits for floating point precision [default: 3]
  --header                     Print header
  --verbose                    Print verbose output
  --debug                      Enable diagnostics
  -h, --help                   Show help
  """ % ["version", version, "env_fasta", env_fasta, "default_flags", $DEFAULT_EXCLUDE_FLAGS])

  let
    args = docopt(doc, version=version, argv=argv)
    digitsPrecision = parseInt($args["--precision"])

  # Parse command-line arguments into structured options
  var opts = ProcessingOptions(
    mapq: uint8(parse_int($args["--mapq"])),
    eflag: uint16(parse_int($args["--flag"])),
    threads: parse_int($args["--threads"]),
    jobs: parse_int($args["--jobs"]),
    fasta: if $args["--fasta"] != "nil": $args["--fasta"] else: "",
    strict: bool(args["--strict"]),
    paired: bool(args["--paired"]),
    stranded: bool(args["--stranded"]),
    do_rpkm: bool(args["--rpkm"]),
    do_norm: bool(args["--norm-len"]),
    do_coords: bool(args["--coords"]),
    debug: bool(args["--debug"]),
    verbose: bool(args["--verbose"])
  )

  if opts.debug:
    stderr.writeLine("args:", args)

  var
    targetCoords = Table[string, feature_coords_t]()
    featureOrder: seq[string] = @[]
    sampleResults: seq[SampleResult] = @[]
    bamPaths = @(args["<BAM-or-CRAM>"])
    targetOpts = TargetOptions(
      parseAsGff: bool(args["--gff"]),
      gffIdentifier: $args["--id"],
      gffSeparator: ";",
      gffField: $args["--type"]
    )

  if not fileExists($args["<Target>"]):
    echo "ERROR: Target file does not exist: ", $args["<Target>"]
    quit(1)

  for bamPath in bamPaths:
    if not fileExists(bamPath):
      stderr.writeLine("ERROR: BAM/CRAM file does not exist: ", bamPath)
      quit(1)


  if ($args["<Target>"]).contains("gff") or ($args["<Target>"]).contains(".gtf"):
    if opts.debug:
      stderr.writeLine("Setting GFF/GTF format for: ", $args["<Target>"])
    targetOpts.parseAsGff = true

  let targetTable =
    try:
      loadTargetTable($args["<Target>"], targetOpts)
    except CatchableError as e:
      stderr.writeLine("ERROR: Unable to parse target file: ", $args["<Target>"], ": ", e.msg)
      quit(1)

  let frozenTarget = freezeTarget(targetTable)

  if opts.debug:
    stderr.writeLine("[OK] Target table loaded")

  if len(targetTable) == 0:
    stderr.writeLine("ERROR: No target regions found (try changing --id and --type): see an example line below")
    for line in lines($args["<Target>"]):
      if line.startsWith("#"):
        continue
      stderr.writeLine(line)
      quit(1)
    quit(1)
  if opts.debug or opts.verbose:
    stderr.writeLine("[OK] Target loaded: ", len(targetTable), " reference sequences")

  (featureOrder, targetCoords) = collectFeatureMetadata($args["<Target>"], targetTable, targetOpts)

  if args["--header"] or len(bamPaths) > 1:
    # Multi-sample output is ambiguous without column labels, so always emit a
    # header when more than one BAM/CRAM is requested.
    let coords = if opts.do_coords: "\tChrom\tStart\tEnd"
                 else: ""
    var header = "#Feature" & coords
    if len(bamPaths) == 1:
      if opts.stranded:
        header &= "\tFor\tRev"
      else:
        header &= "\tCounts"
      if opts.do_rpkm:
        header &= "\tRPKM"
      if opts.do_norm:
        header &= "\tCounts/Length"
    else:
      for bamPath in bamPaths:
        let sample = sampleName(bamPath)
        if opts.stranded:
          header &= "\t" & sample & "_For\t" & sample & "_Rev"
        else:
          header &= "\t" & sample
        if opts.do_rpkm:
          header &= "\t" & sample & "_RPKM"
        if opts.do_norm:
          header &= "\t" & sample & "_Counts/Length"
    echo header

  if opts.debug:
    stderr.writeLine("\\/ Target regions: ", len(featureOrder))

  when defined(threads):
    let requestedJobs =
      if opts.jobs > 0: opts.jobs
      else: min(len(bamPaths), countProcessors())
    let poolSize = max(1, min(len(bamPaths), requestedJobs))
    if poolSize == 1:
      for bamPath in bamPaths:
        sampleResults.add(processSample(bamPath, frozenTarget, opts))
    else:
      # BAMs are independent, so per-file workers are the least risky first
      # layer of parallelism and avoid shared counting state. The target has
      # already been frozen into immutable values, so workers only do BAM-
      # specific header mapping instead of reparsing the same file.
      var tp = Taskpool.new(poolSize)
      defer: tp.shutdown()
      var jobs = newSeq[FlowVar[SampleResult]]()
      for bamPath in bamPaths:
        if opts.debug:
          stderr.writeLine("[threads] Spawning sample worker for ", bamPath)
        jobs.add(tp.spawn processSample(bamPath, frozenTarget, opts))
      for job in jobs:
        sampleResults.add(sync(job))
  else:
    if len(bamPaths) > 1 and (opts.debug or opts.verbose):
      stderr.writeLine("[threads] Build has no thread support, processing BAMs sequentially")
    for bamPath in bamPaths:
      sampleResults.add(processSample(bamPath, frozenTarget, opts))

  for sample in sampleResults:
    if sample.error.len > 0:
      stderr.writeLine("ERROR: ", sample.error)
      quit(1)

  for feature in featureOrder:
    let
      coords = if opts.do_coords: targetCoords[feature].chrom & "\t" & targetCoords[feature].starts & "\t" & targetCoords[feature].stops & "\t"
                 else: ""
    var row = feature & "\t" & coords
    for sample in sampleResults:
      let rawcounts = sample.counts.getOrDefault(feature, (fwd: 0, rev: 0))
      row &= countsToString(rawcounts, opts.stranded)
      if opts.do_rpkm:
        row &= "\t" & rpkm(targetCoords[feature], sample.perMillion, rawcounts).formatBiggestFloat(ffDecimal, digitsPrecision)
      if opts.do_norm:
        row &= "\t" & (counts(rawcounts) / targetCoords[feature].length).formatBiggestFloat(ffDecimal, digitsPrecision)
      row &= "\t"
    echo row[0 .. ^2]


 
when isMainModule:
  var args = commandLineParams()
  try:
    discard main(args)
  except EKeyboardInterrupt:
    stderr.writeLine( "Quitting.")
  except:
    stderr.writeLine( getCurrentExceptionMsg() )
    quit(1)   
