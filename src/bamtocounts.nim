# Standard library
import std/[os, strutils, tables]

# External dependencies
import docopt, hts

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
    fasta: string
    strict: bool        # Require full containment within feature
    paired: bool        # Count pairs instead of individual reads
    stranded: bool      # Track forward/reverse strand counts separately
    do_rpkm: bool       # Calculate RPKM
    do_norm: bool       # Calculate counts/length
    do_coords: bool     # Output feature coordinates
    debug: bool
    verbose: bool

proc handler() {.noconv.} =
  raise newException(EKeyboardInterrupt, "Keyboard Interrupt")

setControlCHook(handler)

# GFF parsing parameters (could be moved to a config object if needed)
var
  gffIdentifier = "ID"
  gffSeparator  = ";"
  gffField      = "CDS"

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

  # GFF parsing options
  var prokkaGff: bool = args["--gff"]
  gffIdentifier = $args["--id"]
  gffField = $args["--type"]

  var
    targetCoords = Table[string, feature_coords_t]()
    bam: Bam
    fasta: cstring = nil

  if opts.fasta.len > 0:
    fasta = cstring(opts.fasta)

  if len(args["<BAM-or-CRAM>"]) > 1:
    echo "Multiple BAM/CRAM files not supported in the current version."
    quit(1)

  if not fileExists($args["<Target>"]):
    echo "ERROR: Target file does not exist: ", $args["<Target>"]
    quit(1)

  try:
    open(bam, cstring($args["<BAM-or-CRAM>"]), threads=opts.threads, fai=fasta)
    if opts.debug:
      stderr.writeLine("Opening BAM/CRAM file: ", $args["<BAM-or-CRAM>"])
  except:
    stderr.writeLine("Unable to open BAM file: ", $args["<BAM-or-CRAM>"])
    quit(1)


  if ($args["<Target>"]).contains("gff") or ($args["<Target>"]).contains(".gtf"):
    if opts.debug:
      stderr.writeLine("Setting GFF/GTF format for: ", $args["<Target>"])
    prokkaGff = true

  var targetTable: TableRef[string, seq[region_t]]

  if prokkaGff == true:
    try:
      targetTable = gff_to_table($args["<Target>"], gffField, gffSeparator, gffIdentifier)
    except Exception as e:
      stderr.writeLine("ERROR: Unable to parse GFF file: ", $args["<Target>"], ": ", e.msg)
      quit(1)
  else: 
    try:
      targetTable = bed_to_table($args["<Target>"])
    except Exception as e:
      stderr.writeLine("ERROR: Unable to parse BED file: ", $args["<Target>"], ": ", e.msg)
      quit(1)  


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

  let cookedTarget = cookTarget(targetTable, bam)
  var targetCounts = OrderedTable[string, stranded_counts_t]()
  for index, chrName in bam.hdr.targets:
    if opts.debug:
      stderr.writeLine(" > BAM targets: ", chrName, " - index:", index)
    if index in cookedTarget:
      if opts.debug:
        stderr.writeLine(" + Coocked targets: ", chrName, " - index:", index)
      for interval in cookedTarget[index]:
        let
          c: feature_coords_t = (chrom: chrName.name, starts: $interval.start, stops: $interval.stop, name: interval.label, length: int(interval.stop - interval.start))
        if opts.debug:
          stderr.writeLine("    > Interval: ", interval.label, "-", interval.start, "-", interval.stop)
        if interval.label notin targetCoords:
          if opts.debug:
            stderr.writeLine("      - Adding")
          targetCoords[interval.label] = c
          targetCounts[interval.label] = (fwd: 0, rev: 0)
        else:
          if opts.debug:
            stderr.writeLine("      - Extending")
          targetCoords[interval.label].add(c)
    else:
      if opts.debug:
        stderr.writeLine("No coocked targets: ", chrName, "-", index)
        
        
  if args["--header"]:
    let coords = if opts.do_coords: "Chrom\tStart\tEnd\t"
                 else: ""
    let header = if opts.stranded: "#Feature\t" & coords & "For\tRev"
                else:   "#Feature\t" & coords & "Counts"
    if opts.do_rpkm and opts.do_norm:
      echo header & "\tRPKM\tCounts/Length"
    elif opts.do_rpkm:
      echo header & "\tRPKM"
    elif opts.do_norm:
      echo header & "\tCounts/Length"
    else:
      echo header

  if opts.debug:
    stderr.writeLine("\\/ Target regions: ", len(targetCounts))

  ## GATHER THE COUNTS
  let perMillion = targetCounts.alignments_count(bam, opts.mapq, opts.eflag, cookedTarget, opts.strict, opts.paired)
  if opts.debug:
    stderr.writeLine("/\\ Counts done: ", perMillion) 
  
   
  for feature, rawcounts in targetCounts:
    let
      coords = if opts.do_coords: targetCoords[feature].chrom & "\t" & targetCoords[feature].starts & "\t" & targetCoords[feature].stops & "\t"
                 else: ""
      counts = countsToString(rawcounts, opts.stranded)
      rpkm   = if opts.do_rpkm: "\t" & rpkm(targetCoords[feature], perMillion, rawcounts).formatBiggestFloat(ffDecimal, digitsPrecision)
               else: ""
      norm   = if opts.do_norm:  "\t" & ( counts(rawcounts) / targetCoords[feature].length ).formatBiggestFloat(ffDecimal, digitsPrecision)
               else: ""
    echo feature, "\t", coords, counts, rpkm, norm


 
when isMainModule:
  var args = commandLineParams()
  try:
    discard main(args)
  except EKeyboardInterrupt:
    stderr.writeLine( "Quitting.")
  except:
    stderr.writeLine( getCurrentExceptionMsg() )
    quit(1)   