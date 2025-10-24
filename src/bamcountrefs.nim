# Standard library
import std/[os, strutils, tables]

# External dependencies
import docopt, hts

const NimblePkgVersion {.strdefine.} = "prerelease"

let
  version = NimblePkgVersion

var
  tableCounts = initTable[string, seq[int]]()
  tableValues = initTable[string, seq[float]]()

type EKeyboardInterrupt = object of CatchableError
 
proc handler() {.noconv.} =
  raise newException(EKeyboardInterrupt, "Keyboard Interrupt")
 
setControlCHook(handler)


var
  do_norm = false
  debug = false
  do_rpkm = false
  gffIdentifier = "ID"
  gffSeparator  = ";"
  gffField      = "CDS"

type
  referenceCounts = tuple[refName: string, order: int, length: int, counts: int, value: float]

proc get_alignments_per_million(bam:Bam): float =
  # Calculate total mapped reads and normalize to "per million"
  for i in bam.hdr.targets:
    result += float(stats(bam.idx,i.tid).mapped)
  result /= 1000000

proc count_alignments_per_ref(bam:Bam, mapq:uint8, eflag:uint16, factor: float): seq[referenceCounts] =
  # Performance optimization: check if we can use BAM index statistics
  # Index stats give us mapped read counts (reads without flag 4 = unmapped)
  # We can only use them when no additional filtering is required:
  #   - mapq == 0: no mapping quality filter
  #   - eflag == 4: only exclude unmapped reads (already excluded by index)
  # For default eflag=1796 (excludes secondary, qc-fail, duplicate), we must iterate
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
    chromCounts.value  = float(rawCounts) / factor

    # Initialize table entries on first access, then append counts for this sample
    discard tableCounts.hasKeyOrPut(chromosome.name, newSeq[int]())
    discard tableValues.hasKeyOrPut(chromosome.name, newSeq[float]())
    tableCounts[chromosome.name].add( rawCounts )
    tableValues[chromosome.name].add( float(rawCounts) / factor )

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
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]

Annotation options:
  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]
  -n, --rpkm                   Add a RPKM column
  -l, --norm-len               Add a counts/length column (after RPKM when both used)

Other options;
  --tag STR                    First column name [default: ViralSequence]
  --multiqc                    Print output as MultiQC table
  --header                     Print header
  --debug                      Enable diagnostics    
  -h, --help                   Show help
  """ % ["version", version, "env_fasta", env_fasta])

  let args = docopt(doc, version=version, argv=argv)
  let
    mapq = parse_int($args["--mapq"])
    columnName = $args["--tag"]
  var prokkaGff : bool = args["--gff"]
  do_rpkm = args["--rpkm"]
  do_norm = args["--norm-len"]
  debug = args["--debug"]
  gffIdentifier = $args["--id"]
  gffField      = $args["--type"]

  var fasta: cstring 
  if $args["--fasta"] != "nil":
    fasta = cstring($args["--fasta"])

  var
    eflag = uint16(parse_int($args["--flag"]))
    threads = parse_int($args["--threads"])
    bam:Bam

  var
    samples = @[columnName]

  for bamFile in @(args["<BAM-or-CRAM>"]):
    var sampleName = extractFilename(bamFile)
    samples.add(sampleName.split('.')[0])
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
      stderr.writeLine("[debug] Sample: ", sampleName.split('.')[0],
                       " - alignments per million: ", formatFloat(currentAlignmentsPerMillion, ffDecimal, 3))

    let sampleCounts = count_alignments_per_ref(bam, uint8(mapq), eflag, currentAlignmentsPerMillion)
    
  
  if args["--multiqc"]:
    echo "# plot_type: 'table'"
    echo "# section_name: 'CovTools count'"
    echo "# description: 'Feature table: counts of mapped reads against predicted viral sequences'"

  echo samples.join("\t")
  for reference in tableCounts.keys:
    if do_rpkm:
      echo reference, "\t", tableValues[reference].join("\t")
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