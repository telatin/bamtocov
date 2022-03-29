import os
import hts
import docopt
 
import strutils
import tables


import ./covutils
 



type EKeyboardInterrupt = object of CatchableError
 
proc handler() {.noconv.} =
  raise newException(EKeyboardInterrupt, "Keyboard Interrupt")
 
setControlCHook(handler)

var
  do_strict = false
  do_paired = false
  do_norm = false
  debug = false
  verbose = false
  do_rpkm = false
  do_strand = false
  do_coords = false
  gffIdentifier = "ID"
  gffSeparator  = ";"
  gffField      = "CDS"
  alignmentsPerMillion : float = 0



 
  

# Expand 'toString' with normalized counts (no longer BED)


#[ 
proc legacytoline(r: region_t, s: var string) {.inline.} =
  r.renderString(s)
  let countFloat = float(1)


 
]#
type
  target_feature  = tuple[chrom: string, cid: int, start: int, stop: int, feature: string]
  stranded_counts = tuple[fwd, rev: int]
  feature_coords  = tuple[chrom, starts, stops, name: string, length: int]
  
 
proc inc(c: var stranded_counts, reverse=false) =
  if reverse == false:
    c.fwd += 1
  else:
    c.rev += 1

proc counts(c: stranded_counts): int =
  c.fwd + c.rev

proc countsToString(c: stranded_counts, stranded: bool): string =
  if stranded:
    $(c.fwd) & "\t" & $(c.rev)
  else:
    $(c.fwd + c.rev)

type
  counts_t*       = TableRef[interval_t[string], int]

proc makeCountsTable(table: var OrderedTable[string, stranded_counts], bam:Bam, mapq:uint8, eflag:uint16, regions: target_t, strict = false): float =
  var total: float = 0
  for read in bam:
    total += 1
    if read.tid notin regions:
      continue
    if read.mapping_quality == 0 or ( (read.flag and 1796) != 0):
      continue

    for region in regions[read.tid]:
      
      if (read.start < region.start   and read.stop > region.stop) or (read.stop > region.start   and read.stop < region.stop) or (read.start > region.start  and read.start < region.stop):
        if strict and ( read.start < region.start or read.stop > region.stop ):
            continue
        table[region.label].inc(read.flag.reverse)
  return total / 1000000


proc alignments_count(table: var OrderedTable[string, stranded_counts], bam:Bam, mapq:uint8, eflag:uint16, regions: target_t, strict = false, paired = false): float =
  var total: float = 0
  for read in bam:
    total += 1
    if read.tid notin regions:
      continue
    if read.mapping_quality == 0 or ( (read.flag and eflag) != 0):
      continue
    if paired and (read.flag.proper_pair == false or read.flag.read2 == true):
      continue


    for region in regions[read.tid]:
      
      if (read.start < region.start   and read.stop > region.stop) or (read.stop > region.start   and read.stop < region.stop) or (read.start > region.start  and read.start < region.stop):
        if strict and ( read.start < region.start or read.stop > region.stop ):
            continue
        table[region.label].inc(read.flag.reverse)
  return total / 1000000

proc add(s: var feature_coords, z: feature_coords) = 
  # feature_coords  = tuple[chrom, starts, stops, name: string]
  s.chrom &= ";" & z.chrom
  s.starts &= ";" & z.starts
  s.stops &= ";" & z.stops
  s.length += z.length
  
proc rpkm(f: feature_coords, alignmentsPerMillion: float, c: stranded_counts ): float =
  let
    kb   : float = f.length / 1000 
  
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
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]
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
  """ % ["version", version, "env_fasta", env_fasta])

  let
    args = docopt(doc, version=version, argv=argv)
    digitsPrecision = parseInt($args["--precision"])

  if args["--debug"]:
    stderr.write_line("args:", args)
  let mapq = parse_int($args["--mapq"])
  var prokkaGff : bool = args["--gff"]

  do_rpkm = bool(args["--rpkm"])
  do_norm = bool(args["--norm-len"])
  do_strand = bool(args["--stranded"])
  do_coords = bool(args["--coords"])
  do_strict = bool(args["--strict"])
  do_paired = bool(args["--paired"])
  debug = bool(args["--debug"])
  verbose = bool(args["--verbose"])
  gffIdentifier = $args["--id"]
  gffField      = $args["--type"]

  var fasta: cstring 
  if $args["--fasta"] != "nil":
    fasta = cstring($args["--fasta"])

  var
    eflag = uint16(parse_int($args["--flag"]))
    threads = parse_int($args["--threads"])
    #targetNames = Table[int, string]()
    targetCoords = Table[string, feature_coords]()
    #targetCounts = Table[string, stranded_counts]()
    bam:Bam

  if len(args["<BAM-or-CRAM>"]) > 1:
    echo "Multiple BAM/CRAM files not supported in the current version."
    quit(1)

  if not fileExists($args["<Target>"]):
    echo "ERROR: Target file does not exist: ", $args["<Target>"]
    quit(1)

  try:                                    #index=true,
    open(bam, cstring($args["<BAM-or-CRAM>"]), threads=threads,  fai=fasta)
    if debug:
      stderr.writeLine("Opening BAM/CRAM file: ", $args["<BAM-or-CRAM>"])
  except:
    stderr.writeLine("Unable to open BAM file: ", $args["<BAM-or-CRAM>"] )
    quit(1)


  if ($args["<Target>"]).contains("gff") or ($args["<Target>"]).contains(".gtf"):
    if debug:
      stderr.writeLine("Setting GFF/GTF format for: ", $args["<Target>"])
    prokkaGff = true

  var targetTable : TableRef[string, seq[region_t]]
  
 
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


  if debug:
    stderr.writeLine("[OK] Target table loaded")
    
    #stderr.writeLine("Target table: ", targetTable)

  if len(targetTable) == 0:
    stderr.writeLine("ERROR: No target regions found (try changing --id and --type): see an example line below")
    for line in lines($args["<Target>"]):
      if line.startsWith("#"):
        continue
      stderr.writeLine(line)
      quit(1)
    quit(1)
  if debug or verbose:
    stderr.writeLine("[OK] Target loaded: ", len(targetTable), " reference sequences")
    
  let cookedTarget = cookTarget(targetTable, bam)
  #let countsTable  = alignments_count(bam, uint8(mapq), eflag, cookedTarget)
  var targetCounts = OrderedTable[string, stranded_counts]()
  for index, chrName in bam.hdr.targets:
    #feature_coords  = tuple[chrom, starts, stops, name]
    if debug:
      stderr.writeLine(" > BAM targets: ", chrName, " - index:", index)
    if index in cookedTarget:
      if debug:
        stderr.writeLine(" + Coocked targets: ", chrName, " - index:", index)
      for interval in cookedTarget[index]:
        let
          c : feature_coords = (chrom: chrName.name, starts: $interval.start, stops: $interval.stop, name: interval.label, length: int(interval.stop - interval.start))
        if debug:
          stderr.writeLine("    > Interval: ",interval.label, "-", interval.start, "-", interval.stop)
        if interval.label notin targetCoords:
          if debug:
            stderr.writeLine("      - Adding")
          targetCoords[interval.label] = c
          targetCounts[interval.label] = (fwd: 0, rev: 0)
          
        else:
          if debug:
            stderr.writeLine("      - Extending")
          targetCoords[interval.label].add(c)
    else:
      if debug:
        stderr.writeLine("No coocked targets: ", chrName, "-", index)
        
        
  if args["--header"]:
    let coords = if do_coords: "Chrom\tStart\tEnd\t"
                 else: ""
    let header = if do_strand: "#Feature\t" & coords & "For\tRev"
                else:   "#Feature\t" & coords & "Counts"
    if do_rpkm and do_norm:
      echo header & "\tRPKM\tCounts/Length"
    elif do_rpkm:
      echo header & "\tRPKM"
    elif do_norm:
      echo header & "\tCounts/Length"
    else:
      echo header

  if debug:
    stderr.writeLine("\\/ Target regions: ", len(targetCounts))

  ## GATHER THE COUNTS
  let perMillion = targetCounts.alignments_count(bam, uint8(mapq), eflag, cookedTarget, do_strict, do_paired)  
  if debug:
    stderr.writeLine("/\\ Counts done: ", perMillion) 
  
   
  for feature, rawcounts in targetCounts:
    let
      coords = if do_coords: targetCoords[feature].chrom & "\t" & targetCoords[feature].starts & "\t" & targetCoords[feature].stops & "\t"
                 else: ""
      counts = countsToString(rawcounts, do_strand)
      rpkm   = if do_rpkm: "\t" & rpkm(targetCoords[feature], perMillion, rawcounts).formatBiggestFloat(ffDecimal, digitsPrecision)
               else: "" 
      norm   = if do_norm:  "\t" & ( counts(rawcounts) / targetCoords[feature].length ).formatBiggestFloat(ffDecimal, digitsPrecision)
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