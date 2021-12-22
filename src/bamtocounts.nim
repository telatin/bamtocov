import os
import hts
import docopt
 
import strutils
import tables
import algorithm

import ./covutils
 
#[
  **bamToCounts**, part of MAGENTA Flow
  based on count-reads in the "hts-nim-tools" suite by Brent Pedersen
  see: "https://github.com/brentp/hts-nim-tools"
  Static binary thanks to  "https://github.com/brentp/hts-nim"

 
  0.4.1   Adding autoguess of GFF also from column count [experimental]

  0.4.0   Fix: start coordinate of GFF; 
          fix: buried total alignemnts print only when debug is ON
          Auto --gff if file ends with gff.
          Added --header; 
  0.3.2   Code refactoring
  0.3.1   Fix warnings while parsing GFF
  0.3.0   Added RPKM calculastion (requires: sum of total alignments)
          Added normalizaation by gene length
  0.2.0   Added support for GFF files
  0.1.0   Added debug mode
]#


type EKeyboardInterrupt = object of CatchableError
 
proc handler() {.noconv.} =
  raise newException(EKeyboardInterrupt, "Keyboard Interrupt")
 
setControlCHook(handler)

var
  do_norm = false
  debug = false
  do_rpkm = false
  do_strand = false
  do_coords = false
  gffIdentifier = "ID"
  gffSeparator  = ";"
  gffField      = "CDS"
  alignmentsPerMillion : float = 0



 
  

# Expand 'toString' with normalized counts (no longer BED)

proc get_alignments_per_million(bam:Bam): float =
  for i in bam.hdr.targets:
    result += float(stats(bam.idx,i.tid).mapped)
  result /= 1000000
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

proc alignments_count(table: var OrderedTable[string, stranded_counts], bam:Bam, mapq:uint8, eflag:uint16, regions: target_t) =
 
  for chrom in regions.keys():
    if debug:
      stderr.writeLine("[alignments_count] Got chrom: ", chrom)
    for aln in bam.query(chrom):
     
      if not regions.contains(chrom) or regions[chrom].len == 0:
        continue

      var
        target_idx: target_index_t

      for aln in bam.query(chrom):
        if aln.mapping_quality < mapq: continue
        if (aln.flag and eflag) != 0: continue
        
        let readAsInterval = (aln.tid, pos_t(aln.start), pos_t(aln.stop), aln.flag.reverse)
        
        for interval in intersections(readAsInterval, regions, target_idx):
          # Returns: genomic_interval_t[tuple[l1: T, l2: string]
          
          #let feature : target_feature = (chrom: "", cid: interval.chrom.int, start: interval.start.int, stop: interval.stop.int, feature: interval.label.l2)
          if interval.label.l2 notin table:
            stderr.writeLine("[alignments_count] Warning: unknown feature: ", interval.label.l2)
            table[interval.label.l2] = (fwd: 0, rev: 0)
          table[interval.label.l2].inc(aln.flag.reverse)
  

proc targetSort(x, y: target_feature): int =
  if x.cid < y.cid:
    -1
  elif x.cid > y.cid:
    1
  else:
    if x.start < y.start:
      -1
    elif x.start > y.start:
      1
    else:
      if x.stop < y.stop:
        -1
      elif x.stop > y.stop:
        1
      else:
        0

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
  -r, --fasta <fasta>          FASTA file for use with CRAM files [default: $env_fasta].
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]
  --stranded                   Print strand-specific counts
  --coords                     Also print coordinates of each feature

  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]

  -n, --rpkm                   Add a RPKM column
  -l, --norm-len               Add a counts/length column (after RPKM when both used)
  -p, --precision INT          Digits for floating point precision [default: 3]
  --header                     Print header
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
  debug = bool(args["--debug"])
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
    open(bam, cstring($args["<BAM-or-CRAM>"]), threads=threads,index=true,  fai=fasta)
    if debug:
      stderr.writeLine("Opening BAM/CRAM file: ", $args["<BAM-or-CRAM>"])
  except:
    stderr.writeLine("Unable to open BAM file: ", $args["<BAM-or-CRAM>"] )
    quit(1)

  if do_rpkm:
    alignmentsPerMillion = bam.get_alignments_per_million()
    if debug:
      stderr.writeLine("Total: ", 1000000 * bam.get_alignments_per_million())

 
  if bam.idx == nil:
    stderr.write_line("ERROR: requires BAM/CRAM index")
    quit(1) 

  if ($args["<Target>"]).contains("gff") or ($args["<Target>"]).contains(".gtf"):
    prokkaGff = true

  var targetTable = if prokkaGff == true: gff_to_table($args["<Target>"], gffField, gffSeparator, gffIdentifier)
                 else: bed_to_table($args["<Target>"])



  if len(targetTable) == 0:
    stderr.writeLine("ERROR: No target regions found (try changing --id and --type): see an example line below")
    for line in lines($args["<Target>"]):
      if line.startsWith("#"):
        continue
      stderr.writeLine(line)
      quit(1)
    quit(1)
  if debug:
    stderr.writeLine("Target loaded: ", len(targetTable), " reference sequences")
    
  let cookedTarget = cookTarget(targetTable, bam)
  #let countsTable  = alignments_count(bam, uint8(mapq), eflag, cookedTarget)
  var targetCounts = OrderedTable[string, stranded_counts]()
  for index, chrName in bam.hdr.targets:
    #feature_coords  = tuple[chrom, starts, stops, name]
    if debug:
      stderr.writeLine("BAM targets: ", chrName, "-", index)
    if index in cookedTarget:
      if debug:
        stderr.writeLine("Coocked targets: ", chrName, "-", index)
      for interval in cookedTarget[index]:
        let
          c : feature_coords = (chrom: chrName.name, starts: $interval.start, stops: $interval.stop, name: interval.label, length: int(interval.stop - interval.start))
        if debug:
          stderr.writeLine("\tInterval: ",interval.label, "-", interval.start, "-", interval.stop)
        if interval.label notin targetCoords:
          if debug:
            stderr.writeLine("\t - Adding")
          targetCoords[interval.label] = c
          targetCounts[interval.label] = (fwd: 0, rev: 0)
          
        else:
          if debug:
            stderr.writeLine("\t - Extending")
          targetCoords[interval.label].add(c)
        
        
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
    stderr.writeLine("Target regions: ", len(targetCounts))
  targetCounts.alignments_count(bam, uint8(mapq), eflag, cookedTarget)  
  if debug:
    stderr.writeLine("Counts done") 
  
   
  for feature, rawcounts in targetCounts:
    let
      coords = if do_coords: targetCoords[feature].chrom & "\t" & targetCoords[feature].starts & "\t" & targetCoords[feature].stops & "\t"
                 else: ""
      counts = countsToString(rawcounts, do_strand)
      rpkm   = if do_rpkm: "\t" & rpkm(targetCoords[feature], alignmentsPerMillion, rawcounts).formatBiggestFloat(ffDecimal, digitsPrecision)
               else: "" 
      norm   = if do_norm:  "\t" & ( counts(rawcounts) / targetCoords[feature].length ).formatBiggestFloat(ffDecimal, digitsPrecision)
               else: ""
    echo feature, "\t", coords, counts, rpkm, norm

#[ 
  # Print table
  var tableKeys = newSeq[target_feature]()
  for j in keys(countsTable):
    tableKeys.add(j)
  #tableKeys.sort(targetSort)
  for feat in tableKeys: 
    var fields = ""
    let featCoords = if do_coords: #[feat.cid] & "\t" &  $feat.start & "\t" &  $feat.stop & "\t"
                     else: ""
    if do_rpkm#:
      let 
        kb   : float = (feat.stop - feat.start ) / 1000 
        RPKM : float = float(counts(countsTable[feat])) / alignmentsPerMillion / kb
      fields &= "\t" & $RPKM.formatBiggestFloat(ffDecimal, digitsPrecision)
    if do_norm:
      let normal = counts(countsTable[feat]) / (feat.stop - feat.start)
      fields &= "\t" & $normal.formatBiggestFloat(ffDecimal, digitsPrecision) 
    echo  featCoords,  "\t", feat.feature, "\t", counts(countsTable[feat]), fields ]#
  return 0
]#
 
when isMainModule:
  var args = commandLineParams()
  try:
    discard main(args)
  except EKeyboardInterrupt:
    stderr.writeLine( "Quitting.")
  except:
    stderr.writeLine( getCurrentExceptionMsg() )
    quit(1)   