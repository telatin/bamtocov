import os
import hts
import docopt
import lapper
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
  gffIdentifier = "ID"
  gffSeparator  = ";"
  gffField      = "CDS"
  alignmentsPerMillion : float = 0




type
  region_t = ref object
    chrom: string
    start: int
    stop: int
    name: string
    count: int

proc inc_count(r:region_t) = inc(r.count)
proc start(r: region_t): int {.inline.} = return r.start
proc stop(r: region_t): int {.inline.} = return r.stop


proc tostring(r: region_t, s:var string) {.inline.} =
  # Print a 'region' to string (BED)
  s.set_len(0)
  s.add(r.chrom & "\t" & $r.start & "\t" & $r.stop & "\t")
  if r.name != "":
    s.add(r.name & "\t")
  s.add($r.count)

# Expand 'toString' with normalized counts (no longer BED)
proc toline(r: region_t, s: var string) {.inline.} =
  r.tostring(s)

  if do_rpkm:  
    let kb : float =(r.stop - r.start ) / 1000 
    let RPKM : float = float(r.count) / alignmentsPerMillion / kb
    s.add("\t" & $RPKM)

  if do_norm:
    let normLen = r.count / (r.stop - r.start)
    s.add("\t" & $normLen)


# Converts a GFF line to region object
proc gff_line_to_region(line: string): region_t =
  var
   cse = line.strip().split('\t')

  if len(cse) < 5:
    stderr.write_line("[warning] skipping GFF line (fields not found):", line.strip())
    return nil

  # Skip non CDS fields (or user provided)
  if cse[2] != gffField:
    return nil

  var
    s = parse_int(cse[3])  - 1
    e = parse_int(cse[4])
    reg = region_t(chrom: cse[0], start: s, stop: e, count:0)
  
  # In the future, 8th field could be requireed [TODO]
  if len(cse) == 9:
    for gffAnnotPart in cse[8].split(gffSeparator):
      if gffAnnotPart.startsWith(gffIdentifier):
        reg.name = gffAnnotPart.split("=")[1] 
        break
  return reg

# Convert a BED line to region object
proc bed_line_to_region(line: string): region_t =
  var
   cse = line.strip().split('\t', 5)

  if len(cse) == 9:
    stderr.writeLine("[warning] GFF format detected.")
    return gff_line_to_region(line)
    

  if len(cse) < 3:
    stderr.write_line("[warning] skipping bad bed line:", line.strip())
    return nil
  var
    s = parse_int(cse[1])
    e = parse_int(cse[2])
    reg = region_t(chrom: cse[0], start: s, stop: e, count:0)
  if len(cse) > 3:
   reg.name = cse[3]
  return reg

# Convert BED file to table
proc bed_to_table(bed: string): TableRef[string, seq[region_t]] =
  var bed_regions = newTable[string, seq[region_t]]()
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr: hts.kstring_t
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("track "):
      continue
    if $kstr.s[0] == "#":
      continue
    var v = bed_line_to_region($kstr.s)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)
  
  for chrom, ivs in bed_regions.mpairs:     # since it is read into mem, can also well sort. (BP)
    sort(ivs, proc (a, b: region_t): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions



proc gff_to_table(bed: string): TableRef[string, seq[region_t]] =
  var bed_regions = newTable[string, seq[region_t]]()
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr: hts.kstring_t
  kstr.l = 0
  kstr.m = 0
  kstr.s = nil
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if ($kstr.s).startswith("##FASTA"):
      break
    if $kstr.s[0] == "#":
      continue

    var v = gff_line_to_region($kstr.s)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, can also well sort.
  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region_t): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions

proc get_alignments_per_million(bam:Bam): float =
  for i in bam.hdr.targets:
    result += float(stats(bam.idx,i.tid).mapped)
  result /= 1000000

proc print_alignments_count(bam:Bam, mapq:uint8, eflag:uint16, regions:TableRef[string, seq[region_t]]) =
  for chrom in regions.keys():
    if not regions.contains(chrom) or regions[chrom].len == 0:
      continue
    var lap: Lapper[region_t] = lapify(regions[chrom])

    for aln in bam.query(chrom):
      if aln.mapping_quality < mapq: continue
      if (aln.flag and eflag) != 0: continue

      lap.each_seek(aln.start.int, aln.stop.int, inc_count)
    var s = new_string_of_cap(1000)         # Returns a new string of length 0 but with capacity cap.
    for region in regions[chrom]:
      region.toline(s)
      echo s
 
#[
  proc each_seek[T: Interval](
    L: var Lapper[T]; 
    start: int; 
    stop: int; f
    n: proc (v: T)) {..}
call fn(x) for each interval x in L that overlaps start..stop this assumes that subsequent calls to this function will be in sorted order
]#

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
  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]
  -n, --rpkm                   Add a RPKM column
  -l, --norm-len               Add a counts/length column (after RPKM when both used)
  --header                     Print header
  --debug                      Enable diagnostics    
  -h, --help                   Show help
  """ % ["version", version, "env_fasta", env_fasta])

  let args = docopt(doc, version=version, argv=argv)
  let mapq = parse_int($args["--mapq"])
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

  if len($args["<BAM-or-CRAM>"]) > 1:
    echo "Multiple BAM/CRAM files not supported in the current version."
    quit(1)

  try:
    open(bam, $args["<BAM-or-CRAM>"], threads=threads, index=true, fai=fasta)
    if debug:
      stderr.writeLine("Opening BAM/CRAM file: ", $args["<BAM-or-CRAM>"])
  except:
    stderr.writeLine("Unable to open BAM file: ", $args["<BAM-or-CRAM>"] )

  if do_rpkm:
    alignmentsPerMillion = bam.get_alignments_per_million()
    if debug:
      stderr.writeLine("Total: ", 1000000 * bam.get_alignments_per_million())

  if bam.idx == nil:
    stderr.write_line("ERROR: requires BAM/CRAM index")
    quit(1)

  if args["--header"]:
    if do_rpkm and do_norm:
      echo "#Chrom\tstart\tend\tcounts\tRPKM\tCounts/Length"
    elif do_rpkm:
      echo "#Chrom\tstart\tend\tcounts\tRPKM"
    elif do_norm:
      echo "#Chrom\tstart\tend\tcounts\tCounts/Length"
    else:
      echo "#Chrom\tstart\tend\tcounts"

  if ($args["<Target>"]).endsWith("gff"):
    prokkaGff = true

  var regions = if prokkaGff == true: gff_to_table($args["<Target>"])
                 else: bed_to_table($args["<Target>"])
  
  if debug:
    stderr.writeLine("Target loaded: ", len(regions), " reference sequences")

  
  print_alignments_count(bam, uint8(mapq), eflag, regions)
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