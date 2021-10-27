import os
import hts
import docopt
import lapper
import strutils
import tables
import algorithm
import ./covutils

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
  alignmentsPerMillion : float = 0

type
  region_t = ref object
    chrom: string
    start: int
    stop: int
    name: string
    count: int

type
  referenceCounts = tuple[refName: string, order: int, length: int, counts: int, value: float]

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
    
    tableCounts[i.name] = @[]
    tableValues[i.name] = @[]
    
  result /= 1000000

proc count_alignments_per_ref(bam:Bam, mapq:uint8, eflag:uint16, factor: float): seq[referenceCounts] =
  for chromosome in bam.hdr.targets:
    var 
      chromCounts : referenceCounts
      rawCounts = 0
    chromCounts.refName = chromosome.name
    chromCounts.order   = chromosome.tid
    chromCounts.length   = int(chromosome.length)

    for aln in bam.query(chromosome.name):
      if aln.mapping_quality < mapq: continue
      if (aln.flag and eflag) != 0: continue
      rawCounts += 1
    chromCounts.counts = rawCounts
    chromCounts.value  = float(rawCounts) / factor
 
    if chromosome.name in tableCounts:
      tableCounts[chromosome.name].add( rawCounts )
      tableValues[chromosome.name].add( float(rawCounts) / factor )
    result.add(chromCounts)



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
      open(bam, bamFile, threads=threads, index=true, fai=fasta)
      if debug:
        stderr.writeLine("Opening BAM/CRAM file: ", bamFile)
    except:
      stderr.writeLine("Unable to open BAM file: ", bamFile )
         

    if bam.idx == nil:
      stderr.write_line("ERROR: requires BAM/CRAM index")
      quit(1)


    if alignmentsPerMillion <= 0:
      alignmentsPerMillion = bam.get_alignments_per_million()

    let sampleCounts = count_alignments_per_ref(bam, uint8(mapq), eflag, alignmentsPerMillion)
    
  
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