import strutils
import tables
import hts
import algorithm
import posix

signal(SIG_PIPE, SIG_IGN)
const NimblePkgVersion {.strdefine.} = "prerelease"

let
  version* = NimblePkgVersion

type
  region_t* = ref object
    chrom: string
    start: int
    stop: int
    name: string
    count: int


type
  covopt* = ref object
    verbose, debug: bool
    outputFmt: string # bed, wig
    targetFmt: string # bed, gff
    gffSep: char
    gffId, gffType: string

proc inc_count*(r:region_t) = inc(r.count)
proc start*(r: region_t): int {.inline.} = return r.start
proc stop*(r: region_t): int {.inline.} = return r.stop
proc chrom*(r: region_t): string {.inline.} = return r.chrom
proc name*(r: region_t): string {.inline.} = return r.name


proc tostring*(r: region_t, s:var string) {.inline.} =
  # Print a 'region' to string (BED)
  s.set_len(0)
  s.add(r.chrom & "\t" & $r.start & "\t" & $r.stop & "\t")
  if r.name != "":
    s.add(r.name & "\t")
  s.add($r.count)

# Converts a GTF line to region object
proc gtf_line_to_region*(line: string, gffField = "exon", gffSeparator = ";", gffIdentifier = "gene_id"): region_t =
  #NC_001422.1     Prodigal:002006 gene    51      221     .       +       0       gene_id "nbis-gene-1"; ID "nbis-gene-1"; inference "ab initio prediction:Prodigal:002006"; locus_tag "PhiX_01"; product "hypothetical protein";
  var
   cse = line.strip().split('\t')

  if len(cse) < 8:
    stderr.write_line("[warning] skipping GTF line (fields not found):", line.strip())
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
    for gffAnnotPartRaw in cse[8].split(gffSeparator):
      let gffAnnotPart = gffAnnotPartRaw.strip(chars = {'"', '\'', ' '})
      if gffAnnotPart.startsWith(gffIdentifier):
        try:
          reg.name = gffAnnotPart.split("=")[1].strip(chars = {'"', '\'', ' '}) 
        except:
          reg.name = gffAnnotPart.split(" ")[1].strip(chars = {'"', '\'', ' '})
        break
  return reg

# Converts a GFF line to region object
proc gff_line_to_region*(line: string, gffField = "CDS", gffSeparator = ";", gffIdentifier = "ID"): region_t =
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
    for gffAnnotPartRaw in cse[8].split(gffSeparator):
      let gffAnnotPart = gffAnnotPartRaw.strip(chars = {'"', '\'', ' '})
      if gffAnnotPart.startsWith(gffIdentifier):
        try:
          reg.name = gffAnnotPart.split("=")[1].strip(chars = {'"', '\'', ' '}) 
        except:
          reg.name = gffAnnotPart.split(" ")[1].strip(chars = {'"', '\'', ' '})
        break
  return reg

# Convert a BED line to region object
proc bed_line_to_region*(line: string): region_t =
  var
   cse = line.strip().split('\t', 5)

  if len(cse) == 9:
    stderr.writeLine("[warning] GFF format detected, attempting detection with defaults.")
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
proc bed_to_table*(bed: string): TableRef[string, seq[region_t]] =
  var bed_regions = newTable[string, seq[region_t]]()
  if bed == "nil":
    return bed_regions

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


proc gtf_to_table*(bed: string, gffField, gffSeparator, gffIdentifier: string): TableRef[string, seq[region_t]] =
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

    var v = gtf_line_to_region($kstr.s, gffField, gffSeparator, gffIdentifier)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, can also well sort.
  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region_t): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions

proc gff_to_table*(bed: string, gffField, gffSeparator, gffIdentifier: string): TableRef[string, seq[region_t]] =
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

    var v = gff_line_to_region($kstr.s, gffField, gffSeparator, gffIdentifier)
    if v == nil: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, can also well sort.
  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region_t): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions