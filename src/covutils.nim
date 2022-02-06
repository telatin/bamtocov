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
#    counts_for: int
#    counts_rev: int
     


type
  chrom_t* = int # reference id
  pos_t* = int64
  # here intervals have a "label", which contains additional information besides the location
  # there is one interval without explicit chromosome to be used in the target table, where intervals are already grouped by chromosome
  interval_t*[T]  = tuple[start, stop: pos_t, label: T]
  genomic_interval_t*[T] = tuple[chrom: chrom_t, start, stop: pos_t, label: T]
  target_index_t* = tuple[chrom: chrom_t, interval: int]
  target_t*       = TableRef[chrom_t, seq[interval_t[string]]]
  raw_target_t*   = TableRef[string, seq[region_t]]

################################
# INTERVAL TYPES AND FUNCTIONS #
################################


proc is_null*(c: chrom_t): bool = c == -1

proc intersection_both*[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): genomic_interval_t[tuple[l1: T1, l2: T2]] =
  (i1.chrom, max(i1.start, i2.start), min(i1.stop, i2.stop), (i1.label, i2.label))

proc intersection_first*[T1, T2](i1: genomic_interval_t[T1], i2: genomic_interval_t[T2]): genomic_interval_t[T1] =
  if i1.chrom == i2.chrom:
    (i1.chrom, max(i1.start, i2.start), min(i1.stop, i2.stop), i1.label)
  else:
    (i1.chrom, pos_t(0), pos_t(0), i1.label)



proc cookTarget*(orig: raw_target_t, bam: Bam): target_t =
  var chrom_map = newTable[string, chrom_t]()
  for t in bam.hdr.targets:
    chrom_map[t.name] = t.tid
  var cooked = newTable[chrom_t, seq[interval_t[string]]]()
  for chrom_str, intervals in orig:
    doAssert(not (":" in chrom_str), "bad target")
    let chrom = chrom_map[chrom_str]
    cooked[chrom] = @[]
    var last_start = 0
    for i in intervals:
      doAssert(i.chrom == chrom_str, "bad target")
      doAssert(i.start >= last_start)
      let name = if i.name == "": ($i.chrom & ":" & $i.start & "-" & $i.stop) else: i.name
      cooked[chrom].add((pos_t(i.start), pos_t(i.stop), name))
      last_start = i.start
  cooked

proc is_empty*[T](i: genomic_interval_t[T]): bool = i.start >= i.stop

# true when interval i1 is before (without intersection) interval i2
proc `<<`*[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): bool = i1.stop < i2.start
proc `<<`*[T1, T2](i1: interval_t[T1], i2: genomic_interval_t[T2]): bool = i1.stop < i2.start
# true when interval i1 is after (without intersection) interval i2
proc `>>`*[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): bool = i2 << i1
proc `>>`*[T1, T2](i1: interval_t[T1], i2: genomic_interval_t[T2]): bool = i2 << i1
# interval length
#proc len[T](i: interval_t[T]): pos_t = max(pos_t(0), i.stop - i.start)
proc len*[T](i: genomic_interval_t[T]): pos_t = max(pos_t(0), i.stop - i.start)
# true when interval is not empty
proc to_bool*[T](i: genomic_interval_t[T]): bool = i.start < i.stop

# seek target to the right, return true when there is an intersection and move index to the leftmost intersercting interval of the target
iterator intersections*[T](query: genomic_interval_t[T], target: target_t, idx: var target_index_t): genomic_interval_t[tuple[l1: T, l2: string]] =
  #dbEcho("target intersection for query interval:", query, ", starting from:", idx)
  let chrom = query.chrom
  if chrom in target:
    let
      intervals = target[chrom]
      max_idx = len(intervals) - 1
    # if the chrom changes, start from leftmost interval
    if idx.chrom != chrom:
      #dbEcho("target seeking: change chrom to:", idx)
      idx = (chrom, 0)
    # advance target interval until we reach the query interval
    while intervals[idx.interval] << query and idx.interval < max_idx:
      idx.interval += 1
      #dbEcho("target seeking: advance target index to:", idx.interval, "at target interval:", intervals[idx.interval])
    # yield all intersections
    var i = idx.interval
    while not (query << intervals[i]) and i <= max_idx:
      yield intersection_both(query, intervals[i])
      i += 1

iterator intersections2*[T](query: genomic_interval_t[T], target: target_t, idx: var target_index_t): genomic_interval_t[tuple[l1: T, l2: string]] =
  #dbEcho("target intersection for query interval:", query, ", starting from:", idx)
  let chrom = query.chrom
  echo "IS2.1 ", chrom
  if chrom in target:
    let
      intervals = target[chrom]
      max_idx = len(intervals) - 1
    # if the chrom changes, start from leftmost interval
    if idx.chrom != chrom:
      #dbEcho("target seeking: change chrom to:", idx)
      idx = (chrom, 0)
    # advance target interval until we reach the query interval
    while intervals[idx.interval] << query and idx.interval < max_idx:
      idx.interval += 1
      echo "IS2.2 idx=", idx.interval
      #dbEcho("target seeking: advance target index to:", idx.interval, "at target interval:", intervals[idx.interval])
    # yield all intersections
    var i = idx.interval
    try:
      discard intervals[i]
    except Exception as e:
      stderr.writeLine("INTERSECT2:", e.msg)
    while not (query << intervals[i]) and i <= max_idx: 
      echo "IS2.3 idx=", idx.interval
      try:
        let r= intersection_both(query, intervals[i])
        yield r
        i += 1
      except Exception as e:
        stderr.writeLine("-------- INTERSECT2:", e.msg)
        
proc intersects*[T](query: genomic_interval_t[T], target: target_t, idx: var target_index_t): bool =
  for i in intersections(query, target, idx):
    return to_bool(i)
  return false

type
  covopt* = ref object
    verbose, debug: bool
    outputFmt: string # bed, wig
    targetFmt: string # bed, gff
    gffSep: char
    gffId, gffType: string

 
#proc inc_for*(r:region_t) = inc(r.counts_for)
#proc inc_rev*(r:region_t) = inc(r.counts_rev)
proc start*(r: region_t): int {.inline.} = return r.start
proc stop*(r: region_t): int {.inline.} = return r.stop
proc chrom*(r: region_t): string {.inline.} = return r.chrom
proc name*(r: region_t): string {.inline.} = return r.name


#[ proc tostring*(r: region_t, s:var string) {.inline.} =
  # Print a 'region' to string (BED)
  s.set_len(0)
  s.add(r.chrom & "\t" & $r.start & "\t" & $r.stop & "\t")
  if r.name != "":
    s.add(r.name & "\t")
  s.add($r.counts) 
]#

#[ 
proc renderString*(r: region_t, alignmentsPerMillion: float, rpkm, norm, stranded: bool): string  =
  let
    total = r.counts_for + r.counts_rev
    counts = if stranded == true: $(r.counts_for) & "\t" & $(r.counts_rev)
            else: $(total)
  result &= r.name & "\t" & r.chrom & "\t" & $(r.start) & "\t" & $(r.stop) & "\t" & counts
  if rpkm:  
    let kb : float = (r.stop - r.start ) / 1000  
    let RPKM : float = float(total) / alignmentsPerMillion / kb
    result &= "\t" & $RPKM

  if norm:
    let normLen =  float(total) / float(r.stop - r.start)
    result &= "\t" & $normLen
 ]#


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
    s = parse_int(cse[3]) - 1
    e = parse_int(cse[4])
    reg = region_t(chrom: cse[0], start: s, stop: e)#counts_for: 0, counts_rev: 0)
  
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

  # Skip unexpectedly short lines
  if len(cse) < 8:
    return nil

  # Skip non CDS fields (or user provided)
  if cse[2] != gffField:
    return nil

  var
    s, e: int
  
  try:
    s = parse_int(cse[3])  - 1
    e = parse_int(cse[4])
  except:
    stderr.write_line("[warning] fields 4 and 5 are not integers):",  cse[3], ", ", cse[4])
    return nil
  var
    reg = region_t(chrom: cse[0], start: s, stop: e)

  # In the future, 8th field could be requireed [TODO]
  if len(cse) == 9:
    try:
      for gffAnnotPartRaw in cse[8].split(gffSeparator):
        
        let gffAnnotPart = gffAnnotPartRaw.strip(chars = {'"', '\'', ' '})
        
        if gffAnnotPart.startsWith(gffIdentifier):
          let splittedField = gffAnnotPart.split("=")

          # Try splitting on "="
          if len(splittedField) == 2:
            reg.name = splittedField[1].strip(chars = {'"', '\'', ' '})
            break
          else:
            let resplittedField = gffAnnotPart.split(" ")
            if len(resplittedField) == 2:
              reg.name = resplittedField[0].strip(chars = {'"', '\'', ' '})
              break
            else:
              reg.name = "Error"
              break
          
    except Exception as e:
      stderr.write_line("[warning] fields 8 is not a string):",  cse[8], "\n  ", e.msg)
      return nil
 
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
    reg = region_t(chrom: cse[0], start: s, stop: e)
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
    
    var v: region_t
    
    try:
      v = gff_line_to_region($kstr.s, gffField, gffSeparator, gffIdentifier)
    except Exception as e:
      stderr.write_line("[GFF/GTF error]:", e.msg)
      continue

    
    if v == nil:
      continue
    
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region_t]())
    
    bed_regions[v.chrom].add(v)
    
  
  # since it is read into mem, can also well sort.
  for chrom, ivs in bed_regions.mpairs:
    sort(ivs, proc (a, b: region_t): int = a.start - b.start)

  hts.free(kstr.s)
  return bed_regions