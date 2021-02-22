import os
import hts
import docopt
import heapqueue
import strutils
import tables
import algorithm 
import ./covutils

# FEATURES paper
# TODO multi-bam report
# TODO quantized output
# TODO max-min coverage?
# TODO number of bases under X per target
# TODO no output option, if one only wants the report
# TODO WIG

# FEATURES piu' tardi
# TODO multi-bam coverage?

# NON FEATURES
# TODO output_t fa un po' un pasticcio con intervalli target sovrapposti perche' cerca di appiccicarli, non mi e' chiaro cosa vogliamo ottenere comunque
# TODO usare Record tid (numeric id of contig) instead of chromosome string, may be faster?
# TODO add explicit check for sortedness
# TODO use array for coverage?


################################
# INTERVAL TYPES AND FUNCTIONS #
################################
type
  chrom_t = string # reference name
  #chrom_t = int # reference id
  pos_t = int64
  # here intervals have a "label", which contains additional information besides the location
  # there is one interval without explicit chromosome to be used in the target table, where intervals are already grouped by chromosome
  interval_t[T] = tuple[start, stop: pos_t, label: T]
  genomic_interval_t[T] = tuple[chrom: chrom_t, start, stop: pos_t, label: T]


#proc intersection[T1, T2](i1: interval_t[T1], i2: interval_t[T2]): interval_t[T1] =
#  (max(i1.start, i2.start), min(i1.stop, i2.stop), i1.label)

# return the intersection of two intervals with the first label
#proc intersection[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): genomic_interval_t[T1] =
#  (i1.chrom, max(i1.start, i2.start), min(i1.stop, i2.stop), i1.label)
# return the intersection of two intervals with both labels
proc intersection_both[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): genomic_interval_t[tuple[l1: T1, l2: T2]] =
  (i1.chrom, max(i1.start, i2.start), min(i1.stop, i2.stop), (i1.label, i2.label))
#proc intersection[T1, T2](i1: interval_t[T1], i2: genomic_interval_t[T2]): interval_t[T1] =
#  (max(i1.start, i2.start), min(i1.stop, i2.stop), i1.label)
#proc intersection[T1, T2](i1: genomic_interval_t[T1], i2: genomic_interval_t[T2]): interval_t[T1] =
#  if i1.chrom == i2.chrom:
#    (i1.chrom, max(i1.start, i2.start), min(i1.stop, i2.stop), i1.label)
#  else:
#    (i1.chrom, 0, 0, i1.label)

#proc isEmpty(i: genomic_interval): bool = 
#  i.start >= i.stop 

# true when interval i1 is before (without intersection) interval i2
proc `<<`[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): bool = i1.stop < i2.start
proc `<<`[T1, T2](i1: interval_t[T1], i2: genomic_interval_t[T2]): bool = i1.stop < i2.start
# true when interval i1 is after (without intersection) interval i2
proc `>>`[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): bool = i2 << i1
proc `>>`[T1, T2](i1: interval_t[T1], i2: genomic_interval_t[T2]): bool = i2 << i1
# interval length
proc len[T](i: genomic_interval_t[T]): pos_t = max(pos_t(0), i.stop - i.start)
# true when interval is not empty
proc to_bool[T](i: genomic_interval_t[T]): bool = i.start < i.stop


type
  coverage_t = object
    forward: int
    reverse: int

var
  debug = false

proc db(things: varargs[string, `$`]) =
  if debug:
    stderr.write("debug:")
    for t in things:
      if t[0] == ',':
        stderr.write(t)
      else:
        stderr.write(" " & t)
    stderr.write("\n")



type
  target_t = TableRef[chrom_t, seq[interval_t[string]]]


proc convertTarget(orig: TableRef[string, seq[region_t]]): target_t =
  var conv = newTable[chrom_t, seq[interval_t[string]]]()
  for chrom, intervals in orig:
    doAssert(not (":" in chrom), "bad target")
    conv[chrom] = @[]
    var last_start = 0
    for i in intervals:
      doAssert(i.chrom == chrom, "bad target")
      doAssert(i.start >= last_start)
      let name = if i.name == "": ($i.chrom & ":" & $i.start & "-" & $i.stop) else: i.name
      conv[chrom].add((pos_t(i.start), pos_t(i.stop), i.name))
      last_start = i.start
  conv


proc `$`[T](i: genomic_interval_t[T]): string =
  $i.chrom & ":" & $i.start & "-" & $i.stop & $i.label

type
  target_index_t = tuple[chrom: chrom_t, interval: int]

# seek target to the right, return true when there is an intersection and move index to the leftmost intersercting interval of the target
iterator intersections[T](query: genomic_interval_t[T], target: target_t, idx: var target_index_t): genomic_interval_t[tuple[l1: T, l2: string]] =
  db("target intersection for query interval:", query, ", starting from:", idx)
  let chrom = query.chrom
  if chrom in target:
    let
      intervals = target[chrom]
      max_idx = len(intervals) - 1
    # if the chrom changes, start from leftmost interval
    if idx.chrom != chrom:
      db("target seeking: change chrom to:", idx)
      idx = (chrom, 0)
    # advance target interval until we reach the query interval
    while intervals[idx.interval] << query and idx.interval < max_idx:
      idx.interval += 1
      db("target seeking: advance target index to:", idx.interval, "at target interval:", intervals[idx.interval])
    # yield all intersections
    var i = idx.interval
    while not (query << intervals[i]) and i < max_idx:
      yield intersection_both(query, intervals[i])
      i += 1
  else:
    db("target seeking: chrom", chrom, "not in target")
    
proc intersects[T](query: genomic_interval_t[T], target: target_t, idx: var target_index_t): bool =
  for i in intersections(query, target, idx):
    return to_bool(i)
  return false

type
  input_option_t = tuple[min_mapping_quality: uint8, eflag: uint16, physical: bool, target: target_t]
proc alignment_stream(bam: Bam, opts: input_option_t): iterator (): genomic_interval_t[bool] =
  result = iterator(): genomic_interval_t[bool] {.closure.} =
    var
      b = bam
      o = opts
      target_idx: target_index_t
    for r in b:
      # alignment filter
      if r.mapping_quality < o.min_mapping_quality or (r.flag and o.eflag) != 0:
        continue

      # alignment processing
      var stop: pos_t = 0
      if o.physical:
        if r.isize > 0: 
          stop = r.start + r.isize
        else:
          continue # skip the mate with negative insert size
      else:
        stop = r.stop
      let i = (chrom_t(r.chrom), pos_t(r.start), stop, r.flag.reverse)

      # return alignment if there is any intersection with target (or if there is no target)
      if len(o.target) == 0 or i.intersects(o.target, target_idx):
        yield i


proc doAssert(condition: bool, message: string) =
  if condition == false:
    stderr.writeLine("ERROR: ", message)
    quit(1)

proc newCov(f = 0, r = 0): coverage_t =
  coverage_t(forward: f, reverse: r)

proc inc(c: var coverage_t, reverse=false) =
  if reverse == false:
    c.forward += 1
  else:
    c.reverse += 1

proc dec(c: var coverage_t, reverse=false) =
  if reverse == false:
    c.forward -= 1
  else:
    c.reverse -= 1

proc tot(c: coverage_t): int =
  c.forward + c.reverse

proc topStop(q: HeapQueue): int64 =
  if not q[0].isNil:
    return q[0].stop
  return -1


proc topReverse(q: HeapQueue): bool =
  if not q[0].isNil:
    return q[0].reverse
  return false


proc empty(q: HeapQueue): bool = len(q) == 0

# Class that stores info about the end position of alignments, used in the alignment queue
type 
  covEnd = ref object
    stop: pos_t
    reverse: bool

proc `<`(a, b: covEnd): bool = a.stop < b.stop

proc `$`(c: coverage_t): string =
  "c=" & $(c.forward + c.reverse) & "(" & $c.forward & "+/" & $c.reverse & "-)"

#type 
#  output_format_t = enum
#    BedFormat
#    WigFormat
#  output_handler_t = object
#    target: target_t,
#    output_format: output_format_t,
#    output_strand: bool,
#    last_chrom: chrom_t,
#    last_start, last_end: pos_t,
#    last_coverage: coverage


type
  coverage_interval_t = genomic_interval_t[tuple[l1: coverage_t, l2: string]] # l2 is the target interval or the chromosome
  #coverage_end_t = tuple[stop: pos_t, rev: bool]
  #coverage_t = tuple[forward, reverse: int]

proc coverage_iter(bam: Bam, opts: input_option_t): iterator(): coverage_interval_t =
  result = iterator(): coverage_interval_t {.closure.} =
    var
      next_alignment                                     = alignment_stream(bam, opts)
      next_change             : pos_t                    = 0
      aln                     : genomic_interval_t[bool] = next_alignment()
      more_alignments         : bool                     = not finished(next_alignment)
      target_idx              : target_index_t

    
    for reference in bam.hdr.targets():
      let
        reflen = pos_t(reference.length)
        refname = chrom_t(reference.name)
      
      db("new reference start:", refname, ",", reflen, "bp")

      var 
        last_pos : pos_t = 0
        coverage_ends = initHeapQueue[covEnd]()
        cov           = newCov()
        more_alignments_for_ref = more_alignments and aln.chrom == refname

      while true:
        db("Aln:", if more_alignments: $aln else: "no more", "in", refname)

        
        # calculate next change
        if more_alignments_for_ref:
          if coverage_ends.empty():
            next_change = aln.start
          else:
            next_change = min(aln.start, coverage_ends.topStop())
        else:
          if coverage_ends.empty():
            next_change = reflen
          else:
            next_change = coverage_ends.topStop()

        if debug: stderr.writeLine("Last pos: " & $last_pos & ", next pos: " & $next_change)
        doAssert(last_pos < next_change, "coverage went backwards from " & $last_pos & " to " & $next_change)   
        # output coverage ...
        doAssert(coverage_ends.len() == cov.tot(), "coverage not equal to queue size")
        #let cov_inter: genomic_interval_t[coverage_t] = 
        if len(opts.target) == 0:
          #yield cov_inter
          yield (refname, last_pos, next_change, (cov, refname))
        else:
          for i in intersections((refname, last_pos, next_change, cov), opts.target, target_idx):
            yield i

        if next_change == reflen:
          break  

        if debug:
          stderr.writeLine("-+-  next=", next_change, "\tMoreAln=", more_alignments, "|", more_alignments_for_ref,";Cov=", cov.tot(), ";Size=", len(coverage_ends))
          if more_alignments:
            stderr.writeLine( " +-> more aln @ chr=", aln.chrom, ",pos=", aln.start)

        # increment coverage with aln starting here
        while more_alignments_for_ref and (next_change == aln.start):
          coverage_ends.push(covEnd(stop: aln.stop, reverse: aln.label))
          cov.inc(aln.label)
          if debug: stderr.writeLine("Added aln: " & $aln)

          aln = next_alignment()
          more_alignments = not finished(next_alignment) # we need to check this after each next_alignment
          more_alignments_for_ref = more_alignments and aln.chrom == refname
          
        #decrement coverage with alignments that end here
        while not coverage_ends.empty() and next_change == coverage_ends.topStop():
          cov.dec(coverage_ends.topReverse())
          discard coverage_ends.pop()
        
        

        # End chromosome loop
        if last_pos == reflen:
          doAssert(cov.tot()==0, "coverage not null at the end of chromosome " & refname & ": cov.tot=" & $cov.tot() & " = For:" & $cov.forward & "+Rev:" & $cov.reverse )
          doAssert(coverage_ends.len() == 0, "coverage queue not null at the end of chromosome "  & refname & ": " & $coverage_ends.len())
          break

        last_pos = next_change
        
      # end while ----
      if not coverage_ends.len() == 0:
        stderr.writeLine("Coverage not zero when expected. Try samtools fixmate.")
        raise
    doAssert(not more_alignments, "Is the BAM sorted?") # FIXME put more explicit check for sortednees


type
  output_t = object
    queued: genomic_interval_t[coverage_t]
    output_strand: bool

proc write_output(o: var output_t, i: genomic_interval_t[coverage_t]) =
  if i.start < i.stop: # skip empty intervals
    if o.output_strand:
      echo $i.chrom & "\t" & $i.start & "\t" & $i.stop & "\t" & $i.label.forward & "\t" & $i.label.reverse
    else:
      echo $i.chrom & "\t" & $i.start & "\t" & $i.stop & "\t" & $(i.label.forward + i.label.reverse)


proc push_interval(o: var output_t, i: coverage_interval_t) =
  let
    q = o.queued
    c: coverage_t = i.label.l1
  #let debug=true
  if debug: stderr.writeLine("push_interval: " & $i)
  if q.label == c and q.stop == i.start and q.chrom == i.chrom: # extend previous interval
    if debug: stderr.writeLine("push_inteval: extend " & $q)
    o.queued.stop = i.stop # FIXME se setto direttamente q.stop la modifica viene persa?
  else: # output previous interval and queue the new one
    o.write_output(q)
    o.queued = (i.chrom, i.start, i.stop, c)

proc `=destroy`(o: var output_t) =
  o.write_output(o.queued)

proc newOutput(output_strand: bool): output_t = 
  output_t(
    queued: (chrom_t(""), pos_t(0), pos_t(0), newCov()),
    output_strand: output_strand
  )

type cov_t = array[0..1, int64]
proc to_array(c: coverage_t): cov_t = [int64(c.forward), int64(c.reverse)]

type #FIXME statistics for rev and for
  stats_t = tuple[tot_bases, min_cov, max_cov: cov_t, tot_length: pos_t]
  target_stat_t = TableRef[string, stats_t]

proc `*`(c: cov_t, l: pos_t): cov_t = [c[0]*l, c[1]*l]
proc `/`(c: cov_t, l: pos_t): array[2, float] = [float(c[0])/float(l), float(c[1])/float(l)]
proc `+`(c1, c2: cov_t): cov_t = [c1[0] + c2[0], c1[1] + c2[1]]
proc min(c1, c2: cov_t): cov_t = [min(c1[0], c2[0]), min(c1[1], c2[1])]
proc max(c1, c2: cov_t): cov_t = [max(c1[0], c2[0]), max(c1[1], c2[1])]
proc tot(c: cov_t): int64 = c[0] + c[1]
proc tot(c: array[2, float]): float = c[0] + c[1]

proc push_interval(self: var target_stat_t, i: coverage_interval_t) =
  let
    l = len(i)
    name: string = i.label.l2 # target interval name
    cov: cov_t = to_array(i.label.l1)
  #let o = if i.label in self: self[i.label] else: (0, c, c, 0)
  if name in self:
    let o = self[name]
    self[name] = (o.tot_bases + cov*l, min(o.min_cov, cov), max(o.max_cov, cov), o.tot_length + l)
  else:
    self[name] = (cov*l, cov, cov, l)
proc `$`(c: cov_t): string =
  $c[0] & "+/" & $c[1] & "-"

proc `$`(s: stats_t): string =
  #$s.min_cov & " " & $(float(s.tot_bases)/float(s.tot_length)) & " " & $s.max_cov & " " & $s.tot_length
  $s.min_cov & " " & $(s.tot_bases/s.tot_length) & " " & $s.max_cov & " " & $s.tot_length
proc to_string(s: stats_t, strand: bool): string =
  if strand:
    $s.min_cov & " " & $(s.tot_bases/s.tot_length) & " " & $s.max_cov & " " & $s.tot_length
  else:
    $tot(s.min_cov) & " " & $tot(s.tot_bases/s.tot_length) & " " & $tot(s.max_cov) & " " & $s.tot_length

proc main(argv: var seq[string]): int =
  let doc = format("""
  BamToCov $version

  Usage: bamtocov [options] [<BAM>]

Arguments:                                                                                                                                                 
  <BAM>          the alignment file for which to calculate depth (default: STDIN)

Core options:
  -p, --physical               Calculate physical coverage
  -s, --stranded               Report coverage separate by strand
  -w, --wig <SPAN>             Output in wig format (using fixed <SPAN>)

Target files:
  -r, --regions <bed>          Target file in BED or GFF format (detected with the extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]

BAM reading options:
  -T, --threads <threads>      BAM decompression threads [default: 0]
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]

Other options:
  --debug                      Enable diagnostics    
  -h, --help                   Show help
  """ % ["version", version])

  let args = docopt(doc, version="covtobed " & $version, argv=argv)


  debug = args["--debug"]
  let output_strand: bool     = args["--stranded"] # FIXME e' giusto covertirlo cosi'?

  let
    threads = parse_int($args["--threads"])
    target_file       = $args["--regions"]
  var
    bam:Bam
    format_gff : bool # FIXME metterci un valore sensato?

  let
    input_opts: input_option_t = (
      min_mapping_quality: uint8(parse_int($args["--mapq"])),
      eflag: uint16(parse_int($args["--flag"])),
      physical: bool(args["--physical"]), # FIXME e' giusto convertirlo cosi'?
      target: convertTarget(
        if format_gff: gff_to_table(target_file)
        else: bed_to_table(target_file)
      )
    )
    
  if $args["<BAM>"] != "nil":
    # Read from FILE
    try:
      if not fileExists($args["<BAM>"]):
        stderr.writeLine("FATAL ERROR: File <", $args["<BAM>"], "> not found")
        quit(1)
      open(bam, $args["<BAM>"], threads=threads)
    except:
      stderr.writeLine("FATAL ERROR: Unable to read input file: ", $args["<BAM>"]) 
      quit(1)
  else:
    # Read STDIN
    stderr.writeLine("Reading from STDIN [Ctrl-C to break]")
    open(bam, "-", threads=threads)

  # TODO copiato da mosdepth, serve a sveltire il parsing per i CRAM
#  var opts = SamField.SAM_FLAG.int or SamField.SAM_RNAME.int or SamField.SAM_POS.int or SamField.SAM_MAPQ.int or SamField.SAM_CIGAR.int
#  if not fast_mode:
#      opts = opts or SamField.SAM_QNAME.int or SamField.SAM_RNEXT.int or SamField.SAM_PNEXT.int #or SamField.SAM_TLEN.int
#  discard bam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, opts)
#  discard bam.set_option(FormatOption.CRAM_OPT_DECODE_MD, 0)



  try:
    var
      output: output_t = newOutput(output_strand) 
      target_stats: target_stat_t = newTable[string, tuple[tot_bases, min_cov, max_cov: cov_t, tot_length: pos_t]]()
      cov_iter = coverage_iter(bam, input_opts)
    db("target:", input_opts.target)
    for cov_inter in cov_iter():
      db("coverage:", cov_inter)
      target_stats.push_interval(cov_inter)
      output.push_interval(cov_inter)

    #var order = newOrderedSet()
    #for chrom, intervals:
    #  order.incl()
    for name, stats in target_stats:
      stderr.writeLine( name & "\t" & to_string(stats, true) )
 

  except:
    stderr.writeLine("FATAL ERROR: Unable to read input input: ", $args["<BAM>"] )
    stderr.writeLine( () )
    quit(1)
  return 0

type EKeyboardInterrupt = object of CatchableError
 
proc handler() {.noconv.} =
  raise newException(EKeyboardInterrupt, "Keyboard Interrupt")
 
setControlCHook(handler)


when isMainModule:
  var args = commandLineParams()
  try:
    discard main(args)
  except EKeyboardInterrupt:
    stderr.writeLine( "Quitting.")
  except:
    stderr.writeLine( getCurrentExceptionMsg() )
    quit(1)   