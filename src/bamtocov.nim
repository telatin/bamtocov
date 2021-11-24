import os
import hts
import docopt
import heapqueue
import strutils
import tables
import algorithm 
import ./covutils
import sets
#import nimprof

# ✅ Fixed: total min/max coverage cannot be computed from the two strands! - DONE
# ✅ Fixed: per la cronaca ho trovato un mini “bachetto”, se il BED non ha nomi, giustamente, ficca tutto in un mega intervallo immaginario. Ora,  potrebbe essere la cosa giusta da fare, l’alternativa è che se il nome non c’è lo creiamo noi tipo “chr2:100-200" e cosi li manteniamo forzatamente separati e se uno vuole il megatarget specifica lo stesso nome in tutto il file
# ✅ Fixed: mini2.bam coverage went backwards error - DONE era il check troppo zelante!
# ✅ Fixed: report without target using chromosomes

# FEATURES paper
# ✅ TODO multi-bam report 
# ✅ TODO quantized output 
# TODO max-min coverage?
# TODO number of bases under X per target
# ✅ TODO no output option, if one only wants the report
# ✅ TODO WIG 

# FEATURES  
# TODO multi-bam coverage?

# NON FEATURES
# TODO output_t fa un po' un pasticcio con intervalli target sovrapposti perche' cerca di appiccicarli, non mi e' chiaro cosa vogliamo ottenere comunque
# TODO usare Record tid (numeric id of contig) instead of chromosome string, may be faster?
# TODO add explicit check for sortedness


################################
# INTERVAL TYPES AND FUNCTIONS #
################################
type
  chrom_t = int # reference id
  pos_t = int64
  # here intervals have a "label", which contains additional information besides the location
  # there is one interval without explicit chromosome to be used in the target table, where intervals are already grouped by chromosome
  interval_t[T] = tuple[start, stop: pos_t, label: T]
  genomic_interval_t[T] = tuple[chrom: chrom_t, start, stop: pos_t, label: T]


proc is_null(c: chrom_t): bool = c == -1
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
proc intersection_first[T1, T2](i1: genomic_interval_t[T1], i2: genomic_interval_t[T2]): genomic_interval_t[T1] =
  if i1.chrom == i2.chrom:
    (i1.chrom, max(i1.start, i2.start), min(i1.stop, i2.stop), i1.label)
  else:
    (i1.chrom, pos_t(0), pos_t(0), i1.label)

proc is_empty[T](i: genomic_interval_t[T]): bool = i.start >= i.stop

# true when interval i1 is before (without intersection) interval i2
proc `<<`[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): bool = i1.stop < i2.start
proc `<<`[T1, T2](i1: interval_t[T1], i2: genomic_interval_t[T2]): bool = i1.stop < i2.start
# true when interval i1 is after (without intersection) interval i2
proc `>>`[T1, T2](i1: genomic_interval_t[T1], i2: interval_t[T2]): bool = i2 << i1
proc `>>`[T1, T2](i1: interval_t[T1], i2: genomic_interval_t[T2]): bool = i2 << i1
# interval length
#proc len[T](i: interval_t[T]): pos_t = max(pos_t(0), i.stop - i.start)
proc len[T](i: genomic_interval_t[T]): pos_t = max(pos_t(0), i.stop - i.start)
# true when interval is not empty
proc to_bool[T](i: genomic_interval_t[T]): bool = i.start < i.stop


type
  coverage_t = object
    forward: int
    reverse: int

var
  debug = false
  developer = false

proc db(things: varargs[string, `$`]) =
  #if debug:
    stderr.write("debug:")
    for t in things:
      if t[0] == ',':
        stderr.write(t)
      else:
        stderr.write(" " & t)
    stderr.write("\n")

proc dev(things: varargs[string, `$`]) =
  stderr.write("dev:")
  for t in things:
    if t[0] == ',':
      stderr.write(t)
    else:
      stderr.write(" " & t)
  stderr.write("\n")

template devEcho(things: varargs[string, `$`]) =
  if developer:
    dev(things)

template dbEcho(things: varargs[string, `$`]) =
  if debug:
    db(things)
type
  target_t = TableRef[chrom_t, seq[interval_t[string]]]
  raw_target_t = TableRef[string, seq[region_t]]


proc cookTarget(orig: raw_target_t, bam: Bam): target_t =
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


proc `$`[T](i: genomic_interval_t[T]): string =
  $i.chrom & ":" & $i.start & "-" & $i.stop & $i.label

type
  target_index_t = tuple[chrom: chrom_t, interval: int]

# seek target to the right, return true when there is an intersection and move index to the leftmost intersercting interval of the target
iterator intersections[T](query: genomic_interval_t[T], target: target_t, idx: var target_index_t): genomic_interval_t[tuple[l1: T, l2: string]] =
  dbEcho("target intersection for query interval:", query, ", starting from:", idx)
  let chrom = query.chrom
  if chrom in target:
    let
      intervals = target[chrom]
      max_idx = len(intervals) - 1
    # if the chrom changes, start from leftmost interval
    if idx.chrom != chrom:
      dbEcho("target seeking: change chrom to:", idx)
      idx = (chrom, 0)
    # advance target interval until we reach the query interval
    while intervals[idx.interval] << query and idx.interval < max_idx:
      idx.interval += 1
      dbEcho("target seeking: advance target index to:", idx.interval, "at target interval:", intervals[idx.interval])
    # yield all intersections
    var i = idx.interval
    while not (query << intervals[i]) and i < max_idx:
      yield intersection_both(query, intervals[i])
      i += 1
  else:
    dbEcho("target seeking: chrom", chrom, "not in target")
    
proc intersects[T](query: genomic_interval_t[T], target: target_t, idx: var target_index_t): bool =
  for i in intersections(query, target, idx):
    return to_bool(i)
  return false

type
  input_option_t = tuple[min_mapping_quality: uint8, eflag: uint16, physical: bool, target: raw_target_t]

proc alignment_stream(bam: Bam, opts: input_option_t, target: target_t): iterator (): genomic_interval_t[bool] =
  result = iterator(): genomic_interval_t[bool] {.closure.} =
    var
      o = opts
      target_idx: target_index_t
    for r in bam:
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
      let i = (r.tid, pos_t(r.start), stop, r.flag.reverse)

      # return alignment if there is any intersection with target (or if there is no target)
      if len(target) == 0 or i.intersects(target, target_idx):
        yield i


#[ 
proc doAssert(condition: bool, message: string) =
  if condition == false:
    stderr.writeLine("ERROR: ", message)
    quit(1) 
]#
template doAssert(condition: bool, message: string) = # FIXME is this already in the standard library?
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

type
  coverage_interval_t = genomic_interval_t[tuple[l1: coverage_t, l2: string]] # l2 is the target interval or the chromosome
  #coverage_end_t = tuple[stop: pos_t, rev: bool]
  #coverage_t = tuple[forward, reverse: int]

proc coverage_iter(bam: Bam, opts: input_option_t, target: target_t): iterator(): coverage_interval_t =
  result = iterator(): coverage_interval_t {.closure.} =
    var
      next_alignment                                     = alignment_stream(bam, opts, target)
      next_change             : pos_t                    = 0
      aln                     : genomic_interval_t[bool] = next_alignment()
      more_alignments         : bool                     = not finished(next_alignment)
      target_idx              : target_index_t

    
    for reference in bam.hdr.targets():
      let
        reflen: pos_t = pos_t(reference.length)
        refname = reference.name
        refid: chrom_t = reference.tid
      
      dbEcho("new reference start:", refname, ",", reflen, "bp")

      var 
        last_pos : pos_t = 0
        coverage_ends = initHeapQueue[covEnd]()
        cov           = newCov()
        more_alignments_for_ref = more_alignments and aln.chrom == refid

      while true:
        dbEcho("Aln:", if more_alignments: $aln else: "no more", "in", refname)

        
        # calculate the position of the next coverage change
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
        # check that we are advancing, this should always be the case if the bam is sorted (except when there is an alignment right at the beginning of a chromosome, hence the alternative condition)
        doAssert(last_pos < next_change or (last_pos == 0 and next_change == 0), 
          "coverage went backwards from " & $last_pos & " to " & $next_change & ", at " & refname & ":" & $aln.start)  
        # output coverage ...
        doAssert(coverage_ends.len() == cov.tot(), "coverage not equal to queue size")
        #let cov_inter: genomic_interval_t[coverage_t] = 
        if len(target) == 0:
          yield (refid, last_pos, next_change, (cov, refname))
        else:
          dbEcho("pre-intersection coverage:", (refname, last_pos, next_change, cov))
          for i in intersections((refid, last_pos, next_change, cov), target, target_idx):
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
          more_alignments_for_ref = more_alignments and aln.chrom == refid
          
        # decrement coverage with alignments that end here
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


#output_coverage_types = enum
#  oc_int_strand,    # integer strand-specific coverage output
#  oc_quant_strand,  # quantized strand-specific coverage output
#  oc_int_tot,       # integer total coverage output
#  oc_quant_tot      # quantized total coverage output
#type
#  output_coverage_t = object
#  case output_strand: bool
#  of true: tuple[forward, reverse: int64]
#  of false: total: int64]

#  output_coverage_t = tuple[total, forward, reverse: int64]

type
  output_format_t = enum # supported output formats
    of_bed, # BED format
    of_wig_fixstep, # fixed step WIG format
  span_func_t = enum sf_max, sf_min, sf_mean # which function use to summarize coverage in WIG span
  output_option_t = tuple[
    strand: bool, # output strand-specific coverage and stats
    # coverage specific options
    no_coverage: bool, # do not output coverage, only stats
    quantization: string, # 
    output_format: output_format_t, # output format
    span_length: pos_t, # span for wig output format
    span_func: span_func_t,
    # stats specific options
    low_cov: int64 # report length of regions under low_cov in stats
  ]

  output_t = object
    queued: genomic_interval_t[coverage_t]
    # data fields needed for spanned output (wig)
    current_span: genomic_interval_t[coverage_t] 

    opts: output_option_t
    quantization_index2label: seq[string]
    quantization_coverage2index: seq[int]
    chrom2str: TableRef[chrom_t, string]

proc write_output(o: var output_t, i: genomic_interval_t[coverage_t]) =
  if i.start < i.stop: # skip empty intervals
    case o.opts.output_format:
      of of_bed: # bed format output
        let interval_str = o.chrom2str[i.chrom] & "\t" & $i.start & "\t" & $i.stop & "\t"
        let coverage_str =
          if o.opts.strand:
            if len(o.quantization_index2label) > 0:
              "\t" & o.quantization_index2label[i.label.forward] & "\t" & o.quantization_index2label[i.label.reverse]
            else:
              $i.label.forward & "\t" & $i.label.reverse
          else:
            if len(o.quantization_index2label) > 0:
              o.quantization_index2label[i.label.forward]
            else:
              $i.label.forward
        echo interval_str & coverage_str
      of of_wig_fixstep: 
        if len(o.quantization_index2label) > 0:
          stderr.writeLine("wig output does not support quantized coverage")
          raise
        if o.opts.strand:
          stderr.writeLine("wig output does not support stranded coverage")
          raise
        let span_length = o.opts.span_length
        if o.current_span.chrom != i.chrom: # start new contig
          if o.current_span.chrom != -1: # output last possibly incomplete span from previous chrom
            let span_value = case o.opts.span_func:
              of sf_max, sf_min: $o.current_span.label.forward
              of sf_mean: $(float(o.current_span.label.forward)/float(span_length)) # FIXME the actual span is less than span_length!
            echo $o.current_span.start & "\t" & span_value
            
          o.current_span.chrom = i.chrom
          o.current_span.start = 0
          o.current_span.stop = o.current_span.start + span_length
          o.current_span.label.forward = 0
          echo "fixedStep chrom=" & o.chrom2str[o.current_span.chrom] & " start=1 step=" & $span_length & " span=" & $span_length
        
        while o.current_span.start <= i.stop:
          let inter = intersection_first(o.current_span, i)
          if not is_empty(inter): # update the current span value
            o.current_span.label.forward = case o.opts.span_func:
              of sf_max: max(o.current_span.label.forward, i.label.forward)
              of sf_min: min(o.current_span.label.forward, i.label.forward)
              of sf_mean: o.current_span.label.forward + i.label.forward*int(len(inter))
          if inter.stop == o.current_span.stop: # span is concluded
            # output span
            let span_value = case o.opts.span_func:
              of sf_max, sf_min: $o.current_span.label.forward
              of sf_mean: $(float(o.current_span.label.forward)/float(span_length))
            echo $o.current_span.start & "\t" & span_value
            # next span
            o.current_span.start += span_length
            o.current_span.stop = o.current_span.start + span_length
            o.current_span.label.forward = case o.opts.span_func
              of sf_max, sf_mean: 0
              of sf_min: high(int)
          else: # span extends beyond the interval, we are done
            break

proc push_interval(o: var output_t, i: coverage_interval_t) =
  let q = o.queued
  var c: coverage_t = i.label.l1

  dbEcho("push_interval: ", i)
  # handle stranded output
  if not o.opts.strand:
    c.forward = c.forward + c.reverse
    c.reverse = 0
  # handle quantized output
  let qmax_cov = len(o.quantization_coverage2index)
  if qmax_cov > 0:
    c.forward = o.quantization_coverage2index[min(c.forward, qmax_cov - 1)]
    if o.opts.strand:
      c.reverse = o.quantization_coverage2index[min(c.reverse, qmax_cov - 1)]

  if q.label == c and q.stop == i.start and q.chrom == i.chrom: # extend previous interval
    if debug: stderr.writeLine("push_inteval: extend " & $q)
    o.queued.stop = i.stop # FIXME se setto direttamente q.stop la modifica viene persa?
  else: # output previous interval and queue the new one
    o.write_output(q)
    o.queued = (i.chrom, i.start, i.stop, c)

proc `=destroy`(o: var output_t) =
  case o.opts.output_format:
    of of_bed: o.write_output(o.queued)
    of of_wig_fixstep: o.write_output((chrom_t(-1), pos_t(0), pos_t(0), newCov()))

import sequtils
# output quantization
proc parse_quantization(o: var output_t, breaks: string) =
  var bb: seq[int] = @[0]
  for b in split(breaks, ','):
    bb.add(parse_int(b))
  bb = deduplicate(sorted(bb), isSorted=true)
  let
    max_right = max(bb)
  for i in 0..(len(bb) - 2):
    let
      left = bb[i]
      right = bb[i + 1] - 1
    o.quantization_index2label.add($left & "-" & $right)
    for c in left..right:
      o.quantization_coverage2index.add(i)
  o.quantization_index2label.add($max_right & "-")
  o.quantization_coverage2index.add(len(o.quantization_index2label) - 1)
  #dev("index2label:", $o.quantization_index2label)
  #dev("coverage2index:", $o.quantization_coverage2index)

proc newOutput(opts: output_option_t, bam: Bam): output_t = 
  var o = output_t(
    opts: opts,
    queued: (chrom_t(-1), pos_t(0), pos_t(0), newCov()),
    quantization_index2label: @[],
    quantization_coverage2index: @[],
    chrom2str: newTable[chrom_t, string](),
    current_span: (chrom_t(-1), pos_t(0), pos_t(0), newCov())
  )
  for t in bam.hdr.targets:
    o.chrom2str[t.tid] = t.name
  if opts.quantization != "nil":
    o.parse_quantization(opts.quantization)
  o

type
  coverage_stats_t[T] = tuple[total, forward, reverse: T]
  interval_stats_t = tuple[bases, min_cov, max_cov: coverage_stats_t[int64], low_length: coverage_stats_t[pos_t], length: pos_t]
  target_stat_t = ref object
    opts: output_option_t
    stats: TableRef[string, interval_stats_t]

#proc bin_apply[T](f: proc (T, T) -> T, s1, s2: coverage_stats_t[T]): coverage_stats_t[T] =
#  (f(s1.total, s2.total), f(s1.forward, s2.forward), f(s1.reverse, s2.reverse))
#proc `*`[T](s: coverage_stats_t[T], x: T): coverage_stats_t[T] =
#  (s.total*x, s.forward*x, s.reverse*x)

proc `max`[T](s1, s2: coverage_stats_t[T]): coverage_stats_t[T] =
  (max(s1.total, s2.total), max(s1.forward, s2.forward), max(s1.reverse, s2.reverse))
proc `min`[T](s1, s2: coverage_stats_t[T]): coverage_stats_t[T] =
  (min(s1.total, s2.total), min(s1.forward, s2.forward), min(s1.reverse, s2.reverse))
proc `+`[T](s1, s2: coverage_stats_t[T]): coverage_stats_t[T] =
  (s1.total + s2.total, s1.forward + s2.forward, s1.reverse + s2.reverse)
proc `+`[T](s: coverage_stats_t[T], x: T): coverage_stats_t[T] =
  (s.total + x, s.forward + x, s.reverse + x)
proc `*`[T](s: coverage_stats_t[T], x: T): coverage_stats_t[T] =
  (s.total*x, s.forward*x, s.reverse*x)
 
proc new_stats(opts: output_option_t): target_stat_t =
  target_stat_t(opts: opts, stats: newTable[string, interval_stats_t]())

proc push_interval(self: var target_stat_t, i: coverage_interval_t) =
  # update coverage statistics
  let
    l = len(i)
    name: string = i.label.l2 # target interval name
    cov: coverage_stats_t[int64] = (int64(i.label.l1.forward + i.label.l1.reverse), int64(i.label.l1.forward), int64(i.label.l1.reverse))
    low_cov = self.opts.low_cov
    low_length: coverage_stats_t[int64] = (
      (if cov.total   < low_cov: int64(l) else: 0),
      (if cov.forward < low_cov: int64(l) else: 0),
      (if cov.reverse < low_cov: int64(l) else: 0)
    )
  #let prev: interval_stats_t =
  #  if name in self.stats:
  #    self.stats[name]
  #  else:
  #    (bases: )
  self.stats[name] = 
    if name in self.stats:
      let o = self.stats[name]
      (
        bases: o.bases + (cov*l), 
        min_cov: min(o.min_cov, cov), 
        max_cov: max(o.max_cov, cov), 
        low_length: o.low_length + low_length,
        length: o.length + l
      )
    else:
      (
        bases: cov*l, 
        min_cov: cov,
        max_cov: cov,
        low_length: low_length,
        length: l
      )

proc mean(s: interval_stats_t): coverage_stats_t[float] = 
  let l = float(s.length)
  (float(s.bases.total)/l, float(s.bases.forward)/l, float(s.bases.reverse)/l)

proc to_string[T](s: coverage_stats_t[T], strand: bool = true, sep: string=" "): string =
  if strand:
    $s.total & sep & $s.forward & sep & $s.reverse
  else:
    $s.total

proc stat_columns(self: target_stat_t, sep: string = " ", prefix: string=""): string =
  let
    cov_cols = @["bases", "mean", "min", "max"] & (if self.opts.low_cov > 0: @["low" & $self.opts.low_cov] else: @[])
    strand_suffixes = if self.opts.strand: @["", "_forward", "_reverse"] else: @[""]
  var r: string = ""
  for mid in cov_cols:
    for suf in strand_suffixes:
      r = r & sep & prefix & mid & suf
      #dev("report: header cols:", r)
  r & sep & prefix & "length"

proc to_string(self: target_stat_t, name: string, sep: string = " "): string =
  let strand = self.opts.strand
  if name in self.stats:
    let s = self.stats[name]
    var r =
      to_string(s.bases, strand, sep) & sep & 
      to_string(s.mean, strand, sep) & sep & 
      to_string(s.min_cov, strand, sep) & sep & 
      to_string(s.max_cov, strand, sep) & sep
    if self.opts.low_cov > 0:
      r = r & to_string(s.low_length, strand, sep) & sep
    r & $s.length
  else: 
    # FIXME this should not be happening if we do things correctly in main
    # this is when an interval has no alignments
    let nullstring = "0\t"
    var r = nullstring.repeat( (1 + (if self.opts.low_cov > 0: 5 else: 4)*(if strand: 3 else: 1)) )
    # remove last tab from string
    r[0 .. ^2] 
    
    #for i in 1..(1 + (if self.opts.low_cov > 0: 5 else: 4)*(if strand: 3 else: 1)):
    #  r = r & sep & "/" # TODO vedere se c'e' un modo di ripetere una string invece di stamparla tutte queste volte
    #r

# process coverage from a single file
# open bam, compute coverage, print coverage output (based on outopts) and return coverage stats
proc bam2stats(bam_path: string, inopts: input_option_t, outopts: output_option_t, bam_threads: int = 0): tuple[stats: target_stat_t, target: target_t] =
  var
    bam: Bam
    target_stats: target_stat_t = new_stats(outopts)

  if bam_path == "-":
    stderr.writeLine("Reading from STDIN [Ctrl-C to break]")
  open(bam, bam_path, threads=bam_threads)
  let target = cookTarget(inopts.target, bam)
  var output: output_t = newOutput(outopts, bam)

# TODO copiato da mosdepth, serve a sveltire il parsing per i CRAM
#  var opts = SamField.SAM_FLAG.int or SamField.SAM_RNAME.int or SamField.SAM_POS.int or SamField.SAM_MAPQ.int or SamField.SAM_CIGAR.int
#  if not fast_mode:
#      opts = opts or SamField.SAM_QNAME.int or SamField.SAM_RNEXT.int or SamField.SAM_PNEXT.int #or SamField.SAM_TLEN.int
#  discard bam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, opts)
#  discard bam.set_option(FormatOption.CRAM_OPT_DECODE_MD, 0)

  var cov_iter = coverage_iter(bam, inopts, target)
  for cov_inter in cov_iter():
    dbEcho("coverage:", cov_inter)
    target_stats.push_interval(cov_inter)
    if not outopts.no_coverage:
      output.push_interval(cov_inter)
  (target_stats, target)


proc main(argv: var seq[string]): int =
  let doc = format("""
  BamToCov $version

  Usage: bamtocov [options] [<BAM>]...

Arguments:                                                                                                                                                 
  <BAM>         the alignment file for which to calculate depth (default: STDIN)

Core options:
  -p, --physical               Calculate physical coverage
  -s, --stranded               Report coverage separate by strand
  -q, --quantize <breaks>      Comma separated list of breaks for quantized output
  -w, --wig <SPAN>             Output in WIG format (using fixed <SPAN>), 0 will print in BED format [default: 0]
  --op <func>                  How to summarize coverage for each WIG span (mean/min/max) [default: max]
  -o, --report <TXT>           Output coverage report
  --skip-output                Do not output per-base coverage
  --report-low <min>           Report coverage for bases with coverage < min [default: 0]

Target files:
  -r, --regions <bed>          Target file in BED or GFF3/GTF format (detected with the extension)
  -t, --gff-type <feat>        GFF feature type to parse [default: CDS]
  -i, --gff-id <ID>            GFF identifier [default: ID]
  --gff-separator <sep>        GFF attributes separator [default: ;]
  --gff                        Force GFF input (otherwise assumed by extension .gff)

BAM reading options:
  -T, --threads <threads>      BAM decompression threads [default: 0]
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]

Other options:
  --debug                      Enable diagnostics
  -h, --help                   Show help
  """ % ["version", version])

  let args = docopt(doc, version=version, argv=argv)

  #

  debug = args["--debug"]

  if debug:
    dbEcho("args:", args)


  let
    threads = parse_int($args["--threads"])
    target_file       = $args["--regions"]
  var
    #bams: seq[BAM]
    format_gff = false 
  
  # Set target format (GFF/BED) using extension or forced by the user
  if ($args["--regions"]).toLower().contains(".gff") or ($args["--regions"]).toLower().contains(".gtf")  or args["--gff"]:
    dbEcho("Parsing target as GFF")
    format_gff = true
  else:
    dbEcho("Parsing target as BED")
  
  
  assert( $args["--op"] in  @["mean", "min", "max"], "--op must be one of mean, min, max, got: " & $args["--op"])

  let 
    gffField = $args["--gff-type"]
    gffSeparator = $args["--gff-separator"]
    gffIdentifier = $args["--gff-id"]

  let
    input_paths =
      if len(@(args["<BAM>"])) > 0:
        @(args["<BAM>"])
      else:
        @["-"]
    input_opts: input_option_t = (
      min_mapping_quality: uint8(parse_int($args["--mapq"])),
      eflag: uint16(parse_int($args["--flag"])),
      physical: bool(args["--physical"]),
      target: if format_gff: gff_to_table(target_file, gffField, gffSeparator, gffIdentifier) else: bed_to_table(target_file)
    )
    output_opts: output_option_t = (
      strand: bool(args["--stranded"]),
      no_coverage: bool(args["--skip-output"]) or len(input_paths) > 1,
      quantization: $args["--quantize"],
      output_format: if (parseInt($args["--wig"]) > 0): of_wig_fixstep else: of_bed,
      span_length: pos_t(parse_int($args["--wig"])),
      span_func: if $args["--op"] == "max": sf_max
                elif $args["--op"] == "min": sf_min
                else: sf_mean,
      low_cov: int64(parse_int($args["--report-low"]))
    )



  #Preflight check input files
  var
    missing_files = 0
  # FIXME warn about not output for multiple bams (if --no-ouput is not given)
  for inputBam in input_paths:
    if inputBam == "-":
      dbEcho("Will read STDIN")
    elif not fileExists(inputBam):
      missing_files += 1
      stderr.writeLine("ERROR: Input file <", inputBam, "> not found.")
  if missing_files > 0:
    quit(1)

  # Multiple BAMs
  if len(input_paths) > 1:
    if not bool(args["--skip-output"]):
      stderr.writeLine("WARNING: coverage output for multiple input files is not implemented, so it will not be produced; use --skip-output to suppress this warning")
    if not args["--regions"]:
      stderr.writeLine("ERROR: Multiple BAMs are handled via target file (--regions). Supply a target.")
      quit(1)
 
  # CIAO ho condensato qui tutta la ciccia dei calcoli, in bam2stats, cosi' dovrebbe essere piu' semplice da mettere in threads
  #interval_t[T] = tuple[start, stop: pos_t, label: T]
  var bam_stats: seq[tuple [stats: target_stat_t, target: target_t]]
  for p in input_paths:
    dbEcho("running", p)
    bam_stats.add(bam2stats(p, input_opts, output_opts, bam_threads=1))
     
  
  if args["--report"]: # print report table
    dbEcho("stats reporting")
    # assemble table index
    var index = initOrderedSet[string]()
    let
      sample_names = input_paths # TODO qui forse possiamo fare un po' meglio che mettere tutto il path ;-)
      target = bam_stats[0].target # get cooked target from the first bam
    
    if len(input_opts.target) > 0:
      # get interval names from target in the order they appear
      for chrom, intervals in target:
        for t in intervals:
          index.incl(t.label)
    else:
      # if there is no target, use chromosomes
      for s in bam_stats:
        for t in s.stats.stats.keys():
          index.incl(t)

    # check that each interval in stats has been put in index
    for s in bam_stats:
      for t in s.stats.stats.keys():
        if not (t in index):
          doAssert t in index, "t not in index: " & t
    dbEcho("target:", target)
    dbEcho("index:", index)

    # print header
    dbEcho("report: header")
    let
      report = open($args["--report"], fmWrite)
      sep = "\t"
    report.write("interval")
    for x in zip(sample_names, bam_stats):
      report.write(stat_columns(x[1].stats, sep=sep, prefix=x[0] & "_"))
    report.write("\n")

    # print body
    dbEcho("report: body")
    for t in index:
      report.write(t)
      for s in bam_stats:
        report.write(sep & to_string(s.stats, t, sep=sep))
      report.write("\n")

  #except:
  #  stderr.writeLine("FATAL ERROR: Unable to read input input: ", $args["<BAM>"] ) # FIXME handle multiple bams here!
  #  stderr.writeLine( () )
  #  quit(1)
  dbEcho("exiting successfully!")
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