import os
import hts
import docopt
import heapqueue
import strutils
import tables
import algorithm
 
import ./covutils



#[
  
]# 

type EKeyboardInterrupt = object of CatchableError
 
proc handler() {.noconv.} =
  raise newException(EKeyboardInterrupt, "Keyboard Interrupt")
 
setControlCHook(handler)

type
  coverage = ref object
    forward: int
    reverse: int

type
  interval = ref object
    chr:     string
    start:   int
    stop:    int
    name:    string
    cov:     int
    forward: int
    reverse: int
 
    
var
  debug = false
  physical_coverage = false
  output_strand = false
  outputQueue = newSeq[interval]()


template initClosure(closureName, thisIterator:untyped, minq: uint8, eflag: uint16) =
  let closureName = iterator():auto{.closure.} =
    for thisElement in thisIterator:
      if thisElement.mapping_quality < mapq: continue
      if (thisElement.flag and eflag) != 0: continue
      yield thisElement

proc doAssert(condition: bool, message: string) =
  if condition == false:
    stderr.writeLine("ERROR: ", message)
    quit(1)

proc newCov(f = 0, r = 0): coverage =
  coverage(forward: f, reverse: r)

proc inc(c: coverage, reverse=false) =
  if reverse == false:
    c.forward += 1
  else:
    c.reverse += 1

proc dec(c: coverage, reverse=false) =
  if reverse == false:
    c.forward -= 1
  else:
    c.reverse -= 1

proc tot(c: coverage): int =
  c.forward + c.reverse

proc topStop(q: HeapQueue): int =
  if not q[0].isNil:
    return q[0].stop
  return -1


proc topReverse(q: HeapQueue): bool =
  if not q[0].isNil:
    return q[0].reverse
  return false


proc empty(q: HeapQueue): bool =
  if len(q) == 0:
    return true
  else:
    return false

# Class that stores info about the end position of alignments, used in the alignment queue
type 
  covEnd = ref object
    stop: int
    reverse: bool

proc `<`(a, b: covEnd): bool = a.stop < b.stop

#[ 
proc `$`(i: interval): string =
  i.chr & ":" & $i.start & "-" & $i.stop & " (" & $i.forward & "+" & $i.reverse & "=" & $i.cov & ")X"
]#

proc printBed(q: seq[interval]) =
  var
    outputLine: string

  if len(q) > 0:
    if output_strand == true:
      outputLine =  q[0].chr & "\t" & $q[0].start & "\t" & $q[^1].stop & "\t" & $q[0].forward & "\t" & $q[0].reverse
    else:
      outputLine =  q[0].chr & "\t" & $q[0].start & "\t" & $q[^1].stop & "\t" & $q[0].cov

    outputQueue.setLen(0)
    try:
      stdout.write(outputLine & '\n')
    except IOError:
      quit(0)
    

proc addIntervalToQueue(chrName: string, last_pos, next_change: int, c: coverage) =
  let
    i = interval(chr: chrName, start: last_pos, stop: next_change, cov: c.tot(), forward: c.forward, reverse: c.reverse)
 
  if len(outputQueue) > 0 and (i.chr != outputQueue[^1].chr or i.cov != outputQueue[^1].cov):
    printBed(outputQueue)
  
  if last_pos < next_change:
    outputQueue.add(i)
  
  #echo chrName, "\t", last_pos, "\t", next_change, "\t", c.tot()


proc covtobed(bam:Bam, mapq:uint8, eflag:uint16) = 
  var
    next_change   = 0
    chrSize       = initTable[string, int]()
    aln           : Record
    referenceEnd  : bool
    more_alignments         : bool
    more_alignments_for_ref : bool

  initClosure(nextAlignment,bam.items(), mapq, eflag)

  aln = nextAlignment()

  for reference in bam.hdr.targets():
    referenceEnd = false
    if debug == true:
      stderr.writeLine("# ", reference.name, '\t', reference.length, "\t(newReference)")

    chrSize[reference.name] = int(reference.length)
    
    var 
      last_pos = 0
      coverage_ends = initHeapQueue[covEnd]()
      cov           = newCov()

    while true:
      
      more_alignments         = not aln.isNil
      more_alignments_for_ref = more_alignments and aln.chrom == reference.name
      
      # calculate next change
      if more_alignments_for_ref:
        if coverage_ends.empty():
          next_change = int(aln.start)
        else:
          next_change = min(int(aln.start), coverage_ends.topStop())
      else:
        if coverage_ends.empty():
          next_change = int(reference.length)
        else:
          next_change = coverage_ends.topStop()
         
      # output coverage ...
      doAssert(coverage_ends.len() == cov.tot(), "coverage not equal to queue size")
      addIntervalToQueue(reference.name, last_pos, next_change, cov)

      #if next_change == int(reference.length):
      #  break  

      if debug == true:
        stderr.writeLine("-+-  next=", next_change, "\tMoreAln=", more_alignments, "|", more_alignments_for_ref,";Cov=", cov.tot(), ";Size=", len(coverage_ends))
        if more_alignments:
          stderr.writeLine( " +-> more aln @ chr=", aln.chrom, ",pos=", aln.start)


      
      # increment coverage with aln starting here
      while more_alignments_for_ref and (next_change == aln.start):
        if physical_coverage == true:
          if aln.isize > 0:
            coverage_ends.push(covEnd(stop: int(aln.start + aln.isize), reverse: aln.flag.reverse) )
            cov.inc(aln.flag.reverse)

          #doAssert(false, "coding fiddling doo")
          #[						
            if (alignment.InsertSize > 0) {
						        debug cerr << "   [phy] pos:" << alignment.Position << " size:" << alignment.InsertSize << endl;
							coverage_ends.push({alignment.Position + alignment.InsertSize, alignment.IsReverseStrand()});
							coverage.inc(alignment.IsReverseStrand());
						}
            ]#
        else:
          coverage_ends.push( covEnd(stop: aln.stop, reverse: aln.flag.reverse) )
          cov.inc(aln.flag.reverse)

        aln = nextAlignment()
        more_alignments = not aln.isNil
        more_alignments_for_ref = more_alignments and aln.chrom == reference.name
        
      #decrement coverage with alignments that end here
      while ( not coverage_ends.empty()  and next_change == coverage_ends.topStop()):
        cov.dec(coverage_ends.topReverse())
        discard coverage_ends.pop()
      
      

      # End chromosome loop
#[       if referenceEnd == true:
        doAssert(cov.tot()==0, "coverage not null at the end of chromosome " & reference.name & ": cov.tot=" & $cov.tot() & " = For:" & $cov.forward & "+Rev:" & $cov.reverse )
        doAssert(coverage_ends.len() == 0, "coverage queue not null at the end of chromosome "  & reference.name & ": " & $coverage_ends.len())
        break ]#

      

      if last_pos == int(reference.length):
#        referenceEnd = true
        doAssert(cov.tot()==0, "coverage not null at the end of chromosome " & reference.name & ": cov.tot=" & $cov.tot() & " = For:" & $cov.forward & "+Rev:" & $cov.reverse )
        doAssert(coverage_ends.len() == 0, "coverage queue not null at the end of chromosome "  & reference.name & ": " & $coverage_ends.len())
        break

      last_pos = next_change
      
    # end while ----
    if not coverage_ends.len() == 0:
      stderr.writeLine("Coverage not zero when expected. Try samtools fixmate.")
      raise
  
  # Flush?
  printBed(outputQueue)

  #if more_alignments:
  #  stderr.writeLine("Is the BAM sorted?")
  #ß  raise


proc main(argv: var seq[string]): int =
  let doc = format("""
  covToBed $version

  Usage: covtobed [options] [<BAM>]

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
  let mapq = parse_int($args["--mapq"])

  debug = args["--debug"]
  physical_coverage = args["--physical"]
  output_strand     = args["--stranded"]

  var
    eflag = uint16(parse_int($args["--flag"]))
    threads = parse_int($args["--threads"])
    bam:Bam
    format_gff : bool
    target_file       = $args["--regions"]


  var
    regions = if format_gff == true: gff_to_table($args["--regions"])
                 else: bed_to_table($args["--regions"])

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

  try:
    covtobed(bam, uint8(mapq), eflag)
  except:
    stderr.writeLine("FATAL ERROR: Unable to read input input: ", $args["<BAM>"] )
    stderr.writeLine( getCurrentExceptionMsg() )
    quit(1)
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