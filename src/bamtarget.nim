import os
import hts
import docopt
import heapqueue
import strutils
import tables
import algorithm 
import ./covutils
import sets 

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



proc `$`[T](i: genomic_interval_t[T]): string =
  $i.chrom & ":" & $i.start & "-" & $i.stop & $i.label





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


 
template doAssert(condition: bool, message: string) =  
  if condition == false:
    stderr.writeLine("ERROR: ", message)
    quit(1)

 
proc main(argv: var seq[string]): int =
  let doc = format("""
  BamTarget $version

  Usage: BamTarget [options] <Target> [<BAM>]

Arguments:                                                                                                                                                 
  <Target>         A target to be parsed (GFF, GTF, BED)
  <BAM>            A BAM file linked to the target

Target files:
  -r, --regions <bed>          Target file in BED or GFF3/GTF format (detected with the extension)
  -t, --gff-type <feat>        GFF feature type to parse [default: CDS]
  -i, --gff-id <ID>            GFF identifier [default: ID]
  --gff-separator <sep>        GFF attributes separator [default: ;]
  --gff                        Force GFF input (otherwise assumed by extension .gff)
  --gtf                        Force GTF input (otherwise assumed by extension .gtf)

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
    gffField = $args["--gff-type"]
    gffSeparator = $args["--gff-separator"]
    gffIdentifier = $args["--gff-id"]

  var
    format_gff = false 
    format_gtf = false
  
  # Set target format (GFF/BED) using extension or forced by the user
  if ($args["<Target>"]).toLower().contains(".gff"):
    dbEcho("Parsing target as GFF")
    format_gff = true
  elif ($args["<Target>"]).toLower().contains(".gtf"):
    format_gtf = true
  else:
    dbEcho("Parsing target as BED")

  if debug and (format_gff or format_gtf):
    stderr.writeLine "GFF field: ", gffField, "\n",
                      "GFF separator: ", gffSeparator, "\n",
                      "GFF identifier: ", gffIdentifier, "\n" 

  if $args["<Target>"] != "nil":
    if not fileExists($args["<Target>"]):
      stderr.writeLine("ERROR: Target file not found:", $args["<Target>"])
      quit(1)

    if (args["--gtf"] and args["--gff"]) or (format_gff and format_gtf):
      echo "ERROR: Target format is ambiguous: specify a GFF or GTF target (ideally autoinferred from the extension)"
      quit(1)
    elif args["--gtf"]:
      format_gtf = true
    elif args["--gff"]:
      format_gff = true
  else:
    echo "ERROR: Target file not specified."
    quit(1)

  var target: TableRef[system.string, seq[region_t]]
  if format_gff: 
    if debug:
      stderr.writeLine "Parsing target as GFF"
      target = gff_to_table($args["<Target>"], gffField, gffSeparator, gffIdentifier) 
  elif format_gtf: 
    if debug:
      stderr.writeLine "Parsing target as GTF"
      
    target = gtf_to_table($args["<Target>"], gffField, gffSeparator, gffIdentifier) 
  else: 
    if debug:
      stderr.writeLine "Parsing target as BED"
    target = bed_to_table($args["<Target>"])
   
  if debug:
    stderr.writeLine "Target parsed: ", len(target), " reference sequences"

  for t in target.keys():
    if debug:
      echo "# Chromosome:" , t
    for region in target[t]:
      echo region.chrom, "\t", region.start, "\t", region.stop, "\t", region.name
  
  if fileExists($args["<BAM>"]):
    if debug:
      stderr.writeLine "BAM file: ", $args["<BAM>"]
    var
      bam: Bam
    open(bam, cstring($args["<BAM>"]), threads=1)
    let cooked = cookTarget(target, bam)
    for t in cooked.keys():
      if debug:
        echo "# Chromosome:" , t
      for region in cooked[t]:
        echo "[",t,"]\t", region.start, "\t", region.stop, "\t", region.label
  
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