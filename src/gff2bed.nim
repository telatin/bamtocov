import os 
import docopt 
import strutils
import tables
import ./covutils


proc main(argv: var seq[string]): int =
  let doc = format("""
  gff2bed $version

  Usage: gff2bed [options] <GFF> 

Input:
  <GFF>                        GFF3 input file (or GTF)

Output options:
  -a, --antitarget             Include antitarget in the output
  -x, --exclude-target         Do not include target intervals
  -o, --output=FILE            Output file (default: stdout)

Input parsing options:
  -t, --gff-type <feat>        GFF feature type to parse [default: CDS]
  -i, --gff-id <ID>            GFF identifier [default: ID]
  --gff-separator <sep>        GFF attributes separator [default: ;]

Other options:
  --debug                      Enable diagnostics
  --verbose                    Enable verbose output
  -h, --help                   Show help
  """ % ["version", version])

  let args = docopt(doc, version=version, argv=argv)

  let
    verbose = args["--verbose"]
    outputfile = $args["--output"]
    antitarget = args["--antitarget"]
    exclude_target = args["--exclude-target"]

    gffType = $args["--gff-type"]
    gffId = $args["--gff-id"]
    gffSeparator = $args["--gff-separator"]

  if not fileExists($args["<GFF>"]):
    stderr.writeLine "Error: file not found: ", $args["<GFF>"]
    quit(1)

  if ( ".gff" notin ($args["<GFF>"]).toLower() ) and ( ".gtf" notin ($args["<GFF>"]).toLower() ):
    stderr.writeLine "WARNING: file does not have .gff extension: ", $args["<GFF>"]

  if verbose:
    stderr.writeLine "Reading GFF file: ", $args["<GFF>"]
    stderr.writeLine "GFF parameters: type=", gffType, " id=", gffId, " separator=", gffSeparator
  let
    annotation = gff_to_table($args["<GFF>"], gffType, gffSeparator, gffId)

  var
    outputText = ""
 

  # gff_to_table*(bed: string): TableRef[string, seq[region_t]] =
  for chrom, intervals in annotation:
      for j in intervals:
          if outputfile == "nil":
            echo chrom, "\t", j.start, "\t", j.stop, "\t", j.name
          else:
            outputText &= chrom & "\t" & $(j.start) & "\t" & $(j.stop) & "\t" & j.name & "\n"


  if outputfile != "nil":
    try:
      writeFile(outputfile, outputText)
    except Exception as e:
        echo outputText
        stderr.writeLine "Error writing output file: ", e.msg

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