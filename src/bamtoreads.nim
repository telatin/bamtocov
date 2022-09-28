import os
import hts
import docopt
import strutils
import tables
import algorithm

type EKeyboardInterrupt = object of CatchableError

const NimblePkgVersion {.strdefine.} = "prerelease"

let
  version = NimblePkgVersion
proc handler() {.noconv.} =
  raise newException(EKeyboardInterrupt, "Keyboard Interrupt")
 
setControlCHook(handler)

proc charToQual*(c: char, offset = 33): int =
  ## returns Illumina quality score for a given character
  c.ord - offset

proc qualToChar*(q: int, offset = 33): char =
  ## returns character for a given Illumina quality score
  (q+offset).char


proc qualToChar*(q: seq[uint8], offset = 33): string =
  ## returns string for a given Illumina quality score
  for i in q:
    result &= qualToChar(int(i), offset)

proc fastq(id, comment, sequence, quality: string): string =
  ## returns a FASTQ string for a given sequence and quality
  let spacer = if len(comment) > 0: " "
           else: ""
  "@" & id & spacer & comment &  "\n" & sequence & "\n+\n" & quality & "\n"

proc main(argv: var seq[string]): int =
  let env_fasta = getEnv("REF_PATH")
  let doc = format("""
  BamCountRefs $version

  Usage: bamcountrefs [options]  <BAM-or-CRAM>...

Arguments:                                                                                                                                                 
 
  <BAM-or-CRAM>  the alignment file for which to calculate depth

Main options:
  -l, --list <contigs>         List of contigs to rescue

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
  
  var fasta: cstring 
  if $args["--fasta"] != "nil":
    fasta = cstring($args["--fasta"])

  var
    eflag = uint16(parse_int($args["--flag"]))
    threads = parse_int($args["--threads"])
    bam:Bam
    chromosomes = initTable[int, string]()

  var
    samples = @[columnName]

  for bamFile in @(args["<BAM-or-CRAM>"]):
    var sampleName = extractFilename(bamFile)
    samples.add(sampleName.split('.')[0])
    try:
      open(bam, cstring(bamFile), threads=threads, index=true, fai=fasta)
    except:
      stderr.writeLine("Unable to open BAM file: ", bamFile )
         

    if bam.idx == nil:
      stderr.write_line("ERROR: requires BAM/CRAM index")
      quit(1)


    # Preload chrnames
    for i in bam.hdr.targets:
        chromosomes[i.tid] = i.name

    for chrom in bam.hdr.targets:
      var 
        s: string
        q: seq[uint8]
      #stderr.writeLine ">", chrom.name, " (", chrom.tid, ")"

      # Unsorted
      for aln in bam:

      #for aln in bam.query(chrom.tid):
        let
          pair = if aln.flag.read1: "/1 "
                  elif aln.flag.read2: "/2 "
                  else: " /nopair "
          comment = pair & chromosomes[aln.tid] & ":" & $(aln.start) & "-" & $(aln.stop)
        echo fastq(aln.qname, comment, aln.sequence(s), aln.base_qualities(q).qualToChar(33))
        echo aln.qname, "\t", chromosomes[aln.mate_tid], ":", aln.mate_pos, "\t", aln.mapping_quality, "\t", aln.flag
        #

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