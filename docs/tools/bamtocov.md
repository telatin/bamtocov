---
sort: 1
---

# BamToCov

The main program of the suite that parses a sorted BAM file and
generates a coverage file in BED format.

## Splash screen

```text
BamToCov 2.2.0

  Usage: bamtocov [options] [<BAM>]...

Arguments:                                                                                                                                                 
  <BAM>         the alignment file for which to calculate depth (default: STDIN)

Core options:
  -p, --physical               Calculate physical coverage
  -s, --stranded               Report coverage separate by strand
  -q, --quantize <breaks>      Comma separated list of breaks for quantized output
  -w, --wig <SPAN>             Output in wig format (using fixed <SPAN>)
  -o, --report <TXT>           Output coverage report
  --skip-output                Do not output per-base coverage
  --report-low <min>

Target files:
  -r, --regions <bed>          Target file in BED or GFF format (detected with the extension)
  -t, --gff-type <feat>        GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]
  --gff                        Force GFF input (otherwise assumed by extension .gff)

BAM reading options:
  -T, --threads <threads>      BAM decompression threads [default: 0]
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]

Other options:
  --debug                      Enable diagnostics
  -h, --help                   Show help
```