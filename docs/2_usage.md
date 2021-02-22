---
sort: 2
permalink: /usage
---

# bamtocov

```text
BamToCov 2.0.2

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
```