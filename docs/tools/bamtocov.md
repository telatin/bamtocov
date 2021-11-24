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
```

## Bam To BedGraph

To produce a [BED graph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html)
file with the coverage of an input BAM file, the syntax can be as minimal as:

```bash
bamtocov input.bam > coverage.bed
```

The output format will be in BED graph, having as columns: chromosome name, start position, end position, coverage.

```text
seq1    0       9       0
seq1    9       109     5
seq1    109     189     0
```

To filter the reads in input:

* `--flag INT` to exclude reads with any of the bits in the FLAG set
* `--mapq INT` to exclude reads with a mapping quality lower than INT

To manipulate the output:

* `--quantize STRING` to produce a quantized output, where the coverage is divided in bins specified in comma separated list as INT,INT...
* `--physical-coverage` to calculate the physical coverage of paired libraries
* and `--stranded` to print a separate column for the forward and negative strands

### Stranded output

When `--stranded` is used, the output is a five column BED-like file, where the columns are:

1. Chromosome name
2. Interval start
3. Interval end
4. Coverage in the forward strand
5. Coverage in the reverse strand

```text
seq1    0       9       0       0
seq1    9       109     5       0
seq1    109     189     0       0
seq1    189     200     2       0
seq1    200     250     4       0
```

## Bam To Wiggle

The [Wiggle](https://genome.ucsc.edu/goldenPath/help/wiggle.html) format is commonly used to prepare
track for genome browsers.
BamToCov can produce a wiggle file with a fixed step specified via `--wig INT`. 
By default the maximum coverage will be printed, but mean and min can be specified via `--op mean|max|min`.

```bash
bamtocov --wig 200 input.bam > coverage.wig
```

The output is similar to:

```text
fixedStep chrom=seq1 span=100
0       5
100     5
200     6
300     3
400     2
500     0
```

An alternative program to perform this analysis is [bam2wig](http://lindenb.github.io/jvarkit/Bam2Wig.html), that will
take into account the CIGAR operations, if required.