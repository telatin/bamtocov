# bamtocov

Tools to extract coverage informations from BAM (and CRAM) files, based on the
[covtobed](https://github.com/telatin/covtobed) algorithm that supports stranded
coverage and physical coverage, input from streams 
and uses a memory-efficient algorithm.

## covtobed
```
covToBed 2.0.0

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
```

## covtotarget
will count the _total nucleotide coverage_ per feature in a BED or GFF file using as input the output of [covtobed](https://github.com/telatin/covtobed) (also from STDIN).
```
covToTarget  

  Usage: covtotarget [options] <Target> [<covtobed-output>]

Arguments:                                                                                                                                                 

  <Target>           the BED (or GFF) file containing regions in which to count reads
  <covtobed-output>  covtobed output, or STDIN if not provided

Options:

  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]
  -l, --norm-len               Normalize by gene length
  -b, --bed-output             Output format is BED-like (default is feature_name [tab] counts)
  -h, --help                   Show help
```

Example, can be used in a stream from the BAM emitter to covtobed:
```bash
cat input/mini.bam | covtobed | covtotarget input/mini.gff
```

Where _covtobed_ output is:
```text
seq1    0       9       0
seq1    9       109     5
[...]
seq2    499     599     10
seq2    599     1000    0
```

and `covtocounts` output is (extracts):
```text
MGLILCEK_00002  0
MGLILCEK_00003  51
MGLILCEK_00010  1000
```
with `--norm-len` and `--bed-output`:
```
seq0    299     400     ZERO_COV_CHR_2  0.0
seq1    199     400     MGLILCEK_00001  3.487562189054726
seq1    599     650     MGLILCEK_00002  0.0
```


## covtocounts
will count the _number of alignments_ in a BAM file per feature of a target BED or GFF file (basically, adds GFF support to `read-count` found in [nim-hts-tools](https://github.com/brentp/hts-nim-tools))
```
covToCounts 0.4.1

  Usage: covtocounts [options] <Target> <BAM-or-CRAM>

Arguments:                                                                                                                                                 

  <Target>       the BED (or GFF) file containing regions in which to count reads
  <BAM-or-CRAM>  the alignment file for which to calculate depth

Options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -r, --fasta <fasta>          FASTA file for use with CRAM files [default: ].
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]
  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]
  -n, --rpkm                   Add a RPKM column
  -l, --norm-len               Add a counts/length column (after RPKM when both used)
  --header                     Print header
  --debug                      Enable diagnostics    
  -h, --help                   Show help
```

## References
 * Brent Pedersen,  Aaron Quinlan, [hts-nim: scripting high-performance genomic analyses](https://academic.oup.com/bioinformatics/article/34/19/3387/4990493) (Bioinformatics)
 * Giovanni Birolo, Andrea Telatin, [covtobed: a simple and fast tool to extract coverage tracks from BAM files](https://joss.theoj.org/papers/10.21105/joss.02119) (JOSS)
