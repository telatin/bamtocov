---
sort: 3
---

# BamCountsRefs

A program to build a count table from multiple BAM
files (having the same reference sequence).

```text
BamCountRefs 2.2.0

  Usage: bamcountrefs [options]  <BAM-or-CRAM>...

Arguments:                                                                                                                                                 
 
  <BAM-or-CRAM>  the alignment file for which to calculate depth

BAM/CRAM processing options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -r, --fasta <fasta>          FASTA file for use with CRAM files [default: ].
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
```

## Example

```bash
bin/bamcountrefs --tag "Chrom" input/mini.bam input/mini2.bam  
```

Output:

```text
Chrom   mini    mini2
seq0    0       1
seq1    15      15
seq2    10      10
```