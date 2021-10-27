---
sort: 2
---

# BamToCounts

A program that counts the number of reads per target in
a BAM file.

## Help screen

```text
  BamToCounts $version

  Usage: bamtocounts [options] <Target> <BAM-or-CRAM>...

Arguments:                                                                                                                                                 

  <Target>       the BED (or GFF) file containing regions in which to count reads
  <BAM-or-CRAM>  the alignment file for which to calculate depth

Options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -r, --fasta <fasta>          FASTA file for use with CRAM files [default: $env_fasta].
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

## Example output

Given a BED file with the annotation and a BAM file, the
output produced by the program is a table with the following
columns:

1. chromosome
2. start
3. end
4. feature name
5. read count
6. Optional columns (RPKM, counts/length)

```text
seq1    200     400     target1_8X      8
seq1    600     650     target2_0X      0
seq1    700     900     target3_1X      1
seq2    530     540     for_rev_10Xa    10
seq2    570     580     for_rev_10Xb    10
seq2    590     610     for_rev_10Xc    10
```
