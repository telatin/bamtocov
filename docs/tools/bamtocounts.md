---
sort: 2
---

# BamToCounts

A program that counts the number of reads per target in
a BAM file.

## Help screen

```text
BamToCounts 2.5.0

  Usage: bamtocounts [options] <Target> <BAM-or-CRAM>...

Arguments:                                                                                                                                                 

  <Target>       the BED (or GFF) file containing regions in which to count reads
  <BAM-or-CRAM>  the alignment file for which to calculate depth

Options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -r, --fasta <fasta>          FASTA file for use with CRAM files [default: ].
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 0]
  --stranded                   Print strand-specific counts
  --coords                     Also print coordinates of each feature

  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]

  -n, --rpkm                   Add a RPKM column
  -l, --norm-len               Add a counts/length column (after RPKM when both used)
  -p, --precision INT          Digits for floating point precision [default: 3]
  --header                     Print header
  --debug                      Enable diagnostics    
  -h, --help                   Show help
```

## Example output

Given a BED file with the annotation and a BAM file, the
output produced by the program is a table with the following
columns:

* Feature identifier (always present)

If `--coords` is enabled:

* Chromosome(s)
* Start(s)
* Ends(s)

If `--stranded` is enabled:
* Counts Forward
* Counts Reverse
Otherwise 
* Counts

If `--rpkm`:
* RPKM (Reads Per Kilobase per Million)


If `--norm`:
* Normalized counts (per target length)

### Examples

Using only `--header`:
```text
#Feature        counts
feature_1       7
feature_2       7
second_half     0
```

Adding `--stranded`:
```text
#Feature        for     rev
feature_1       7       0
feature_2       0       7
second_half     0       0
```


Adding `--stranded` and `--coords`:
```text
#Feature        Chrom   start   end     for     rev
feature_1       seq1    100     200     7       0
feature_2       seq1    300     400     0       7
second_half     seq1    500     1000    0       0
```