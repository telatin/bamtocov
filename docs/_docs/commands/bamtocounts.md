---
title: BamToCounts
parent: Tools
nav_order: 2
---
 

A program that counts the number of reads per target in
a BAM file. The BAM file does not need to be sorted.

By default, overlapping features will be reported each with all the
counts, so a read mapped in a shared interval will be counted
more than once. This can be changed with `--exclude-ambiguous`
(drop reads/fragments overlapping more than one feature, as
featureCounts does by default) or `--fraction` (split the credit
1/n across the n overlapping features, as `featureCounts -O --fraction`).
The two options are mutually exclusive.

Other options to reproduce featureCounts filters are
`--exclude-multimapping` (drop reads with `NH` tag > 1) and,
with `--paired`, `-B/--both-mapped` (require the mate to be mapped).

## Help screen

```text
BamToCounts 

  Usage: bamtocounts [options] <Target> <BAM-or-CRAM>...

Arguments:

  <Target>       the BED (or GFF) file containing regions in which to count reads
  <BAM-or-CRAM>  the alignment file for which to calculate depth

Options:

  -T, --threads <threads>      BAM decompression threads [default: 0]
  -j, --jobs <jobs>            BAM files to process concurrently; 0 uses available CPUs [default: 0]
  -r, --fasta <fasta>          FASTA file for use with CRAM files [default: ]
  -F, --flag <FLAG>            Exclude reads with any of the bits in FLAG set [default: 1796]
  -Q, --mapq <mapq>            Mapping quality threshold [default: 1]
  --exclude-multimapping       Exclude reads with NH tag > 1 (matches featureCounts'
                               default multi-mapping filter; no effect if the aligner
                               does not write an NH tag)
  --paired                     Count fragments not reads: count each pair once over
                               the full fragment span (either mate overlapping counts)
  -B, --both-mapped            Require the mate to be mapped (matches featureCounts' -B;
                               only excludes mate_unmapped reads, unlike --proper-pairs
                               which also excludes discordant-but-both-mapped pairs)
  --proper-pairs               Require reads from paired experiments to be properly paired
  --fraction                   Split an ambiguous fragment's credit 1/n across the n
                               regions it overlaps, instead of crediting each with a
                               full 1 (matches featureCounts -O --fraction). Mutually
                               exclusive with --exclude-ambiguous.
  --exclude-ambiguous          Drop a fragment entirely if it overlaps more than one
                               region, instead of crediting every region it overlaps
                               (matches featureCounts' own default, without -O).
                               Mutually exclusive with --fraction.
  --strict                     Read must be contained, not just overlap, with feature
  --stranded                   Print strand-specific counts
  --coords                     Also print coordinates of each feature

  -g, --gff                    Force GFF for input (otherwise autodetected by .gff extension)
  -t, --type <feat>            GFF feature type to parse [default: CDS]
  -i, --id <ID>                GFF identifier [default: ID]

  -n, --rpkm                   Add a RPKM column
  -m, --tpm                    Add a TPM (Transcripts Per Million) column
  -l, --norm-len               Add a counts/length column (after RPKM when both used)
  -p, --precision INT          Digits for floating point precision [default: 3]
  --header                     Print header
  --verbose                    Print verbose output
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

If `--tpm`:
* TPM (Transcripts Per Million)

If `--norm-len`:
* Normalized counts (per target length)

Counts are integers, unless `--fraction` is used: fractional
credits are then printed with `--precision` decimal digits.

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
