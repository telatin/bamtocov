---
sort: 4
---

# feat-counts.py

Counts feature across multiple samples producing a count matrix.

```text
usage: feat-counts.py [-h] -t TARGET [-o OUTPUT] [--feat FEAT] [--id ID]
                      [--binary BINARY] [--verbose]
                      BAM [BAM ...]

Count reads in multiple BAM files using a target file (BED or GFF or GTF)

positional arguments:
  BAM                   BAM file

optional arguments:
  -h, --help            show this help message and exit
  -t TARGET, --target TARGET
                        Target file
  -o OUTPUT, --output OUTPUT
                        Output file
  --feat FEAT           Feature type [default: exon]
  --id ID               ID attribute [default: gene_id]
  --binary BINARY       Binary to bamtocounts [default: bamtocounts]
  --verbose             Verbose output
```

## Example

:warning: All the input files must be sorted and indexed.

Using the files in the repository it is possible to test the tool
using a minimal dataset:

```bash
feat-counts.py -t input/regions.bed input/min{i,i2,i3}.bam 
```

## Example output

By default the matrix is produced with four columns (BED coordinates and name)
and then the counts for each sample. The sample name is the filename, removed the `.bam` extension.

```text
Chr     Start   Stop    Name             mini3  mini2   mini
seq0    0       210     empty0_0         0      1       0
seq1    5       112     include_5        5      5       5
seq1    410     532     overlap_2        1      2       2
seq1    800     950     empty1_0         0      0       0
seq2    300     532     shared1_10       0      10      10
seq2    566     769     shared2_10       0      10      10
```