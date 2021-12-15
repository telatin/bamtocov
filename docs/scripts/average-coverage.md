---
sort: 1
---

# average-coverage.py

Calculates the average coverage per contig over multiple samples

## Usage

```text
usage: average-coverage.py [-h] [-o OUTPUT] [-t] [-b BIN] bam [bam ...]

Return statistics of coverage of a BAM file using "bamtocov". Will print a
three columns report: chromosome, length, average coverage

positional arguments:
  bam                   BAM file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output file
  -t, --total           Print total coverage
  -b BIN, --bin BIN     bamtocov binary [default: bamtocov]

```

## Example

```bash
average-coverage.py input/{mini,mini2,mini3}.bam -b bin/bamtocov -t
```

Produces:

```text
Chromosome       mini   mini2   mini3
seq1    1.5     1.5     2.3
seq2    1.0     1.0     0.0
seq0    0.0     0.1     0.0
#Total  0.83    0.87    0.77
```