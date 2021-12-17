---
sort: 3
---

# make-target-from-bam.py

Generate a BED target comprising the whole length of all 
the chromosomes.

## Usage

```text
usage: make-target-from-bam.py [-h] [-o OUT] bam

Generate a target covering the full chromosomes from a BAM file

positional arguments:
  bam                BAM file

optional arguments:
  -h, --help         show this help message and exit
  -o OUT, --out OUT  Output file
```

## Example

```bash
make-target-from-bam.py input/mini.bam  
```

Produces a BED having as feature name the chromosome names.

```text
seq1    0       1000    seq1
seq2    0       1000    seq2
seq0    0       1000    seq0
```
