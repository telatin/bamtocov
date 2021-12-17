---
sort: 2
---

# low-cov-multisample.py

Identifies low coverage regions common to multiple samples.

Can be configured to use bedtools, bamtocov or megadepth, but
**requires** bedtools for the intersection.

## Usage

```text
usage: low-cov-multisample.py [-h] [--min MIN] [--max MAX] [--shared SHARED]
                              [-t {bamtocov,bedtools,megadepth}]
                              [--bamtocov BAMTOCOV] [--bedtools BEDTOOLS]
                              [--megadepth MEGADEPTH] [--verbose]
                              BAM [BAM ...]

Detect low coverage regions in multiple BAM files using bamtocov, megadepth or
bedtools

positional arguments:
  BAM                   BAM file(s) to process

optional arguments:
  -h, --help            show this help message and exit
  --verbose             Print verbose output

Coverage filters:
  --min MIN             Minimum coverage to consider a region as low coverage
                        [default: 10]
  --max MAX             Maximum coverage to consider a region as low coverage
                        [default: 10]
  --shared SHARED       Ratio of input samples sharing the interval [default:
                        1.0]

External tools:
  -t {bamtocov,bedtools,megadepth}, --tool {bamtocov,bedtools,megadepth}
                        Tool to use for coverage detection
  --bamtocov BAMTOCOV   Path to bamtocov executable
  --bedtools BEDTOOLS   Path to bedtools executable
  --megadepth MEGADEPTH
                        Path to megadepth executable
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
