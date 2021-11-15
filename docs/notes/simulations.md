# Simulated BAMs

A script to generate simulated BAM files for a single-chromosome genome ships
with the repository, to measure the effects of read length and coverage on the
performance of the program.

## Using the script

```text
usage: simulate-long-bam.py [-h] -o OUTPUT [-l LENGTH] [-c COVERAGE]
                            [-m MIN_LEN] [-M MAX_LEN] [-s SEED] [-r]
                            [--progress PROGRESS]

Simulate a BAM file with long reads mapped against a hypothetical genome

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output BAM file
  -l LENGTH, --length LENGTH
                        Length of the genome [default: 10M]
  -c COVERAGE, --coverage COVERAGE
                        Approximate coverage [default: 50]
  -m MIN_LEN, --min-len MIN_LEN
                        Minimum read length [default: 1000]
  -M MAX_LEN, --max-len MAX_LEN
                        Minimum read length [default: 10000]
  -s SEED, --seed SEED  Random seed [default: 42]
  -r, --randomcigar     Generate a random CIGAR
  --progress PROGRESS   Print progress every INT reads [default: 10000]
```

In normal mode all reads have no INDELs, while in `--randomcigar` mode
they have random insertion and deletions. Some program will check CIGAR
operations (BamToCov does not).

## Example datasets

The basic dataset is a 10 Mbp genome, 10X coverage:

```bash
# Long reads (genome size: 10 Mbp) with perfect matches and random CIGARs
simulate-long-bam.py -m 1000 -M 10000  -o test/sim/short_match.bam 
simulate-long-bam.py -m 1000 -M 10000 -r -o test/sim/short_cigar.bam 

# Short reads (genome size: 10 Mbp) with perfect matches and random CIGARs
simulate-long-bam.py -m 100 -M 300 -o test/sim/short_match.bam 
simulate-long-bam.py -m 100 -M 300 -r  -o test/sim/short_match.bam 
```

Effect of coverage and genome size have been checked varying the
appropriate parameters.

## Results

### Memory 

Peak memory usage in synthetic BAM files with a single chromosome **10 Mbp long**
and **10X coverage** (long and short refers to read lengths)

| Tool          | long match | long INDELs | short match | short INDELs |
| ------------- | ---------: | ----------: | ----------: | -----------: |
| bamtocov      |      2,172 |       2,176 |       2,176 |        2,236 |
| covtobed      |      3,784 |       3,768 |       3,856 |        3,840 |
| megadepth     |     45,516 |      45,488 |      45,612 |       45,552 |
| mosdepth      |     82,268 |      82,148 |      82,272 |       82,284 |
| mosdepth fast |     82,180 |      82,204 |      82,468 |       82,244 |

### Speed

Long reads (1000 bp to 10,000 bp):

| Command                                    |      Mean (s) | Min (s) | Max (s) |    Relative |
| :----------------------------------------- | ------------: | ------: | ------: | ----------: |
| `megadepth --coverage long_100M_match.bam` | 1.954 ± 0.169 |   1.838 |   2.148 | 1.21 ± 0.13 |
| `bamtocov long_100M_match.bam`             | 1.617 ± 0.106 |   1.510 |   1.722 |        1.00 |
| `covtobed long_100M_match.bam`             | 4.703 ± 0.066 |   4.628 |   4.748 | 2.91 ± 0.19 |
| `mosdepth -x prefix long_100M_match.bam`  | 3.873 ± 0.437 |   3.568 |   4.373 | 2.40 ± 0.31 |

Short reads (100 bp to 300 bp):

| Command                                     |       Mean (s) | Min (s) | Max (s) |    Relative |
| :------------------------------------------ | -------------: | ------: | ------: | ----------: |
| `megadepth --coverage short_100M_match.bam` |  5.552 ± 0.334 |   5.288 |   5.928 |        1.00 |
| `bamtocov short_100M_match.bam`             | 21.147 ± 0.362 |  20.742 |  21.439 | 3.81 ± 0.24 |
| `covtobed short_100M_match.bam`             | 25.036 ± 0.788 |  24.579 |  25.946 | 4.51 ± 0.31 |
| `mosdepth -x prefix short_100M_match.bam`  | 10.778 ± 0.894 |   9.949 |  11.725 | 1.94 ± 0.20 |

## Large files

```bash
simulate-long-bam.py -l 400M -c 100 -m 1000 -M 10000 -x 3 -r --progress 100 -o local/400Mb_100X3_1000-10000.bam
simulate-long-bam.py -l 400M -c 100 -m 100 -M 300    -x 3 -r --progress 100 -o local/400Mb_100X3_100-300.bam
```