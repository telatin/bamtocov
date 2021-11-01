# Speed tests

## Panel, human, 16 genes

### Comparison with bedtools

```bash
$ hyperfine "covtobed $FILE" "./bin/bamtocov $FILE" "bedtools genomecov -bga -ibam $FILE" --max-runs 6
```

```text
Benchmark #1: covtobed ../covtobed/example_data/panel_02.bam
  Time (mean ± σ):      1.110 s ±  0.034 s    [User: 1.056 s, System: 0.039 s]
  Range (min … max):    1.068 s …  1.149 s    6 runs
 
Benchmark #2: ./bin/bamtocov ../covtobed/example_data/panel_02.bam
  Time (mean ± σ):     811.6 ms ±  27.8 ms    [User: 1.117 s, System: 0.125 s]
  Range (min … max):   783.0 ms … 847.2 ms    6 runs
 
Benchmark #3: bedtools genomecov -bga -ibam ../covtobed/example_data/panel_02.bam
  Time (mean ± σ):     42.297 s ±  4.590 s    [User: 14.759 s, System: 27.363 s]
  Range (min … max):   39.610 s … 51.546 s    6 runs
 
Summary
  './bin/bamtocov ../covtobed/example_data/panel_02.bam' ran
    1.37 ± 0.06 times faster than 'covtobed ../covtobed/example_data/panel_02.bam'
   52.12 ± 5.93 times faster than 'bedtools genomecov -bga -ibam ../covtobed/example_data/panel_02.bam'
```

### Comparison with mosdepth

```text
Benchmark #1: covtobed ../covtobed/example_data/panel_02.bam
  Time (mean ± σ):      1.088 s ±  0.028 s    [User: 1.034 s, System: 0.035 s]
  Range (min … max):    1.056 s …  1.125 s    6 runs
 
Benchmark #2: ./bin/bamtocov ../covtobed/example_data/panel_02.bam
  Time (mean ± σ):     789.6 ms ±  31.9 ms    [User: 1.088 s, System: 0.135 s]
  Range (min … max):   760.1 ms … 844.7 ms    6 runs
 
Benchmark #3: mosdepth /tmp/prefix ../covtobed/example_data/panel_02.bam
  Time (mean ± σ):     62.022 s ±  4.574 s    [User: 50.410 s, System: 11.236 s]
  Range (min … max):   57.778 s … 67.262 s    6 runs
 
Summary
  './bin/bamtocov ../covtobed/example_data/panel_02.bam' ran
    1.38 ± 0.07 times faster than 'covtobed ../covtobed/example_data/panel_02.bam'
   78.55 ± 6.61 times faster than 'mosdepth /tmp/prefix ../covtobed/example_data/panel_02.bam'
```

## Exome chr2 at 40% (300 Mb)

```bash
hyperfine "covtobed $FILE" "./bin/bamtocov $FILE" "mosdepth /tmp/prefix $FILE" "megadepth --coverage $FILE" --max-runs 6
```

```text
Benchmark #1: covtobed /local/giovanni/exome/HG00258.bam.04-chr2.bam
  Time (mean ± σ):     28.692 s ±  1.064 s    [User: 23.684 s, System: 4.780 s]
  Range (min … max):   27.852 s … 30.552 s    6 runs
 
Benchmark #2: ./bin/bamtocov /local/giovanni/exome/HG00258.bam.04-chr2.bam
  Time (mean ± σ):     18.735 s ±  0.360 s    [User: 17.998 s, System: 5.463 s]
  Range (min … max):   18.076 s … 19.072 s    6 runs
 
Benchmark #3: mosdepth /tmp/prefix /local/giovanni/exome/HG00258.bam.04-chr2.bam
  Time (mean ± σ):     19.486 s ±  1.535 s    [User: 15.571 s, System: 3.470 s]
  Range (min … max):   18.098 s … 22.342 s    6 runs
 
Benchmark #4: megadepth --coverage /local/giovanni/exome/HG00258.bam.04-chr2.bam
  Time (mean ± σ):      9.097 s ±  0.054 s    [User: 7.614 s, System: 1.379 s]
  Range (min … max):    9.025 s …  9.169 s    6 runs


  'megadepth --coverage /local/giovanni/exome/HG00258.bam.04-chr2.bam' ran
    2.06 ± 0.04 times faster than './bin/bamtocov /local/giovanni/exome/HG00258.bam.04-chr2.bam'
    2.14 ± 0.17 times faster than 'mosdepth /tmp/prefix /local/giovanni/exome/HG00258.bam.04-chr2.bam'
    3.15 ± 0.12 times faster than 'covtobed /local/giovanni/exome/HG00258.bam.04-chr2.bam'
```

## Exome chr1 full

## Exome chr22 full

```

Benchmark #1: covtobed /local/giovanni/exome/HG00258.chr21.bam
  Time (mean ± σ):     11.188 s ±  0.314 s    [User: 9.169 s, System: 1.867 s]
  Range (min … max):   10.637 s … 11.570 s    6 runs
 
Benchmark #2: ./bin/bamtocov /local/giovanni/exome/HG00258.chr21.bam
  Time (mean ± σ):      7.045 s ±  0.194 s    [User: 6.906 s, System: 2.092 s]
  Range (min … max):    6.830 s …  7.281 s    6 runs
 
Benchmark #3: mosdepth /tmp/prefix /local/giovanni/exome/HG00258.chr21.bam
  Time (mean ± σ):      7.541 s ±  0.654 s    [User: 6.536 s, System: 0.867 s]
  Range (min … max):    7.206 s …  8.873 s    6 runs
 
  Warning: The first benchmarking run for this command was significantly slower than the rest (8.873 s). This could be caused by (filesystem) caches that were not filled until after the first run. You should consider using the '--warmup' option to fill those caches before the actual benchmark. Alternatively, use the '--prepare' option to clear the caches before each timing run.
 
Benchmark #4: megadepth --coverage /local/giovanni/exome/HG00258.chr21.bam
  Time (mean ± σ):      3.199 s ±  0.018 s    [User: 2.849 s, System: 0.303 s]
  Range (min … max):    3.185 s …  3.235 s    6 runs
 
Summary
  'megadepth --coverage /local/giovanni/exome/HG00258.chr21.bam' ran
    2.20 ± 0.06 times faster than './bin/bamtocov /local/giovanni/exome/HG00258.chr21.bam'
    2.36 ± 0.20 times faster than 'mosdepth /tmp/prefix /local/giovanni/exome/HG00258.chr21.bam'
    3.50 ± 0.10 times faster than 'covtobed /local/giovanni/exome/HG00258.chr21.bam'
```