# Automatic benchmark

This directory contains a Dockerfile to automatically benchmark bamtocov, covtobed,
mosdepth and megadepth using a set of BAM files.

The image should be build from the root directory.


* [Example log](log.txt)

## Output

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `covtobed 'HG00258.bam' #1.3.5` | 563.326 ± 8.611 | 555.070 | 574.771 | 3.25 ± 0.05 |
| `bamtocov 'HG00258.bam' #2.5.0` | 363.811 ± 8.651 | 350.075 | 373.599 | 2.10 ± 0.05 |
| `megadepth --coverage 'HG00258.bam' #1.1.2` | 173.520 ± 0.228 | 173.290 | 173.810 | 1.00 ± 0.00 |
| `megadepth --coverage --longreads 'HG00258.bam' #1.1.2` | 173.281 ± 0.478 | 172.890 | 173.919 | 1.00 |
| `mosdepth mosdepth.pMD4Cm 'HG00258.bam' #0.3.3` | 469.604 ± 0.839 | 468.725 | 471.000 | 2.71 ± 0.01 |
| `mosdepth --fast-mode mosdepth.pMD4Cm 'HG00258.bam' #0.3.3` | 459.380 ± 1.134 | 458.283 | 461.134 | 2.65 ± 0.01 |

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `covtobed 'cpara-illumina-noseq.bam' #1.3.5` | 14.553 ± 0.246 | 14.289 | 14.895 | 3.70 ± 0.10 |
| `bamtocov 'cpara-illumina-noseq.bam' #2.5.0` | 10.774 ± 0.144 | 10.550 | 10.939 | 2.74 ± 0.07 |
| `megadepth --coverage 'cpara-illumina-noseq.bam' #1.1.2` | 3.932 ± 0.077 | 3.842 | 4.030 | 1.00 |
| `megadepth --coverage --longreads 'cpara-illumina-noseq.bam' #1.1.2` | 4.129 ± 0.112 | 4.012 | 4.313 | 1.05 ± 0.04 |
| `mosdepth mosdepth.X0pr36 'cpara-illumina-noseq.bam' #0.3.3` | 11.090 ± 0.135 | 10.908 | 11.289 | 2.82 ± 0.07 |
| `mosdepth --fast-mode mosdepth.X0pr36 'cpara-illumina-noseq.bam' #0.3.3` | 9.387 ± 0.074 | 9.291 | 9.468 | 2.39 ± 0.05 |

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `covtobed 'cpara-ont-noseq.bam' #1.3.5` | 494.1 ± 14.3 | 479.1 | 512.8 | 2.09 ± 0.11 |
| `bamtocov 'cpara-ont-noseq.bam' #2.5.0` | 236.3 ± 10.7 | 224.2 | 253.5 | 1.00 |
| `megadepth --coverage 'cpara-ont-noseq.bam' #1.1.2` | 400.3 ± 11.5 | 387.4 | 419.2 | 1.69 ± 0.09 |
| `megadepth --coverage --longreads 'cpara-ont-noseq.bam' #1.1.2` | 406.4 ± 11.4 | 392.1 | 424.8 | 1.72 ± 0.09 |
| `mosdepth mosdepth.XMG5Wi 'cpara-ont-noseq.bam' #0.3.3` | 664.9 ± 29.3 | 629.0 | 687.7 | 2.81 ± 0.18 |
| `mosdepth --fast-mode mosdepth.XMG5Wi 'cpara-ont-noseq.bam' #0.3.3` | 542.8 ± 10.0 | 535.5 | 560.0 | 2.30 ± 0.11 |

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `covtobed 'panel_01.bam' #1.3.5` | 422.6 ± 3.7 | 417.5 | 427.5 | 2.21 ± 0.06 |
| `bamtocov 'panel_01.bam' #2.5.0` | 191.3 ± 4.9 | 185.7 | 204.4 | 1.00 |
| `megadepth --coverage 'panel_01.bam' #1.1.2` | 6068.9 ± 36.9 | 6021.0 | 6121.8 | 31.72 ± 0.83 |
| `megadepth --coverage --longreads 'panel_01.bam' #1.1.2` | 6050.0 ± 19.7 | 6023.0 | 6078.5 | 31.63 ± 0.81 |
| `mosdepth mosdepth.zCQZhW 'panel_01.bam' #0.3.3` | 42686.6 ± 7201.3 | 35745.8 | 51457.9 | 223.14 ± 38.07 |
| `mosdepth --fast-mode mosdepth.zCQZhW 'panel_01.bam' #0.3.3` | 25916.5 ± 166.2 | 25730.0 | 26144.3 | 135.47 ± 3.55 |

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `covtobed 'w1118_f_hd_R1.junc.bam' #1.3.5` | 1.530 ± 0.021 | 1.502 | 1.556 | 2.50 ± 0.06 |
| `bamtocov 'w1118_f_hd_R1.junc.bam' #2.5.0` | 1.432 ± 0.074 | 1.375 | 1.523 | 2.34 ± 0.13 |
| `megadepth --coverage 'w1118_f_hd_R1.junc.bam' #1.1.2` | 0.613 ± 0.013 | 0.604 | 0.634 | 1.00 |
| `megadepth --coverage --longreads 'w1118_f_hd_R1.junc.bam' #1.1.2` | 0.614 ± 0.003 | 0.609 | 0.618 | 1.00 ± 0.02 |
| `mosdepth mosdepth.7dbdyA 'w1118_f_hd_R1.junc.bam' #0.3.3` | 2.939 ± 0.017 | 2.913 | 2.958 | 4.80 ± 0.10 |
| 
`mosdepth --fast-mode mosdepth.7dbdyA 'w1118_f_hd_R1.junc.bam' #0.3.3` | 2.907 ± 0.030 | 2.881 | 2.953 | 4.74 ± 0.11 |
