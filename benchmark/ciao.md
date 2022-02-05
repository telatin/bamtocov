| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `covtobed ../_test/cpara-ont-noseq.bam` | 689.3 ± 102.4 | 629.1 | 871.0 | 2.21 ± 0.37 |
| `bamtocov ../_test/cpara-ont-noseq.bam` | 312.0 ± 24.5 | 283.2 | 372.8 | 1.00 |
| `mosdepth _ciao ../_test/cpara-ont-noseq.bam` | 820.5 ± 37.1 | 780.1 | 853.6 | 2.63 ± 0.24 |
| `mosdepth --fast-mode _ciao ../_test/cpara-ont-noseq.bam` | 673.4 ± 31.5 | 637.7 | 724.4 | 2.16 ± 0.20 |
| `megadepth --coverage ../_test/cpara-ont-noseq.bam` | 588.7 ± 13.1 | 571.2 | 603.4 | 1.89 ± 0.15 |
| `megadepth --coverage --longreads ../_test/cpara-ont-noseq.bam` | 553.9 ± 15.7 | 540.4 | 579.3 | 1.78 ± 0.15 |
