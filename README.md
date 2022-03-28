# bamtocov

[![Build Bamtocov](https://github.com/telatin/bamtocov/actions/workflows/build.yml/badge.svg)](https://github.com/telatin/bamtocov/actions/workflows/build.yml)
[![Downloads](https://img.shields.io/conda/dn/bioconda/bamtocov)](https://anaconda.org/bioconda/bamtocov)
[![Platforms](https://anaconda.org/bioconda/bamtocov/badges/platforms.svg)](https://bioconda.github.io/recipes/bamtocov/README.html)
[![Paper](https://img.shields.io/badge/paper-bioinformatics-yellow)](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac125/6535233)

[![bamtocov logo](docs/bamtocov-banner.png)](https://telatin.github.io/bamtocov/)

Tools to extract coverage informations from BAM (and CRAM) files, based on the
[covtobed](https://github.com/telatin/covtobed) algorithm that supports 
_stranded coverage_ and _physical coverage_, input from _streams_ 
and uses a _memory-efficient algorithm_. 

## :book: Documentation

Full documentation is available online at the :book: **[ dedicated website](https://telatin.github.io/bamtocov/)**, or in
this repository under `docs`.

## Installation

The BamToCov package is available from :package: [BioConda](https://bioconda.github.io/recipes/bamtocov/README.html)

```bash
conda install -y -c bioconda bamtocov
```

## Benchmarks

Bamtocov has the smallest memory footprint of any other coverage tool, while maintaining reasonable speeds

* Four BAM files used for the benchmarks are [available from Zenodo](https://zenodo.org/record/5636944#.Yf_cMe7P36Y)
* A [Dockerfile](benchmark/) with an automatic speed test is available in this repository
  
## See also

* Giovanni Birolo, Andrea Telatin [BamToCov: an efficient toolkit for sequence coverage calculations](https://doi.org/10.1101/2021.11.12.466787) (BioRχiv)
* Brent Pedersen,  Aaron Quinlan,
[hts-nim: scripting high-performance genomic analyses](https://academic.oup.com/bioinformatics/article/34/19/3387/4990493) (Bioinformatics)
* Giovanni Birolo, Andrea Telatin,
[covtobed: a simple and fast tool to extract coverage tracks from BAM files](https://joss.theoj.org/papers/10.21105/joss.02119) (JOSS)

Initially we also used [Lapper](https://brentp.github.io/nim-lapper/index.html), and I recommend checking out this library for fast interval operations.

## Cite

Giovanni Birolo, Andrea Telatin, 
[BamToCov: an efficient toolkit for sequence coverage calculations](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac125/6535233), Bioinformatics, 2022
