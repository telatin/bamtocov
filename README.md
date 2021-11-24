# bamtocov

[![Build Bamtocov](https://github.com/telatin/bamtocov/actions/workflows/build.yml/badge.svg)](https://github.com/telatin/bamtocov/actions/workflows/build.yml)
[![Downloads](https://img.shields.io/conda/dn/bioconda/bamtocov)](https://anaconda.org/bioconda/bamtocov)
[![Platforms](https://anaconda.org/bioconda/bamtocov/badges/platforms.svg)](https://bioconda.github.io/recipes/bamtocov/README.html)

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

## References

* Giovanni Birolo, Andrea Telatin [BamToCov: an efficient toolkit for sequence coverage calculations](https://doi.org/10.1101/2021.11.12.466787) (BioRχiv)
* Brent Pedersen,  Aaron Quinlan,
[hts-nim: scripting high-performance genomic analyses](https://academic.oup.com/bioinformatics/article/34/19/3387/4990493) (Bioinformatics)
* Giovanni Birolo, Andrea Telatin,
[covtobed: a simple and fast tool to extract coverage tracks from BAM files](https://joss.theoj.org/papers/10.21105/joss.02119) (JOSS)
* Brent Pedersen, [nim-lapper: fast, simple interval overlapping (Version 0.1.7)](https://brentp.github.io/nim-lapper/index.html)
