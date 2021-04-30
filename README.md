# bamtocov

Tools to extract coverage informations from BAM (and CRAM) files, based on the
[covtobed](https://github.com/telatin/covtobed) algorithm that supports stranded
coverage and physical coverage, input from streams 
and uses a memory-efficient algorithm.

## :book: Documentation

Full documentation is available online at the **[dedicated website](https://telatin.github.io/bamtocov/)**, or in
this repository under `docs`.

## Installation

The BamToCov package is available from [BioConda](https://bioconda.github.io/recipes/bamtocov/README.html)

```bash
conda install -y -c bioconda bamtocov
```

## References
 * Brent Pedersen,  Aaron Quinlan, 
 [hts-nim: scripting high-performance genomic analyses](https://academic.oup.com/bioinformatics/article/34/19/3387/4990493) (Bioinformatics)
 * Giovanni Birolo, Andrea Telatin, 
 [covtobed: a simple and fast tool to extract coverage tracks from BAM files](https://joss.theoj.org/papers/10.21105/joss.02119) (JOSS)
