# BamToCov

![CovToBam logo](bamtocov.png)

- :book: Documentation: <https://telatin.github.io/bamtocov/>
- :package: Github: <https://github.com/telatin/bamtocov>

**BamToCov** is a suite of tools for coverage analysis written in Nim and built upon the
memory efficient algorithm of [**covtobed**](https://github.com/telatin/covtobed).

The program uses [_htslib_](https://www.htslib.org) to parse BAM/CRAM files, and specifically the Nim-wrapper
[_hts-nim_](https://github.com/brentp/hts-nim).

We designed **BamToCov** to fill some gaps in coverage analysis:

1. Accepting input streams (i.e. not requiring the BAM index, at the expense of speed)
1. Enabling _per strand_ coverage analysis
1. Enabling the _physical coverage_ analysis
1. Using a minimum amount of memory (minimum memory usage)

This makes **BamToCov** a useful companion especially when testing pipelines on small datasets.

For coverage analysis of large datasets it can be useful to consider
[Mosdepth](https://github.com/brentp/mosdepth), that uses indexed
BAM/CRAM files as input. Another interesting alternative is
[Megadepth](https://github.com/ChristopherWilks/megadepth), a fast and memory
efficient tool.