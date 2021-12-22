# Tiny target dataset

![Screenshot of IGV](igv.png)

## Source files

* `target.bed` is the original target
* `reads.bed` is the list of reads

## Produced files

* `ref.fa` is a random reference (generated for IGV)
* `reads.bam` is the alignment file
* `target.gtf` is the target for featureCounts
* `fc` is FeatureCounts output (counts)
* `mos-target.regions.bed` is Mosdepth output (coverage)

A **Makefile** produces all file but mosdepth
