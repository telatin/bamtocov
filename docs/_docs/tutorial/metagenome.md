---
title: Counting reads in a metagenome assembly
description: "Map reads to a metagenome assembly with minimap2 and summarise contig abundance with bamcountrefs"
date created: "2026-10-14"
last updated: "2026-10-14"
---

This tutorial shows how to go from raw reads and a metagenome assembly to
per-contig tables with [bamcountrefs]({% link _docs/commands/bamcountrefs.md %}).
You can produce either a simple count table or a set of metrics that are
useful in metagenomics, such as abundance, coverage depth and breadth.

## What you need

- A metagenome **assembly** in FASTA format, e.g. `assembly.fasta`. It can be
  the contigs from MEGAHIT or metaSPAdes, a set of MAGs merged into one file,
  or a reference catalogue.
- One or more **samples** of reads (Illumina paired-end in this example).
- [minimap2](https://github.com/lh3/minimap2), [samtools](https://www.htslib.org/)
  and BamToCov. All three are available from BioConda:

```bash
conda create -n metacov -c conda-forge -c bioconda minimap2 samtools bamtocov
conda activate metacov
```

Here we assume a layout like this:

```text
assembly.fasta
reads/
  sample1_R1.fastq.gz  sample1_R2.fastq.gz
  sample2_R1.fastq.gz  sample2_R2.fastq.gz
  ...
```

## 1. Index the assembly

minimap2 can map against the FASTA directly, but building the index once
saves time when you map many samples:

```bash
minimap2 -x sr -d assembly.mmi assembly.fasta
```

The preset (`-x sr`) is stored in the index, so it must match the preset used
when mapping.

## 2. Map the reads and make sorted BAM files

Map each sample, then sort and index the output. `bamcountrefs` needs
**coordinate-sorted, indexed** BAM files, and all the samples must be mapped
against **the same assembly**.

```bash
mkdir -p bam

for R1 in reads/*_R1.fastq.gz; do
  R2=${R1/_R1/_R2}
  SAMPLE=$(basename "$R1" _R1.fastq.gz)

  minimap2 -t 8 -a -x sr assembly.mmi "$R1" "$R2" \
    | samtools sort -@ 4 -o bam/${SAMPLE}.bam -

  samtools index bam/${SAMPLE}.bam
done
```

Which preset to use depends on the reads:

| Reads | minimap2 preset |
|---|---|
| Illumina short reads | `-x sr` |
| Oxford Nanopore | `-x map-ont` |
| PacBio HiFi | `-x map-hifi` |

For long reads, map a single FASTQ file per sample instead of a pair.

:::tip
Each sample's column in the output is named after its BAM file, up to the
first dot. `bam/sample1.bam` becomes `sample1`, so give each BAM file a unique,
meaningful name.
:::

## 3a. Just a count table

The simplest output is the number of reads mapped to each contig in each
sample, printed to standard output:

```bash
bamcountrefs bam/*.bam > counts.tsv
```

```text
Contig      sample1  sample2  sample3
contig_1    10234    8812     15020
contig_2    0        12       3
contig_3    5521     4980     6107
...
```

This is a good input for differential abundance tools that expect raw counts,
such as DESeq2 or edgeR.

You can change the first column header with `--tag` (e.g. `--tag MAG` if each
reference is a genome bin) and speed things up with threads:

```bash
bamcountrefs --tag Contig -T 4 bam/*.bam > counts.tsv
```

### Which reads are counted

By default, unmapped, secondary, supplementary, QC-failed and duplicate
alignments are skipped (`-F 3844`), so each read is counted at most once. You
can filter further:

```bash
# keep only confidently placed reads, and only proper pairs
bamcountrefs -Q 10 -P bam/*.bam > counts.tsv
```

:::caution
Metagenome assemblies often contain closely related strains or repeated
elements. A read that maps equally well to several contigs gets MAPQ 0 and is
assigned to one of them at random. A MAPQ filter (`-Q`) removes these reads,
which makes the counts more specific but lowers them for contigs that share
sequence.
:::

## 3b. Metrics for metagenomics

Raw counts don't account for contig length or sequencing depth, and they don't
tell you whether a contig is truly present or just collecting reads in one
small region (e.g. a conserved gene or a mobile element). To get more metrics,
pass an output basename with `-o`. `bamcountrefs` then writes one table per
metric:

```bash
mkdir -p results
bamcountrefs -o results/meta -T 4 \
  --length --tpm --mean --trimmed-mean --covered-ratio --discov \
  bam/*.bam
```

This creates:

| File | Metric | Use in metagenomics |
|---|---|---|
| `results/meta_counts.tsv` | raw read counts | input for differential abundance |
| `results/meta_length.tsv` | contig length | filter out short contigs, normalise |
| `results/meta_tpm.tsv` | TPM | relative abundance, comparable across samples |
| `results/meta_mean.tsv` | mean depth | contig coverage, e.g. for binning |
| `results/meta_trimmed_mean.tsv` | trimmed mean depth | depth that is robust to spikes from repeats |
| `results/meta_covered_fraction.tsv` | breadth of coverage | presence/absence of a contig or genome |
| `results/meta_discov.tsv` | [DisCov]({% link _docs/discov.md %}) score | how evenly reads are spread along the contig |

Or ask for everything at once:

```bash
bamcountrefs -o results/meta -T 4 --all-metrics bam/*.bam
```

### Reading the metrics together

- **Abundance:** TPM corrects for both contig length and sample depth, so it's
  the best choice for comparing relative abundance across samples. RPKM
  (`--rpkm`) is also available.
- **Presence:** a high count with a **low breadth** (covered fraction) usually
  means reads are piling up on a shared region, not that the organism is
  present. A common rule of thumb is to treat a genome as present only above a
  breadth threshold (e.g. 10–50%, depending on depth).
- **Depth:** the trimmed mean removes the lowest and highest coverage
  positions (5% each by default, see `--trim-min`/`--trim-max`). Like the
  `trimmed_mean` in CoverM, it's a more stable estimate of a
  genome's depth than the plain mean.
- **Evenness:** a low DisCov score flags patchy coverage, which can point to
  misassemblies, chimeric bins, or a related strain that only shares part of
  the sequence.

:::warning
Breadth, trimmed mean, variance and DisCov track coverage at every base. They
need memory proportional to the total length of the contigs that receive
reads. With large assemblies, run fewer samples in parallel (`-W`) if memory
is tight. Counts, TPM, RPKM, mean and length need no extra memory.
:::

### MultiQC report

Add `--multiqc` to also write `results/meta.mqc.txt`, which
[MultiQC](https://multiqc.info) can include in its report:

```bash
bamcountrefs -o results/meta --tpm --covered-ratio --multiqc bam/*.bam
```

## Summary

```bash
# index once
minimap2 -x sr -d assembly.mmi assembly.fasta

# map, sort and index each sample
minimap2 -t 8 -a -x sr assembly.mmi s1_R1.fq.gz s1_R2.fq.gz \
  | samtools sort -o bam/s1.bam - && samtools index bam/s1.bam

# 1) count table
bamcountrefs bam/*.bam > counts.tsv

# 2) metagenomics metrics
bamcountrefs -o results/meta --length --tpm --trimmed-mean --covered-ratio --discov bam/*.bam
```
