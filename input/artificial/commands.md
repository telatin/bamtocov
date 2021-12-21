# Commands

```bash
# Make BAM from list of reads
../../scripts/benchmarking/make-bam-from-bed.py  -i reads.bed -o reads.bam

# Make FASTA and GTF from BED
../../scripts/benchmarking/bed-to-gtf.py -i annotation.bed -o annotation.gtf --min-len 500
../../scripts/benchmarking/bed-to-fasta.py -i annotation.bed -o reference.fasta --min-len 500

# Indexes
samtools faidx reference.fasta
samtools index reads.bam
```

## Tests

```bash
# Generate featureCounts output
featureCounts -a annotation.gtf -o featureCounts reads.bam

# BamToCounts output 
../../bin/bamtocounts annotation.bed reads.bam > counts.tsv
```
