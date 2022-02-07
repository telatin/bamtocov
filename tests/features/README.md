# Compare featureCounts and bamtocounts

1) Generate a target in bed format (**target.bed**), like:
   
```text
seq1    1000    2000    feature1
seq1    5000    6000    feature2
seq1    10000   11000   feature3
```

2) Generate a list of reads mapping position in bed format (**seqs.bed**), like:

```text
seq1    1100    1200    +
seq1    1000    1200    +
seq1    900 1000        -
```

3) Run `make` to generate the bam file, the counts, and print a comparison. Example output:

```text
==> bamtocounts.tsv <==
feature1        3
feature2        7
feature3        2

==> fc.tsv <==
Geneid  seqs.bam
feature1        3
feature2        7
feature3        2
```