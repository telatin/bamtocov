# Compare featureCounts and bamtocounts

1) Generate a target in bed format, like:
   
```text
seq1    1000    2000    feature1
seq1    5000    6000    feature2
seq1    10000   11000   feature3
```

2) Generate a list of reads mapping position in bed format, like:

```text
seq1    1100    1200    +
seq1    1000    1200    +
seq1    900 1000        -
```

3) Run `make` to generate the bam file, the counts, and print a comparison