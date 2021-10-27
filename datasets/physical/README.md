# Physical coverage

## Single fragment: `paired.bam`

### Description

_paired.bam_ contains a single fragment, composed by two paired-end reads covering 
the beginning (R1) and the end (R2) of the reference sequence (250 bp).
Each fragment is 74 bp long.

### Generation

```
seqfu cat reference.fa --fastq --truncate 47 > file_R1.fq
seqfu rc reference.fa  | seqfu cat --fastq --truncate 47 > file_R2.fq
bwa mem -t 4 reference.fa file_R1.fq file_R2.fq | samtools view -bS - > paired.bam
```

### Output
* Sequence coverage: `bamtocov paired.bam`

```text
filt.1  0       46      1
filt.1  46      203     0
filt.1  203     249     1
```

* Physical: `bamtocov --physical paired.bam`

```text
filt.1  0       249     1
```



