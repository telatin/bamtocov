# Physical coverage

## Single fragment: `paired.bam`

### Description

_paired.bam_ contains a single fragment, composed by two paired-end reads covering 
the beginning (R1) and the end (R2) of the reference sequence (250 bp).
Each fragment is 74 bp long.

### Generation

```text
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


## Single fragment: `overlap-paired.bam`

### Description

_overlap-paired.bam_ contains a single fragment, composed by two paired-end reads covering 
the beginning (R1) and the end (R2) of the reference sequence (250 bp).
Each fragment is 300 bp long.

### Generation

```text
seqfu cat ref-250.fa --fastq --truncate 200 > overlap_R1.fq
seqfu rc ref-250.fa  | seqfu cat --fastq --truncate 200 > overlap_R2.fq
bwa mem -t 4 ref-250.fa overlap_R1.fq overlap_R2.fq | samtools view -bS - > overlap-paired.bam
```

### Output

* Sequence coverage: `bamtocov overlap-paired.bam`

```text
filt.1  0       50      1
filt.1  50      200     2
filt.1  200     250     1
```

* Physical: `bamtocov --physical overlap-paired.bam`

```text
filt.1  0       250     1
```
