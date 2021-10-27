# Stranded coverage

## Single pair: `stranded.bam`

### Description

_stranded.bam_ contains two reads, one of which is reverse-complemented.
Together, they cover the whole 250bp genome.

### Generation

```
seqfu cat --fastq --truncate 125 --prefix for --strip-name ref-250.fa > for.fq
seqfu rc  ref-250.fa  | seqfu cat --fastq --truncate 125  --prefix rev --strip-name > rev.fq
cat for.fq rev.fq > stranded.fq && rm {for,rev}.fq
bwa mem -t 4 ref-250.fa stranded.fq | samtools view -bS - > stranded.bam
```

this generates the following SAM:
```
for_1  0       filt.1  1       60      125M    *       0       0       CTTGGTCATTTA...GCTTTGGTAGG  BBBBBBBBB...BBBBBBBB   NM:i:0  MD:Z:125     AS:i:125 XS:i:0
rev_1  16      filt.1  126     60      124M1S  *       0       0       CCTTCTATATGG...CATCGATGAAGA  BBBBBBB...BBBBBBBBBB   NM:i:0  MD:Z:124     AS:i:124 XS:i:0
```

### Output
* Sequence coverage: `bamtocov stranded.bam`

```text
filt.1  0       249     1
```

* Physical: `bamtocov --stranded stranded.bam`

```text
filt.1  0       125     1       0
filt.1  125     249     0       1
```



