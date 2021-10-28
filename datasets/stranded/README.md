# Stranded coverage

## Full stranded: `stranded.bam`

### Description

_stranded.bam_ contains two reads, one of which is reverse-complemented.
Together, they cover the whole 250bp genome.

### Generation

```text
seqfu cat --fastq --truncate 125 --prefix for --strip-name ref-250.fa > for.fq
seqfu rc  ref-250.fa  | seqfu cat --fastq --truncate 125  --prefix rev --strip-name > rev.fq
cat for.fq rev.fq > stranded.fq && rm {for,rev}.fq
bwa mem -t 4 ref-250.fa stranded.fq | samtools view -bS - > stranded.bam
```

this generates the following SAM:

```text
for_1  0       filt.1  1       60      125M    *       0       0       CTTGGTCATTTA...GCTTTGGTAGG  BBBBBBBBB...BBBBBBBB   NM:i:0  MD:Z:125     AS:i:125 XS:i:0
rev_1  16      filt.1  126     60      124M1S  *       0       0       CCTTCTATATGG...CATCGATGAAGA  BBBBBBB...BBBBBBBBBB   NM:i:0  MD:Z:124     AS:i:124 XS:i:0
```

### Output

* Sequence coverage: `bamtocov stranded.bam`

```text
filt.1  0       249     1
```

* Stranded: `bamtocov --stranded stranded.bam`

```text
filt.1  0       125     1       0
filt.1  125     249     0       1
```

## Small stranded with gap: `gap-stranded.bam`

### Description

_gap-stranded.bam_ contains two reads, one of which is reverse-complemented.
Together, they physically cover the whole genome, but having a gap between
them their sequence coverage is smaller.

### Generation

```text
seqfu cat --fastq --truncate 100 --prefix for --strip-name ref-250.fa > for.fq
seqfu rc  ref-250.fa  | seqfu cat --fastq --truncate 100  --prefix rev --strip-name > rev.fq
cat for.fq rev.fq > gap-stranded.fq && rm {for,rev}.fq
bwa mem -t 4 ref-250.fa gap-stranded.fq | samtools view -bS - > gap-stranded.bam
```

### Output

* Sequence coverage: `bamtocov gap-stranded.bam`

```text
filt.1  0       100     1
filt.1  100     150     0
filt.1  150     250     1
```

* Stranded: `bamtocov --stranded gap-stranded.bam`

```text
filt.1  0       100     1       0
filt.1  100     150     0       0
filt.1  150     250     0       1
```


## Small stranded with overlap: `overlap-stranded.bam`

### Description

_overlap-stranded.bam_ contains two reads, one of which is reverse-complemented.
Together, they physically cover the whole genome, but having a gap between
them their sequence coverage is smaller.

### Generation

```text
seqfu cat --fastq --truncate 200 --prefix for --strip-name ref-250.fa > for.fq
seqfu rc  ref-250.fa  | seqfu cat --fastq --truncate 200  --prefix rev --strip-name > rev.fq
cat for.fq rev.fq > overlap-stranded.fq && rm {for,rev}.fq
bwa mem -t 4 ref-250.fa overlap-stranded.fq | samtools view -bS - > overlap-stranded.bam
```

### Output

* Sequence coverage: `bamtocov overlap-stranded.bam`

```text
filt.1  0       50      1
filt.1  50      200     2
filt.1  200     250     1
```

* Stranded: `bamtocov --stranded overlap-stranded.bam`

```text
filt.1  0       50      1       0
filt.1  50      200     1       1
filt.1  200     250     0       1
```
