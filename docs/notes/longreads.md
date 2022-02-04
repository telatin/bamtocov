# Testing on long reads


## Benchmark

The speed test is performed via hyperfine: 

```bash
hyperfine \
  "bin/bamtocov $FILE" \
  "megadepth --coverage --longreads $FILE" \
  "mosdepth -x /tmp/x $FILE" \
  --max-runs 3 --export-markdown $FILE.md
```

## Datasets 

### Fungal dataset

The alignment of _whole genome shotgun_ reads coming from a unicellular fungus sequences with MinION
in our laboratory is available in the [Zenodo dataset](https://zenodo.org/record/5636944#.Yf1AQPXP36M)
accompanying the paper.

```bash
wget -O cpara-ont.bam "https://zenodo.org/record/5636944/files/cpara-ont-noseq.bam?download=1"
``
### Human sample

[SRR13615770](https://www.ncbi.nlm.nih.gov/sra/?term=SRR13615770), from the study
_Cas9 targeted enrichment of mobile elements using nanopore sequencing_ is composed by
123,254 reads.

```bash
# Download reference
wget -O db/hg19.fa.gz http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.masked.gz
# Download sample
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR136/070/SRR13615770/SRR13615770_1.fastq.gz \
  -o SRR13615770_L1Hs_MinION_1.fastq.gz

# Select reads > 2500bp
seqfu cat -m 2500 SRR13615770_L1Hs_MinION_1.fastq.gz | gzip -c > SRR13615770_L1Hs_MinION_1.min2500.fastq.gz

# Map 
minimap2 -ax map-ont -t 4 db/hg19.fa.gz SRR13615770_L1Hs_MinION_1.min2500.fastq.gz | \
  samtools view -bS | samtools sort -o SRR13615770.bam -
```


## Benchmark

```bash
hyperfine \
  "bin/bamtocov $FILE" \
  "megadepth --coverage --longreads $FILE" \
  "mosdepth -x /tmp/x $FILE" \
  --max-runs 2 --export-markdown $FILE.md
```

## Results

* **SRR13615770_L1Hs_MinION_1** (full sample)

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `bin/bamtocov SRR13615770_L1Hs_MinION_1.bam` | 32.357 ± 0.261 | 32.172 | 32.541 | 1.00 |
| `megadepth --coverage --longreads SRR13615770_L1Hs_MinION_1.bam` | 46.115 ± 0.954 | 45.440 | 46.789 | 1.43 ± 0.03 |
| `mosdepth -x /tmp/x SRR13615770_L1Hs_MinION_1.bam` | 92.359 ± 5.448 | 88.506 | 96.211 | 2.85 ± 0.17 |

* **SRR13615770_2500bp** (size selection: > 2,500 bp)

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `bin/bamtocov SRR13615770.bam` | 29.113 ± 0.121 | 28.981 | 29.219 | 1.00 |
| `megadepth --coverage --longreads SRR13615770.bam` | 43.366 ± 0.327 | 42.994 | 43.607 | 1.49 ± 0.01 |
| `mosdepth -x /tmp/x SRR13615770.bam` | 86.483 ± 2.312 | 83.878 | 88.293 | 2.97 ± 0.08 |

* **Fungus LR**

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `bin/bamtocov _test/cpara-ont-noseq.bam` | 310.5 ± 23.0 | 275.5 | 339.0 | 1.00 |
| `megadepth --coverage --longreads _test/cpara-ont-noseq.bam` | 557.4 ± 13.5 | 540.5 | 572.7 | 1.80 ± 0.14 |
| `mosdepth -x /tmp/x _test/cpara-ont-noseq.bam` | 792.5 ± 176.9 | 654.5 | 1092.7 | 2.55 ± 0.60 |