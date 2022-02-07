| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `../../bin/bamtocounts --id gene_id --type exon target.gtf seqs.bam > bamtocounts.tsv` | 5.9 ± 1.6 | 4.7 | 19.1 | 1.00 |
| `featureCounts -a target.gtf -o fc seqs.bam` | 5.9 ± 1.3 | 4.9 | 18.0 | 1.01 ± 0.35 |
