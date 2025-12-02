set -euxo pipefail

coverm contig   \
  -m mean count rpkm tpm covered_bases length variance trimmed_mean reads_per_base \
  -b input/metagenome/*bam > temp/coverm.raw
votuderep splitcoverm -i temp/coverm.raw -o temp/coverm

bin/bamcountrefs --all-metrics -o temp/this input/metagenome/*bam 
