LS="
.rw-r--r--@  296 telatina 12 Mar 15:51 sim_counts.tsv
.rw-r--r--@  378 telatina 12 Mar 15:51 sim_covered_bases.tsv
.rw-r--r--@  534 telatina 12 Mar 15:51 sim_covered_fraction.tsv
.rw-r--r--@  422 telatina 12 Mar 15:51 sim_length.tsv
.rw-r--r--@  538 telatina 12 Mar 15:51 sim_mean.tsv
.rw-r--r--@  534 telatina 12 Mar 15:51 sim_reads_per_base.tsv
.rw-r--r--@  606 telatina 12 Mar 15:51 sim_rpkm.tsv
.rw-r--r--@  698 telatina 12 Mar 15:51 sim_tpm.tsv
.rw-r--r--@  538 telatina 12 Mar 15:51 sim_trimmed_mean.tsv
.rw-r--r--@  539 telatina 12 Mar 15:51 sim_variance.tsv
"
hyperfine -w 2 -m 13 \
  "coverm contig -m tpm rpkm  length count mean variance covered_bases covered_fraction trimmed_mean reads_per_base -b input/sim-chr/*bam > tmp/coverm.txt" \
  "bin/bamcountrefs -a -W 8 -P input/sim-chr/*bam -o tmp/sim"
