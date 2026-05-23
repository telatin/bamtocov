echo "Benchmark 1: featureCounts"
if [[ -e input/sim_out/sample1.bam ]]; then
  echo "Run:"
  hyperfine -w 1 -m 5 \
  "featureCounts -F GTF -p -a input/genome_10.gtf -o input/sim_out/fc_counts.txt \
    input/sim_out/sample1.bam input/sim_out/sample2.bam input/sim_out/sample3.bam input/sim_out/sample4.bam \
    -t gene \
    --primary --ignoreDup -Q 1 --countReadPairs >/dev/null 2>/dev/null" \
  "bin/bamtocounts input/genome_10.gtf input/sim_out/*bam --type gene --id gene_id \
    --paired \
    > input/sim_out/btc_counts.txt"
else
  echo Missing input: input/sim_out/sample1.bam
fi

echo "Benchmark 2: featureCounts?"
if [[ -e input/sim-chr/sample1.bam ]];
then
  echo "Run:"
  hyperfine -w 1 -m 5 \
   "featureCounts -a input/genome_10.gtf  input/sim-chr/*bam -o tmp/fc -t gene" \
   "bin/bamtocounts input/genome_10.gtf  input/sim-chr/*bam  --type gene --id gene_id"
else
  echo Missing input: input/sim-chr/sample1.bam
fi


hyperfine -w 1 -m 5 \
  "coverm contig -m tpm rpkm count mean variance covered_bases covered_fraction trimmed_mean reads_per_base -b input/sim-chr/*bam > tmp/coverm.txt" \
  "bin/bamcountrefs -a -W 8 -P input/sim-chr/*bam -o tmp/sim"

