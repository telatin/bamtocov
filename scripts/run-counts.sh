INPUT=input/genome_10.gtf
OUT=input/sim_out/
mkdir -p "$OUT"
python scripts/sim-counts.py -t "$INPUT" -n 4 -o "$OUT" \
  --reads-per-feature 100 --read-length 75 \
  --suppl-ratio 0.05 \
  --secondary-ratio 0.03 \
  --duplicate-ratio 0.04 \
  --unmapped-ratio 0.02 \
  --edge-ratio 0.15 \
  --low-mapq-ratio 0.05 \
  --paired \
  --seed 42

featureCounts -F GTF -p -a input/genome_10.gtf -o input/sim_out/fc_counts.txt \
  input/sim_out/sample1.bam input/sim_out/sample2.bam input/sim_out/sample3.bam input/sim_out/sample4.bam \
  -t gene \
  --primary --ignoreDup -Q 1 --countReadPairs >/dev/null 2>/dev/null

bin/bamtocounts input/genome_10.gtf input/sim_out/*bam --type gene --id gene_id \
  --paired \
  > input/sim_out/btc_counts.txt

FC=$(head -n 44 input/sim_out/fc_counts.txt | cut -f 1,7- | grep ^gene  | md5sum)
BTC=$(cat input/sim_out/btc_counts.txt | grep gene | md5sum)

if [[ "$FC" == "$BTC" ]]; then
  echo OK
else
  echo FAIL
  echo FC:
  head -n 44 input/sim_out/fc_counts.txt | cut -f 1,7- | grep ^gene | head -n 3
  echo ---
  echo BTC:
  cat input/sim_out/btc_counts.txt | grep gene | head -n 3
fi
