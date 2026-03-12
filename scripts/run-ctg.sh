set -exuo pipefail
OUT=input/sim-chr/
mkdir -p "$OUT"

python scripts/sim-chrs.py \
  -c 1222 -n 4 -o "$OUT" \
  --mean-reads 300 --read-length 100 \
  --suppl-ratio 0.05 \
  --secondary-ratio 0.03 \
  --duplicate-ratio 0.04 \
  --unmapped-ratio 0.02 \
  --low-mapq-ratio 0.05 \
  --seed 42

BAMS=("$OUT"*.bam)
echo "bams: ${BAMS[*]}"
echo "------ Run Coverm"
for M in count mean rpkm tpm covered_bases covered_fraction length variance reads_per_base;
do
coverm contig -m $M \
  -b "${BAMS[@]}" \
  > "$OUT"/coverm_${M}.txt
done
echo "------ Run BamCountRefs"
bin/bamcountrefs --all-metrics -o "${OUT}btc" "${BAMS[@]}"

# Compare: read counts for all contigs, first sample only
# coverm output: "Contig \t sample1.bam Count \t ..."  (col 2)
# bamcountrefs:  "Contig \t sample1 \t ..."            (col 2)
COVERM=$(tail -n +2 "$OUT"coverm_counts.txt | sort -k1 | cut -f1,2)
BTC=$(tail -n +2 "$OUT"btc_counts.tsv     | sort -k1 | cut -f1,2)

COVERM_MD5=$(echo "$COVERM" | md5sum)
BTC_MD5=$(echo "$BTC"     | md5sum)

if [[ "$COVERM_MD5" == "$BTC_MD5" ]]; then
  echo "OK — coverm and bamcountrefs counts match"
else
  echo "DIFF — counts differ (expected: flag filtering may vary)"
  echo "  coverm    : $COVERM_MD5"
  echo "  bamcountrefs: $BTC_MD5"
  # Show first few discrepancies
  diff <(echo "$COVERM") <(echo "$BTC") | head -20
fi
