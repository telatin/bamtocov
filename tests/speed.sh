FILE=$1

if [[ ! -e "$FILE.bai" ]];
then
  echo "Index missing"; exit
fi

hyperfine --max-runs 4 --export-markdown $FILE.md --export-csv $FILE.csv  \
 "/local/miniconda3/envs/mega/bin/megadepth --coverage $FILE"  \
 "bin/bamtocov $FILE"  \
 "covtobed $FILE"  \
 " mosdepth -x /tmp/prefix $FILE"  