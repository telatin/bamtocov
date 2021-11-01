FILE=$1

if [[ ! -e "$FILE.bai" ]];
then
  echo "Index missing"; exit
fi
echo megadepth
memusg /local/miniconda3/envs/mega/bin/megadepth --coverage $1 > /dev/null

sleep 3; echo ""

echo bamtocov

echo covtobed
memusg covtobed $1 > /dev/null

sleep 3; echo ""

echo mosdepth
memusg mosdepth /tmp/prefix $1

sleep 3; echo ""
echo mosdepth.fast
memusg mosdepth -x /tmp/prefix $1

sleep 3; echo ""

