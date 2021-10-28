SKIP=""
SKIP2=$SKIP
for i in *bam;do
  FLAG=$(basename $i| cut -f2 -d_|cut -f1 -d.)
  ../../bin/bamtocov $i  | grep  $FLAG || SKIP="$SKIP $FLAG" ;
  ../../bin/bamtocov $i  --flag 4  | grep  $FLAG || SKIP2="$SKIP2 $FLAG" ;
done > /dev/null

echo Default $SKIP
echo
echo Relaxed $SKIP2

