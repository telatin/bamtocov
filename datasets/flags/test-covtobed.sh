SKIP="Skipped: "
SKIP2=$SKIP
for i in *bam;do
  covtobed $i    | grep  $(basename $i| cut -f2 -d_|cut -f1 -d.) || SKIP="$SKIP $(basename $i| cut -f2 -d_|cut -f1 -d.)" ;
  covtobed --discard-invalid-alignments $i | grep  $(basename $i| cut -f2 -d_|cut -f1 -d.) || SKIP2="$SKIP2 $(basename $i| cut -f2 -d_|cut -f1 -d.)" ;
done > /dev/null

echo Default $SKIP
echo Discard $SKIP

