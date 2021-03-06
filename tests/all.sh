#!/bin/bash
set -euxo pipefail
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
BIN="$DIR/../bin/"
DATA="$DIR/../input"

TMP=$(mktemp -d)

BamToCov="$BIN/bamtocov"

# Check binaries
for bin in bamtocov bamcountrefs bamtocounts covtotarget; do
  if [ ! -e "$BIN/$bin" ]; then
    echo "ERROR: Missing required binary <$bin> in $BIN/"
    exit 1
  fi
done

# Check test files
for testFile in mini.bam mini.bed; do
  if [ ! -e "$DATA/$testFile" ]; then
    echo "ERROR: Missing required binary <$testFile> in $DATA/"
    exit 1
  fi
done

$BamToCov $DATA/mini.bam > $TMP/mini.BamToCov
tail -n 4 $TMP/mini.BamToCov
$BamToCov -r $DATA/mini.bed -o $TMP/report.txt $DATA/mini.bam > $TMP/mini.BamToCov
tail -n 4 $TMP/mini.BamToCov $TMP/report.txt

rm -rf $TMP
VER=$(grep version bamtocov.nimble  | cut -f2 -d\")
for i in bin/*; do 
  $i --version | grep $VER; 
done