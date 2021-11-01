#!/bin/bash
set -euo pipefail
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
BIN="$DIR/../bin/"
DATA="$DIR/../input"

TMP=$(mktemp -d)

BamToCov="$BIN/bamtocov"

PASS=0
FAIL=0

# Get last release from github
LAST_RELEASE=$(curl -s https://api.github.com/repos/telatin/bamtocov/releases/latest | grep tag_name | cut -d '"' -f 4)
CURRENT_RELEASE=$(grep ver "$DIR"/../bamtocov.nimble  | perl -ne 'if ($_=~/"([0-9.]+)"/) {print $1}')
# Check if the "fu-tabcheck" command is available in the systsm
if command -v fu-tabcheck >/dev/null 2>&1; then
    TABCHECK=1
fi
# Check binaries
for bin in bamtocov bamcountrefs bamtocounts covtotarget; do
  if [ ! -e "$BIN/$bin" ]; then
    echo "ERROR: Missing required binary <$bin> in $BIN/"
    exit 1
  else
    PASS=$((PASS+1))
  fi
done

# Check test files
for testFile in mini.bam mini.bed; do
  if [ ! -e "$DATA/$testFile" ]; then
    echo "ERROR: Missing required binary <$testFile> in $DATA/"
    exit 1
  else
    PASS=$((PASS+1))
  fi
done

# Check covtobed style
if [[ $("$BamToCov" "$DATA"/mini.bam | wc -l) -eq 21 ]]; then
  echo "PASS: covtobed style output, lines"
  PASS=$((PASS+1))
else
  echo "FAIL: covtobed style output $($BamToCov $DATA/mini.bam | wc -l) lines, but 21 expected"
  FAIL=$((FAIL+1))
fi 

MD5="e09d11db350851b41b97b3ea3c7c41c0"
if [[ $("$BamToCov" "$DATA"/mini.bam | md5sum | cut -f 1 -d " ") == $MD5 ]]; then
  echo "PASS: covtobed style output, MD5"
  PASS=$((PASS+1))
else
  echo "FAIL: covtobed style output MD5 $($BamToCov $DATA/mini.bam | md5sum) lines, but e09d11db350851b41b97b3ea3c7c41c0 expected"
  FAIL=$((FAIL+1))
fi 

"$BamToCov" --regions "$DATA"/mini.bed --report "$TMP"/report.tsv "$DATA"/mini.bam --skip-output

## Check that the tabular report is a proper table (will likely work locally where seqfu is installed)
if [[ $TABCHECK ]]; then
  set +e
  fu-tabcheck "$TMP"/report.tsv >/dev/null
  if [[ $? == 0 ]]; then
    echo "PASS: target report, tabcheck $(fu-tabcheck "$TMP"/report.tsv >/dev/null)"
    PASS=$((PASS+1))
  else
    echo "FAIL: target report tabcheck"
    FAIL=$((FAIL+1))
  fi
  set -e
else
  echo "SKIP: tab check [requires bioconda:seqfu]"
fi

## Check release github vs nimble

if [[ $LAST_RELEASE == $CURRENT_RELEASE ]]; then
  echo "FAIL: current release $CURRENT_RELEASE"
  FAIL=$((FAIL+1))
else
  echo "PASS: Current nimble version is different from GitHub release (should be newer)"
  PASS=$((PASS+1))
fi

if [[ $($BIN/bamtocov --version | sed 's/covtobed //') != $CURRENT_RELEASE ]]; then
  echo "FAIL: Nimble version differs from binary: recompile!"
  FAIL=$((FAIL+1))
else
  echo "PASS: Nimble version matches binary"
  PASS=$((PASS+1))
fi

rm -rf $TMP
VER=$(grep version bamtocov.nimble  | cut -f2 -d\")
for i in "$BIN"/*; do 
  b=$(basename "$i")
  echo "Checking version for $b: $($i --version | grep $VER)"; 
done

echo "--------------------------"
echo "SUMMARY (PASS=$PASS,FAIL=$FAIL)"
echo "Last release:    $LAST_RELEASE"
echo "Current release: v$CURRENT_RELEASE"
echo "Binary release:  $($BIN/bamtocov --version | sed 's/covtobed /v/')"
echo "--------------------------"
if [[ $FAIL -gt 0 ]]; then
  echo "FINAL RESULT: FAIL"
  exit 1
else
  echo "FINAL RESULT: PASS"
fi