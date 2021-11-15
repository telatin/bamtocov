#!/bin/bash


# check if username is runner
if [[ $(whoami) == "runner" ]]; then
  echo "Running tests on the cloud"
  CLOUD=1
  set -euxo pipefail
else
  CLOUD=0
  echo user: $(whoami)
  set -euo pipefail
fi

# Check OS

if [[ "$OSTYPE" == "darwin"* ]]; then
  MD5="md5"
else
  MD5="md5sum"
fi
# Function to print to stderr
echoerr() { echo "$@" 1>&2; }

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
BIN="$DIR/../bin/"
DATA="$DIR/../input"

echoerr "Running tests: $DIR"
echoerr "Binary dir: $BIN"
TMP=$(mktemp -d)

BamToCov="$BIN/bamtocov"

PASS=0
FAIL=0

if [[ $CLOUD == 0 ]]; then
  # Get last release from github
  LAST_RELEASE=$(curl -s https://api.github.com/repos/telatin/bamtocov/releases/latest | grep tag_name | cut -d '"' -f 4)
  CURRENT_RELEASE=$(grep ver "$DIR"/../bamtocov.nimble  | perl -ne 'if ($_=~/"([0-9.]+)"/) {print $1}')
fi
# Check if the "fu-tabcheck" command is available in the systsm
if command -v fu-tabcheck >/dev/null 2>&1; then
    echo " Tabcheck ON"
    TABCHECK=1
else
    echo " Will skip tabular check"
    TABCHECK=0
fi

# Check binaries
for binName in bamtocov bamcountrefs bamtocounts covtotarget; do
  if [ ! -e "$BIN/$binName" ]; then
    echo "ERROR: Missing required binary <$binName> in $BIN/"
    ls -l "$binName"/*
    exit 1
  else
    echo -n "[$binName]  "; "$BIN"/$binName --help | head -n 1 
    
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

ORIG_MD5="e09d11db350851b41b97b3ea3c7c41c0"
CALCULATED=$("$BamToCov" "$DATA"/mini.bam | $MD5 | cut -f 1 -d " ")
if [[ $CALCULATED == $ORIG_MD5 ]]; then
  echo "PASS: covtobed style output, MD5"
  PASS=$((PASS+1))
else
  echo "FAIL: covtobed style output MD5: got $CALCULATED, but $ORIG_MD5 expected"
  "$BamToCov" "$DATA"/mini.bam | head
  echo "---"
  FAIL=$((FAIL+1))
fi 

# Wig line
if [[ $("$BamToCov" "$DATA"/mini.bam --wig 250 --op max | wc -l ) -eq $(cat $DIR/results/mini_wig250_max.wig | wc -l) ]]; then
  echo "PASS: Wiggle output fixed step: line numbers"
  PASS=$((PASS+1))
else
  echo "FAIL: Wiggle output has different lines than $DIR/results/mini_wig250_max.wig"
  FAIL=$((FAIL+1))
fi
# Wig lines (no headers)
if [[ $("$BamToCov" "$DATA"/mini.bam --wig 250 --op max | grep -v "fixed" | $MD5 ) == $(cat $DIR/results/mini_wig250_max.wig | grep -v "fixed" | $MD5) ]]; then
  echo "PASS: Wiggle output fixed step"
  PASS=$((PASS+1))
else
  echo "FAIL: Wiggle output differs from $DIR/results/mini_wig250_max.wig"
  FAIL=$((FAIL+1))
fi

"$BamToCov" --regions "$DATA"/mini.bed --report "$TMP"/report.tsv "$DATA"/mini.bam --skip-output

## Check that the tabular report is a proper table (will likely work locally where seqfu is installed)
if [[ $TABCHECK == 1 ]]; then
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
  if [[ ! -e "$TMP"/report.tsv ]]; then
    echo "FAIL: target report missing $TMP/report.tsv"
    FAIL=$((FAIL+1))
  else
    echo "PASS: target report found"
    PASS=$((PASS+1))
  fi
fi

## Check release github vs nimble
if [[ $CLOUD == 0 ]]; then
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

  echo "Last release:    $LAST_RELEASE"
  echo "Current release: v$CURRENT_RELEASE"
  echo "Binary release:  $($BIN/bamtocov --version | sed 's/covtobed /v/')"
fi
echo "--------------------------"
echo "SUMMARY (PASS=$PASS,FAIL=$FAIL)"

echo "--------------------------"
if [[ $FAIL -gt 0 ]]; then
  echo "FINAL RESULT: FAIL"
  exit 1
else
  echo "FINAL RESULT: PASS"
fi