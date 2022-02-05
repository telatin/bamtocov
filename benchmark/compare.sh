#!/usr/bin/env bash

if [[ ! -z ${2+x} ]]; then
  out=$2
else
  out="benchmark"
fi
COVTOBED_VERSION=$(covtobed --version | head -n 1 | grep -oP \\d\+\.\\d\+\.\\d\+)
BAMTOCOV_VERSION=$(bamtocov --version)
MEGADEPTH_VERSION=$(megadepth --version | grep -oP \\d\+\.\\d\+\.\\d\+)
MOSDEPTH_VERSION=$(mosdepth --version | grep -oP \\d\+\.\\d\+\.\\d\+)

function compare {
 tmpfile=$(mktemp mosdepth.XXXXXX)
 samtools index "$1"
 hyperfine --export-markdown $2.md --min-runs 5 \
  "covtobed '$1' #$COVTOBED_VERSION" \
  "bamtocov '$1' #$BAMTOCOV_VERSION" \
  "megadepth --coverage '$1' #$MEGADEPTH_VERSION" \
  "megadepth --coverage --longreads '$1' #$MEGADEPTH_VERSION" \
  "mosdepth $tmpfile '$1' #$MOSDEPTH_VERSION" \
  "mosdepth --fast-mode $tmpfile '$1' #$MOSDEPTH_VERSION"

  rm "$tmpfile"*
}

echo "covtobed  : $COVTOBED_VERSION"
echo "bamtocov  : $BAMTOCOV_VERSION"
echo "megadepth : $MEGADEPTH_VERSION"
echo "mosdepth  : $MOSDEPTH_VERSION"


if [  -e "$1.bai" ]; then
 compare "$1" "$out"
 cat $out.md

else
  for URL in "https://zenodo.org/record/5636944/files/cpara-illumina-noseq.bam?download=1" \
            "https://zenodo.org/record/5636944/files/cpara-ont-noseq.bam?download=1" \
            "https://zenodo.org/record/5636944/files/HG00258.bam?download=1" \
            "https://zenodo.org/record/5636944/files/panel_01.bam?download=1" \
            "https://zenodo.org/record/1251249/files/w1118_f_hd_R1.junc.bam?download=1" \
            "https://zenodo.org/record/1221361/files/Library_01.bam?download=1" ;

  do
    
    NAME=$(basename "$URL" | cut -f 1 -d "?")
    wget --quiet -O "$NAME" "$URL"
 
    compare "$NAME" "$NAME"

  done
fi

