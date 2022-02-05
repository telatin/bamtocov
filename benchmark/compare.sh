#!/usr/bin/env bash

if [[ ! -z ${2+x} ]]; then
  out=$2
else
  out="benchmark"
fi

function compare {
 tmpfile=$(mktemp mosdepth.XXXXXX)
 hyperfine --export-markdown $2.md --min-runs 5 \
  "covtobed '$1'" \
  "bamtocov '$1'" \
  "mosdepth "$tmpfile" '$1'" \
  "mosdepth --fast-mode _$out '$1'" \
  "megadepth --coverage '$1'" \
  "megadepth --coverage --longreads '$1'"
  rm "$tmpfile"*
}

if [  -e "$1.bai" ]; then
 compare "$1" "$out"
 cat $out.md

else
  for URL in "https://zenodo.org/record/5636944/files/cpara-illumina-noseq.bam?download=1" \
	"https://zenodo.org/record/5636944/files/cpara-ont-noseq.bam?download=1" \
	"https://zenodo.org/record/5636944/files/HG00258.bam?download=1" \
	"https://zenodo.org/record/5636944/files/panel_01.bam?download=1";
  do
    echo pretest
    tmpfile=$(mktemp mosdepth.XXXXXX)
    mosdepth "$tmpfile" '$1' && rm "$tmpfile"*
    NAME=$(basename "$URL" | cut -f 1 -d "?")
    wget --quiet -O "$NAME" "$URL"
    compare "$NAME" "$NAME"

  done
fi

