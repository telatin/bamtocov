#!/bin/bash

for DATA in ../exome/HG00258-100.bam ../exome/HG00258-1000.bam ../pannelli/panel_01.bam _test/cpara-illumina-noseq.bam  _test/cpara-ont-noseq.bam ../exome/HG00258.chr21.bam ../exome/HG00258.bam;
do
    NAME=$(basename "$DATA" | sed 's/.bam//g')
    echo "=========== $NAME ==========="
    
    if [[ ! -e mem.covtobed.$NAME.txt ]]; then
      echo CovToBed 
      memusg2.py --no-humanize -p 1 -t -o mem.covtobed.$NAME.txt covtobed $DATA > /dev/null
    fi

    if [[ ! -e mem.bamtocov.$NAME.txt ]]; then
      echo BamTocov
      memusg2.py --no-humanize -p 1 -t -o mem.bamtocov.$NAME.txt bamtocov $DATA > /dev/null
    fi

    if [[ ! -e mem.mosdepth.$NAME.txt ]]; then
      echo MosDepth
      memusg2.py --no-humanize -p 1 -t -o mem.mosdepth.$NAME.txt  mosdepth -x /tmp/tmpmd $DATA > /dev/null
    fi

    if [[ ! -e mem.megadepth.$NAME.txt ]]; then
      echo MegaDepth
      memusg2.py --no-humanize -p 1 -t -o mem.megadepth.$NAME.txt megadepth --coverage $DATA > /dev/null
    fi
done
