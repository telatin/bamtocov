# Memory usage

*BamToCov* is an extremely memory efficient tools, and the only one not
using a vector (as long as the chromosome) to store the changes in coverage.

## Results

Memory usage (in bytes) has been measured for alignments against:

* Candida parapsilosis, Illumina dataset
* Candida parapsilosis, Nanopore dataset
* Human exome HG00258, Illumina
* Human targeted sequencing, 16 genes, Illumina

|           | Fungus, Illumina | Fungus, ONT | HG00258 Exome | Human panel |
| --------- | ---------------: | ----------: | ------------: | ----------: |
| bamtocov  |        **2,740** |   **4,376** |     **5,700** |   **2,172** |
| covtobed  |            4,080 |       5,008 |         6,588 |       4,052 |
| mosdepth  |           13,952 |      19,140 |     1,983,928 |   6,425,744 |
| megadepth |           11,644 |      11,636 |       995,232 |     980,040 |
| bedtools  |           12,940 |      14,876 |           n/a |   1,951,288 |

For reference, this is the speed:

| Command                               |       Mean (ms) | Min (ms) | Max (ms) |      Relative |
| :------------------------------------ | --------------: | -------: | -------: | ------------: |
| `bamtocov "panel_01.bam"`             |    358.9 ± 18.3 |    341.6 |    387.9 |          1.00 |
| `covtobed "panel_01.bam"`             |     533.2 ± 8.8 |    527.1 |    548.4 |   1.49 ± 0.08 |
| `megadepth --coverage "panel_01.bam"` | 9246.8 ± 1509.7 |   8026.7 |  11072.8 |  25.77 ± 4.41 |
| `mosdepthprefix "panel_01.bam"`       | 53499.6 ± 875.1 |  52284.1 |  54548.5 | 149.08 ± 7.97 |

## Scripts

### memusg 

Evaluation of memory usage has been performed with `memusg` bu Jaeho Sigh,
as reported below:

```bash
#!/usr/bin/env bash
# memusg -- Measure memory usage of processes
# Usage: memusg COMMAND [ARGS]...
#
# Author: Jaeho Shin <netj@sparcs.org>
# Created: 2010-08-16
############################################################################
# Copyright 2010 Jaeho Shin.                                               #
#                                                                          #
# Licensed under the Apache License, Version 2.0 (the "License");          #
# you may not use this file except in compliance with the License.         #
# You may obtain a copy of the License at                                  #
#                                                                          #
#     http://www.apache.org/licenses/LICENSE-2.0                           #
#                                                                          #
# Unless required by applicable law or agreed to in writing, software      #
# distributed under the License is distributed on an "AS IS" BASIS,        #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. #
# See the License for the specific language governing permissions and      #
# limitations under the License.                                           #
############################################################################
set -um

# check input
[[ $# -gt 0 ]] || { sed -n '2,/^#$/ s/^# //p' <"$0"; exit 1; }

# TODO support more options: peak, footprint, sampling rate, etc.

pgid=$(ps -o pgid= $$)
# make sure we're in a separate process group
if [[ "$pgid" == "$(ps -o pgid= $(ps -o ppid= $$))" ]]; then
    cmd=
    set -- "$0" "$@"
    for a; do cmd+="'${a//"'"/"'\\''"}' "; done
    exec bash -i -c "$cmd"
fi

# detect operating system and prepare measurement
case $(uname) in
    Darwin|*BSD) sizes() { /bin/ps -o rss= -g $1; } ;;
    Linux) sizes() { /bin/ps -o rss= -$1; } ;;
    *) echo "$(uname): unsupported operating system" >&2; exit 2 ;;
esac

# monitor the memory usage in the background.
(
peak=0
while sizes=$(sizes $pgid)
do
    set -- $sizes
    sample=$((${@/#/+}))
    let peak="sample > peak ? sample : peak"
    sleep 0.1
done
echo "memusg: peak=$peak" >&2
) &
monpid=$!


# run the given command
exec "$@"
```

### Memory.sh

To compare the memory usage of multiple tools:

```bash
 
echo MegaDepth:
memusg /local/miniconda3/envs/mega/bin/megadepth --coverage $1 > /dev/null
sleep 3; echo ""

echo BamToCov:
memusg bin/bamtocov $1 > /dev/null
sleep 3; echo ""

echo CovToBed:
memusg covtobed $1 > /dev/null
sleep 3; echo ""

echo Mosdepth:
memusg mosdepth /tmp/prefix $1
sleep 3; echo ""
```

