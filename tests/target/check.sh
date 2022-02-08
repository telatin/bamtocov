#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

DATA="${DIR}/../../input"
BIN="${DIR}/../../bin"

# Create array of binaries: $BIN/bamtocov and $BIN/bamtocov-check
declare -a BINARIES=($BIN/bamtocov $BIN/bamtocov-2.4.0 $BIN/bamtocov-2.5)
# Cycle through array of binaries
for BINARY in "${BINARIES[@]}"
do
    # Check if binary exists
    if [ ! -f "${BINARY}" ]
    then
        echo "Binary ${BINARY} does not exist"
        exit 1
    fi
    # Check if binary is executable
    if [ ! -x "${BINARY}" ]
    then
        echo "Binary ${BINARY} is not executable"
        exit 1
    fi
done

## Mini test

BAM="${DATA}/mini.bam"
BED="${DATA}/regions.bed"

for BINARY in "${BINARIES[@]}"
do
    VERSION=$($BINARY --version)
    echo Testing $VERSION
    $BINARY -r $BED $BAM > /dev/null
done
