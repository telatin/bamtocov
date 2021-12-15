# Testing BamToCov

## Test suite

The program default test suite can be invoked by:

```bash
make test
```

An extended (and redundant) set of tests can be run with:

```bash
make testall
```

This runs both the default tests, plus a legacy Bash script

## Shpec tests

The default test suite was written using the [Shpec](https://github.com/rylnd/shpec)
framework, and can be individually run as:

```bash
# Same as "make test"
./tests/bin/shpec ./tests/shpec/bamtocov.sh
```

And will produce an indented output covering the binaries tested,
similar to the example below:

```text
Make test
  BamToCov
    Binary exist
    Version emitted is 2.x
    Mini coverage, verify output line number
    One line with empty chromosome seq0 (empty chromosome)
    Works with sorted file
    Fails with unsorted file
    Produces wig output (check lines)
    Produces wig output (check lines at 750)
    Produces wig output (check lines at 1000 bases, unexpected)
    Produces wig output header
  BamToCounts
    Binary exists
    Version 2.x
    Counts target
    Coverage check(x6)
14 examples, 0 failures
```

## Legacy test suite

A minimal test suite is provided in `tests/all.sh`. The output will
show the passing and failing tests as well as a check on the versions
(ensuring that the current development version is different from the
released one):

```text
PASS: covtobed style output, lines
PASS: covtobed style output, MD5
PASS: target report, tabcheck 
PASS: Current nimble version is different from GitHub release (should be newer)
PASS: Nimble version matches binary
Checking version for bamcountrefs: 2.2.1
Checking version for bamtocounts: 2.2.1
Checking version for bamtocov: covtobed 2.2.1
Checking version for covtotarget: 2.2.1
--------------------------
SUMMARY (PASS=11,FAIL=0)
Last release:    v2.1.0
Current release: v2.2.1
Binary release:  v2.2.1
--------------------------
FINAL RESULT: PASS
```
