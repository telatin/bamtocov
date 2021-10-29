# Testing BamToCov

```note
To do
```

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