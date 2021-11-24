#!/bin/sh
# vim:et:ft=sh:sts=2:sw=2

oneTimeSetUp() {
  outputDir="${SHUNIT_TMPDIR}"/output
  mkdir "${outputDir}"
  stdoutF="${outputDir}/stdout"
  stderrF="${outputDir}/stderr"
  SCRIPT_SELF_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
  mkdirCmd='mkdir'  # save command name in variable to make future changes easy
  testDir="${SHUNIT_TMPDIR}/some_test_dir"

  # Input directories
  binDir="$SCRIPT_SELF_DIR"/../../bin/
  dataDir="$SCRIPT_SELF_DIR"/../../input/
  phiDir="$SCRIPT_SELF_DIR"/../../input/phi/

  # Binaries
  bamtocov="$binDir"/bamtocov
  bamcountrefs="$binDir"/bamcountrefs
  covtotarget="$binDir"/covtotarget
  bamtocounts="$binDir"/bamtocounts
}

testBinDirectoriesFound() {
  assertTrue "Missing binaries dir: $binDir"  "[ -d '${binDir}' ]"
  assertTrue "Missing data dir: $dataDir"     "[ -d '${dataDir}' ]"
  assertTrue "Missing phiX dir: $phiDir"      "[ -d '${phiDir}' ]"
}

testBinariesFound() {
  assertTrue "Missing bamtocov binary: $bamtocov"          "[ -x '${bamtocov}' ]"
  assertTrue "Missing bamcountrefs binary: $bamcountrefs"  "[ -x '${bamcountrefs}' ]"
  assertTrue "Missing bamtocounts binary: $bamtocounts"    "[ -x '${bamtocounts}' ]"
  assertTrue "Missing covtotarget binary: $covtotarget"    "[ -x '${covtotarget}' ]"
}
 
testCreateTemporaryDir() {
  ${mkdirCmd} "${testDir}" #>${stdoutF} 2>${stderrF}
  rtrn=$?
  th_assertTrueWithNoOutput ${rtrn} "${stdoutF}" "${stderrF}"
  assertTrue 'directory missing' "[ -d '${testDir}' ]"
}

testVersion() {
  ${bamtocov} --version > "${outputDir}"/version.txt
  rtrn=$?
  assertTrue 'version command failed' ${rtrn}
  assertTrue 'version output empty' "[ -s '${outputDir}/version.txt' ]"
  assertTrue 'version is missing dots' "[ $(grep -c \\. "${outputDir}"/version.txt) -eq 1 ]"
}

testBamToCov() {
  ${bamtocov} "${dataDir}"/mini.bam > "${outputDir}"/mini.bed
  rtrn=$?
  
  # Program was executed successfully
  assertTrue "bamtocov failed: ${bamtocov} ${dataDir}/mini.bam > ${testDir}/mini.bed" ${rtrn}
  # Output file was created
  assertTrue "Expecing 21 lines from mini.bam"            "[ $(wc -l < ${outputDir}/mini.bed) -eq 21 ]"
  
  # Output is consistent with expectation
  assertTrue "Non empty chromosome seq1 returns 17 line"  "[ $(grep -c 'seq1' ${outputDir}/mini.bed) -eq 17 ]"
  assertTrue "Maximum coverage in seq1 is 6"              "[ $(grep  seq1 ${outputDir}/mini.bed | cut -f 4 | sort -n | tail -n 1 ) -eq 6 ]"
  assertTrue "Empty chromosome seq0 returns 1 line"       "[ $(grep -c 'seq0' ${outputDir}/mini.bed) -eq 1 ]"
}

testBamToCov() {
  ${bamtocov} "${dataDir}"/mini.bam > "${outputDir}"/mini.bed
  rtrn=$?
  
  # Program was executed successfully
  assertTrue "bamtocov failed: ${bamtocov} ${dataDir}/mini.bam > ${testDir}/mini.bed" ${rtrn}
  # Output file was created
  assertTrue "Expecing 21 lines from mini.bam"            "[ $(wc -l < ${outputDir}/mini.bed) -eq 21 ]"
  
  # Output is consistent with expectation
  assertTrue "Non empty chromosome seq1 returns 17 line"  "[ $(grep -c 'seq1' ${outputDir}/mini.bed) -eq 17 ]"
  assertTrue "Maximum coverage in seq1 is 6"              "[ $(grep  seq1 ${outputDir}/mini.bed | cut -f 4 | sort -n | tail -n 1 ) -eq 6 ]"
  assertTrue "Empty chromosome seq0 returns 1 line"       "[ $(grep -c 'seq0' ${outputDir}/mini.bed) -eq 1 ]"
}

testBamToCovStream() {
  cat "${dataDir}"/mini.bam | ${bamtocov} 2>/dev/null > "${outputDir}"/mini-stream.bed
  # Program was executed successfully
  assertTrue "bamtocov failed: ${bamtocov} ${dataDir}/mini.bam > ${testDir}/mini-stream.bed" ${rtrn}
  # Output file was created
  assertTrue "Expecing 21 lines from mini.bam"            "[ $(wc -l < ${outputDir}/mini-stream.bed) -eq 21 ]"
}

testBamToCovQuantized() {
  ${bamtocov} --quantize 2,4 "${dataDir}"/mini.bam > "${outputDir}"/mini-quantized.bed
  rtrn=$?
  
  # Program was executed successfully
  assertTrue "bamtocov quantized failed: ${bamtocov} --quantize 2,4 " ${rtrn}
  # Output file was created
  assertTrue "Expecing output file"                            "[ -f '${outputDir}/mini-quantized.bed' ]"
  assertTrue "Expecting three classes of coverage "            "[ $(cat ${outputDir}/mini-quantized.bed | cut -f 4 | sort | uniq | wc -l ) -eq 3 ]"  
  assertTrue "First class is open interval 0-1"                "[ $(cat ${outputDir}/mini-quantized.bed | cut -f 4 | sort | uniq | head -n 1) = 0-1 ]"  
  assertTrue "Last class is open interval 4-"                  "[ $(cat ${outputDir}/mini-quantized.bed | cut -f 4 | sort | uniq | tail -n 1) = 4- ]"  
}

testBamToCovPhiX() {
  ${bamtocov} "${phiDir}"/shotgun.bam > "${outputDir}"/phi.bed
  # Program was executed successfully
  assertTrue "bamtocov failed on shotgun.bam" ${rtrn}
  # Output file was created
  assertTrue "Expecing 1375 lines from mini.bam"         "[ $(wc -l < ${outputDir}/phi.bed) -eq 1375 ]"
  
  # Output is consistent with expectation
  assertTrue "Maximum coverage is 36X"                   "[ $(cat ${outputDir}/phi.bed | cut -f 4 | sort -n | tail -n 1 ) -eq 36 ]"
  assertTrue "Maximum coverage is 36X at pos424"         "[ $(cat ${outputDir}/phi.bed |  sort -n -k 4 | tail -n 1 | cut -f 2) -eq 424 ]"
  
}

testStrandedCoverage() {
  ${bamtocov} --stranded "${dataDir}"/mini.bam > "${outputDir}"/mini-stranded.bed
  assertTrue "bamtocov stranded in mini.bam" ${rtrn}
  assertTrue "Maximum REVERSE coverage is 5X"  "[ $(cat "${outputDir}"/mini-stranded.bed | cut -f 5 | sort -n | tail -n 1) -eq 5 ]"
}

testOutputStillExists() {
  assertTrue "bamtocov failed on shotgun.bam" "[ -e ${outputDir}/phi.bed ]"
}

th_assertTrueWithNoOutput() {
  th_return_=$1
  th_stdout_=$2
  th_stderr_=$3

  assertFalse 'unexpected output to STDOUT' "[ -s '${th_stdout_}' ]"
  assertFalse 'unexpected output to STDERR' "[ -s '${th_stderr_}' ]"

  unset th_return_ th_stdout_ th_stderr_
}


tearDown() {
  # Executed at the end of the tests
  rm -fr "${testDir}"
}

# Load and run shUnit2.
[ -n "${ZSH_VERSION:-}" ] && SHUNIT_PARENT=$0
SCRIPT_SELF_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
. "$SCRIPT_SELF_DIR"/../shunit2
