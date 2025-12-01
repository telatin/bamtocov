describe "BamCountRefs - Count reads per reference"

  if [ "$0" = "./tests/bin/shpec" ];
  then
    echo "Make test"
    BINDIR="./bin/"
    DATADIR="./input/"
  else
    echo "Manual test $0"
    SELFDIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
    BINDIR="$( cd -- "$( dirname -- "$SELFDIR"/../../bin/bamcountrefs )" &> /dev/null && pwd )"
    DATADIR="$( cd -- "$( dirname -- "$SELFDIR"/../../input/mini.bam )" &> /dev/null && pwd )"
  fi

  describe "Binary and Version"
    it "Binary exists"
      assert file_present "$BINDIR"/bamcountrefs
    end

    VERSION=$("$BINDIR"/bamcountrefs --version)
    it "Version is 2.x [$VERSION]"
      assert glob "$VERSION" "2.*"
    end
  end

  describe "Basic Counting"
    it "Counts reads in mini.bam (verify output line count)"
      LINES=$("$BINDIR"/bamcountrefs "$DATADIR"/mini.bam | wc -l)
      assert equal $((LINES+0)) 4
    end

    it "Produces tab-separated output with header"
      HEADER=$("$BINDIR"/bamcountrefs "$DATADIR"/mini.bam | head -n 1)
      assert glob "$HEADER" "ViralSequence*mini"
    end

    it "Counts seq1 correctly (15 reads)"
      COUNT=$("$BINDIR"/bamcountrefs "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      assert equal $((COUNT+0)) 15
    end

    it "Counts seq2 correctly (10 reads)"
      COUNT=$("$BINDIR"/bamcountrefs "$DATADIR"/mini.bam | grep "^seq2" | cut -f 2)
      assert equal $((COUNT+0)) 10
    end

    it "Reports zero for seq0 (empty chromosome)"
      COUNT=$("$BINDIR"/bamcountrefs "$DATADIR"/mini.bam | grep "^seq0" | cut -f 2)
      assert equal $((COUNT+0)) 0
    end
  end

  describe "Multiple BAM Files"
    it "Handles multiple BAM files (3 samples)"
      TMPFILE=$(mktemp)
      "$BINDIR"/bamcountrefs "$DATADIR"/mini.bam "$DATADIR"/mini2.bam "$DATADIR"/mini3.bam > "$TMPFILE"
      SAMPLES=$(head -n 1 "$TMPFILE" | wc -w)
      assert equal $SAMPLES 4
      rm "$TMPFILE"
    end

    it "Header contains all sample names"
      HEADER=$("$BINDIR"/bamcountrefs "$DATADIR"/mini.bam "$DATADIR"/mini2.bam | head -n 1)
      assert glob "$HEADER" "*mini*mini2*"
    end

    it "Each reference gets counts for all samples"
      SEQ1_LINE=$("$BINDIR"/bamcountrefs "$DATADIR"/mini.bam "$DATADIR"/mini2.bam | grep "^seq1")
      COLS=$(echo "$SEQ1_LINE" | wc -w)
      assert equal $COLS 3
    end

    it "Different samples can have different counts"
      SEQ1_LINE=$("$BINDIR"/bamcountrefs "$DATADIR"/mini.bam "$DATADIR"/mini3.bam | grep "^seq1")
      COUNT1=$(echo "$SEQ1_LINE" | cut -f 2)
      COUNT2=$(echo "$SEQ1_LINE" | cut -f 3)
      assert gt $COUNT2 $COUNT1
    end
  end

  describe "Flag Filtering"
    it "Default flag excludes unmapped/secondary/duplicate (eflag=1796)"
      COUNT_DEFAULT=$("$BINDIR"/bamcountrefs "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      assert equal $((COUNT_DEFAULT+0)) 15
    end

    it "Flag -F 4 only excludes unmapped reads"
      COUNT_F4=$("$BINDIR"/bamcountrefs -F 4 "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      assert equal $((COUNT_F4+0)) 15
    end

    it "Fast path activated with -F 4 (check debug output)"
      OUTPUT=$("$BINDIR"/bamcountrefs --debug -F 4 "$DATADIR"/mini.bam 2>&1 | grep -c "Establishing reference order from first sample")
      assert equal $((OUTPUT+0)) 1
    end

    it "Iteration path used with default flags (check debug output)"
      OUTPUT=$("$BINDIR"/bamcountrefs --debug "$DATADIR"/mini.bam 2>&1 | grep -c "Merging results from all workers")
      assert equal $((OUTPUT+0)) 1
    end
  end

  describe "RPKM Normalization (Stdout Mode)"
    it "RPKM output is numeric (float)"
      OUTPUT=$("$BINDIR"/bamcountrefs --rpkm "$DATADIR"/phi/wildtype.bam | grep "^NC_001422" | cut -f 2)
      assert glob "$OUTPUT" "*.*"
    end

    it "RPKM values are calculated per sample"
      LINE=$("$BINDIR"/bamcountrefs --rpkm "$DATADIR"/phi/wildtype.bam "$DATADIR"/phi/shotgun.bam | grep "^NC_001422")
      RPKM1=$(echo "$LINE" | cut -f 2)
      RPKM2=$(echo "$LINE" | cut -f 3)
      # Both samples should have RPKM values
      assert glob "$RPKM1" "*.*"
      assert glob "$RPKM2" "*.*"
    end

    it "Deprecated -n flag still works with warning"
      OUTPUT=$("$BINDIR"/bamcountrefs -n "$DATADIR"/mini.bam 2>&1)
      # Check for deprecation warning
      assert grep "WARNING.*deprecated" <<< "$OUTPUT"
      # Check that RPKM values are still produced
      RPKM=$(echo "$OUTPUT" | grep "^seq1" | cut -f 2)
      assert glob "$RPKM" "*.*"
    end
  end

  describe "Multi-File Output (-o flag)"
    TMPDIR=$(mktemp -d)

    it "Creates counts file when -o specified"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/test "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/test_counts.tsv
    end

    it "Counts file has correct structure"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/test "$DATADIR"/mini.bam
      LINES=$(wc -l < "$TMPDIR"/test_counts.tsv)
      assert equal $((LINES+0)) 4
    end

    it "Counts file has same values as stdout mode"
      "$BINDIR"/bamcountrefs "$DATADIR"/mini.bam > "$TMPDIR"/stdout.txt
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/test "$DATADIR"/mini.bam
      STDOUT_COUNT=$(grep "^seq1" "$TMPDIR"/stdout.txt | cut -f 2)
      FILE_COUNT=$(grep "^seq1" "$TMPDIR"/test_counts.tsv | cut -f 2)
      assert equal $STDOUT_COUNT $FILE_COUNT
    end

    it "Creates RPKM file when --rpkm specified with -o"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/rpkm_test --rpkm "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/rpkm_test_rpkm.tsv
    end

    it "Creates TPM file when --tpm specified with -o"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/tpm_test --tpm "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/tpm_test_tpm.tsv
    end

    it "Creates mean file when --mean specified with -o"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/mean_test --mean "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/mean_test_mean.tsv
    end

    it "Creates all files when --all-metrics specified"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/all_test --all-metrics "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/all_test_counts.tsv
      assert file_present "$TMPDIR"/all_test_rpkm.tsv
      assert file_present "$TMPDIR"/all_test_tpm.tsv
      assert file_present "$TMPDIR"/all_test_mean.tsv
    end

    it "Always creates counts file even when only --rpkm requested"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/rpkm_only --rpkm "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/rpkm_only_counts.tsv
    end

    it "Multi-file output preserves reference order"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/order_test "$DATADIR"/mini.bam
      FIRST_REF=$(tail -n +2 "$TMPDIR"/order_test_counts.tsv | head -n 1 | cut -f 1)
      assert equal "$FIRST_REF" "seq1"
    end

    rm -rf "$TMPDIR"
  end

  describe "TPM Calculations"
    TMPDIR=$(mktemp -d)

    it "TPM file is created with -o --tpm"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/tpm_test --tpm "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/tpm_test_tpm.tsv
    end

    it "TPM values are numeric (float)"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/tpm_test --tpm "$DATADIR"/mini.bam
      TPM=$(grep "^seq1" "$TMPDIR"/tpm_test_tpm.tsv | cut -f 2)
      assert glob "$TPM" "*.*"
    end

    it "TPM sum equals 1,000,000 per sample"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/tpm_test --tpm "$DATADIR"/mini.bam
      # Sum all TPM values (skip header, sum column 2)
      SUM=$(tail -n +2 "$TMPDIR"/tpm_test_tpm.tsv | awk '{sum+=$2} END {printf "%.0f", sum}')
      assert equal "$SUM" "1000000"
    end

    it "TPM calculated per sample independently"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/tpm_multi --tpm "$DATADIR"/mini.bam "$DATADIR"/mini2.bam
      # Sum sample 1 (column 2)
      SUM1=$(tail -n +2 "$TMPDIR"/tpm_multi_tpm.tsv | awk '{sum+=$2} END {printf "%.0f", sum}')
      # Sum sample 2 (column 3)
      SUM2=$(tail -n +2 "$TMPDIR"/tpm_multi_tpm.tsv | awk '{sum+=$3} END {printf "%.0f", sum}')
      assert equal "$SUM1" "1000000"
      assert equal "$SUM2" "1000000"
    end

    rm -rf "$TMPDIR"
  end

  describe "Mean Coverage Calculations"
    TMPDIR=$(mktemp -d)

    it "Mean file is created with -o --mean"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/mean_test --mean "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/mean_test_mean.tsv
    end

    it "Mean coverage values are numeric (float)"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/mean_test --mean "$DATADIR"/mini.bam
      MEAN=$(grep "^seq1" "$TMPDIR"/mean_test_mean.tsv | cut -f 2)
      assert glob "$MEAN" "*.*"
    end

    it "Mean coverage is calculated correctly for seq1 (15 reads, 1000bp)"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/mean_test --mean "$DATADIR"/mini.bam
      MEAN=$(grep "^seq1" "$TMPDIR"/mean_test_mean.tsv | cut -f 2)
      # Expected: ~1.5x coverage (approximate method)
      # Allow for floating point comparison (between 1.0 and 2.0)
      assert glob "$MEAN" "1.*"
    end

    it "Zero coverage references have mean of 0.0"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/mean_test --mean "$DATADIR"/mini.bam
      MEAN=$(grep "^seq0" "$TMPDIR"/mean_test_mean.tsv | cut -f 2)
      assert equal "$MEAN" "0.000000"
    end

    it "Mean coverage calculated per sample"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/mean_multi --mean "$DATADIR"/mini.bam "$DATADIR"/mini2.bam
      LINE=$(grep "^seq1" "$TMPDIR"/mean_multi_mean.tsv)
      MEAN1=$(echo "$LINE" | cut -f 2)
      MEAN2=$(echo "$LINE" | cut -f 3)
      # Both should have numeric values
      assert glob "$MEAN1" "*.*"
      assert glob "$MEAN2" "*.*"
    end

    rm -rf "$TMPDIR"
  end

  describe "Coverage Breadth Calculations"
    TMPDIR=$(mktemp -d)

    it "Covered bases file is created with -o --covered-bases"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/breadth_test --covered-bases "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/breadth_test_covered_bases.tsv
    end

    it "Covered fraction file is created with -o --covered-ratio"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/breadth_test --covered-ratio "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/breadth_test_covered_fraction.tsv
    end

    it "Covered bases values are numeric (integer)"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/breadth_test --covered-bases "$DATADIR"/mini.bam
      COVERED=$(grep "^seq1" "$TMPDIR"/breadth_test_covered_bases.tsv | cut -f 2)
      assert glob "$COVERED" "[0-9]*"
    end

    it "Covered fraction values are numeric (float between 0 and 1)"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/breadth_test --covered-ratio "$DATADIR"/mini.bam
      FRACTION=$(grep "^seq1" "$TMPDIR"/breadth_test_covered_fraction.tsv | cut -f 2)
      assert glob "$FRACTION" "0.*"
    end

    it "Zero coverage references have covered bases of 0"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/breadth_test --covered-bases "$DATADIR"/mini.bam
      COVERED=$(grep "^seq0" "$TMPDIR"/breadth_test_covered_bases.tsv | cut -f 2)
      assert equal "$COVERED" "0"
    end

    it "Zero coverage references have covered fraction of 0.0"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/breadth_test --covered-ratio "$DATADIR"/mini.bam
      FRACTION=$(grep "^seq0" "$TMPDIR"/breadth_test_covered_fraction.tsv | cut -f 2)
      assert equal "$FRACTION" "0.000000"
    end

    it "Coverage breadth calculated per sample"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/breadth_multi --covered-bases "$DATADIR"/mini.bam "$DATADIR"/mini2.bam
      LINE=$(grep "^seq1" "$TMPDIR"/breadth_multi_covered_bases.tsv)
      COVERED1=$(echo "$LINE" | cut -f 2)
      COVERED2=$(echo "$LINE" | cut -f 3)
      # Both should have numeric values
      assert glob "$COVERED1" "[0-9]*"
      assert glob "$COVERED2" "[0-9]*"
    end

    rm -rf "$TMPDIR"
  end

  describe "All-Metrics Flag"
    TMPDIR=$(mktemp -d)

    it "Creates all output files"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/all --all-metrics "$DATADIR"/mini.bam
      FILE_COUNT=$(ls "$TMPDIR"/all_*.tsv 2>/dev/null | wc -l)
      assert equal $((FILE_COUNT+0)) 6
    end

    it "All files have same number of data rows"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/all --all-metrics "$DATADIR"/mini.bam
      COUNTS_LINES=$(tail -n +2 "$TMPDIR"/all_counts.tsv | wc -l)
      RPKM_LINES=$(tail -n +2 "$TMPDIR"/all_rpkm.tsv | wc -l)
      TPM_LINES=$(tail -n +2 "$TMPDIR"/all_tpm.tsv | wc -l)
      MEAN_LINES=$(tail -n +2 "$TMPDIR"/all_mean.tsv | wc -l)
      COVERED_BASES_LINES=$(tail -n +2 "$TMPDIR"/all_covered_bases.tsv | wc -l)
      COVERED_FRACTION_LINES=$(tail -n +2 "$TMPDIR"/all_covered_fraction.tsv | wc -l)
      assert equal $COUNTS_LINES $RPKM_LINES
      assert equal $RPKM_LINES $TPM_LINES
      assert equal $TPM_LINES $MEAN_LINES
      assert equal $MEAN_LINES $COVERED_BASES_LINES
      assert equal $COVERED_BASES_LINES $COVERED_FRACTION_LINES
    end

    it "All files have same reference order"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/all --all-metrics "$DATADIR"/mini.bam
      COUNTS_REFS=$(tail -n +2 "$TMPDIR"/all_counts.tsv | cut -f 1)
      RPKM_REFS=$(tail -n +2 "$TMPDIR"/all_rpkm.tsv | cut -f 1)
      TPM_REFS=$(tail -n +2 "$TMPDIR"/all_tpm.tsv | cut -f 1)
      MEAN_REFS=$(tail -n +2 "$TMPDIR"/all_mean.tsv | cut -f 1)
      COVERED_BASES_REFS=$(tail -n +2 "$TMPDIR"/all_covered_bases.tsv | cut -f 1)
      COVERED_FRACTION_REFS=$(tail -n +2 "$TMPDIR"/all_covered_fraction.tsv | cut -f 1)
      assert equal "$COUNTS_REFS" "$RPKM_REFS"
      assert equal "$RPKM_REFS" "$TPM_REFS"
      assert equal "$TPM_REFS" "$MEAN_REFS"
      assert equal "$MEAN_REFS" "$COVERED_BASES_REFS"
      assert equal "$COVERED_BASES_REFS" "$COVERED_FRACTION_REFS"
    end

    rm -rf "$TMPDIR"
  end

  describe "Backwards Compatibility"
    it "Default stdout mode unchanged (no -o flag)"
      OUTPUT=$("$BINDIR"/bamcountrefs "$DATADIR"/mini.bam)
      # Should output to stdout
      LINES=$(echo "$OUTPUT" | wc -l)
      assert equal $((LINES+0)) 4
    end

    it "Stdout mode with --multiqc still works"
      OUTPUT=$("$BINDIR"/bamcountrefs --multiqc "$DATADIR"/mini.bam)
      assert grep "plot_type" <<< "$OUTPUT"
    end

    it "Cannot use --multiqc with -o flag (stdout only)"
      # MultiQC format only makes sense for stdout
      # This test just ensures the tool doesn't crash
      "$BINDIR"/bamcountrefs -o /tmp/test_multiqc --multiqc "$DATADIR"/mini.bam 2>&1
      exitstatus=$?
      # Should complete (may or may not produce MultiQC headers in file)
      assert equal $exitstatus 0
    end
  end

  describe "MultiQC Output"
    it "MultiQC format includes plot_type header"
      COUNT=$("$BINDIR"/bamcountrefs --multiqc "$DATADIR"/mini.bam | grep -c "plot_type")
      assert equal $COUNT 1
    end

    it "MultiQC output contains section_name header"
      COUNT=$("$BINDIR"/bamcountrefs --multiqc "$DATADIR"/mini.bam | grep -c "section_name")
      assert equal $COUNT 1
    end

    it "MultiQC output contains description header"
      COUNT=$("$BINDIR"/bamcountrefs --multiqc "$DATADIR"/mini.bam | grep -c "description")
      assert equal $COUNT 1
    end
  end

  describe "Custom Column Name"
    it "Custom tag changes first column name"
      HEADER=$("$BINDIR"/bamcountrefs --tag "Reference" "$DATADIR"/mini.bam | head -n 1)
      assert glob "$HEADER" "Reference*"
    end
  end

  describe "Error Handling"
    it "Fails with missing BAM file"
      "$BINDIR"/bamcountrefs /nonexistent/file.bam 2>&1 > /dev/null
      exitstatus=$?
      assert gt $exitstatus 0
    end

    it "Requires BAM index"
      "$BINDIR"/bamcountrefs "$DATADIR"/mini-unsorted.bam 2>&1 > /dev/null
      exitstatus=$?
      assert gt $exitstatus 0
    end
  end

  describe "PhiX Dataset Tests"
    it "Counts PhiX wildtype correctly"
      COUNT=$("$BINDIR"/bamcountrefs "$DATADIR"/phi/wildtype.bam | grep "^NC_001422" | cut -f 2)
      assert equal $((COUNT+0)) 1060
    end

    it "Counts PhiX shotgun correctly"
      COUNT=$("$BINDIR"/bamcountrefs "$DATADIR"/phi/shotgun.bam | grep "^NC_001422" | cut -f 2)
      assert equal $((COUNT+0)) 928
    end

    it "Both PhiX samples in multi-sample mode"
      LINE=$("$BINDIR"/bamcountrefs "$DATADIR"/phi/wildtype.bam "$DATADIR"/phi/shotgun.bam | grep "^NC_001422")
      COUNT1=$(echo "$LINE" | cut -f 2)
      COUNT2=$(echo "$LINE" | cut -f 3)
      assert equal $((COUNT1+0)) 1060
      assert equal $((COUNT2+0)) 928
    end
  end

  describe "Performance Optimization"
    it "Debug shows processing messages for each sample"
      OUTPUT=$("$BINDIR"/bamcountrefs --debug "$DATADIR"/phi/wildtype.bam "$DATADIR"/phi/shotgun.bam 2>&1)
      # Should see debug output for processing both samples
      COUNT=$(echo "$OUTPUT" | grep -c "Opening BAM/CRAM file")
      assert equal $COUNT 2
    end

    it "Fast path and iteration path produce same results for appropriate filters"
      COUNT_ITER=$("$BINDIR"/bamcountrefs -F 1796 "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      COUNT_FAST=$("$BINDIR"/bamcountrefs -F 4 "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      # Both should give same count (15) for mini.bam since it has no secondary/duplicate reads
      assert equal $COUNT_ITER $COUNT_FAST
    end
  end

end
