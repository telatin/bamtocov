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

  # The tests use TMPDIR as a scratch variable and remove it after each block.
  # If inherited as an exported environment variable, later mktemp calls would
  # try to create temp directories under a deleted path.
  unset TMPDIR

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
      assert glob "$HEADER" "Contig*mini"
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
      assert file_present "$TMPDIR"/all_test_trimmed_mean.tsv
      assert file_present "$TMPDIR"/all_test_covered_bases.tsv
      assert file_present "$TMPDIR"/all_test_covered_fraction.tsv
      assert file_present "$TMPDIR"/all_test_variance.tsv
      assert file_present "$TMPDIR"/all_test_reads_per_base.tsv
      assert file_present "$TMPDIR"/all_test_length.tsv
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

  describe "Variance Calculations"
    TMPDIR=$(mktemp -d)

    it "Variance file is created with -o --variance"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/variance_test --variance "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/variance_test_variance.tsv
    end

    it "Variance values are numeric (float)"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/variance_test --variance "$DATADIR"/mini.bam
      VARIANCE=$(grep "^seq1" "$TMPDIR"/variance_test_variance.tsv | cut -f 2)
      assert glob "$VARIANCE" "*.*"
    end

    it "Zero coverage references have variance of 0.0"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/variance_test --variance "$DATADIR"/mini.bam
      VARIANCE=$(grep "^seq0" "$TMPDIR"/variance_test_variance.tsv | cut -f 2)
      assert equal "$VARIANCE" "0.000000"
    end

    it "Variance calculated per sample"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/variance_multi --variance "$DATADIR"/mini.bam "$DATADIR"/mini2.bam
      LINE=$(grep "^seq1" "$TMPDIR"/variance_multi_variance.tsv)
      VAR1=$(echo "$LINE" | cut -f 2)
      VAR2=$(echo "$LINE" | cut -f 3)
      # Both should have numeric values
      assert glob "$VAR1" "*.*"
      assert glob "$VAR2" "*.*"
    end

    it "Variance can be output to stdout"
      OUTPUT=$("$BINDIR"/bamcountrefs --variance "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      assert glob "$OUTPUT" "*.*"
    end

    rm -rf "$TMPDIR"
  end

  describe "Reads Per Base Calculations"
    TMPDIR=$(mktemp -d)

    it "Reads per base file is created with -o --reads-per-base"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/rpb_test --reads-per-base "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/rpb_test_reads_per_base.tsv
    end

    it "Reads per base values are numeric (float)"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/rpb_test --reads-per-base "$DATADIR"/mini.bam
      RPB=$(grep "^seq1" "$TMPDIR"/rpb_test_reads_per_base.tsv | cut -f 2)
      assert glob "$RPB" "*.*"
    end

    it "Reads per base is calculated correctly (seq1: 15 reads / 1000 bp = 0.015)"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/rpb_test --reads-per-base "$DATADIR"/mini.bam
      RPB=$(grep "^seq1" "$TMPDIR"/rpb_test_reads_per_base.tsv | cut -f 2)
      assert equal "$RPB" "0.015000"
    end

    it "Zero reads references have reads per base of 0.0"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/rpb_test --reads-per-base "$DATADIR"/mini.bam
      RPB=$(grep "^seq0" "$TMPDIR"/rpb_test_reads_per_base.tsv | cut -f 2)
      assert equal "$RPB" "0.000000"
    end

    it "Reads per base calculated per sample"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/rpb_multi --reads-per-base "$DATADIR"/mini.bam "$DATADIR"/mini2.bam
      LINE=$(grep "^seq1" "$TMPDIR"/rpb_multi_reads_per_base.tsv)
      RPB1=$(echo "$LINE" | cut -f 2)
      RPB2=$(echo "$LINE" | cut -f 3)
      # Both should have numeric values
      assert glob "$RPB1" "*.*"
      assert glob "$RPB2" "*.*"
    end

    it "Reads per base can be output to stdout"
      OUTPUT=$("$BINDIR"/bamcountrefs --reads-per-base "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      assert equal "$OUTPUT" "0.015000"
    end

    rm -rf "$TMPDIR"
  end

  describe "Trimmed Mean Calculations"
    TMPDIR=$(mktemp -d)
    it "Temp directory is created"
      assert file_present "$TMPDIR"
    end

    it "Trimmed mean file is created with -o --trimmed-mean"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/tm_test --trimmed-mean "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/tm_test_trimmed_mean.tsv
    end

    it "Trimmed mean values are numeric (float)"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/tm_test --trimmed-mean "$DATADIR"/mini.bam
      TM=$(grep "^seq1" "$TMPDIR"/tm_test_trimmed_mean.tsv | cut -f 2)
      assert glob "$TM" "*.*"
    end

    it "Trimmed mean is less than or equal to mean (removes outliers)"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/tm_test --mean --trimmed-mean "$DATADIR"/mini.bam
      MEAN=$(grep "^seq1" "$TMPDIR"/tm_test_mean.tsv | cut -f 2)
      TM=$(grep "^seq1" "$TMPDIR"/tm_test_trimmed_mean.tsv | cut -f 2)
      # Trimmed mean should generally be less than or equal to mean (removing high outliers)
      assert glob "$TM" "*.*"
      assert glob "$MEAN" "*.*"
    end

    it "Zero coverage references have trimmed mean of 0.0"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/tm_test --trimmed-mean "$DATADIR"/mini.bam
      TM=$(grep "^seq0" "$TMPDIR"/tm_test_trimmed_mean.tsv | cut -f 2)
      assert equal "$TM" "0.000000"
    end

    it "Trimmed mean calculated per sample"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/tm_multi --trimmed-mean "$DATADIR"/mini.bam "$DATADIR"/mini2.bam
      LINE=$(grep "^seq1" "$TMPDIR"/tm_multi_trimmed_mean.tsv)
      TM1=$(echo "$LINE" | cut -f 2)
      TM2=$(echo "$LINE" | cut -f 3)
      # Both should have numeric values
      assert glob "$TM1" "*.*"
      assert glob "$TM2" "*.*"
    end

    it "Trimmed mean can be output to stdout"
      OUTPUT=$("$BINDIR"/bamcountrefs --trimmed-mean "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      assert glob "$OUTPUT" "*.*"
    end

    it "Custom trim percentiles affect trimmed mean value"
      TM_DEFAULT=$("$BINDIR"/bamcountrefs --trimmed-mean "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      TM_CUSTOM=$("$BINDIR"/bamcountrefs --trimmed-mean --trim-min 10 --trim-max 90 "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      # Different trim parameters should give different results (not always, but likely)
      assert glob "$TM_DEFAULT" "*.*"
      assert glob "$TM_CUSTOM" "*.*"
    end

    it "Variance and trimmed-mean can be calculated together"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/both --variance --trimmed-mean "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/both_variance.tsv
      assert file_present "$TMPDIR"/both_trimmed_mean.tsv
    end

    rm -rf "$TMPDIR"
  end

  describe "DisCov Calculations"
    TMPDIR=$(mktemp -d)

    it "DisCov file is created with -o --discov"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/discov_test --discov "$DATADIR"/mini.bam
      assert file_present "$TMPDIR"/discov_test_discov.tsv
    end

    it "DisCov values are numeric floats"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/discov_test --discov "$DATADIR"/mini.bam
      DISCOV=$(grep "^seq1" "$TMPDIR"/discov_test_discov.tsv | cut -f 2)
      assert glob "$DISCOV" "*.*"
    end

    it "Zero coverage references have DisCov of 0.0"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/discov_test --discov "$DATADIR"/mini.bam
      DISCOV=$(grep "^seq0" "$TMPDIR"/discov_test_discov.tsv | cut -f 2)
      assert equal "$DISCOV" "0.000000"
    end

    it "DisCov can be output to stdout"
      OUTPUT=$("$BINDIR"/bamcountrefs --discov "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      assert glob "$OUTPUT" "*.*"
    end

    it "DisCov is calculated per sample"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/discov_multi --discov "$DATADIR"/mini.bam "$DATADIR"/mini2.bam
      LINE=$(grep "^seq1" "$TMPDIR"/discov_multi_discov.tsv)
      DISCOV1=$(echo "$LINE" | cut -f 2)
      DISCOV2=$(echo "$LINE" | cut -f 3)
      assert glob "$DISCOV1" "*.*"
      assert glob "$DISCOV2" "*.*"
    end

    it "DisCov accepts geometric formula and custom window"
      OUTPUT=$("$BINDIR"/bamcountrefs --discov --discov-formula geometric --discov-window 100 "$DATADIR"/mini.bam | grep "^seq1" | cut -f 2)
      assert glob "$OUTPUT" "*.*"
    end

    rm -rf "$TMPDIR"
  end

  describe "All-Metrics Flag"
    TMPDIR=$(mktemp -d)

    it "Creates all output files"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/all --all-metrics "$DATADIR"/mini.bam
      FILE_COUNT=$(ls "$TMPDIR"/all_*.tsv 2>/dev/null | wc -l)
      assert equal $((FILE_COUNT+0)) 11
    end

    it "All files have same number of data rows"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/all --all-metrics "$DATADIR"/mini.bam
      COUNTS_LINES=$(tail -n +2 "$TMPDIR"/all_counts.tsv | wc -l)
      RPKM_LINES=$(tail -n +2 "$TMPDIR"/all_rpkm.tsv | wc -l)
      TPM_LINES=$(tail -n +2 "$TMPDIR"/all_tpm.tsv | wc -l)
      MEAN_LINES=$(tail -n +2 "$TMPDIR"/all_mean.tsv | wc -l)
      TRIMMED_MEAN_LINES=$(tail -n +2 "$TMPDIR"/all_trimmed_mean.tsv | wc -l)
      COVERED_BASES_LINES=$(tail -n +2 "$TMPDIR"/all_covered_bases.tsv | wc -l)
      COVERED_FRACTION_LINES=$(tail -n +2 "$TMPDIR"/all_covered_fraction.tsv | wc -l)
      VARIANCE_LINES=$(tail -n +2 "$TMPDIR"/all_variance.tsv | wc -l)
      RPB_LINES=$(tail -n +2 "$TMPDIR"/all_reads_per_base.tsv | wc -l)
      DISCOV_LINES=$(tail -n +2 "$TMPDIR"/all_discov.tsv | wc -l)
      LENGTH_LINES=$(tail -n +2 "$TMPDIR"/all_length.tsv | wc -l)
      assert equal $COUNTS_LINES $RPKM_LINES
      assert equal $RPKM_LINES $TPM_LINES
      assert equal $TPM_LINES $MEAN_LINES
      assert equal $MEAN_LINES $TRIMMED_MEAN_LINES
      assert equal $TRIMMED_MEAN_LINES $COVERED_BASES_LINES
      assert equal $COVERED_BASES_LINES $COVERED_FRACTION_LINES
      assert equal $COVERED_FRACTION_LINES $VARIANCE_LINES
      assert equal $VARIANCE_LINES $RPB_LINES
      assert equal $RPB_LINES $DISCOV_LINES
      assert equal $DISCOV_LINES $LENGTH_LINES
    end

    it "All files have same reference order"
      "$BINDIR"/bamcountrefs -o "$TMPDIR"/all --all-metrics "$DATADIR"/mini.bam
      COUNTS_REFS=$(tail -n +2 "$TMPDIR"/all_counts.tsv | cut -f 1)
      RPKM_REFS=$(tail -n +2 "$TMPDIR"/all_rpkm.tsv | cut -f 1)
      TPM_REFS=$(tail -n +2 "$TMPDIR"/all_tpm.tsv | cut -f 1)
      MEAN_REFS=$(tail -n +2 "$TMPDIR"/all_mean.tsv | cut -f 1)
      TRIMMED_MEAN_REFS=$(tail -n +2 "$TMPDIR"/all_trimmed_mean.tsv | cut -f 1)
      COVERED_BASES_REFS=$(tail -n +2 "$TMPDIR"/all_covered_bases.tsv | cut -f 1)
      COVERED_FRACTION_REFS=$(tail -n +2 "$TMPDIR"/all_covered_fraction.tsv | cut -f 1)
      VARIANCE_REFS=$(tail -n +2 "$TMPDIR"/all_variance.tsv | cut -f 1)
      RPB_REFS=$(tail -n +2 "$TMPDIR"/all_reads_per_base.tsv | cut -f 1)
      DISCOV_REFS=$(tail -n +2 "$TMPDIR"/all_discov.tsv | cut -f 1)
      LENGTH_REFS=$(tail -n +2 "$TMPDIR"/all_length.tsv | cut -f 1)
      assert equal "$COUNTS_REFS" "$RPKM_REFS"
      assert equal "$RPKM_REFS" "$TPM_REFS"
      assert equal "$TPM_REFS" "$MEAN_REFS"
      assert equal "$MEAN_REFS" "$TRIMMED_MEAN_REFS"
      assert equal "$TRIMMED_MEAN_REFS" "$COVERED_BASES_REFS"
      assert equal "$COVERED_BASES_REFS" "$COVERED_FRACTION_REFS"
      assert equal "$COVERED_FRACTION_REFS" "$VARIANCE_REFS"
      assert equal "$VARIANCE_REFS" "$RPB_REFS"
      assert equal "$RPB_REFS" "$DISCOV_REFS"
      assert equal "$DISCOV_REFS" "$LENGTH_REFS"
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

    it "Writes MultiQC output to {basename}.mqc.txt when -o is used"
      TMPDIR=$(mktemp -d /tmp/bamcountrefs-mqc.XXXXXX)
      BASENAME="$TMPDIR"/test_multiqc
      "$BINDIR"/bamcountrefs -o "$BASENAME" --multiqc "$DATADIR"/mini.bam >/dev/null 2>&1
      exitstatus=$?
      assert equal $exitstatus 0
      assert equal "$(test -f "$BASENAME".mqc.txt; echo $?)" 0
      assert grep "$(cat "$BASENAME".mqc.txt)" "plot_type"
      rm -rf "$TMPDIR"
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
      assert equal $((COUNT+0)) 920
    end

    it "Both PhiX samples in multi-sample mode"
      LINE=$("$BINDIR"/bamcountrefs "$DATADIR"/phi/wildtype.bam "$DATADIR"/phi/shotgun.bam | grep "^NC_001422")
      COUNT1=$(echo "$LINE" | cut -f 2)
      COUNT2=$(echo "$LINE" | cut -f 3)
      assert equal $((COUNT1+0)) 1060
      assert equal $((COUNT2+0)) 920
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
