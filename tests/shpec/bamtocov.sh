describe "Coverage tools tested by Shpec"

  if [ "$0" = "./tests/bin/shpec" ];
  then
    echo "Make test"
    BINDIR="./bin/"
    DATADIR="./input/"
  else
    echo "Manual test $0"
    SELFDIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
    BINDIR="$( cd -- "$( dirname -- "$SELFDIR"/../../bin/bamtocov )" &> /dev/null && pwd )"
    DATADIR="$( cd -- "$( dirname -- "$SELFDIR"/../../input/mini.bam )" &> /dev/null && pwd )"
  fi

  it "Binary found"
    assert file_present "$BINDIR/bamtocov"
    assert file_present "$BINDIR/bamtocounts"
  end
  # PROGRAM: bamtocov
  describe "BamToCov"
    it "Binary exist"
        assert file_present "$BINDIR"/bamtocov
    end
    VERSION=$("$BINDIR"/bamtocov --version)
    it "Version emitted is 2.x [$VERSION]"
      
      assert glob "$VERSION" "2.*"
    end

    it "Mini coverage, verify output line number"
      COV=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam | wc -l)
      assert equal $((COV+0)) 21
    end

    it "One line with empty chromosome seq0 (empty chromosome)"
      COV=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam | grep -c seq0)
      assert equal $((COV+0)) 1
    end

    it "Target inputs"
      assert file_present "$DATADIR"/regions.gtf
      assert file_present "$DATADIR"/regions.bed
      assert file_present "$DATADIR"/mini.bam 

      TMPFILE=$(mktemp)
      "$BINDIR"/bamtocov -o $TMPFILE --regions "$DATADIR"/regions.bed "$DATADIR"/mini.bam 2> /dev/null > /dev/null
      it "Target BED: Execution"
        exitstatus=$?
        assert equal $exitstatus 0
      end     
      it "Target BED: Report"
        assert file_present $TMPFILE
      end    
      
      it "Target BED: Report file ($TMPFILE)"
        LINES=$(cat "$TMPFILE" | wc -l)
        STRING=$(cat "$TMPFILE")
        assert equal $LINES 7
        assert glob  "${STRING}" "interval*"
        rm $TMPFILE
      end

      TMPFILE=$(mktemp)
      "$BINDIR"/bamtocov -o "$TMPFILE" --regions "$DATADIR"/regions.bed "$DATADIR"/mini.bam > /dev/null 2> /dev/null
      REPORT_ORDER=$(cut -f 1 "$TMPFILE" | tail -n +2 | paste -sd "," -)
      rm -f "$TMPFILE"
      it "Target BED: Report preserves target order"
        assert equal "$REPORT_ORDER" "include_5,overlap_2,empty1_0,shared1_10,shared2_10,null_0"
      end

    end

    it "Target GTF"
      TMPFILE=$(mktemp)
      "$BINDIR"/bamtocov -o $TMPFILE --regions "$DATADIR"/regions.gtf "$DATADIR"/mini.bam 2> /dev/null > /dev/null
      it "Target GTF: execution"
        exitstatus=$?
        assert equal $exitstatus 0
      end
      it "Target GTF: report found  ($TMPFILE)"
        assert file_present $TMPFILE
      end
      # it "Target GTF: report file"
      #   LINES=$(cat "$TMPFILE" | wc -l)
      #   assert equal $LINES 4
      #   rm $TMPFILE 
      # end
    end

    it "Fails when report path is not writable as a file"
      TMP_REPORT_DIR=$(mktemp -d)
      "$BINDIR"/bamtocov -o "$TMP_REPORT_DIR" --regions "$DATADIR"/regions.bed "$DATADIR"/mini.bam > /dev/null 2>&1
      exitstatus=$?
      rm -rf "$TMP_REPORT_DIR"
      assert gt $exitstatus 0
    end

    it "Fails with invalid --op value"
      TMPERR=$(mktemp)
      "$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 250 --op median > /dev/null 2> "$TMPERR"
      exitstatus=$?
      OUTPUT=$(cat "$TMPERR")
      rm -f "$TMPERR"
      assert gt $exitstatus 0
      assert equal "$OUTPUT" "ERROR: --op must be one of mean, min, max; got: median"
    end

    it "Fails with target contig missing from BAM header"
      TMPTARGET=$(mktemp)
      TMPERR=$(mktemp)
      printf "missing_chr\t0\t10\n" > "$TMPTARGET"
      "$BINDIR"/bamtocov --regions "$TMPTARGET" "$DATADIR"/mini.bam > /dev/null 2> "$TMPERR"
      exitstatus=$?
      OUTPUT=$(cat "$TMPERR")
      rm -f "$TMPTARGET" "$TMPERR"
      assert equal $exitstatus 1
      assert equal "$OUTPUT" "ERROR: Target contig not found in BAM/CRAM header: missing_chr"
    end

    it "Produces stranded quantized BED without empty columns"
      FIRSTLINE=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --stranded --quantize 1,5 | head -n 1)
      FIELDS=$(printf "%s\n" "$FIRSTLINE" | awk -F '\t' '{print NF}')
      assert equal "$FIRSTLINE" "seq1	0	9	0-0	0-0"
      assert equal $FIELDS 5
    end

    it "Works with sorted file"
      "$BINDIR"/bamtocov "$DATADIR"/mini-sorted.bam 2> /dev/null > /dev/null
      exitstatus=$?
      assert equal $exitstatus 0
    end

    it "Fails with unsorted file"
      "$BINDIR"/bamtocov "$DATADIR"/mini-unsorted.bam 2> /dev/null > /dev/null
      exitstatus=$?
      assert gt $exitstatus 0
    end

    it "Produces wig output (check lines: 15)"
      LINES=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 250 --op max | wc -l)
      assert equal $LINES 15
    end
    it "Produces wig output (check coverage 10X, twice)"
      LINES=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 250 --op max | grep -w 10 | wc -l)
      assert equal $LINES 2
    end
    it "Produces wig output (check coverage 1000X, unexpected)"
      LINES=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 250 --op max |  grep -w 1000 |wc -l)
      assert equal $LINES 0
    end
    it "Produces wig output header"
      LINES=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 250 --op max | grep "fixed" | wc -l)
      assert equal $LINES 3
    end

    it "Does not leak stranded WIG state across contigs"
      SEQ2_LAST=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 550 --op max --stranded | awk '/^fixedStep chrom=seq2 / {getline; getline; print; exit}')
      SEQ0_FIRST=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 550 --op max --stranded | awk '/^fixedStep chrom=seq0 / {getline; print; exit}')
      assert equal "$SEQ2_LAST" "$(printf '5\t5')"
      assert equal "$SEQ0_FIRST" "$(printf '0\t0')"
    end

    it "Uses the actual width for the final WIG mean bin"
      EXPECTED=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam | awk -F '\t' '$1=="seq2" {start=($2<550?550:$2); stop=($3>1000?1000:$3); if (stop>start) sum += (stop-start)*$4} END {printf "%.15f", sum/450}')
      ACTUAL=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 550 --op mean | awk '/^fixedStep chrom=seq2 / {getline; getline; printf "%.15f", $1; exit}')
      assert equal "$ACTUAL" "$EXPECTED"
    end
  end

  # PROGRAM: bamtocount
  describe "BamToCounts"
    it "Binary exists"
        assert file_present "$BINDIR"/bamtocounts
    end
    it "Version 2.x"
      VERSION=$("$BINDIR"/bamtocounts --version)
      assert glob "$VERSION" "2.*"
    end

    it "Counts target"    
      COV=$("$BINDIR"/bamtocounts "$DATADIR"/regions.bed "$DATADIR"/mini.bam | wc -l)
      assert equal $((COV+0)) 6
    end

    it "Coverage check"    
      TMPFILE=$(mktemp)
      "$BINDIR"/bamtocounts "$DATADIR"/regions.bed "$DATADIR"/mini.bam > "$TMPFILE"
      while read -r line; do
        COV=$(echo "$line" | cut -f 2)
        REGION=$(echo "$line" | cut -f 1)
        assert glob "$REGION" "*_$COV"
      done < "$TMPFILE"
      rm "$TMPFILE"
    end 
  end
end
