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

    it "Produces wig output (check lines)"
      LINES=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 250 --op max | wc -l)
      assert equal $LINES 15
    end
    it "Produces wig output (check coverage 10X)"
      LINES=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 250 --op max | grep -w 10 | wc -l)
      assert equal $LINES 2
    end
    it "Produces wig output (check lines at 1000, unexpected)"
      LINES=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 250 --op max |  grep -w 1000 |wc -l)
      assert equal $LINES 0
    end
    it "Produces wig output header"
      LINES=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam --wig 250 --op max | grep "fixed" | wc -l)
      assert equal $LINES 3
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
