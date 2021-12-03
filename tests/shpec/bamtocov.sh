describe "Coverage tools"

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
  describe "BamToCov"
    it "Binary exist"
        assert file_present "$BINDIR"/bamtocov
    end
    it "Version 2.x"
      VERSION=$("$BINDIR"/bamtocov --version)
      assert glob "$VERSION" "2.*"
    end

    it "Mini coverage lines"
      COV=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam | wc -l)
      assert equal $((COV+0)) 21
    end

    it "One line with empty chromosome seq0"
      COV=$("$BINDIR"/bamtocov "$DATADIR"/mini.bam | grep -c seq0)
      assert equal $((COV+0)) 1
    end

  end

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
        COV=$(echo "$line" | cut -f 5)
        REGION=$(echo "$line" | cut -f 4)
        assert glob "$REGION" "*_$COV"
      done < "$TMPFILE"
      rm "$TMPFILE"
    end 
  end
end
