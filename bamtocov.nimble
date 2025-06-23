# Package

version       = "2.8.0"
author        = "Andrea Telatin, Giovanni Birolo"
description   = "BAM to Coverage"
license       = "MIT"

# Dependencies

when defined(macosx):
  --passL:"-L/opt/homebrew/Cellar/htslib/1.21/lib"
  --passL:"-lhts"
  --passL:"-rpath /opt/homebrew/Cellar/htslib/1.21/lib"
  
requires "hts >= 0.3.1", "docopt >= 0.6.8", "nim >= 1.6.6", "lapper >= 0.1.8"
srcDir = "src"
binDir = "bin"
bin = @["covtotarget", "bamtocounts", "bamtocov", "bamcountrefs", "gff2bed", "bamtarget"]
skipDirs = @["tests", "docs"]
skipFiles = @["example.bam"]


