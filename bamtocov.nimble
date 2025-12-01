# Package

version       = "2.9.0"
author        = "Andrea Telatin, Giovanni Birolo"
description   = "BAM to Coverage"
license       = "MIT"

# Dependencies
when defined(macosx):
  let htslib_prefix = gorge("brew --prefix htslib").strip()
  --passL:"-L" & htslib_prefix & "/lib"
  --passL:"-lhts"
  --passL:"-rpath " & htslib_prefix & "/lib"

requires "hts >= 0.3.20", "docopt >= 0.6.8", "nim >= 1.6.6", "lapper >= 0.1.8"

srcDir = "src"
binDir = "bin"
bin = @["covtotarget", "bamtocounts", "bamtocov", "bamcountrefs", "gff2bed", "bamtarget"]
skipDirs = @["tests", "docs"]
skipFiles = @["example.bam"]


