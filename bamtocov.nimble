# Package

version       = "2.3.0"
author        = "Andrea Telatin, Giovanni Birolo"
description   = "BAM to Coverage"
license       = "MIT"

# Dependencies

requires "hts >= 0.3.1", "docopt >= 0.6.8", "nim >= 1.0.0", "lapper"
srcDir = "src"
binDir = "bin"
bin = @["covtotarget", "bamtocounts", "bamtocov", "bamcountrefs", "gff2bed"]
skipDirs = @["tests", "docs"]
skipFiles = @["example.bam"]


