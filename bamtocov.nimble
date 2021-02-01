# Package

version       = "2.0.1"
author        = "Andrea Telatin, Giovanni Birolo"
description   = "BAM to Coverage"
license       = "MIT"

# Dependencies

requires "hts >= 0.3.1", "docopt >= 0.6.8", "nim >= 1.0.0", "lapper"

bin = @["src/covtotarget", "src/covtocounts", "src/covtobed"]
skipDirs = @["tests"]
skipFiles = @["example.bam"]


