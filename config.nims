import std/os

when defined(macosx):
  # Try homebrew locations in order
  const htslibPaths = [
    "/opt/homebrew/lib",
    "/usr/local/lib"
  ]
  for path in htslibPaths:
    if dirExists(path):
      switch("passL", "-L" & path)
      switch("passL", "-lhts")
      switch("passL", "-rpath " & path)
      # Force ARM64 build if running on Apple Silicon with x86_64 Nim
      switch("passC", "-arch arm64")
      switch("passL", "-arch arm64")
      break
elif defined(linux):
  # Use pkg-config when available
  when gorgeEx("pkg-config --exists htslib").exitCode == 0:
    switch("passL", gorge("pkg-config --libs htslib"))
  else:
    switch("passL", "-lhts")
