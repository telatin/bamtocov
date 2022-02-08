set -euxo pipefail
BIN=bin/bamtocov
touch $BIN && rm $BIN
test="bin/bamtocov --version && bin/bamtocov -r input/regions.bed -o ciao input/mini.bam --debug"
nim c  --colors:on --noNimblePath --path:/local/giovanni/bamtocov/nimbledeps/pkgs/hts-0.3.17 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/docopt-0.6.8 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/regex-0.19.0 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/unicodedb-0.9.0 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/lapper-0.1.7   -d:NimblePkgVersion=V1 -d:release  --opt:speed -d:danger --out:bin/bamtocov src/bamtocov.nim 2>&1 > /dev/null

$test

nim c  --path:/local/giovanni/bamtocov/nimbledeps/pkgs/hts-0.3.17 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/docopt-0.6.8 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/regex-0.19.0 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/unicodedb-0.9.0 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/lapper-0.1.7   -d:NimblePkgVersion=V2 -d:release  --opt:speed -d:danger --out:bin/bamtocov src/bamtocov.nim 2>&1 > /dev/null

$test

nim c  -d:NimblePkgVersion=V3 -d:release  --opt:speed -d:danger --out:bin/bamtocov src/bamtocov.nim 2>&1 > /dev/null

$test

nim c  -d:NimblePkgVersion=V4 -d:release  --opt:speed --out:bin/bamtocov src/bamtocov.nim 2>&1 > /dev/null

$test

nim c  -d:NimblePkgVersion=V5 --opt:speed --out:bin/bamtocov src/bamtocov.nim 2>&1 > /dev/null

$test

nim c  -d:NimblePkgVersion=V6 --out:bin/bamtocov src/bamtocov.nim 2>&1 > /dev/null

$test
