# Create "make test"
.PHONY: testshpec testshunit2 test clean build

BIN=./bin
BIN2=./static
SOURCE=./src
LIBPATH=$(shell dirname $(shell which samtools) | sed 's/bin/lib/')
NIM_PATHS= --colors:on --noNimblePath --path:/local/giovanni/bamtocov/nimbledeps/pkgs/hts-0.3.17 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/docopt-0.6.8 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/regex-0.19.0 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/unicodedb-0.9.0 --path:/local/giovanni/bamtocov/nimbledeps/pkgs/lapper-0.1.7 
VERSION := $(shell grep version bamtocov.nimble  | grep  -o "[0-9]\\+\.[0-9]\.[0-9]\\+")
LIST=$(BIN)/bamtocov $(BIN)/bamtocounts $(BIN)/covtotarget $(BIN)/bamcountrefs $(BIN)/gff2bed $(BIN)/bamtocounts_legacy $(BIN)/bamtarget
STATIC=$(BIN2)/bamtocov $(BIN2)/bamtocounts $(BIN2)/covtotarget $(BIN2)/bamcountrefs $(BIN2)/gff2bed $(BIN2)/bamtocounts_legacy $(BIN2)/bamtarget


$(BIN)/%: $(SOURCE)/%.nim $(SOURCE)/covutils.nim bamtocov.nimble
	nim c  $(NIM_PATHS) -d:NimblePkgVersion=$(VERSION) -d:release --opt:speed  --out:$@ $<

all: $(LIST)

static: $(STATIC)

$(BIN2)/%: $(SOURCE)/%.nim bamtocov.nimble
	mkdir -p $(BIN2)
	echo target: $@
	echo base: $(shell basename $@)
	./hts_nim_static_builder  -n bamtocov.nimble -s $< -- -d:NimblePkgVersion=2.7.0
	mv ./$(shell basename $@) $(BIN2)
test: testshpec

testall: testbash testshpec testshunit2

testshpec:
	@echo " --- Test shpec --- "
	./tests/bin/shpec ./tests/shpec/bamtocov.sh

testshunit2:
	@echo " --- Test shunit2 --- "
	./tests/unit/bamtocov-base.sh

testbash:
	@echo " --- Test (legacy) --- "
	bash tests/all.sh


build:
	nimble build


clean:
	@echo "Cleaning..."
	@for i in $(LIST); \
	do \
		if [ -e "$$i" ]; then rm -f $$i; echo "Removing $$i"; else echo "$$i Not found"; fi \
	done
