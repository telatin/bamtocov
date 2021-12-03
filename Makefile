# Create "make test"
.PHONY: testshpec testshunit2 test clean build

BIN=./bin
SOURCE=./src
VERSION := $(shell grep version bamtocov.nimble  | grep  -o "[0-9]\\+\.[0-9]\.[0-9]\\+")
LIST=$(BIN)/bamtocov $(BIN)/bamtocounts $(BIN)/covtotarget $(BIN)/bamcountrefs $(BIN)/gff2bed

$(BIN)/%: $(SOURCE)/%.nim
	nim c -d:NimblePkgVersion=$(VERSION) -d:release -d:danger --opt:speed --gc:arc --out:$@ $<

all: $(LIST)

bin/bamtocov:
	nimble build

test: testshpec testshunit2

testshpec:
	@echo "Test shpec"
	./tests/bin/shpec ./tests/shpec/bamtocov.sh

testshunit2:
	@echo "Test shunit2"
	./tests/unit/bamtocov-base.sh
	
testbash:
	bash tests/all.sh

build:
	nimble build


clean:
	@echo "Cleaning..."
	@for i in $(LIST); \
	do \
		if [ -e "$$i" ]; then rm -f $$i; echo "Removing $$i"; else echo "$$i Not found"; fi \
	done