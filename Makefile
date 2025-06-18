# Create "make test"
.PHONY: testshpec testshunit2 test clean build


# Auto-detect htslib path with proper rpath handling
HTSLIB_PATH := $(shell \
	if [ -f /usr/lib/x86_64-linux-gnu/libhts.so ]; then \
		echo "-L/usr/lib/x86_64-linux-gnu -lhts"; \
	elif [ -f /opt/homebrew/lib/libhts.dylib ]; then \
		echo "-L/opt/homebrew/lib -lhts -rpath /opt/homebrew/lib"; \
	elif [ -f /usr/local/lib/libhts.dylib ]; then \
		echo "-L/usr/local/lib -lhts -rpath /usr/local/lib"; \
	else \
		echo "-lhts"; \
	fi)

BIN=./bin
BIN2=./static
SOURCE=./src
LIBPATH=$(shell dirname $(shell which samtools) | sed 's/bin/lib/')
NIM_PATHS= --colors:on 
VERSION := $(shell grep version bamtocov.nimble  | grep  -o "[0-9]\\+\.[0-9]\.[0-9]\\+")
LIST=$(BIN)/bamtocov $(BIN)/bamtocounts $(BIN)/covtotarget $(BIN)/bamcountrefs $(BIN)/gff2bed $(BIN)/bamtocounts_legacy $(BIN)/bamtarget
STATIC=$(BIN2)/bamtocov $(BIN2)/bamtocounts $(BIN2)/covtotarget $(BIN2)/bamcountrefs $(BIN2)/gff2bed $(BIN2)/bamtocounts_legacy $(BIN2)/bamtarget

$(BIN)/%: $(SOURCE)/%.nim $(SOURCE)/covutils.nim bamtocov.nimble
	nim c  $(NIM_PATHS) -d:NimblePkgVersion=$(VERSION) -d:release \
		--passL:"$(HTSLIB_PATH)" \
		--opt:speed  --out:$@ $<

all: $(LIST)

./hts_nim_static_builder:
	wget "https://github.com/brentp/hts-nim/releases/download/v0.2.8/hts_nim_static_builder"
	chmod +x hts_nim_static_builder

static: $(STATIC)


$(BIN2)/%: $(SOURCE)/%.nim bamtocov.nimble ./hts_nim_static_builder
	mkdir -p $(BIN2)
	echo target: $@
	echo base: $(shell basename $@)
	./hts_nim_static_builder  -n bamtocov.nimble -s $< -- -d:NimblePkgVersion=$(VERSION)
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
