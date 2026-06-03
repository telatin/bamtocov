SHELL := /usr/bin/env bash
.DEFAULT_GOAL := all

# Create "make test"
.PHONY: testnim testshpec testshunit2 test clean build static-docker static

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
VERSION := $(shell grep version bamtocov.nimble  | grep  -o "[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+")
LIST=$(BIN)/bamtocov $(BIN)/bamtocounts $(BIN)/covtotarget $(BIN)/bamcountrefs $(BIN)/gff2bed $(BIN)/bamtocounts_legacy $(BIN)/bamtarget
STATIC=$(BIN2)/bamtocov $(BIN2)/bamtocounts $(BIN2)/covtotarget $(BIN2)/bamcountrefs $(BIN2)/gff2bed $(BIN2)/bamtocounts_legacy $(BIN2)/bamtarget

$(BIN)/%: $(SOURCE)/%.nim $(SOURCE)/covutils.nim bamtocov.nimble
	# Standard release binaries are thread-enabled so tools like bamtocounts
	# can use taskpools for per-BAM parallelism.
	nim c  $(NIM_PATHS) -d:NimblePkgVersion=$(VERSION) -d:release \
		--threads:on --mm:arc \
		--passL:"$(HTSLIB_PATH)" \
		--opt:speed  --out:$@ $<

$(BIN)/bamcountrefs: $(SOURCE)/discov.nim

all: $(LIST)

# Docker-based static build — produces fully static musl binaries in ./static/
# Requires only Docker on the host (no local Nim or htslib needed).
#
#   make static-docker               # build with versions pinned in Dockerfile.static
#   make static-docker NIM_VERSION=2.2.4 HTSLIB_VERSION=1.21
#
NIM_VERSION    ?= 2.2.4
HTSLIB_VERSION ?= 1.21
DOCKER_IMAGE   := bamtocov-static-builder

static-docker:
	docker build -f Dockerfile.static \
	  --build-arg NIM_VERSION=$(NIM_VERSION) \
	  --build-arg HTSLIB_VERSION=$(HTSLIB_VERSION) \
	  --build-arg VERSION=$(VERSION) \
	  -t $(DOCKER_IMAGE) .
	mkdir -p $(BIN2)
	docker run --rm -v "$(CURDIR)/$(BIN2):/output" $(DOCKER_IMAGE)
	@echo "Static binaries written to $(BIN2)/"
	@ls -lh $(BIN2)/

# Alias
static: static-docker

test: testnim testshpec

testall: testbash testnim testshpec testshunit2

testnim:
	@echo " --- Test Nim units --- "
	mkdir -p nimcache/discov_test
	nim c -r --colors:on --nimcache:nimcache/discov_test \
		--out:nimcache/discov_test/discov_test tests/unit/discov_test.nim

testshpec:
	@echo " --- Test shpec --- "
	./tests/bin/shpec ./tests/shpec/bamtocov.sh
	./tests/bin/shpec ./tests/shpec/bamcountrefs.sh

testshunit2:
	@echo " --- Test shunit2 --- "
	./tests/unit/bamtocov-base.sh

testbash:
	@echo " --- Test (legacy) --- "
	bash tests/all.sh


build: all


clean:
	@echo "Cleaning..."
	@for i in $(LIST); \
	do \
		if [ -e "$$i" ]; then rm -f $$i; echo "Removing $$i"; else echo "$$i Not found"; fi \
	done
