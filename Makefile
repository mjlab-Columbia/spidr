# Makefile for the SPIDR pipeline
# Builds all custom C++ binaries and (eventually) custom STAR.
# Run from the repo root: make
#
# Requires: g++ (or clang++) with C++17 support, zlib development headers

CXX      := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra
LDFLAGS  := -lz
BIN_DIR  := bin

# ---------------------------------------------------------------------------
# Platform-specific adjustments
# ---------------------------------------------------------------------------
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Darwin)
    # On Apple Silicon, Homebrew lives under /opt/homebrew.
    # On Intel Macs it lives under /usr/local.
    # Try to find Homebrew's zlib; fall back gracefully if brew isn't present.
    BREW_PREFIX := $(shell brew --prefix 2>/dev/null)
    ifneq ($(BREW_PREFIX),)
        ZLIB_PREFIX := $(shell brew --prefix zlib 2>/dev/null)
        ifneq ($(ZLIB_PREFIX),)
            CXXFLAGS += -I$(ZLIB_PREFIX)/include
            LDFLAGS  += -L$(ZLIB_PREFIX)/lib
        else
            # zlib not installed via brew – the macOS SDK ships it
            CXXFLAGS += -I$(BREW_PREFIX)/include
            LDFLAGS  += -L$(BREW_PREFIX)/lib
        endif
    endif
    # Silence unused-parameter warnings from clang's <zlib.h>
    CXXFLAGS += -Wno-unused-parameter
endif

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------
.PHONY: all clean

TARGETS := $(BIN_DIR)/identify_barcode

all: $(BIN_DIR) $(TARGETS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/identify_barcode: scripts/cpp/identify_barcode.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
	@echo "Built $@ successfully."

clean:
	rm -f $(TARGETS)
