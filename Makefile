PYTHON ?= python3
TECTONIC ?= tectonic
export PYTHONPATH := $(CURDIR)/src
export MPLCONFIGDIR := $(CURDIR)/.mplconfig
export TECTONIC_CACHE_DIR := $(CURDIR)/.tectonic-cache

.PHONY: all check results paper clean
all: check paper

check:
	$(PYTHON) -m pytest -q
	$(PYTHON) scripts/symbolic_checks.py

results:
	$(PYTHON) scripts/reproduce.py

paper: results
	mkdir -p build
	cd paper && $(TECTONIC) -o ../build main.tex
	mkdir -p output/pdf
	cp build/main.pdf output/pdf/geometry-of-spacetime-v0.5.pdf

clean:
	rm -rf build .pytest_cache .mplconfig
