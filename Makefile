PYTHON ?= python3
TECTONIC ?= tectonic
export PYTHONPATH := $(CURDIR)/src
export MPLCONFIGDIR := $(CURDIR)/.mplconfig
export TECTONIC_CACHE_DIR := $(CURDIR)/.tectonic-cache

.PHONY: all check results paper comparison clean
all: check paper comparison

check:
	$(PYTHON) -m pytest -q
	$(PYTHON) scripts/symbolic_checks.py

results:
	$(PYTHON) scripts/reproduce.py
	$(PYTHON) scripts/compare_diagrams.py

paper: results
	mkdir -p build
	cd paper && $(TECTONIC) -o ../build main.tex
	mkdir -p output/pdf
	cp build/main.pdf output/pdf/geometry-of-spacetime-v0.5.pdf

comparison: results
	mkdir -p build output/pdf
	cd paper && $(TECTONIC) -o ../build diagram-comparison.tex
	cp build/diagram-comparison.pdf output/pdf/diagram-comparison.pdf

clean:
	rm -rf build .pytest_cache .mplconfig
