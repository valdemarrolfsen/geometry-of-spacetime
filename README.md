# The Geometry of Spacetime (LaTeX Manuscript Repository)

This repository contains a scientific-manuscript style LaTeX project for the article:

**The Geometry of Spacetime: A Temporal-Radial Hyperspherical Model of Cosmological Expansion**.

## Structure

- `main.tex` — top-level LaTeX file.
- `sections/` — sectioned manuscript content.
- `references.bib` — BibTeX bibliography entries.
- `.gitignore` — common LaTeX build artifacts.

## Build

### Option A: latexmk (recommended)

```bash
latexmk -pdf main.tex
```

### Option B: pdflatex + bibtex

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Notes

The manuscript clearly separates interpretive geometry from established cosmological constraints and keeps speculative extensions explicitly labeled as speculative.
