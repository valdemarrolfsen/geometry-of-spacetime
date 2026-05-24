# The Geometry of Spacetime (LaTeX Manuscript Repository)

This repository contains a scientific LaTeX manuscript project for:

**The Geometry of Spacetime: A Temporal-Radial Hyperspherical Model of Cosmological Expansion**.

## Repository layout

- `main.tex` — manuscript entry point (title, abstract, packages, bibliography setup).
- `sections/manuscript.tex` — main scientific content.
- `figures/*.tikz` — mathematical schematic illustrations (TikZ source).
- `references.bib` — bibliography database.
- `.gitignore` — LaTeX build artifacts.

## Build

### Option A (recommended if available)

```bash
latexmk -pdf main.tex
```

### Option B

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

## Included mathematical illustrations

1. Temporal-radial foliation of hypersurfaces.
2. Schematic past-directed view curve through nested hypersurfaces.
3. Expansion-history and comoving-integrand intuition panel.

All figures are vector-based TikZ and can be edited directly in `figures/`.
