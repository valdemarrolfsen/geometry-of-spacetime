# Geometry of Spacetime

A reproducible LaTeX research project about **expanding spatial layers, local gravity, and the past-directed light paths that form an observer's sky**.

The current manuscript is **A Temporal–Radial View of an Expanding Three-Sphere with a Localized Mass Perturbation**, working draft v0.5. It develops a smooth compensated dust overdensity on a closed FLRW background, derives spatial deformation and clock rates, and calculates light deflection and timing changes. This is a first-order illustrative GR calculation, not a black-hole model or an observational fit. It has not been peer reviewed.

- [Read the compiled paper](output/pdf/geometry-of-spacetime-v0.5.pdf)
- [Diagram comparison and reader-study proposal](output/pdf/diagram-comparison.pdf)
- [LaTeX entry point](paper/main.tex) · [Sections](paper/sections) · [Bibliography](paper/references.bib)
- [Numerical results and diagnostics](results) · [Generated figures](paper/figures)
- [Original v0.4 manuscript and figures](archive/v0.4)

## Build and reproduce

Use Python 3.12 for the recorded environment. The source supports Python 3.11+, subject to dependency availability. Install [Tectonic](https://tectonic-typesetting.github.io/) (tested with 0.17.0) to compile the paper.

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.lock
make PYTHON=.venv/bin/python all
```

`make all` runs numerical and symbolic checks, regenerates six physics figures plus three comparison figures and all result tables, and compiles both PDFs in `output/pdf/`. Tectonic downloads standard TeX packages on first use and caches them locally in `.tectonic-cache/`; subsequent builds can use that cache. No observations, external datasets, API keys, or credentials are required for the calculations.

Individual steps:

```sh
make PYTHON=.venv/bin/python check    # numerical tests and symbolic identities
make PYTHON=.venv/bin/python results  # CSV, JSON, PDF/SVG/PNG figures, LaTeX numbers
make PYTHON=.venv/bin/python paper    # regenerate results and compile the manuscript
make PYTHON=.venv/bin/python comparison # compile the diagram assessment
```

An existing TeX Live installation can compile the generated paper with:

```sh
make PYTHON=.venv/bin/python results
mkdir -p build output/pdf
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../build main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../build diagram-comparison.tex
cd ..
cp build/main.pdf output/pdf/geometry-of-spacetime-v0.5.pdf
cp build/diagram-comparison.pdf output/pdf/diagram-comparison.pdf
```

For Overleaf, upload the `paper/` directory with its generated figures and `generated/numbers.tex`, and set `main.tex` as the main document. Regenerate the numbers locally after changing the calculation. The committed figures allow the paper to compile without running Python in Overleaf.

## Project structure

```text
paper/
  main.tex                 LaTeX entry point
  diagram-comparison.tex   Literature comparison, map limits and reader-study proposal
  sections/                Argument and derivations, split by topic
  references.bib           Active bibliography
  figures/                 Generated vector figures and PNG previews
  generated/numbers.tex    Numerical values imported into the paper
src/spacetime/
  model.py                 Closed dust background, potential and matter constraints
  rays.py                  Null rays, redshift and matched-source shooting
  diagrams.py              Matched event maps and visible spatial subsets
scripts/
  symbolic_checks.py       SymPy algebraic checks
  reproduce.py             Generates all current numerical outputs and figures
  compare_diagrams.py      Generates the three diagram comparisons and diagnostics
tests/                     Analytic limits and independent physics checks
results/                   CSV data, diagnostics and symbolic-check results
output/pdf/                Compiled review copy
archive/v0.4/              Previous draft, notes and figures, preserved verbatim
requirements.lock          Exact recorded Python dependency versions
Makefile                   Check, reproduce and build entry points
```

## Mathematical model

The background is the expanding branch of a closed dust universe: `A/A_max = sin²(eta/2)`, `ct/A_max = (eta - sin eta)/2`, and `0 < eta < pi`. It is an analytic test case, not a fitted cosmology.

The dimensionless scalar potential is `phi = epsilon * g(eta) * F(cos chi)`. Its time dependence solves the dust Bardeen equation. Its smooth spatial profile removes the homogeneous and dipole modes on S³. The linear Einstein constraints reconstruct the density and velocity, including compensation around the positive mass excess. A default potential depth of `1e-4` gives a density perturbation below 0.4%.

The same metric determines:

1. The lapse and accumulated proper time of coordinate-stationary test clocks.
2. A first-order spatial embedding. In this construction the negative potential gives an **outward** radial deformation and a reduced clock rate.
3. Null paths, source-matched observer-frame angular shifts, timing changes, and normal-observer redshift.

The ray solver integrates the displayed truncated metric. Its physical accuracy remains first order even though the numerical integration retains some higher-order dependence. An independent first-order timing integral and an amplitude scan verify an O(epsilon²) remainder. The paper explains the gauge, observer, and source-matching conventions. Dust Doppler terms and black-hole horizons are not included.

## Diagram comparison

The basic polar map has a direct precedent in [Ned Wright's cosmology tutorial](https://www.astro.ucla.edu/~wright/cosmoall.htm). The [comparison brief](output/pdf/diagram-comparison.pdf) evaluates the same events in conventional and nested views, distinguishes whole spatial slices from observed subsets, and demonstrates early-epoch compression and failure of radius as a clock through recollapse. It recommends testing paired views; no human learning study has been run. Its diagnostics are in `results/diagram-comparison.json`.

## Change the example

Default model parameters are in `LocalizedPerturbation` in `src/spacetime/model.py`. Observer/source geometry and plotting choices are in `scripts/reproduce.py`. Re-run `make all` after changes. Narrowing the source can make the density nonlinear even when the potential stays small, so recheck the density and velocity diagnostics as well as the potential.

All reported times are in units of `A_max/c`, masses in `c² A_max/G`, and angular shifts in radians. Choosing a physical `A_max` supplies units; it does not make this toy background an observational fit.

## Development history

v0.4 focused on the homogeneous temporal–radial map and proposed anisotropy/rotation extensions. v0.5 makes the localized-gravity question the active research direction and converts the project to LaTeX plus tested Python calculations. The original v0.4 files were moved without changing their content. See [CHANGELOG.md](CHANGELOG.md) for scope and remaining work.
