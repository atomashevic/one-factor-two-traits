# Fitting One Factor to Two Traits: reproduction materials

These materials accompany the article *Fitting One Factor to Two Traits: Population
Solutions and Fit Indices* by Aleksandar Tomašević
([ORCID 0000-0003-4863-6051](https://orcid.org/0000-0003-4863-6051)). The article
derives the population one-factor solution and its fit indices when the items
measure two correlated traits.

The repository holds:
- the code;
- the archived simulation summaries;
- portable evidence drawn from the raw simulation output;
- the complete proofs (Online Resources 1–3);
- the manuscript sources.

Section, table and appendix numbers refer to version 1 of the preprint.

## Contents

| Folder | Contents | Used in the article |
|---|---|---|
| `manuscript/` | LaTeX sources and the five figures | The article |
| `online-resources/` | Online Resources 1–3 as PDFs, built from `proofs/` | Complete proofs cited in Sections 3 and 5 and Appendices A–D |
| `proofs/` | Proof sources, verification scripts and their saved output | Online Resources 1–3 |
| `checks/` | Checks of displayed numbers, exact algebra and the PHQ-ADS calculations | Tables 1, 2, 4 and 5; Theorem 1; Sections 5–7 |
| `simulation/` | Scripts of the finite-sample study: 1,104 cells, 2,314,000 one-factor fits | Section 6, Appendix G |
| `results/` | Archived simulation summaries, population values and run logs | Section 6, Appendix G |
| `evidence/` | Selections of the raw simulation output, with the scripts that exported them | Figures 3–5, Table 4, Appendix G |
| `figures/` | Data behind the figures, and the scripts that write and plot it | Figures 1–5 |
| `numerics/exploratory/` | Scripts written while the results were derived, kept as they were used | Table 1, Table 3, Appendix H and side checks (see its README) |
| `empirical/kroenke2016/` | Values transcribed from Kroenke et al. (2016), Tables S2 and S3 | Section 7, Table 5 |
| `undeRfactor/` | The R package undeRfactor 0.2.0, identical to its archived release ([10.5281/zenodo.23668208](https://doi.org/10.5281/zenodo.23668208)) | Population calculations throughout |

## Reproducing the checks

The default build needs Python 3.10 or later (standard library only), TeX Live
with latexmk, pdfLaTeX and XeLaTeX, and Pandoc.

```sh
python3 build.py
```

This builds the following PDFs into `build/`:
- the anonymised manuscript;
- the named preprint;
- the title page;
- Online Resources 1–3.

It then runs three checks: the portable evidence against the figure and table data, the manuscript's internal references and numbering, and selected displayed numbers.

`--checks` adds more checks:
- the exact algebra;
- the unequal-block verifier and its interval certificate;
- the package tests;
- the item-order, score-share and PHQ-ADS calculations.

It needs NumPy, SciPy, SymPy and mpmath, and R 4.1 or later with lavaan, pbivnorm, numDeriv, pkgload and testthat. With [uv](https://docs.astral.sh/uv/):

```sh
uv run --no-project --python 3.10 --with numpy==1.26.4 --with scipy==1.10.1 \
  --with sympy==1.14.0 --with mpmath==1.3.0 python build.py --checks
```

**Tested setup.** These materials were last checked on 26 October 2026 with R 4.6.1, lavaan 0.7.2, TeX Live 2026 and Pandoc 3.10.2, where `--checks` takes a few minutes.

**Outputs.** The build writes only to `build/`. Each step's log is in `build/manuscript/` or `build/checks/`.

## What the build does not rerun

The build never reruns the simulation. `simulation/README.md` describes how to run it again.

**Raw output.** The raw output is one file per cell plus the combined fit table, 324 MB in all. It is not included here and is available from the author on request.

**What the checks use instead.** The checks use the summaries in `results/` and the selections in `evidence/`.

**Software versions.** The archived simulation run used R 4.5.1 and lavaan 0.6-20. The controls in Appendix G used R 4.6.1 and lavaan 0.7-2.

## Licence

- **Code** (`build.py` and the `.py` and `.R` files): MIT licence, see [LICENSE](LICENSE). `undeRfactor/` carries its own MIT licence.
- **Text, proofs, data and results:** Creative Commons Attribution 4.0 International licence, see [LICENSE-CC-BY-4.0.md](LICENSE-CC-BY-4.0.md).
- **Kroenke et al. (2016) values:** those in `empirical/kroenke2016/` come from that article and are not covered by these licences. Cite the original when using them.

## Citation

Please cite the article and these materials; the citation details are in [CITATION.cff](CITATION.cff). The R package has its own DOI, [10.5281/zenodo.23668208](https://doi.org/10.5281/zenodo.23668208).
