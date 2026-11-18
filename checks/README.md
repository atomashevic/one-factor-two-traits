# Canonical focused checks

Run `python3 ../build.py --checks` from this directory (with the scientific
Python/R dependencies listed in the root README). New outputs go to `build/`.
`order_ctrl.R` and `score_share.R` use the bundled package source.

Current manuscript locations are Section 6 and Appendix G for item-order
controls, and Sections 5.2 and 7 for score shares. The historical table below
retains the earlier manuscript's section numbers for provenance.

| Script | Paper | What it does |
|---|---|---|
| `order_ctrl.R` | Section 7.3, "Starting values, not identification, drive the order effect when the blocks differ" | Refits population matrices in both item orders, with the factor variance fixed (std.lv) and with marker identification, from lavaan's default start and from a common starting covariance |
| `score_share.R` | Sections 3.1 and 8 | Trait shares of the fitted factor, of the regression and Bartlett scores, and of the unit-weighted sum |

Both use bundled `undeRfactor` functions and need lavaan.

Python checks added on 11 October 2026:

| Script | Paper | What it does |
|---|---|---|
| `table2_ratios.py` | Table 2 | Recomputes the smallest composite-making slope ratios from Proposition 3 and Theorem 1 and checks the rounded-up table cells (standard library only; part of the default build) |
| `theorem1_exact.py` | Theorem 1, Corollary 1, §4, Lemma A.5 | 49 reported checks, of which 45 are computed SymPy identities for the closed forms and the boundary sum of Lemma A.5 (`--checks`) |
| `unequal_example_certificate.py` | §5.2, Online Resource 3 §2.6 | Interval-arithmetic certificate that the longer block is the unique global choice for 0 < φ ≤ .51 in the 2 × .85 / 3 × .75 example (`--checks`) | Prefer the root build
entry point, which sets numerical thread limits and preserves archived outputs.

R check added on 23 October 2026:

| Script | Paper | What it does |
|---|---|---|
| `phq_table5.R` | Table 5, Section 7 | Recomputes every Table 5 cell and each Section 7 number from the published two-factor loadings of Kroenke et al. (2016, Table S3; `empirical/kroenke2016/`) and the common-loading sensitivity range, and compares them with the LaTeX source (`--checks`) |

Python checks added on 17 November 2026 (all `--checks`):

| Script | Paper | What it does |
|---|---|---|
| `share_ordering.py` | Section 5.2, Appendix D.6 | Exact SymPy check of the coefficient-share and score-share identities and of φ − q > 0 in the exchangeable design; reproduces the shares .786 and .750 |
| `least_squares_relaxation.py` | Appendix C | Symbolic identities for the relaxed least-squares objective, and a comparison of the closed-form minimizers with a bounded multistart search on 180 cases |
| `appendix_d_chain.py` | Appendix D | Exact identities used in the proofs of Lemmas D.3–D.5, Proposition 2 and Theorem 3, including the contrast quadratic form at five rational stationary points |

These scripts check the identities and cases they compute. They do not verify the theorems.

