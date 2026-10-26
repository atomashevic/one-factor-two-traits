# Finite-sample study (Section 6 and Appendix G)

R scripts for the simulation study of "Fitting One Factor to Two Traits:
Population Solutions and Fit Indices". Requires R (4.3 or later) and lavaan;
`pbivnorm` and `parallel` come with lavaan and R. The archived run used R 4.5.1
with lavaan 0.6-20 (Appendix G).

## Files

| File | Purpose |
|---|---|
| `R/sim_functions.R` | Design grid, data generation, model fitting, outcomes, population values |
| `R/cli.R` | Command-line parsing (no extra packages) |
| `R/thresholds.R` | Exact asymptotic variances of polychoric correlations and the WLSMV symmetry-breaking threshold |
| `analysis/wlsmv_threshold.R` | Derives and checks the WLSMV thresholds that set the low-φ points of Design A (about 20 s) |
| `checks.R` | Checks the building blocks against independent computations (about 10 s) |
| `run_population.R` | Population values for every population and estimator (about 1 CPU hour) |
| `run_sim.R` | Runs cells of the design, one output file per cell; resumable |
| `summarise_sim.R` | Tables (CSV) and draft figures (PDF) for Questions 1 to 3 |

## Design

92 populations × 4 sample sizes (250, 500, 1,000, 2,000) × 3 response formats
(continuous; four categories from thresholds −.8, 0, .8; four categories from
.3, 1.0, 1.7) = 1,104 cells. Replications: 1,000 per continuous cell and 500 per
ordinal cell (736,000 in total).

| Design | Blocks | Loadings | Trait correlations |
|---|---|---|---|
| A | 4+4 or 8+8 | .5 or .7 | 0 to .9 by .1, plus φ*, 1/p, 1/(2p), and .08 at m = 4, λ = .5 |
| B | 4+4 | .68/.72 and .60/.80 | 0 to .9 by .1, plus .01, .24 and the φ* of the exchangeable design at .7 (.409) |
| C | 9+7 | .75 | 0 to .9 by .1, plus .653, .865, .912 |

φ* = (1 + (m − 1)λ²) / (m(2 − λ²)): .250 and .196 at λ = .5, .409 and .367 at λ = .7
(m = 4 and 8 items per trait). Design C has no symmetric threshold, because the
blocks differ in length.

The low-φ points of Design A locate the thresholds of the other estimators. ULS breaks
symmetry below 1/p (Theorem 2): .125 with 8 items, .0625 with 16. WLSMV fits the polychoric
correlations by least squares with weights equal to their inverse asymptotic variances, w_W
within blocks and w_C across them; with κ = w_C/w_W its symmetric solution is unstable exactly
when φ < 1/2 − (m − 1)/(2mκ), which gives 1/p at κ = 1. Because within-block correlations are
larger and so estimated more precisely, κ < 1 and the WLSMV threshold is lower: .081 (symmetric
categories) and .076 (skewed) at m = 4, λ = .5; .010 and .004 at m = 8, λ = .5; and none at
λ = .7, where κ(0) < (m − 1)/m and the WLSMV fit stays a composite at every φ ≥ 0.
`analysis/wlsmv_threshold.R` derives these and checks them against lavaan.

Every replication fits the one-factor and the correlated two-factor model:

* continuous data: ML on the sample covariance matrix; ULS as well in Design A at N = 500 and 2,000;
* ordinal data: WLSMV (`ordered`, delta parameterization) and ML on the 0 to 3 scores;
* Designs B and C at N = 500 and 2,000: every one-factor ML fit (continuous ML and ML on scores) is
  repeated on the sample correlation matrix from lavaan's default start, once with the items in
  design order and once with the blocks reversed, and from 20 random starts (loadings U(0, 1),
  uniquenesses U(.2, 1)). Converged solutions are clustered into distinct minima (standardized
  loadings within .002). The two default fits matter because lavaan's default starting values use
  the first-listed item as reference: when the weaker-block minimum exists, the default start
  tends to land in the basin of whichever block is listed first.

Seeds: replication r of population k at the j-th sample size uses seed
100000·k + 10000·j + r. All three formats of a population and sample size share the same
latent draws (common random numbers), so format differences are not sampling noise.

## Outcomes recorded per fit

`fx` (lavaan's objective; for ML, F = 2·fx), χ², df, p, CFI, TLI, RMSEA, SRMR; for WLSMV also
the scaled (`_sc`, as Mplus reports them) and robust (`_rob`) versions; block-average standardized
loadings `lA`, `lB`; the share of block A, `shareA` = w_A/(w_A + w_B), from
l_A = λ_A(w_A + φw_B), l_B = λ_B(φw_A + w_B) with the population two-factor loadings and φ on the
estimator's scale; the factor correlation of the two-factor model; the smallest standardized
uniqueness; and, when it is below .05, the Proposition 2 rate on the fitted correlation matrix
(Pearson for ML, polychoric for WLSMV) at that item.

## Population values

* Continuous ML: exact. The one-factor global minimum comes from multi-start fits to Σ;
  CFI = 1 − F/F_B, TLI = 1 − (F/df)/(F_B/df_B), RMSEA = √(F/df); SRMR from lavaan at the population fit.
  All distinct minima are kept, so the fold points of Proposition 1 are checked directly.
* ML on 0 to 3 scores: exact, from the covariance matrix of the scores (bivariate normal probabilities).
* WLSMV: the true polychoric matrix and thresholds, with Γ estimated once from 200,000 cases
  (`--ngamma`) and the correlation weights set to their exact values, at nominal N = 10⁶. Exact
  weights matter: sampling noise in estimated weights acts as a small asymmetry, which near a
  threshold shifts the tilt by its cube root (Theorem 3). Fitting from summary statistics
  reproduces a full-data WLSMV fit exactly.

## Running

```sh
Rscript checks.R                                   # all checks should pass
Rscript run_population.R --cores 8                 # writes population.rds
Rscript run_sim.R --list                           # the design
Rscript run_sim.R --pilot --cores 8                # 72 cells, few reps, prints a cost estimate (about 155 CPU hours)
Rscript run_sim.R --cores 8                        # everything, into results/
Rscript summarise_sim.R                            # tables and figures into summary/
```

On a cluster, use an array job and leave `--cores` at the cores per task:

```sh
#SBATCH --array=1-100
#SBATCH --cpus-per-task=4
Rscript run_sim.R --task $SLURM_ARRAY_TASK_ID --ntasks 100 --cores 4
```

Cells already written are skipped, so a failed or interrupted task can be resubmitted.
