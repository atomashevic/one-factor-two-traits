# Archived results of the finite-sample study

These files are the output of the simulation in `../simulation/` (Section 6 and
Appendix G). The run had 1,104 cells, 2,314,000 one-factor fits and 156,000
multistart replications, and used R 4.5.1 with lavaan 0.6-20. `build.py` reads
these files and never overwrites them. File names that include `stage4` or
`stage5` record the stage of the analysis that produced them.

| File | Contents |
|---|---|
| `population.rds` | Population reference values for every population and estimator; a solution is certified as global only where the article proves it |
| `summary/q1_fit.csv` | Sample fit indices and the rates at which conventional cutoffs are met |
| `summary/q2_diagnostics_stage5.csv` | Signed, magnitude and polarity diagnostics of the fitted solutions, as used in the article |
| `summary/q2_diagnostics_stage4.csv` | The earlier diagnostics (1,946 solution groups) from which `simulation/posthoc/q2fix.R` derives the stage 5 file |
| `summary/q2_tilt.csv` | An earlier magnitude-normalized summary, corrected and kept for comparison; the article uses `q2_diagnostics_stage5.csv` |
| `summary/q3_minima.csv` | Solution groups retained from the multistart fits; they do not certify stationarity, local minimality or global optimality |
| `summary/heywood.csv` | Records of fits with zero or negative uniquenesses |
| `summary/order_ctrl_stage4.csv` | Item-order controls of Appendix G: 48 controlled fits and 68 reference candidates, with all fitted parameters |
| `summary/sample_ctrl_stage4.csv` | Three selected sample refits of Appendix G, with complete parameters and the archived comparisons |
| `logs/` | Logs of the runs: `sim.log` has one line per cell; `progress.log`, `population.log` and `summary.log` cover the main run; the other three cover the diagnostics and controls |

**Not included.**
- The raw output, `raw/cell_*.rds` with one file per cell and `fits1f.rds`, 324 MB together. It is available from the author on request.
- The draft figures that `summarise_sim.R` writes to `summary/`. The article's figures are built from `../figures/`.
