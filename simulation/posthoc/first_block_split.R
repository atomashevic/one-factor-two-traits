## First-listed-block proportions among magnitude-selected fits (|M| > .1) in
## Design A, split by the population regime of each estimator. Read-only.
##   Rscript simulation/posthoc/first_block_split.R [--out FILE.csv]
## Inputs: results/summary/q2_diagnostics_stage5.csv (per-cell counts written by
## q2fix.R) and results/population.rds (population solutions). If the local
## archive results/fits1f.rds exists, the cell behind Appendix G's averaging
## example is also recomputed from its rows.
## Regime rule. "Below" means that the population solution of that estimator
## has unequal block loadings; threshold equality counts as the composite
## regime ("at or above").
##   ML, continuous: phi < phi* (Theorem 1).
##   ULS:            phi < 1/p (Theorem 2).
##   ML on scores and WLSMV: the best population solution found for that
##     estimator and response format has |lA - lB| > .005. For ML and ULS this
##     numerical rule is checked against the analytic one.
argv <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), '../..'))
out <- if (length(argv) == 2L && argv[1] == '--out') argv[2] else ''
d <- read.csv(file.path(root, 'results/summary/q2_diagnostics_stage5.csv'))
pop <- readRDS(file.path(root, 'results/population.rds'))
pop <- pop[pop$design == 'A' & pop$model == '1F',
           c('pop_id', 'format', 'est', 'lA', 'lB')]
pop$pop_imbalance <- abs(pop$lA - pop$lB)
a <- merge(d[d$design == 'A', ], pop[c('pop_id', 'format', 'est', 'pop_imbalance')],
           by = c('pop_id', 'format', 'est'), all.x = TRUE)
stopifnot(!anyNA(a$pop_imbalance))
p <- a$mA + a$mB
numeric_below <- a$pop_imbalance > .005
analytic_below <- ifelse(a$est == 'ML', a$phi < a$phistar - 1e-9,
                  ifelse(a$est == 'ULS', a$phi < 1 / p - 1e-9, NA))
chk <- a$est %in% c('ML', 'ULS')
stopifnot(all(numeric_below[chk] == analytic_below[chk]))
cat('Analytic and numerical regime rules agree in all', sum(chk), 'ML and ULS cell rows.\n')
a$below <- ifelse(chk, analytic_below, numeric_below)

analyses <- list(c('ML', 'cont'), c('MLs', 'ord_sym'), c('ULS', 'cont'), c('WLSMV', 'ord_sym'),
                 c('MLs', 'ord_skew'), c('WLSMV', 'ord_skew'))
ratio <- function(x) {
  s <- sum(x$n_selected_M_gt_10); f <- sum(x$n_selected_first)
  data.frame(cells = nrow(x), attempted = sum(x$n_attempted), selected = s, first = f,
             pct_first = if (s) round(100 * f / s, 1) else NA_real_)
}
rows <- list()
for (N in c(2000, 250)) for (an in analyses) {
  x <- a[a$N == N & a$est == an[1] & a$format == an[2], ]
  if (!nrow(x)) next
  parts <- list(all = x, phi0 = x[x$phi == 0, ], below = x[x$below, ],
                at_or_above = x[!x$below, ])
  for (k in names(parts))
    rows[[length(rows) + 1L]] <- cbind(N = N, est = an[1], format = an[2], subset = k,
                                       ratio(parts[[k]]))
}
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)
cat('\nCells below the population threshold at N = 2,000 (cond: phi values):\n')
for (an in analyses) {
  x <- a[a$N == 2000 & a$est == an[1] & a$format == an[2] & a$below, ]
  cat(sprintf('  %s %s: %s\n', an[1], an[2],
              paste(tapply(x$phi, x$cond, function(v) paste(sort(v), collapse = ', ')),
                    collapse = ' | ')))
}
if (nzchar(out)) { write.csv(tab, out, row.names = FALSE); cat('Wrote', out, '\n') }

## The averaging example of Appendix G: A_m4_l70, phi = 0, N = 2,000, continuous ML.
e <- d[d$cond == 'A_m4_l70' & d$phi == 0 & d$N == 2000 & d$est == 'ML', ]
stopifnot(nrow(e) == 1L)
cat(sprintf('\nAveraging example (summary row): common flips %d, opposite signs %d, mean D %.4f, mean |D| %.4f, mean |M| %.4f, selected %d, first %d\n',
            e$n_common_flips, e$n_opposite_sign, e$mean_D, e$mean_abs_D,
            e$mean_abs_magnitude_D, e$n_selected_M_gt_10, e$n_selected_first))
arc <- file.path(root, 'results/fits1f.rds')
if (file.exists(arc)) {
  f <- readRDS(arc)
  x <- f[f$cell_id == e$cell_id & f$est == 'ML', ]
  flip <- e$mA * x$lA + e$mB * x$lB < 0
  la <- ifelse(flip, -x$lA, x$lA); lb <- ifelse(flip, -x$lB, x$lB)
  D <- la - lb; M <- abs(la) - abs(lb)
  cat(sprintf('Archive rows %d: flips %d; mean D %.4f; mean |D| %.4f; mean |M| %.4f; D < 0 in %d; M > .1 in %d; M < -.1 in %d\n',
              nrow(x), sum(flip), mean(D), mean(abs(D)), mean(abs(M)), sum(D < 0),
              sum(M > .1), sum(M < -.1)))
  cat(sprintf('  block A mean loading: %.3f to %.3f; block B: %.3f to %.3f, negative in %d\n',
              min(la), max(la), min(lb), max(lb), sum(lb < 0)))
} else cat('results/fits1f.rds not present; archive-row recomputation skipped.\n')
