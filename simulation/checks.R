## checks.R
## Sanity checks of the building blocks against independent computations.
## Run before the study: Rscript checks.R

.here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
source(file.path(.here, "R", "sim_functions.R"))
set.seed(20261008)
ok <- function(label, cond, detail = "") {
  cat(sprintf("[%s] %s %s\n", if (isTRUE(cond)) "pass" else "FAIL", label, detail))
  invisible(cond)
}
pops <- make_populations()
P <- function(cond, phi) pops[pops$cond == cond & abs(pops$phi - phi) < 1e-9, ][1, ]

## 1. Exact score covariance against Monte Carlo
pop <- P("C_9_7_l75", 0.653); Sigma <- pop_sigma(pop)
for (fmt in c("ord_sym", "ord_skew")) {
  C <- score_cov(Sigma, THRESHOLDS[[fmt]])
  X <- discretize(gen_latent(4e5, Sigma), THRESHOLDS[[fmt]])
  d <- max(abs(cov(X) - C))
  ok(sprintf("score covariance, %s", fmt), d < 0.01, sprintf("(max abs diff %.4f)", d))
}

## 2. Population identities against lavaan's fit measures at very large N
pop <- P("A_m4_l70", 0.6); Sigma <- pop_sigma(pop)
f <- fit_ml_cov(model_1f(8), Sigma, 1e9, rescale = FALSE)
pi <- pop_indices(2 * f@optim$fx, Sigma, 8 * 9 / 2 - 16)
fm <- fitMeasures(f, c("cfi", "tli", "rmsea"))
d <- max(abs(pi[c("cfi", "tli", "rmsea")] - fm))
ok("CFI, TLI, RMSEA closed forms vs lavaan at N = 1e9", d < 1e-6, sprintf("(max abs diff %.2e)", d))

## 3. Proposition 2 rate against a finite difference on the boundary
F_ml <- function(Sig, R) {
  as.numeric(determinant(Sig)$modulus) + sum(diag(solve(Sig, R))) -
    as.numeric(determinant(R)$modulus) - nrow(R)
}
worst <- 0
for (k in 1:50) {
  p <- 6
  A <- matrix(rnorm(p * 3), p); R <- cov2cor(A %*% t(A) + diag(runif(p, .2, 1)))
  j <- sample(p, 1)
  l <- R[, j]; psi <- 1 - l^2; psi[j] <- 0
  star <- outer(l, l) + diag(psi)
  e <- 1e-6; E <- matrix(0, p, p); E[j, j] <- e
  fd <- (F_ml(star + E, R) - F_ml(star, R)) / e
  worst <- max(worst, abs(fd - prop2_rate(R, j)) / max(1, abs(fd)))
}
ok("Proposition 2 rate vs finite difference (50 random R)", worst < 1e-3, sprintf("(max rel diff %.1e)", worst))

## 4. Symmetry breaking at phi* (Theorem 1), Design A populations
for (cd in c("A_m4_l50", "A_m4_l70", "A_m8_l50", "A_m8_l70")) {
  pp <- pops[pops$cond == cd, ][1, ]
  ps <- pp$phistar
  sh <- sapply(c(ps - 0.03, ps + 0.03), function(ph) {
    q <- pp; q$phi <- ph
    r <- pop_ml(pop_sigma(q), q$mA, q$mB, ref_latent(q), NA, nstart = 5)
    max(r$best$shareA, 1 - r$best$shareA)
  })
  ok(sprintf("%s: tilted below phi* = %.3f, composite above", cd, ps),
     sh[1] > 0.6 && abs(sh[2] - 0.5) < 1e-3, sprintf("(share %.3f / %.3f)", sh[1], sh[2]))
}

## 5. Fold points of Proposition 1, Design B
for (cd in names(PHI_FOLD)) {
  pp <- pops[pops$cond == cd, ][1, ]
  ex <- sapply(PHI_FOLD[[cd]] + c(-0.01, 0.01), function(ph) {
    q <- pp; q$phi <- max(ph, 0)
    pop_ml(pop_sigma(q), q$mA, q$mB, ref_latent(q), q$weak, nstart = 10)$weak_exists
  })
  ok(sprintf("%s: weaker-block minimum just below phi_f = %.3f, none just above", cd, PHI_FOLD[[cd]]),
     ex[1] && !ex[2])
}

## 6. Share mapping: a population composite returns the composite weights
q <- P("B_l68_l72", 0.5)
lA <- q$lamA * (0.3 + q$phi * 0.7); lB <- q$lamB * (q$phi * 0.3 + 0.7)
ok("share mapping inverts l = lambda (w_A + phi w_B)", abs(share_A(lA, lB, ref_latent(q)) - 0.3) < 1e-12)
