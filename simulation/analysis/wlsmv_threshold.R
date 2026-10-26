## wlsmv_threshold.R
## Population symmetry-breaking threshold of the one-factor WLSMV (DWLS) fit in
## the exchangeable two-trait design.
##
## Argument. With a diagonal weight matrix and free thresholds, the DWLS
## discrepancy of the one-factor model reduces to a weighted least-squares fit
## of the off-diagonal polychoric correlations,
##     G(l) = sum_{i<j} w_ij (rho_ij - l_i l_j)^2,
## with w_ij the inverse asymptotic variance of the polychoric estimate of
## rho_ij. In the exchangeable design there are two weights, w_W for the
## within-block pairs (rho = L = lambda^2) and w_C for the cross-block pairs
## (rho = L phi). Write kappa = w_C / w_W and n = m(m-1)/2. On block-uniform
## loadings (a, b), the symmetric stationary point a = b = c has
##     c^2 = L (2n + kappa m^2 phi) / (2n + kappa m^2),
## and the second derivative of G along a = c + t, b = c - t is 8 n w_W (2c^2 - L).
## The symmetric solution is therefore unstable exactly when c^2 < L/2, i.e.
##     phi < phi_LS = 1/2 - (m - 1) / (2 m kappa).
## kappa = 1 (ULS) gives phi_LS = 1/(2m) = 1/p, Theorem 2. For WLSMV kappa
## depends on phi through the asymptotic variances, so the threshold solves
## phi = 1/2 - (m - 1)/(2 m kappa(phi)); if kappa(0) <= (m - 1)/m there is no
## threshold and the symmetric solution is stable for every phi >= 0.
##
## The functions are in R/thresholds.R. This script checks the exact
## asymptotic variances against lavaan's Gamma from a large sample, solves for
## the threshold in every Design A condition, and checks it with direct DWLS fits.

.here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
source(file.path(.here, "..", "R", "sim_functions.R"))
source(file.path(.here, "..", "R", "thresholds.R"))

## ------------------------------------------- direct weighted LS fit ----

## Global one-factor minimum of G with exact weights; returns the share of the
## followed trait computed from the block-average loadings.
dwls_share <- function(m, lam, phi, thr, nstart = 8) {
  L <- lam^2; p <- 2 * m
  wW <- 1 / avar_polychoric(L, thr, thr); wC <- 1 / avar_polychoric(L * phi, thr, thr)
  pop <- data.frame(mA = m, mB = m, lamA = lam, lamB = lam, phi = phi)
  R <- pop_sigma(pop)
  blk <- rep(1:2, each = m)
  W <- ifelse(outer(blk, blk, "=="), wW, wC); diag(W) <- 0
  G <- function(l) sum(W * (R - outer(l, l))^2) / 2
  gr <- function(l) -2 * as.vector((W * (R - outer(l, l))) %*% l)
  best <- NULL
  starts <- c(list(c(rep(0.9, m), rep(0.1, m)), rep(lam, p)), replicate(nstart, runif(p), simplify = FALSE))
  for (s in starts) {
    o <- optim(s, G, gr, method = "BFGS", control = list(reltol = 1e-14, maxit = 2000))
    if (is.null(best) || o$value < best$value - 1e-12) best <- o
  }
  l <- abs(best$par)
  share_A(mean(l[1:m]), mean(l[m + 1:m]), list(lamA = lam, lamB = lam, phi = phi)) -> sA
  max(sA, 1 - sA)
}

## ---------------------------------------------------------------- run ----

if (sys.nframe() == 0) {
  set.seed(1)
  cat("Check 1: exact polychoric asymptotic variances vs lavaan Gamma (N = 200,000)\n")
  pop <- data.frame(mA = 4, mB = 4, lamA = 0.7, lamB = 0.7, phi = 0.3)
  Sigma <- pop_sigma(pop)
  for (fmt in c("ord_sym", "ord_skew")) {
    thr <- THRESHOLDS[[fmt]]
    X <- discretize(gen_latent(2e5, Sigma), thr)
    f0 <- quiet(cfa(model_1f(8), data = as.data.frame(X), ordered = colnames(X),
                    estimator = "WLSMV", do.fit = FALSE))
    Gd <- diag(lavInspect(f0, "gamma"))[-(1:(8 * 3))]
    pr <- which(lower.tri(Sigma), arr.ind = TRUE)
    within <- (pr[, 1] <= 4) == (pr[, 2] <= 4)
    cat(sprintf("  %-8s within: lavaan %.4f exact %.4f | cross: lavaan %.4f exact %.4f\n", fmt,
                mean(Gd[within]), avar_polychoric(0.49, thr, thr),
                mean(Gd[!within]), avar_polychoric(0.49 * 0.3, thr, thr)))
  }

  cat("\nThresholds (phi* = ML, Theorem 1; 1/p = ULS, Theorem 2; WLSMV from the formula)\n")
  out <- list()
  for (m in c(4, 8)) for (lam in c(0.5, 0.7)) for (fmt in c("ord_sym", "ord_skew")) {
    thr <- THRESHOLDS[[fmt]]
    k0 <- kappa_wlsmv(lam^2, 0, thr)
    pw <- wlsmv_threshold(m, lam, thr)
    ## numerical check: global DWLS share just below and just above
    chk <- if (is.na(pw)) c(dwls_share(m, lam, 0, thr), NA) else
      c(dwls_share(m, lam, max(pw - 0.02, 0), thr), dwls_share(m, lam, pw + 0.02, thr))
    out[[length(out) + 1]] <- data.frame(m = m, lambda = lam, format = fmt,
      phistar_ML = round(phistar(m, lam), 3), uls = round(1 / (2 * m), 3),
      kappa0 = round(k0, 3), kappa_crit = round((m - 1) / m, 3),
      phistar_WLSMV = round(pw, 3), share_below = round(chk[1], 3), share_above = round(chk[2], 3))
  }
  res <- do.call(rbind, out)
  print(res, row.names = FALSE)
  saveRDS(res, file.path(.here, "wlsmv_threshold.rds"))
}
