## thresholds.R
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
## kappa(phi) is computed exactly, from the influence function of the two-step
## polychoric estimator. analysis/wlsmv_threshold.R checks it against lavaan.

## ---------------------------------------- exact polychoric asymptotics ----

F2 <- function(x, y, rho) {
  if (x == -Inf || y == -Inf) return(0)
  if (x == Inf && y == Inf) return(1)
  if (x == Inf) return(pnorm(y))
  if (y == Inf) return(pnorm(x))
  if (abs(rho) < 1e-14) return(pnorm(x) * pnorm(y))
  pbivnorm::pbivnorm(x, y, rho = rho)
}

cell_probs <- function(rho, t1, t2) {
  a <- c(-Inf, t1, Inf); b <- c(-Inf, t2, Inf)
  P <- matrix(0, length(a) - 1, length(b) - 1)
  for (k in seq_len(nrow(P))) for (l in seq_len(ncol(P)))
    P[k, l] <- F2(a[k + 1], b[l + 1], rho) - F2(a[k], b[l + 1], rho) -
               F2(a[k + 1], b[l], rho) + F2(a[k], b[l], rho)
  P
}

## Asymptotic variance (times N) of the two-step polychoric estimator:
## thresholds from the margins, then rho by bivariate ML given the thresholds.
avar_polychoric <- function(rho, t1, t2, h = 1e-5) {
  P <- cell_probs(rho, t1, t2)
  d_rho <- (cell_probs(rho + h, t1, t2) - cell_probs(rho - h, t1, t2)) / (2 * h)
  s_rho <- d_rho / P
  I_rr <- sum(P * s_rho^2)
  K1 <- length(t1); K2 <- length(t2)
  corr <- 0 * P
  for (a in seq_len(K1)) {               # thresholds of item 1
    e <- replace(0 * t1, a, h)
    s_t <- ((cell_probs(rho, t1 + e, t2) - cell_probs(rho, t1 - e, t2)) / (2 * h)) / P
    I_rt <- sum(P * s_rho * s_t)
    IF_t <- outer((seq_len(K1 + 1) <= a) - pnorm(t1[a]), rep(1, K2 + 1)) / dnorm(t1[a])
    corr <- corr + I_rt * IF_t
  }
  for (b in seq_len(K2)) {               # thresholds of item 2
    e <- replace(0 * t2, b, h)
    s_t <- ((cell_probs(rho, t1, t2 + e) - cell_probs(rho, t1, t2 - e)) / (2 * h)) / P
    I_rt <- sum(P * s_rho * s_t)
    IF_t <- outer(rep(1, K1 + 1), (seq_len(K2 + 1) <= b) - pnorm(t2[b])) / dnorm(t2[b])
    corr <- corr + I_rt * IF_t
  }
  IF <- (s_rho - corr) / I_rr
  sum(P * IF^2)
}

kappa_wlsmv <- function(L, phi, thr) avar_polychoric(L, thr, thr) / avar_polychoric(L * phi, thr, thr)
phi_ls <- function(m, kappa) 0.5 - (m - 1) / (2 * m * kappa)

wlsmv_threshold <- function(m, lam, thr) {
  L <- lam^2
  g <- function(phi) phi - phi_ls(m, kappa_wlsmv(L, phi, thr))
  if (g(0) >= 0) return(NA_real_)       # stable at every phi >= 0
  uniroot(g, c(0, 0.5), tol = 1e-8)$root
}
