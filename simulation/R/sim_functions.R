## sim_functions.R
## Functions for the finite-sample study of Section 7.3 of
## "When one factor represents two traits: symmetry breaking and population
## comparative fit".
##
## Conventions
##   * Two standard normal traits correlate phi. Block A (items 1..mA) loads
##     lamA on trait 1, block B (items mA+1..mA+mB) loads lamB on trait 2.
##     Items are standardized: Sigma = Lambda Phi Lambda' + Psi, diag(Sigma) = 1.
##   * F is the ML discrepancy log|Sigma_hat| + tr(S Sigma_hat^-1) - log|S| - p,
##     equal to 2 * lavaan's fmin for ML. For other estimators fx is lavaan's fmin.
##   * The followed-trait share (Section 3.1): write the fitted factor as
##     w_A theta_1 + w_B theta_2 whose covariances with the items reproduce the
##     block-average standardized loadings, l_A = lamA (w_A + phi w_B) and
##     l_B = lamB (phi w_A + w_B). shareA = w_A / (w_A + w_B); the followed-trait
##     share is max(shareA, 1 - shareA); the signed tilt is shareA - .5.
##     The reference loadings and phi are the population two-factor values on the
##     scale of the estimator (latent for ML, ULS and WLSMV; score scale for ML
##     on 0-3 scores).

suppressPackageStartupMessages(library(lavaan))

## ---------------------------------------------------------------- design ----

NS <- c(250, 500, 1000, 2000)
FORMATS <- c("cont", "ord_sym", "ord_skew")
THRESHOLDS <- list(cont = NULL, ord_sym = c(-0.8, 0, 0.8), ord_skew = c(0.3, 1.0, 1.7))
REPS <- c(cont = 1000, ord_sym = 500, ord_skew = 500)
ULS_N <- c(500, 2000)          # ULS arm: Design A, continuous, these N
MS_N <- c(500, 2000)           # multi-start arm: Designs B and C, these N
NSTART <- 20                   # random starts per multi-start fit

## Symmetry-breaking threshold for the exchangeable design (Theorem 1)
phistar <- function(m, lambda) {
  L <- lambda^2
  (1 + (m - 1) * L) / (m * (2 - L))
}

## Population fold points of Proposition 1 (Design B): the weaker-block
## minimum exists for phi < phi_f. Checked independently by run_population.R.
PHI_FOLD <- c(B_l68_l72 = 0.242, B_l60_l80 = 0.017)

make_populations <- function() {
  rows <- list()
  add <- function(design, cond, mA, mB, lamA, lamB, phis, ps, weak, fold) {
    for (ph in phis) {
      rows[[length(rows) + 1]] <<- data.frame(
        design = design, cond = cond, mA = mA, mB = mB, lamA = lamA, lamB = lamB,
        phi = ph, phistar = ps, weak = weak, phi_fold = fold,
        stringsAsFactors = FALSE)
    }
  }
  base <- round(seq(0, 0.9, by = 0.1), 2)
  ## Design A: exchangeable. Besides phi*, the grid holds the ULS threshold
  ## 1/p of Theorem 2 and 1/(2p) below it, and, where it exists above .02, the
  ## WLSMV threshold of R/thresholds.R (.081 symmetric, .076 skewed at m = 4,
  ## lambda = .5; .010 and .004 at m = 8, lambda = .5; none at lambda = .7).
  for (m in c(4, 8)) for (lam in c(0.5, 0.7)) {
    ps <- phistar(m, lam)
    extra <- c(round(ps, 3), 1 / (2 * m), 1 / (4 * m))
    if (m == 4 && lam == 0.5) extra <- c(extra, 0.08)
    add("A", sprintf("A_m%d_l%02d", m, round(100 * lam)), m, m, lam, lam,
        sort(unique(round(c(base, extra), 4))), ps, NA_character_, NA_real_)
  }
  ## Design B: unequal loadings, 4 and 4 items; block A is the weaker block.
  ## phi* of the exchangeable design at the mean loading .7 is the reference.
  for (pr in list(c(0.68, 0.72), c(0.60, 0.80))) {
    cond <- sprintf("B_l%02d_l%02d", round(100 * pr[1]), round(100 * pr[2]))
    add("B", cond, 4, 4, pr[1], pr[2],
        sort(c(base, 0.01, 0.24, round(phistar(4, 0.7), 3))),
        phistar(4, 0.7), "A", PHI_FOLD[[cond]])
  }
  ## Design C: PHQ-ADS layout, 9 and 7 items at .75; block B (7 items) is the
  ## weaker block. No symmetric threshold exists for unequal lengths, so the
  ## grid is 0 to .9 plus the three trial correlations of Kroenke et al. (2016).
  add("C", "C_9_7_l75", 9, 7, 0.75, 0.75, sort(c(base, 0.653, 0.865, 0.912)),
      NA_real_, "B", NA_real_)
  pops <- do.call(rbind, rows)
  pops$pop_id <- seq_len(nrow(pops))
  pops
}

make_cells <- function(pops = make_populations()) {
  g <- expand.grid(format = FORMATS, N = NS, pop_id = pops$pop_id,
                   stringsAsFactors = FALSE)[, c("pop_id", "N", "format")]
  cells <- merge(g, pops, by = "pop_id", sort = FALSE)
  cells <- cells[order(cells$pop_id, cells$N, match(cells$format, FORMATS)), ]
  cells$cell_id <- seq_len(nrow(cells))
  ## One seed per population x N. All formats of that population and N share the
  ## same latent draws (common random numbers); replication r uses seed + r.
  cells$seed <- cells$pop_id * 100000 + match(cells$N, NS) * 10000
  cells$reps <- unname(REPS[cells$format])
  cells$uls <- cells$design == "A" & cells$format == "cont" & cells$N %in% ULS_N
  cells$ms <- cells$design %in% c("B", "C") & cells$N %in% MS_N
  rownames(cells) <- NULL
  cells
}

## ------------------------------------------------------------ population ----

item_names <- function(p) paste0("y", seq_len(p))

pop_sigma <- function(pop) {
  mA <- pop$mA; mB <- pop$mB; p <- mA + mB
  Lam <- matrix(0, p, 2)
  Lam[seq_len(mA), 1] <- pop$lamA
  Lam[mA + seq_len(mB), 2] <- pop$lamB
  S <- Lam %*% matrix(c(1, pop$phi, pop$phi, 1), 2) %*% t(Lam)
  diag(S) <- 1
  dimnames(S) <- list(item_names(p), item_names(p))
  S
}

gen_latent <- function(N, Sigma) {
  p <- nrow(Sigma)
  Y <- matrix(rnorm(N * p), N, p) %*% chol(Sigma)
  colnames(Y) <- colnames(Sigma)
  Y
}

discretize <- function(Y, thr) {
  X <- matrix(findInterval(Y, thr), nrow(Y), ncol(Y))
  colnames(X) <- colnames(Y)
  X
}

## Exact covariance matrix of 0..K scores s = sum_a 1(y > tau_a) when the
## latent responses are standard bivariate normal with correlation Sigma[i, k].
score_cov <- function(Sigma, thr) {
  p <- nrow(Sigma)
  mu <- sum(pnorm(-thr))
  ex2 <- sum(pnorm(-outer(thr, thr, pmax)))
  ta <- rep(thr, times = length(thr))
  tb <- rep(thr, each = length(thr))
  C <- matrix(NA_real_, p, p, dimnames = dimnames(Sigma))
  for (i in seq_len(p - 1)) for (k in (i + 1):p) {
    rho <- Sigma[i, k]
    e <- if (abs(rho) < 1e-12) sum(pnorm(-ta) * pnorm(-tb)) else
      sum(pbivnorm::pbivnorm(-ta, -tb, rho = rep(rho, length(ta))))
    C[i, k] <- C[k, i] <- e - mu^2
  }
  diag(C) <- ex2 - mu^2
  C
}

## ---------------------------------------------------------------- models ----

model_1f <- function(p) paste("f =~", paste(item_names(p), collapse = " + "))

model_2f <- function(mA, mB) {
  nm <- item_names(mA + mB)
  paste0("fA =~ ", paste(nm[seq_len(mA)], collapse = " + "), "\n",
         "fB =~ ", paste(nm[mA + seq_len(mB)], collapse = " + "))
}

## One-factor model with user start values (loadings l0, uniquenesses psi0)
model_1f_start <- function(l0, psi0) {
  nm <- item_names(length(l0))
  paste0("f =~ ", paste(sprintf("start(%.6f)*%s", l0, nm), collapse = " + "), "\n",
         paste(sprintf("%s ~~ start(%.6f)*%s", nm, psi0, nm), collapse = "\n"))
}

## --------------------------------------------------------------- fitting ----

quiet <- function(expr) withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))
try_fit <- function(expr) tryCatch(quiet(expr), error = function(e) NULL)

fit_ml_cov <- function(model, S, N, test = TRUE, rescale = TRUE) {
  cfa(model, sample.cov = S, sample.nobs = N, std.lv = TRUE, estimator = "ML",
      se = "none", test = if (test) "standard" else "none",
      h1 = test, baseline = test, sample.cov.rescale = rescale)
}

fit_uls_cov <- function(model, S, N, rescale = TRUE) {
  cfa(model, sample.cov = S, sample.nobs = N, std.lv = TRUE, estimator = "ULS",
      se = "none", test = "none", h1 = FALSE, baseline = FALSE,
      sample.cov.rescale = rescale)
}

fit_wlsmv <- function(model, X) {
  X <- as.data.frame(X)
  cfa(model, data = X, ordered = names(X), estimator = "WLSMV", std.lv = TRUE,
      se = "none")
}

## WLSMV from summary statistics (used for population values): polychoric
## matrix R, thresholds th (length p * K), and Gamma / W estimated elsewhere.
fit_wlsmv_summary <- function(model, R, th, K, Gamma, W, N) {
  p <- nrow(R); nm <- item_names(p)
  attr(th, "th.idx") <- rep(seq_len(p), each = K)
  names(th) <- paste0(rep(nm, each = K), "|t", rep(seq_len(K), p))
  thl <- paste(vapply(seq_len(p), function(i) {
    paste0(nm[i], " | ", paste(sprintf("start(%.6f)*t%d", th[(i - 1) * K + seq_len(K)],
                                       seq_len(K)), collapse = " + "))
  }, ""), collapse = "\n")
  cfa(paste(model, thl, sep = "\n"), sample.cov = R, sample.th = th,
      sample.mean = setNames(rep(0, p), nm), sample.nobs = N, NACOV = Gamma,
      WLS.V = W, estimator = "WLSMV", std.lv = TRUE, se = "none", ordered = nm)
}

## ------------------------------------------------------------- outcomes ----

## Proposition 2: rate at which F changes when the zero uniqueness of item j is
## released, computed on correlation matrix R. Negative: releasing pays, so a
## Heywood case on j cannot be a local minimum of the bounded problem.
prop2_rate <- function(R, j) {
  r <- R[-j, j]
  w <- r / (1 - r^2)
  M <- R[-j, -j] - outer(r, r)
  diag(M) <- 0
  -sum(outer(w, w) * M)
}

share_A <- function(lA, lB, ref) {
  if (is.null(ref) || anyNA(c(lA, lB))) return(NA_real_)
  w <- solve(matrix(c(1, ref$phi, ref$phi, 1), 2), c(lA / ref$lamA, lB / ref$lamB))
  w[1] / sum(w)
}

FIT_COLS <- c("est", "model", "ok", "converged", "fx", "chisq", "df", "pvalue",
              "cfi", "tli", "rmsea", "srmr", "chisq_sc", "pvalue_sc", "cfi_sc",
              "tli_sc", "rmsea_sc", "cfi_rob", "tli_rob", "rmsea_rob", "lA", "lB",
              "shareA", "phi_hat", "min_psi", "heywood_item", "p2_rate")

empty_row <- function(est, model) {
  row <- setNames(as.list(rep(NA_real_, length(FIT_COLS))), FIT_COLS)
  row$est <- est; row$model <- model; row$ok <- FALSE; row$converged <- FALSE
  row
}

FM_ML <- c(chisq = "chisq", df = "df", pvalue = "pvalue", cfi = "cfi",
           tli = "tli", rmsea = "rmsea", srmr = "srmr")
FM_CAT <- c(FM_ML, chisq_sc = "chisq.scaled", pvalue_sc = "pvalue.scaled",
            cfi_sc = "cfi.scaled", tli_sc = "tli.scaled", rmsea_sc = "rmsea.scaled",
            cfi_rob = "cfi.robust", tli_rob = "tli.robust", rmsea_rob = "rmsea.robust")

## Sign-normalized standardized loadings of a one-factor fit
std_loadings_1f <- function(fit) {
  l <- as.numeric(lavInspect(fit, "std")$lambda[, 1])
  if (sum(l) < 0) l <- -l
  l
}

extract_fit <- function(fit, est, model, mA, mB, R = NULL, ref = NULL,
                        measures = NULL) {
  row <- empty_row(est, model)
  if (is.null(fit)) return(row)
  row$ok <- TRUE
  row$converged <- isTRUE(lavInspect(fit, "converged"))
  row$fx <- fit@optim$fx
  if (!is.null(measures)) {
    fm <- tryCatch(quiet(fitMeasures(fit, unname(measures))), error = function(e) NULL)
    if (!is.null(fm)) for (k in names(measures)) row[[k]] <- unname(fm[measures[[k]]])
  }
  std <- tryCatch(lavInspect(fit, "std"), error = function(e) NULL)
  if (is.null(std)) return(row)
  ia <- seq_len(mA); ib <- mA + seq_len(mB)
  if (model == "1F") {
    l <- as.numeric(std$lambda[, 1])
    if (sum(l) < 0) l <- -l
    ## Block-wise sign normalization: when a sample's cross-block correlations
    ## are negative the fitted factor can load with opposite signs on the two
    ## blocks (a reflected composite), which a global sign flip would misread.
    row$lA <- abs(mean(l[ia])); row$lB <- abs(mean(l[ib]))
    row$shareA <- share_A(row$lA, row$lB, ref)
    psi <- 1 - l^2
    row$min_psi <- min(psi)
    row$heywood_item <- which.min(psi)
    if (!is.null(R) && row$min_psi < 0.05) row$p2_rate <- prop2_rate(R, which.min(psi))
  } else {
    L <- std$lambda
    sa <- sign(sum(L[ia, 1])); sb <- sign(sum(L[ib, 2]))
    row$lA <- sa * mean(L[ia, 1]); row$lB <- sb * mean(L[ib, 2])
    row$phi_hat <- sa * sb * std$psi[1, 2]
    psi <- 1 - rowSums(L^2)
    row$min_psi <- min(psi)
    row$heywood_item <- which.min(psi)
  }
  row
}

rows_to_df <- function(rows) {
  do.call(rbind, lapply(rows, function(r) as.data.frame(r[FIT_COLS], stringsAsFactors = FALSE)))
}

## ---------------------------------------------------------- multi-start ----

## Cluster converged one-factor solutions into distinct stationary points.
distinct_solutions <- function(sols, tol_l = 2e-3) {
  sols <- Filter(function(s) isTRUE(s$converged), sols)
  if (!length(sols)) return(list())
  sols <- sols[order(vapply(sols, `[[`, 0, "fx"))]
  out <- list()
  for (s in sols) {
    hit <- FALSE
    for (k in seq_along(out)) {
      if (max(abs(out[[k]]$l - s$l)) < tol_l) {
        out[[k]]$count <- out[[k]]$count + 1
        out[[k]]$starts <- c(out[[k]]$starts, s$start)
        hit <- TRUE; break
      }
    }
    if (!hit) { s$count <- 1; s$starts <- s$start; out[[length(out) + 1]] <- s }
  }
  out
}

## perm: item order used in the fit (fit item k = original item perm[k]);
## loadings are mapped back to the original order.
one_solution <- function(fit, start, mA, mB, ref, perm = NULL) {
  if (is.null(fit)) return(list(converged = FALSE, start = start))
  l <- std_loadings_1f(fit)
  if (!is.null(perm)) l <- l[order(perm)]
  list(converged = isTRUE(lavInspect(fit, "converged")), fx = fit@optim$fx, l = l,
       shareA = share_A(abs(mean(l[seq_len(mA)])), abs(mean(l[mA + seq_len(mB)])), ref),
       min_psi = min(1 - l^2), start = start)
}

## Fit the one-factor ML model on correlation matrix R from lavaan's default
## start twice, with the items in design order (block A listed first, start
## id 0) and with the blocks reversed (block B first, start id -1), and from
## nstart random starts (loadings U(0,1), uniquenesses U(.2,1)). lavaan's
## default starting values use the first-listed item as reference, so the two
## default fits can land in different basins.
## weak: "A" or "B", the block whose minimum is the spurious one in the
## population. Returns one row of Question 3 outcomes.
multistart_1f <- function(R, N, mA, mB, ref, weak, nstart = NSTART, tol_fx = 1e-7) {
  p <- nrow(R)
  perm <- c(mA + seq_len(mB), seq_len(mA))
  Rr <- R[perm, perm]; dimnames(Rr) <- dimnames(R)
  sols <- vector("list", nstart + 2)
  sols[[1]] <- one_solution(try_fit(fit_ml_cov(model_1f(p), R, N, test = FALSE)), 0L, mA, mB, ref)
  sols[[2]] <- one_solution(try_fit(fit_ml_cov(model_1f(p), Rr, N, test = FALSE)), -1L, mA, mB, ref, perm)
  for (s in seq_len(nstart)) {
    fs <- try_fit(fit_ml_cov(model_1f_start(runif(p), runif(p, 0.2, 1)), R, N, test = FALSE))
    sols[[s + 2]] <- one_solution(fs, s, mA, mB, ref)
  }
  d <- distinct_solutions(sols)
  n_conv <- sum(vapply(sols[-(1:2)], function(s) isTRUE(s$converged), TRUE))
  out <- list(n_conv = n_conv, n_distinct = length(d), n_admissible = NA_integer_,
              best_fx = NA_real_, best_shareA = NA_real_, def_conv = isTRUE(sols[[1]]$converged),
              def_nonglobal = NA, def_on_weak = NA, def_rev_nonglobal = NA, def_rev_on_weak = NA,
              global_on_weak = NA, spurious_weak = NA, weak_fx_gap = NA_real_,
              weak_shareA = NA_real_, frac_starts_weak = NA_real_, any_heywood_min = NA)
  if (!length(d)) return(out)
  adm <- Filter(function(s) s$min_psi > 0, d)
  out$n_admissible <- length(adm)
  out$any_heywood_min <- any(vapply(d, function(s) s$min_psi <= 0, TRUE))
  best <- d[[1]]
  on_weak <- function(s) if (weak == "A") s$shareA > 0.5 else s$shareA < 0.5
  out$best_fx <- best$fx
  out$best_shareA <- best$shareA
  out$global_on_weak <- on_weak(best)
  nonglobal_weak <- Filter(function(s) s$fx > best$fx + tol_fx && s$min_psi > 0 && on_weak(s), d)
  out$spurious_weak <- length(nonglobal_weak) > 0
  if (out$spurious_weak) {
    wk <- nonglobal_weak[[1]]
    out$weak_fx_gap <- wk$fx - best$fx
    out$weak_shareA <- wk$shareA
    out$frac_starts_weak <- sum(wk$starts > 0) / max(n_conv, 1)
  } else out$frac_starts_weak <- 0
  for (st in c(0L, -1L)) {
    k0 <- which(vapply(d, function(s) st %in% s$starts, TRUE))
    if (length(k0)) {
      s0 <- d[[k0]]
      ng <- s0$fx > best$fx + tol_fx
      if (st == 0L) { out$def_nonglobal <- ng; out$def_on_weak <- ng && on_weak(s0) }
      else { out$def_rev_nonglobal <- ng; out$def_rev_on_weak <- ng && on_weak(s0) }
    }
  }
  out
}

## --------------------------------------------------- reference mappings ----

ref_latent <- function(cell) list(lamA = cell$lamA, lamB = cell$lamB, phi = cell$phi)

## Population two-factor ML solution on the 0..K score scale
ref_scores <- function(Sigma, thr, mA, mB) {
  C <- score_cov(Sigma, thr)
  f2 <- try_fit(fit_ml_cov(model_2f(mA, mB), C, 1e6, test = FALSE, rescale = FALSE))
  r <- extract_fit(f2, "MLs", "2F", mA, mB)
  list(lamA = r$lA, lamB = r$lB, phi = r$phi_hat)
}

## ------------------------------------------------------------- one cell ----

run_cell <- function(cell, reps = cell$reps, cores = 1, nstart = NSTART) {
  Sigma <- pop_sigma(cell)
  p <- nrow(Sigma); mA <- cell$mA; mB <- cell$mB; N <- cell$N
  thr <- THRESHOLDS[[cell$format]]
  m1 <- model_1f(p); m2 <- model_2f(mA, mB)
  rl <- ref_latent(cell)
  rs <- if (!is.null(thr)) ref_scores(Sigma, thr, mA, mB) else NULL

  one_rep <- function(r) {
    t0 <- proc.time()[["elapsed"]]
    set.seed(cell$seed + r)
    Y <- gen_latent(N, Sigma)
    rows <- list(); ms <- list()
    if (cell$format == "cont") {
      S <- cov(Y); R <- cov2cor(S)
      rows[[1]] <- extract_fit(try_fit(fit_ml_cov(m1, S, N)), "ML", "1F", mA, mB, R, rl, FM_ML)
      rows[[2]] <- extract_fit(try_fit(fit_ml_cov(m2, S, N)), "ML", "2F", mA, mB, R, rl, FM_ML)
      if (cell$uls) rows[[3]] <- extract_fit(try_fit(fit_uls_cov(m1, S, N)), "ULS", "1F", mA, mB, R, rl)
      if (cell$ms) ms[[1]] <- c(est = "ML", multistart_1f(R, N, mA, mB, rl, cell$weak, nstart))
    } else {
      X <- discretize(Y, thr)
      w1 <- try_fit(fit_wlsmv(m1, X))
      Rp <- if (!is.null(w1)) tryCatch(lavInspect(w1, "sampstat")$cov, error = function(e) NULL)
      rows[[1]] <- extract_fit(w1, "WLSMV", "1F", mA, mB, Rp, rl, FM_CAT)
      rows[[2]] <- extract_fit(try_fit(fit_wlsmv(m2, X)), "WLSMV", "2F", mA, mB, Rp, rl, FM_CAT)
      Sx <- cov(X); Rx <- cov2cor(Sx)
      rows[[3]] <- extract_fit(try_fit(fit_ml_cov(m1, Sx, N)), "MLs", "1F", mA, mB, Rx, rs, FM_ML)
      rows[[4]] <- extract_fit(try_fit(fit_ml_cov(m2, Sx, N)), "MLs", "2F", mA, mB, Rx, rs, FM_ML)
      if (cell$ms) ms[[1]] <- c(est = "MLs", multistart_1f(Rx, N, mA, mB, rs, cell$weak, nstart))
    }
    fits <- rows_to_df(rows)
    fits$rep <- r
    msdf <- if (length(ms)) {
      d <- as.data.frame(ms[[1]], stringsAsFactors = FALSE); d$rep <- r; d
    } else NULL
    list(fits = fits, ms = msdf, time = proc.time()[["elapsed"]] - t0)
  }

  t0 <- proc.time()[["elapsed"]]
  res <- if (cores > 1) parallel::mclapply(seq_len(reps), one_rep, mc.cores = cores,
                                           mc.preschedule = TRUE)
         else lapply(seq_len(reps), one_rep)
  failed <- vapply(res, function(x) inherits(x, "try-error") || is.null(x$fits), TRUE)
  res <- res[!failed]
  id <- cell[c("cell_id", "pop_id", "design", "cond", "phi", "N", "format")]
  fits <- do.call(rbind, lapply(res, `[[`, "fits"))
  fits <- cbind(id[rep(1, nrow(fits)), ], fits, row.names = NULL)
  ms <- do.call(rbind, lapply(res, `[[`, "ms"))
  if (!is.null(ms)) ms <- cbind(id[rep(1, nrow(ms)), ], ms, row.names = NULL)
  list(cell = cell, fits = fits, ms = ms, ref_scores = rs,
       rep_time = vapply(res, `[[`, 0, "time"), n_failed = sum(failed),
       wall = proc.time()[["elapsed"]] - t0, cores = cores,
       lavaan = as.character(packageVersion("lavaan")), R = R.version.string)
}

## ----------------------------------------------------- population values ----

## ML-type population fit of covariance S: global one-factor minimum from the
## default start, two tilted starts and nstart random starts, plus the
## two-factor model. Returns population indices and the distinct minima.
pop_indices <- function(Fm, S, df_M) {
  p <- nrow(S)
  FB <- sum(log(diag(S))) - as.numeric(determinant(S, logarithm = TRUE)$modulus)
  dfB <- p * (p - 1) / 2
  c(F = Fm, FB = FB, cfi = 1 - Fm / FB, tli = 1 - (Fm / df_M) / (FB / dfB),
    rmsea = sqrt(max(Fm, 0) / df_M))
}

pop_ml <- function(S, mA, mB, ref, weak, est = "ML", nstart = 20, Nnom = 1e6) {
  p <- nrow(S)
  R <- cov2cor(S)
  fitter <- if (est == "ULS") function(m) fit_uls_cov(m, R, Nnom, rescale = FALSE)
            else function(m) fit_ml_cov(m, R, Nnom, test = TRUE, rescale = FALSE)
  starts <- list(NULL,
                 list(l = c(rep(0.9, mA), rep(0.1, mB)), psi = rep(0.5, p)),
                 list(l = c(rep(0.1, mA), rep(0.9, mB)), psi = rep(0.5, p)))
  for (s in seq_len(nstart)) starts[[length(starts) + 1]] <- list(l = runif(p), psi = runif(p, 0.2, 1))
  fits <- lapply(starts, function(st) try_fit(fitter(
    if (is.null(st)) model_1f(p) else model_1f_start(st$l, st$psi))))
  sols <- Map(function(f, k) one_solution(f, k, mA, mB, ref), fits, seq_along(fits) - 1L)
  d <- distinct_solutions(sols)
  kbest <- which(vapply(sols, function(s) isTRUE(s$converged) && isTRUE(all.equal(s$fx, d[[1]]$fx)), TRUE))[1]
  fbest <- fits[[kbest]]
  df1 <- p * (p + 1) / 2 - 2 * p
  out <- list(minima = d, best = d[[1]])
  if (est == "ULS") {
    out$row1 <- list(F = 2 * d[[1]]$fx)
  } else {
    pi1 <- pop_indices(2 * d[[1]]$fx, R, df1)
    out$row1 <- c(as.list(pi1), srmr = unname(quiet(fitMeasures(fbest, "srmr"))))
    f2 <- try_fit(fit_ml_cov(model_2f(mA, mB), R, Nnom, test = TRUE, rescale = FALSE))
    pi2 <- pop_indices(2 * f2@optim$fx, R, df1 - 1)
    out$row2 <- c(as.list(pi2), srmr = unname(quiet(fitMeasures(f2, "srmr"))))
    out$fit2 <- f2
  }
  if (!is.na(weak)) {
    on_weak <- function(s) if (weak == "A") s$shareA > 0.5 else s$shareA < 0.5
    wk <- Filter(function(s) s$fx > d[[1]]$fx + 1e-9 && on_weak(s), d)
    out$weak_exists <- length(wk) > 0
    out$weak_shareA <- if (length(wk)) wk[[1]]$shareA else NA_real_
    out$weak_gap <- if (length(wk)) 2 * (wk[[1]]$fx - d[[1]]$fx) else NA_real_
  }
  out
}
