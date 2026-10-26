## run_population.R
## Population values for every population and estimator of Section 7.3.
##
##   Continuous, ML:   exact. One-factor global minimum by multi-start on Sigma;
##                     CFI, TLI and RMSEA from the closed-form identities
##                     CFI = 1 - F/F_B, TLI = 1 - (F/df)/(F_B/df_B),
##                     RMSEA = sqrt(F/df); SRMR from lavaan at the population fit.
##                     All distinct one-factor minima are recorded, which checks
##                     the fold points of Proposition 1.
##   Continuous, ULS:  exact loadings and share (Theorem 2 reference).
##   Ordinal, ML on 0-3 scores: exact, from the covariance matrix of the scores
##                     (bivariate normal orthant probabilities).
##   Ordinal, WLSMV:   true polychoric matrix and thresholds; the asymptotic
##                     covariance matrix Gamma and weight matrix are estimated
##                     once from a sample of size --ngamma; nominal N = 1e6.
##
## Usage:
##   Rscript run_population.R [--out population.rds] [--ngamma 200000]
##                            [--pops 1-83] [--cores 1]

.here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
source(file.path(.here, "R", "sim_functions.R"))
source(file.path(.here, "R", "cli.R"))
source(file.path(.here, "R", "thresholds.R"))

args <- parse_args(commandArgs(TRUE), list(out = "population.rds", ngamma = "200000",
                                            pops = "", cores = "1", nstart = "20"))
pops <- make_populations()
ids <- if (nzchar(args$pops)) parse_ids(args$pops) else pops$pop_id
ngamma <- as.numeric(args$ngamma)
cores <- as.integer(args$cores)
nstart <- as.integer(args$nstart)

model_1f_lstart <- function(l0) {
  nm <- item_names(length(l0))
  paste("f =~", paste(sprintf("start(%.6f)*%s", l0, nm), collapse = " + "))
}

row_of <- function(pop, format, est, model, vals) {
  base <- list(pop_id = pop$pop_id, design = pop$design, cond = pop$cond, phi = pop$phi,
               format = format, est = est, model = model)
  cols <- c("F", "FB", "cfi", "tli", "rmsea", "srmr", "cfi_sc", "tli_sc", "rmsea_sc",
            "cfi_rob", "tli_rob", "rmsea_rob", "lA", "lB", "shareA", "n_minima",
            "weak_exists", "weak_shareA", "weak_gap", "refA", "refB", "refphi")
  v <- setNames(as.list(rep(NA_real_, length(cols))), cols)
  for (k in intersect(names(vals), cols)) v[[k]] <- vals[[k]]
  as.data.frame(c(base, v), stringsAsFactors = FALSE)
}

ml_rows <- function(pop, format, est, S, ref) {
  mA <- pop$mA; mB <- pop$mB
  r <- pop_ml(S, mA, mB, ref, pop$weak, est = est, nstart = nstart)
  b <- r$best
  common <- list(lA = mean(b$l[seq_len(mA)]), lB = mean(b$l[mA + seq_len(mB)]),
                 shareA = b$shareA, n_minima = length(r$minima),
                 weak_exists = if (is.null(r$weak_exists)) NA else r$weak_exists,
                 weak_shareA = if (is.null(r$weak_shareA)) NA else r$weak_shareA,
                 weak_gap = if (is.null(r$weak_gap)) NA else r$weak_gap,
                 refA = ref$lamA, refB = ref$lamB, refphi = ref$phi)
  out <- row_of(pop, format, est, "1F", c(r$row1, common))
  if (!is.null(r$row2)) {
    e2 <- extract_fit(r$fit2, est, "2F", mA, mB)
    out <- rbind(out, row_of(pop, format, est, "2F", c(r$row2, lA = e2$lA, lB = e2$lB)))
  }
  out
}

wlsmv_rows <- function(pop, format, thr) {
  mA <- pop$mA; mB <- pop$mB; p <- mA + mB; K <- length(thr)
  Sigma <- pop_sigma(pop)
  set.seed(9e6 + 10 * pop$pop_id + match(format, FORMATS))
  X <- discretize(gen_latent(ngamma, Sigma), thr)
  f0 <- quiet(cfa(model_1f(p), data = as.data.frame(X), ordered = colnames(X),
                  estimator = "WLSMV", std.lv = TRUE, se = "none", do.fit = FALSE))
  G <- lavInspect(f0, "gamma"); W <- lavInspect(f0, "wls.v")
  rm(X, f0); gc(verbose = FALSE)
  ## Replace the correlation weights by their exact values (inverse asymptotic
  ## variances, R/thresholds.R). Sampling noise in the estimated weights acts as
  ## a small asymmetry, which near a threshold moves the tilt by its cube root
  ## (Theorem 3). Threshold weights do not affect the estimates (saturated).
  pr <- which(lower.tri(Sigma), arr.ind = TRUE)        # lavaan's vech order
  rho <- round(Sigma[pr], 12)
  av <- setNames(vapply(unique(rho), avar_polychoric, 0, t1 = thr, t2 = thr), unique(rho))
  dW <- diag(W); dW[p * K + seq_along(rho)] <- 1 / av[as.character(rho)]
  W <- diag(dW)
  th <- rep(thr, p)
  ref <- ref_latent(pop)
  fitter <- function(model) try_fit(fit_wlsmv_summary(model, Sigma, th, K, G, W, 1e6))
  starts <- list(NULL, c(rep(0.9, mA), rep(0.1, mB)), c(rep(0.1, mA), rep(0.9, mB)))
  for (s in seq_len(5)) starts[[length(starts) + 1]] <- runif(p)
  fits <- lapply(starts, function(l0) fitter(if (is.null(l0)) model_1f(p) else model_1f_lstart(l0)))
  sols <- Map(function(f, k) one_solution(f, k, mA, mB, ref), fits, seq_along(fits) - 1L)
  d <- distinct_solutions(sols)
  kbest <- which(vapply(sols, function(s) isTRUE(s$converged) && isTRUE(all.equal(s$fx, d[[1]]$fx)), TRUE))[1]
  e1 <- extract_fit(fits[[kbest]], "WLSMV", "1F", mA, mB, NULL, ref, FM_CAT)
  f2 <- fitter(model_2f(mA, mB))
  e2 <- extract_fit(f2, "WLSMV", "2F", mA, mB, NULL, ref, FM_CAT)
  tidy <- function(e) list(F = 2 * e$fx, cfi = e$cfi, tli = e$tli, rmsea = e$rmsea, srmr = e$srmr,
                           cfi_sc = e$cfi_sc, tli_sc = e$tli_sc, rmsea_sc = e$rmsea_sc,
                           cfi_rob = e$cfi_rob, tli_rob = e$tli_rob, rmsea_rob = e$rmsea_rob,
                           lA = e$lA, lB = e$lB, shareA = e$shareA)
  wk <- list()
  if (!is.na(pop$weak)) {
    on_weak <- function(s) if (pop$weak == "A") s$shareA > 0.5 else s$shareA < 0.5
    wk <- Filter(function(s) s$fx > d[[1]]$fx + 1e-9 && on_weak(s), d)
  }
  wv <- if (is.na(pop$weak)) list() else
    list(weak_exists = length(wk) > 0, weak_shareA = if (length(wk)) wk[[1]]$shareA else NA,
         weak_gap = if (length(wk)) 2 * (wk[[1]]$fx - d[[1]]$fx) else NA)
  rbind(row_of(pop, format, "WLSMV", "1F", c(tidy(e1), n_minima = length(d), wv,
                                              refA = ref$lamA, refB = ref$lamB, refphi = ref$phi)),
        row_of(pop, format, "WLSMV", "2F", tidy(e2)))
}

one_pop <- function(id) {
  pop <- pops[pops$pop_id == id, ]
  t0 <- proc.time()[["elapsed"]]
  set.seed(8e6 + id)
  Sigma <- pop_sigma(pop)
  out <- list(ml_rows(pop, "cont", "ML", Sigma, ref_latent(pop)),
              ml_rows(pop, "cont", "ULS", Sigma, ref_latent(pop)))
  for (fmt in c("ord_sym", "ord_skew")) {
    thr <- THRESHOLDS[[fmt]]
    C <- score_cov(Sigma, thr)
    rs <- ref_scores(Sigma, thr, pop$mA, pop$mB)
    out[[length(out) + 1]] <- ml_rows(pop, fmt, "MLs", C, rs)
    out[[length(out) + 1]] <- wlsmv_rows(pop, fmt, thr)
  }
  res <- do.call(rbind, out)
  res$seconds <- proc.time()[["elapsed"]] - t0
  message(sprintf("population %d (%s, phi = %.3f) done in %.1f s", id, pop$cond, pop$phi, res$seconds[1]))
  res
}

res <- if (cores > 1) parallel::mclapply(ids, one_pop, mc.cores = cores) else lapply(ids, one_pop)
bad <- vapply(res, inherits, TRUE, "try-error")
if (any(bad)) warning("failed populations: ", paste(ids[bad], collapse = ", "))
pop_values <- do.call(rbind, res[!bad])
rownames(pop_values) <- NULL
saveRDS(pop_values, args$out)
message("wrote ", args$out, " (", nrow(pop_values), " rows)")
