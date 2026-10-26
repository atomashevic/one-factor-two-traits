## summarise_sim.R
## Aggregates the cell files into the tables and figures of Section 7.3.
##
## Usage:
##   Rscript summarise_sim.R [--results results] [--population population.rds]
##                           [--out summary]
##
## Writes to --out:
##   q1_fit.csv          Question 1: mean, SD, bias, RMSE and acceptance rate of
##                       every fit index (and F) against its population value
##   q2_tilt.csv         Question 2: tilt criterion, tilt and win rates, mean
##                       signed and absolute loading difference and tilt
##   q3_minima.csv       Question 3: multi-start outcomes (Designs B and C)
##   heywood.csv         nonconvergence, Heywood cases and the sign of the
##                       Proposition 2 rate
##   fig_acceptance.pdf  acceptance curves (draft of figure 1)
##   fig_tilt.pdf        tilt probability (draft of figure 2) and the
##                       averaging artifact
##   fig_spurious.pdf    spurious-minimum frequency (draft of figure 3)

.here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
source(file.path(.here, "R", "sim_functions.R"))
source(file.path(.here, "R", "cli.R"))
source(file.path(.here, "R", "thresholds.R"))

args <- parse_args(commandArgs(TRUE), list(results = "results", population = "population.rds",
                                            out = "summary"))
dir.create(args$out, showWarnings = FALSE, recursive = TRUE)
outf <- function(f) file.path(args$out, f)

files <- list.files(args$results, "^cell_[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no cell files in ", args$results)
res <- lapply(files, readRDS)
cells <- do.call(rbind, lapply(res, `[[`, "cell"))
fits <- do.call(rbind, lapply(res, `[[`, "fits"))
ms <- do.call(rbind, lapply(res, `[[`, "ms"))
pop <- if (file.exists(args$population)) readRDS(args$population) else NULL
if (is.null(pop)) message("no population file: population columns will be NA")
message(sprintf("%d cells, %d fits, %d multi-start rows", nrow(cells), nrow(fits),
                if (is.null(ms)) 0 else nrow(ms)))

cell_info <- cells[, c("cell_id", "mA", "mB", "lamA", "lamB", "phistar", "weak", "phi_fold")]
fits <- merge(fits, cell_info, by = "cell_id", sort = FALSE)
fits$F <- 2 * fits$fx

pop_lookup <- function(d, col) {
  if (is.null(pop) || !col %in% names(pop)) return(rep(NA_real_, nrow(d)))
  fmt <- ifelse(d$est %in% c("ML", "ULS"), "cont", d$format)
  key_d <- paste(d$pop_id, fmt, d$est, d$model)
  key_p <- paste(pop$pop_id, pop$format, pop$est, pop$model)
  pop[[col]][match(key_d, key_p)]
}

grp <- function(d, keys) {
  if (anyNA(d[keys])) stop("grouping keys contain NA: ", paste(keys[vapply(d[keys], anyNA, TRUE)], collapse = ", "))
  interaction(d[keys], drop = TRUE, lex.order = TRUE)
}
first_rows <- function(d, g, keys) {
  i <- match(levels(g), as.character(g))
  out <- d[i, keys, drop = FALSE]; rownames(out) <- NULL; out
}

## ------------------------------------------------------------ Question 1 ----

keys <- c("cell_id", "pop_id", "design", "cond", "phi", "N", "format", "est", "model")
idx <- c("F", "cfi", "tli", "rmsea", "srmr", "cfi_sc", "tli_sc", "rmsea_sc",
         "cfi_rob", "tli_rob", "rmsea_rob")
cutoff <- function(k) switch(sub("_.*", "", k), cfi = c(.95, 1), tli = c(.95, 1),
                             rmsea = c(.06, -1), srmr = c(.08, -1), c(NA, NA))
conv <- fits[fits$converged %in% TRUE & fits$est != "ULS", ]
q1 <- list()
for (k in idx) {
  d <- conv[!is.na(conv[[k]]), ]
  if (!nrow(d)) next
  g <- grp(d, keys)
  base <- first_rows(d, g, keys)
  x <- d[[k]]
  base$index <- k
  base$n <- as.vector(tapply(x, g, length))
  base$mean <- as.vector(tapply(x, g, mean))
  base$sd <- as.vector(tapply(x, g, sd))
  base$pop <- pop_lookup(base, k)
  base$bias <- base$mean - base$pop
  base$rmse <- sqrt(as.vector(tapply(x, g, function(v) mean(v^2))) - 2 * base$mean * base$pop + base$pop^2)
  cu <- cutoff(k)
  base$cutoff <- cu[1]
  base$p_accept <- if (is.na(cu[1])) NA_real_ else
    as.vector(tapply(x, g, function(v) mean(if (cu[2] > 0) v >= cu[1] else v <= cu[1])))
  base$pop_accept <- if (is.na(cu[1])) NA else
    (if (cu[2] > 0) base$pop >= cu[1] else base$pop <= cu[1])
  q1[[k]] <- base
}
q1 <- do.call(rbind, q1)
q1 <- q1[order(q1$cell_id, q1$est, q1$model, match(q1$index, idx)), ]
write.csv(q1, outf("q1_fit.csv"), row.names = FALSE)

## ------------------------------------------------------------ Question 2 ----

## Recompute the share from block-wise sign-normalized loadings (result files
## written before this fix stored globally sign-normalized block means).
ref_tab <- do.call(rbind, lapply(res, function(r) data.frame(cell_id = r$cell$cell_id,
  rA = if (is.null(r$ref_scores)) NA_real_ else r$ref_scores$lamA,
  rB = if (is.null(r$ref_scores)) NA_real_ else r$ref_scores$lamB,
  rphi = if (is.null(r$ref_scores)) NA_real_ else r$ref_scores$phi)))
one <- fits[fits$model == "1F" & fits$converged %in% TRUE & !is.na(fits$lA), ]
one$lA <- abs(one$lA); one$lB <- abs(one$lB)
k <- match(one$cell_id, ref_tab$cell_id)
sc <- one$est == "MLs"
rA <- ifelse(sc, ref_tab$rA[k], one$lamA); rB <- ifelse(sc, ref_tab$rB[k], one$lamB)
rp <- ifelse(sc, ref_tab$rphi[k], one$phi)
u <- one$lA / rA; v <- one$lB / rB
wA <- (u - rp * v) / (1 - rp^2); wB <- (v - rp * u) / (1 - rp^2)
one$shareA <- wA / (wA + wB)
one$share <- pmax(one$shareA, 1 - one$shareA)
one$d <- one$lA - one$lB
one$tiltA <- one$shareA - 0.5
## Tilt criterion: 95th percentile of the share at phi = .9 for the same
## condition, N, format and estimator, where the population is a composite.
tk <- c("cond", "N", "format", "est")
ref9 <- one[abs(one$phi - 0.9) < 1e-9, ]
q95 <- tapply(ref9$share, interaction(ref9[tk], drop = TRUE), quantile, probs = 0.95)
one$q95 <- as.vector(q95[as.character(interaction(one[tk]))])
one$tilted <- one$share > one$q95
k2 <- c("cell_id", "pop_id", "design", "cond", "phi", "N", "format", "est")
g <- grp(one, k2)
q2 <- first_rows(one, g, k2)
q2$phistar <- cell_info$phistar[match(q2$cell_id, cell_info$cell_id)]
q2$phi_minus_phistar <- q2$phi - q2$phistar
## Estimator-specific population threshold in Design A: phi* for ML (Theorem 1),
## 1/p for ULS (Theorem 2), the diagonally weighted threshold for WLSMV
## (R/thresholds.R; NA = symmetric solution stable at every phi). ML on scores
## keeps phi* as its reference.
mm <- cell_info$mA[match(q2$cell_id, cell_info$cell_id)]
ll <- cell_info$lamA[match(q2$cell_id, cell_info$cell_id)]
wth <- local({
  u <- unique(data.frame(m = mm, lam = ll, format = q2$format)[q2$design == "A" & q2$est == "WLSMV", ])
  u$thr <- mapply(function(m, lam, fm) wlsmv_threshold(m, lam, THRESHOLDS[[fm]]), u$m, u$lam, u$format)
  u
})
q2$threshold <- ifelse(q2$design != "A", NA_real_,
  ifelse(q2$est == "ULS", 1 / (2 * mm),
  ifelse(q2$est == "WLSMV", wth$thr[match(paste(mm, ll, q2$format), paste(wth$m, wth$lam, wth$format))],
         q2$phistar)))
q2$phi_minus_threshold <- q2$phi - q2$threshold
q2$n <- as.vector(tapply(one$share, g, length))
q2$q95 <- as.vector(tapply(one$q95, g, `[`, 1))
q2$p_tilted <- as.vector(tapply(one$tilted, g, mean))
q2$p_win_A <- as.vector(tapply(one$tilted & one$shareA > 0.5, g, mean))
q2$p_win_B <- as.vector(tapply(one$tilted & one$shareA < 0.5, g, mean))
q2$p_follow_A <- as.vector(tapply(one$shareA > 0.5, g, mean))
q2$mean_share <- as.vector(tapply(one$share, g, mean))
q2$median_share <- as.vector(tapply(one$share, g, median))
q2$mean_tiltA <- as.vector(tapply(one$tiltA, g, mean))
q2$mean_abs_tilt <- as.vector(tapply(abs(one$tiltA), g, mean))
q2$mean_d <- as.vector(tapply(one$d, g, mean))
q2$mean_abs_d <- as.vector(tapply(abs(one$d), g, mean))
q2$model <- "1F"
q2$pop_shareA <- pop_lookup(q2, "shareA")
q2$pop_share <- pmax(q2$pop_shareA, 1 - q2$pop_shareA)
q2$pop_abs_d <- abs(pop_lookup(q2, "lA") - pop_lookup(q2, "lB"))
q2$model <- NULL
write.csv(q2, outf("q2_tilt.csv"), row.names = FALSE)

## ------------------------------------------------------------ Question 3 ----

q3 <- NULL
if (!is.null(ms) && nrow(ms)) {
  ms <- merge(ms, cell_info[, c("cell_id", "weak", "phi_fold")], by = "cell_id", sort = FALSE)
  k3 <- c("cell_id", "pop_id", "design", "cond", "phi", "N", "format", "est")
  g <- grp(ms, k3)
  q3 <- first_rows(ms, g, k3)
  q3$phi_fold <- cell_info$phi_fold[match(q3$cell_id, cell_info$cell_id)]
  q3$weak <- cell_info$weak[match(q3$cell_id, cell_info$cell_id)]
  m <- function(x) as.vector(tapply(x, g, function(v) mean(v, na.rm = TRUE)))
  q3$n <- as.vector(tapply(ms$rep, g, length))
  q3$mean_n_conv <- m(ms$n_conv)
  q3$p_two_minima <- m(ms$n_admissible >= 2)
  q3$p_spurious_weak <- m(ms$spurious_weak)
  ## Default start with the weaker block listed first or last
  wfirst <- ifelse(ms$weak == "A", ms$def_on_weak, ms$def_rev_on_weak)
  wlast <- ifelse(ms$weak == "A", ms$def_rev_on_weak, ms$def_on_weak)
  q3$p_default_nonglobal <- m(ms$def_nonglobal)
  q3$p_default_on_weak <- m(ms$def_on_weak)
  q3$p_default_on_weak_listed_first <- m(wfirst)
  q3$p_default_on_weak_listed_last <- m(wlast)
  q3$p_global_on_weak <- m(ms$global_on_weak)
  q3$mean_frac_starts_weak <- m(ms$frac_starts_weak)
  q3$p_heywood_min <- m(ms$any_heywood_min)
  q3$model <- "1F"
  q3$pop_weak_exists <- as.logical(pop_lookup(q3, "weak_exists"))
  q3$model <- NULL
  write.csv(q3, outf("q3_minima.csv"), row.names = FALSE)
}

## ------------------------------------------- convergence and Heywood cases --

kh <- c("cell_id", "pop_id", "design", "cond", "phi", "N", "format", "est", "model")
g <- grp(fits, kh)
hw <- first_rows(fits, g, kh)
hw$n <- as.vector(tapply(fits$rep, g, length))
hw$p_error <- as.vector(tapply(!fits$ok, g, mean))
hw$p_nonconverged <- as.vector(tapply(!(fits$converged %in% TRUE), g, mean))
hw$p_heywood <- as.vector(tapply(fits$min_psi <= 0, g, function(v) mean(v, na.rm = TRUE)))
hw$p_near_heywood <- as.vector(tapply(fits$min_psi < 0.01, g, function(v) mean(v, na.rm = TRUE)))
hw$n_rate <- as.vector(tapply(fits$min_psi <= 0 & !is.na(fits$p2_rate), g, sum))
hw$p_rate_nonneg <- as.vector(tapply(ifelse(fits$min_psi <= 0, fits$p2_rate >= 0, NA), g,
                                     function(v) if (all(is.na(v))) NA else mean(v, na.rm = TRUE)))
write.csv(hw, outf("heywood.csv"), row.names = FALSE)

## --------------------------------------------------------------- figures ----

pal <- c("250" = "#9ecae1", "500" = "#4292c6", "1000" = "#2171b5", "2000" = "#08306b")
lines_by_N <- function(d, x, y, lty = 1, pch = 16) {
  for (n in sort(unique(d$N))) {
    s <- d[d$N == n, ]; s <- s[order(s[[x]]), ]
    lines(s[[x]], s[[y]], col = pal[as.character(n)], lty = lty, lwd = 1.6)
    points(s[[x]], s[[y]], col = pal[as.character(n)], pch = pch, cex = 0.6)
  }
}
legend_N <- function(Ns, where = "topleft") legend(where, legend = paste("N =", Ns),
  col = pal[as.character(Ns)], lwd = 1.6, bty = "n", cex = 0.7)
pop_cross <- function(d) {
  ## phi at which the population verdict changes, by linear interpolation
  d <- d[order(d$phi), ]
  s <- which(diff(d$pop_accept) != 0)
  if (!length(s) || anyNA(d$pop[c(s, s + 1)])) return(NA)
  i <- s[1]; c0 <- d$cutoff[1]
  d$phi[i] + (c0 - d$pop[i]) * (d$phi[i + 1] - d$phi[i]) / (d$pop[i + 1] - d$pop[i])
}
conds <- unique(q1$cond)

pdf(outf("fig_acceptance.pdf"), width = 11, height = 8.5)
arms <- list(c("cont", "ML", ""), c("ord_sym", "WLSMV", "_rob"), c("ord_skew", "WLSMV", "_rob"),
             c("ord_sym", "WLSMV", "_sc"), c("ord_skew", "WLSMV", "_sc"),
             c("ord_sym", "MLs", ""), c("ord_skew", "MLs", ""))
for (a in arms) {
  dd <- q1[q1$format == a[1] & q1$est == a[2] & q1$model == "1F", ]
  if (!nrow(dd)) next
  cs <- intersect(conds, unique(dd$cond))
  par(mfrow = c(4, length(cs)), mar = c(3, 3, 2, 0.5), mgp = c(1.8, 0.5, 0), oma = c(0, 0, 2, 0))
  for (k in c("cfi", "tli", "rmsea", "srmr")) {
    kk <- if (k == "srmr") k else paste0(k, a[3])
    for (cd in cs) {
      d <- dd[dd$index == kk & dd$cond == cd, ]
      plot(NA, xlim = c(0, 0.95), ylim = c(0, 1), xlab = "phi", ylab = "P(accept)",
           main = sprintf("%s, %s", toupper(kk), cd), cex.main = 0.8)
      if (nrow(d)) {
        lines_by_N(d, "phi", "p_accept")
        pc <- pop_cross(d[d$N == max(d$N), ])
        if (!is.na(pc)) abline(v = pc, lty = 2, col = "grey40")
      }
    }
  }
  legend_N(sort(unique(dd$N)))
  mtext(sprintf("Acceptance of the one-factor model: %s, %s%s (dashed: population verdict changes)",
                a[1], a[2], c("_rob" = " robust", "_sc" = " scaled", " ")[[if (nzchar(a[3])) a[3] else 3]]),
        outer = TRUE, cex = 0.9)
}
invisible(dev.off())

pdf(outf("fig_tilt.pdf"), width = 11, height = 8.5)
dA <- q2[q2$design == "A", ]
thr_lines <- function(d, which) {
  ## vertical lines at the population thresholds: ML phi* (dashed), ULS 1/p
  ## (dotted), WLSMV (dot-dash; absent when the symmetric solution never breaks)
  for (e in which) {
    v <- unique(d$threshold[d$est == e])
    v <- v[!is.na(v)]
    if (length(v)) abline(v = v, col = "grey40", lty = c(ML = 2, ULS = 3, WLSMV = 4, MLs = 2)[[e]])
  }
}
if (nrow(dA)) {
  cs <- unique(dA$cond)
  par(mfrow = c(2, length(cs)), mar = c(3, 3, 2, 0.5), mgp = c(1.8, 0.5, 0), oma = c(0, 0, 2, 0))
  for (cd in cs) {
    d <- dA[dA$cond == cd & dA$format == "cont", ]
    plot(NA, xlim = c(0, 0.95), ylim = c(0, 1), xlab = "phi", ylab = "P(tilted)", main = cd, cex.main = 0.8)
    thr_lines(d, c("ML", "ULS"))
    lines_by_N(d[d$est == "ML", ], "phi", "p_tilted")
    lines_by_N(d[d$est == "ULS", ], "phi", "p_tilted", lty = 3, pch = 1)
  }
  for (cd in cs) {
    d <- dA[dA$cond == cd & dA$format == "cont" & dA$est == "ML", ]
    plot(NA, xlim = c(0, 0.95), ylim = c(-0.1, 0.8), xlab = "phi",
         ylab = "loading difference", main = paste(cd, ": mean signed (solid), mean |.| (dashed)"),
         cex.main = 0.7)
    abline(h = 0, col = "grey70"); thr_lines(d, "ML")
    lines_by_N(d, "phi", "mean_d")
    lines_by_N(d, "phi", "mean_abs_d", lty = 2, pch = 1)
    s <- d[d$N == max(d$N), ]; s <- s[order(s$phi), ]
    lines(s$phi, s$pop_abs_d, col = "red", lwd = 1)
  }
  legend_N(sort(unique(dA$N)), "topright")
  mtext("Design A, continuous: tilt probability (solid ML, dotted ULS; vertical lines at phi* and 1/p) and the averaging artifact (red: population |l_A - l_B|)",
        outer = TRUE, cex = 0.8)
  for (fm in c("ord_sym", "ord_skew")) {
    par(mfrow = c(2, length(cs)), mar = c(3, 3, 2, 0.5), mgp = c(1.8, 0.5, 0), oma = c(0, 0, 2, 0))
    for (e in c("WLSMV", "MLs")) for (cd in cs) {
      d <- dA[dA$cond == cd & dA$format == fm & dA$est == e, ]
      plot(NA, xlim = c(0, 0.95), ylim = c(0, 1), xlab = "phi", ylab = "P(tilted)",
           main = paste(cd, e), cex.main = 0.8)
      if (nrow(d)) { thr_lines(d, e); lines_by_N(d, "phi", "p_tilted") }
    }
    legend_N(sort(unique(dA$N)), "topright")
    mtext(sprintf("Design A, %s: tilt probability (vertical lines: WLSMV threshold, dot-dash; phi*, dashed)", fm),
          outer = TRUE, cex = 0.85)
  }
}
dBC <- q2[q2$design %in% c("B", "C"), ]
if (nrow(dBC)) {
  cs <- unique(dBC$cond)
  par(mfrow = c(3, length(cs)), mar = c(3, 3, 2, 0.5), mgp = c(1.8, 0.5, 0), oma = c(0, 0, 2, 0))
  for (arm in list(c("cont", "ML"), c("ord_sym", "WLSMV"), c("ord_skew", "WLSMV"))) for (cd in cs) {
    d <- dBC[dBC$cond == cd & dBC$format == arm[1] & dBC$est == arm[2], ]
    plot(NA, xlim = c(0, 0.95), ylim = c(0, 1), xlab = "phi", ylab = "mean share of block A",
         main = paste(cd, arm[1], arm[2]), cex.main = 0.8)
    abline(h = 0.5, col = "grey70")
    if (nrow(d)) {
      d$mean_shareA <- d$mean_tiltA + 0.5
      lines_by_N(d, "phi", "mean_shareA")
      s <- d[d$N == max(d$N), ]; s <- s[order(s$phi), ]
      lines(s$phi, s$pop_shareA, col = "red")
    }
  }
  legend_N(sort(unique(dBC$N)))
  mtext("Designs B and C: mean share of block A (red: population)", outer = TRUE, cex = 0.9)
}
invisible(dev.off())

if (!is.null(q3)) {
  pdf(outf("fig_spurious.pdf"), width = 11, height = 4.8)
  cs <- unique(q3$cond)
  par(mfrow = c(1, length(cs)), mar = c(3, 3, 2, 0.5), mgp = c(1.8, 0.5, 0), oma = c(0, 0, 3, 0))
  for (cd in cs) {
    d <- q3[q3$cond == cd, ]
    plot(NA, xlim = c(0, 0.95), ylim = c(0, 1), xlab = "phi",
         ylab = "P(weaker-block minimum)", main = cd, cex.main = 0.8)
    if (!is.na(d$phi_fold[1])) abline(v = d$phi_fold[1], lty = 2, col = "grey40")
    lines_by_N(d[d$est == "ML" & d$format == "cont", ], "phi", "p_spurious_weak")
    lines_by_N(d[d$est == "MLs" & d$format == "ord_sym", ], "phi", "p_spurious_weak", lty = 3, pch = 1)
    s <- d[d$est == "ML" & d$format == "cont", ]
    if (nrow(s)) {
      points(s$phi, s$p_default_on_weak_listed_first, pch = 4, col = pal[as.character(s$N)], cex = 0.8)
      points(s$phi, s$p_default_on_weak_listed_last, pch = 2, col = pal[as.character(s$N)], cex = 0.6)
    }
  }
  legend_N(sort(unique(q3$N)), "topright")
  mtext(paste0("Weaker-block minima: continuous ML (solid), ML on symmetric 0-3 scores (dotted); dashed: phi_f\n",
               "lavaan's default start lands there with the weaker block listed first (x) or last (triangle)"),
        outer = TRUE, cex = 0.8, line = 0.2)
  invisible(dev.off())
}
message("wrote tables and figures to ", args$out)
