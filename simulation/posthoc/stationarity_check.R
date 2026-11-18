## Stationarity check of multistart solution groups (continuous ML, Designs B
## and C). Regenerates archived replications from their seeds with the
## committed simulation code, repeats the multistart protocol of
## multistart_1f(), and keeps for every start the objective, the loadings and
## the analytic gradient and Hessian of the ML discrepancy.
##
##   generate (machine with the archived runtime, R 4.5.1 and lavaan 0.6-20):
##     Rscript simulation/posthoc/stationarity_check.R --mode generate \
##       --cells 700,706 --cores 2 --out DIR [--reps N]
##   summarise and compare with the archive (machine holding results/raw):
##     Rscript simulation/posthoc/stationarity_check.R --mode summarise --out DIR
##
## Tolerances, fixed before the run and equal to those of checks/order_ctrl.R:
##   stationary:            max |gradient of F| <= 1e-6
##   strict local minimum:  stationary and smallest Hessian eigenvalue > 1e-5
##   saddle:                stationary and smallest Hessian eigenvalue < -1e-5
##   interior:              smallest uniqueness > 1e-8
## F is the ML discrepancy for the matrix lavaan fitted (the sample correlation
## matrix rescaled by (N - 1)/N); coordinates are the p loadings and p
## uniquenesses with factor variance 1. lavaan's objective is fx = F/2.
## Fits that exceed the gradient tolerance are refined by Newton steps from the
## lavaan estimate (at most 25, taken only while the Hessian is positive
## definite). A refined fit counts as a strict local minimum near the lavaan
## estimate when it meets the tolerances above and no standardized loading
## moved by .002 or more, the grouping tolerance of the multistart protocol.
## A positive definite Hessian at a stationary point certifies a local minimum
## of that sample objective only. It says nothing about global optimality.
argv <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), '../..'))
opt <- list(mode = 'generate', cells = '700,706', cores = '1', out = '', reps = '',
            raw = file.path(root, 'results/raw'))
if (length(argv) %% 2L) stop('Options require a value')
if (length(argv)) for (i in seq(1L, length(argv), by = 2L)) {
  key <- sub('^--', '', argv[i])
  if (!key %in% names(opt)) stop('Unknown option: ', argv[i])
  opt[[key]] <- argv[i + 1L]
}
stopifnot(nzchar(opt$out))
grad_tol <- 1e-6; eig_tol <- 1e-5; psi_tol <- 1e-8; move_tol <- 2e-3; fx_tol <- 1e-7

## Analytic discrepancy, gradient and Hessian; same formulas as full_ml() in
## checks/order_ctrl.R, which checks them against numerical derivatives.
full_ml <- function(l, psi, R) {
  p <- length(l); Sigma <- tcrossprod(l) + diag(psi, p)
  U <- chol(Sigma); Ki <- chol2inv(U); B <- Ki %*% R %*% Ki
  A <- Ki - B
  F <- sum(Ki * R) + 2 * sum(log(diag(U))) - as.numeric(determinant(R, logarithm = TRUE)$modulus) - p
  gradient <- c(2 * as.vector(A %*% l), diag(A))
  hessian <- matrix(0, 2 * p, 2 * p)
  for (j in seq_len(2 * p)) {
    ej <- numeric(p); ej[if (j <= p) j else j - p] <- 1
    dS <- if (j <= p) tcrossprod(ej, l) + tcrossprod(l, ej) else tcrossprod(ej)
    dA <- -Ki %*% dS %*% Ki + Ki %*% dS %*% B + B %*% dS %*% Ki
    hessian[, j] <- c(2 * as.vector(dA %*% l + if (j <= p) A %*% ej else numeric(p)), diag(dA))
  }
  list(F = F, gradient = gradient, hessian = (hessian + t(hessian)) / 2)
}
diag_point <- function(l, psi, S) {
  f <- tryCatch(full_ml(l, psi, S), error = function(e) NULL)
  if (is.null(f)) return(c(F = NA, grad = NA, eig = NA))
  c(F = f$F, grad = max(abs(f$gradient)),
    eig = min(eigen(f$hessian, symmetric = TRUE, only.values = TRUE)$values))
}
refine <- function(l, psi, S, maxit = 25L) {
  p <- length(l); z <- c(l, psi); it <- 0L; status <- 'max_iterations'
  repeat {
    f <- tryCatch(full_ml(z[seq_len(p)], z[p + seq_len(p)], S), error = function(e) NULL)
    if (is.null(f)) { status <- 'not_positive_definite_covariance'; break }
    if (max(abs(f$gradient)) <= 1e-10) { status <- 'converged'; break }
    if (it >= maxit) break
    ev <- min(eigen(f$hessian, symmetric = TRUE, only.values = TRUE)$values)
    if (ev <= 0) { status <- 'indefinite_hessian'; break }
    z <- z - solve(f$hessian, f$gradient); it <- it + 1L
  }
  list(l = z[seq_len(p)], psi = z[p + seq_len(p)], iterations = it, status = status)
}
classify <- function(grad, eig, min_psi) {
  if (!is.finite(grad) || !is.finite(eig)) return('not_evaluated')
  if (min_psi <= psi_tol) return('boundary_or_negative_uniqueness')
  if (grad > grad_tol) return('exceeds_gradient_tolerance')
  if (eig > eig_tol) return('strict_minimum')
  if (eig < -eig_tol) return('saddle')
  'flat_or_uncertain'
}

if (opt$mode == 'generate') {
  source(file.path(root, 'simulation/R/sim_functions.R'))
  dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
  cells <- make_cells()
  one_start <- function(fit, start, S, mA, mB, ref, perm = NULL) {
    base <- data.frame(start = start, converged = FALSE, fx = NA_real_, shareA = NA_real_,
                       min_psi = NA_real_, F = NA_real_, grad = NA_real_, eig = NA_real_,
                       ref_status = NA_character_, ref_iter = NA_integer_, ref_F = NA_real_,
                       ref_grad = NA_real_, ref_eig = NA_real_, ref_min_psi = NA_real_,
                       ref_move = NA_real_, stringsAsFactors = FALSE)
    if (is.null(fit)) return(list(row = base, l = NULL))
    s <- one_solution(fit, start, mA, mB, ref, perm)   # archived definitions
    est <- lavInspect(fit, 'est')
    l <- as.numeric(est$lambda[, 1]); psi <- as.numeric(diag(est$theta))
    if (!is.null(perm)) { l <- l[order(perm)]; psi <- psi[order(perm)] }
    if (sum(l) < 0) l <- -l
    d <- diag_point(l, psi, S)
    r <- refine(l, psi, S); dr <- diag_point(r$l, r$psi, S)
    std <- function(l, psi) l / sqrt(l^2 + psi)
    base$converged <- s$converged; base$fx <- s$fx; base$shareA <- s$shareA
    base$min_psi <- min(psi); base$F <- d[['F']]; base$grad <- d[['grad']]; base$eig <- d[['eig']]
    base$ref_status <- r$status; base$ref_iter <- r$iterations; base$ref_F <- dr[['F']]
    base$ref_grad <- dr[['grad']]; base$ref_eig <- dr[['eig']]; base$ref_min_psi <- min(r$psi)
    base$ref_move <- max(abs(std(r$l, r$psi) - std(l, psi)))
    s$std_min_psi <- s$min_psi
    list(row = base, sol = s)
  }
  run_one_cell <- function(id) {
    cell <- cells[cells$cell_id == id, ]
    stopifnot(nrow(cell) == 1L, cell$format == 'cont', isTRUE(cell$ms))
    Sigma <- pop_sigma(cell); p <- nrow(Sigma); mA <- cell$mA; mB <- cell$mB; N <- cell$N
    m1 <- model_1f(p); m2 <- model_2f(mA, mB); rl <- ref_latent(cell); weak <- cell$weak
    reps <- if (nzchar(opt$reps)) as.integer(opt$reps) else cell$reps
    perm <- c(mA + seq_len(mB), seq_len(mA))
    one_rep <- function(r) {
      set.seed(cell$seed + r)
      Y <- gen_latent(N, Sigma); S <- cov(Y); R <- cov2cor(S)
      ## Same calls, in the same order, as run_cell() before the multistart.
      row1 <- extract_fit(try_fit(fit_ml_cov(m1, S, N)), 'ML', '1F', mA, mB, R, rl, FM_ML)
      row2 <- extract_fit(try_fit(fit_ml_cov(m2, S, N)), 'ML', '2F', mA, mB, R, rl, FM_ML)
      Rs <- R * (N - 1) / N                      # matrix lavaan fits (sample.cov.rescale)
      Rr <- R[perm, perm]; dimnames(Rr) <- dimnames(R)
      st <- vector('list', NSTART + 2L)
      st[[1]] <- one_start(try_fit(fit_ml_cov(m1, R, N, test = FALSE)), 0L, Rs, mA, mB, rl)
      st[[2]] <- one_start(try_fit(fit_ml_cov(m1, Rr, N, test = FALSE)), -1L, Rs, mA, mB, rl, perm)
      for (s in seq_len(NSTART)) {
        fs <- try_fit(fit_ml_cov(model_1f_start(runif(p), runif(p, 0.2, 1)), R, N, test = FALSE))
        st[[s + 2L]] <- one_start(fs, s, Rs, mA, mB, rl)
      }
      starts <- do.call(rbind, lapply(st, `[[`, 'row'))
      sols <- lapply(st, function(x) if (is.null(x$sol)) list(converged = FALSE, start = x$row$start) else x$sol)
      d <- distinct_solutions(sols)
      starts$group <- NA_integer_
      for (k in seq_along(d)) starts$group[starts$start %in% d[[k]]$starts] <- k
      on_weak <- function(s) if (weak == 'A') s$shareA > 0.5 else s$shareA < 0.5
      groups <- if (length(d)) do.call(rbind, lapply(seq_along(d), function(k) {
        s <- d[[k]]
        data.frame(group = k, rep_start = s$start, fx = s$fx, count = s$count,
                   n_random = sum(s$starts > 0), has_default = 0L %in% s$starts,
                   has_reversed = -1L %in% s$starts, shareA = s$shareA, std_min_psi = s$min_psi,
                   on_weak = on_weak(s), nonglobal = s$fx > d[[1]]$fx + fx_tol,
                   l = paste(sprintf('%.6f', s$l), collapse = ';'), stringsAsFactors = FALSE)
      })) else NULL
      list(rep = r, default = data.frame(rep = r, fx = row1$fx, cfi = row1$cfi, rmsea = row1$rmsea,
                                         lA = row1$lA, lB = row1$lB),
           starts = cbind(rep = r, starts), groups = if (!is.null(groups)) cbind(rep = r, groups))
    }
    t0 <- proc.time()[['elapsed']]
    res <- parallel::mclapply(seq_len(reps), one_rep, mc.cores = as.integer(opt$cores),
                              mc.preschedule = TRUE)
    bad <- vapply(res, function(x) inherits(x, 'try-error') || is.null(x$starts), TRUE)
    if (any(bad)) stop('Failed replications in cell ', id, ': ', sum(bad))
    ans <- list(cell = cell, default = do.call(rbind, lapply(res, `[[`, 'default')),
                starts = do.call(rbind, lapply(res, `[[`, 'starts')),
                groups = do.call(rbind, lapply(res, `[[`, 'groups')),
                wall = proc.time()[['elapsed']] - t0, R = R.version.string,
                lavaan = as.character(packageVersion('lavaan')),
                code_md5 = unname(tools::md5sum(file.path(root, 'simulation/R/sim_functions.R'))),
                tolerances = c(grad = grad_tol, eig = eig_tol, psi = psi_tol, move = move_tol, fx = fx_tol))
    saveRDS(ans, file.path(opt$out, sprintf('stationarity_cell_%04d.rds', id)))
    cat(sprintf('cell %d done: %d reps, %.0f s\n', id, reps, ans$wall))
  }
  for (id in as.integer(strsplit(opt$cells, ',')[[1]])) {
    if (file.exists(file.path(opt$out, sprintf('stationarity_cell_%04d.rds', id)))) next
    run_one_cell(id)
  }
} else if (opt$mode == 'summarise') {
  fs <- sort(list.files(opt$out, '^stationarity_cell_[0-9]{4}\\.rds$', full.names = TRUE))
  stopifnot(length(fs) > 0)
  all_rows <- list()
  for (f in fs) {
    x <- readRDS(f); cell <- x$cell; n <- nrow(x$default)
    cat(sprintf('\n== cell %d: %s, phi = %s, N = %s; %s, lavaan %s; %d replications ==\n',
                cell$cell_id, cell$cond, cell$phi, cell$N, x$R, x$lavaan, n))
    ## Gate: agreement with the archive.
    arc_file <- file.path(opt$raw, sprintf('cell_%04d.rds', cell$cell_id)); gate <- NA
    g <- x$groups; s <- x$starts
    s$class <- mapply(classify, s$grad, s$eig, s$min_psi)
    s$ref_class <- mapply(classify, s$ref_grad, s$ref_eig, s$ref_min_psi)
    s$ref_near <- s$ref_class == 'strict_minimum' & s$ref_move < move_tol
    key <- paste(s$rep, s$start)
    pick <- match(paste(g$rep, g$rep_start), key)
    g <- cbind(g, s[pick, c('F', 'grad', 'eig', 'min_psi', 'class', 'ref_status', 'ref_grad',
                            'ref_eig', 'ref_move', 'ref_class', 'ref_near')])
    g$admissible <- g$std_min_psi > 0
    weak <- g[g$nonglobal & g$admissible & g$on_weak, ]
    weak <- weak[!duplicated(weak$rep), ]            # first such group, as in multistart_1f()
    best <- g[g$group == 1L, ]
    def <- g[g$has_default, ]
    if (file.exists(arc_file)) {
      a <- readRDS(arc_file)
      af <- a$fits[a$fits$model == '1F' & a$fits$est == 'ML', ]
      af <- af[match(x$default$rep, af$rep), ]; am <- a$ms[match(best$rep, a$ms$rep), ]
      spur <- best$rep %in% weak$rep
      defweak <- def$nonglobal & def$on_weak & def$admissible
      gate <- c(default_fx = max(abs(af$fx - x$default$fx)), default_cfi = max(abs(af$cfi - x$default$cfi)),
                best_fx = max(abs(am$best_fx - best$fx)),
                n_distinct_mismatch = sum(am$n_distinct != as.integer(table(factor(g$rep, best$rep)))),
                spurious_weak_mismatch = sum(am$spurious_weak != spur),
                def_on_weak_mismatch = sum(am$def_on_weak[match(def$rep, best$rep)] != (def$nonglobal & def$on_weak)))
      cat('Gate against archive (max abs differences and mismatch counts):\n'); print(gate)
      cat(sprintf('Archive: def_on_weak %d, spurious_weak %d of %d\n', sum(a$ms$def_on_weak), sum(a$ms$spurious_weak), nrow(a$ms)))
    } else cat('Archive cell not found; gate not evaluated.\n')
    tab <- function(z, label) {
      cat(sprintf('%s: %d\n', label, nrow(z)))
      if (!nrow(z)) return(invisible())
      print(table(lavaan_estimate = z$class))
      print(table(after_refinement = ifelse(z$class == 'strict_minimum', 'strict_minimum (no refinement needed)',
                  ifelse(z$ref_near, 'refined strict minimum, moved < .002', paste0('refined: ', z$ref_class,
                         ifelse(z$ref_class == 'strict_minimum', ', moved >= .002', ''))))))
      cat(sprintf('  max |gradient| at lavaan estimate: median %.2g, max %.2g; smallest Hessian eigenvalue: min %.4g, median %.4g\n',
                  median(z$grad), max(z$grad), min(z$eig), median(z$eig)))
      cat(sprintf('  after refinement: max |gradient| %.2g; smallest eigenvalue min %.4g; largest loading move %.2g\n',
                  max(z$ref_grad, na.rm = TRUE), min(z$ref_eig, na.rm = TRUE), max(z$ref_move, na.rm = TRUE)))
    }
    defweak <- def[def$nonglobal & def$on_weak & def$admissible, ]
    tab(defweak, 'Default-order fits in a higher-objective weaker-block group')
    tab(weak, 'Replications with a higher-objective weaker-block group (group representative)')
    tab(best, 'Best group found (representative)')
    tab(g, 'All retained groups')
    cat(sprintf('All converged starts: %d; classes at the lavaan estimate:\n', sum(s$converged)))
    print(table(s$class[s$converged]))
    ok <- function(z) sum(z$class == 'strict_minimum' | z$ref_near)
    all_rows[[length(all_rows) + 1L]] <- data.frame(
      cell_id = cell$cell_id, cond = cell$cond, phi = cell$phi, N = cell$N, reps = n,
      R = x$R, lavaan = x$lavaan,
      gate_default_fx = gate[['default_fx']], gate_default_cfi = gate[['default_cfi']],
      gate_best_fx = gate[['best_fx']],
      default_weak = nrow(defweak), default_weak_strict_raw = sum(defweak$class == 'strict_minimum'),
      default_weak_strict = ok(defweak), default_weak_saddle = sum(defweak$ref_class == 'saddle' | defweak$class == 'saddle'),
      weak_groups = nrow(weak), weak_strict_raw = sum(weak$class == 'strict_minimum'), weak_strict = ok(weak),
      best_groups = nrow(best), best_strict_raw = sum(best$class == 'strict_minimum'), best_strict = ok(best),
      groups = nrow(g), groups_strict_raw = sum(g$class == 'strict_minimum'), groups_strict = ok(g),
      min_eig_weak = if (nrow(weak)) min(weak$eig) else NA, max_grad_weak = if (nrow(weak)) max(weak$grad) else NA,
      stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, all_rows)
  write.csv(res, file.path(opt$out, 'stationarity_summary.csv'), row.names = FALSE)
  cat('\nSummary written to', file.path(opt$out, 'stationarity_summary.csv'), '\n')
  print(t(res))
} else stop('Unknown mode')
