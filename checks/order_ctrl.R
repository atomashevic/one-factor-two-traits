## Bounded controlled item-order check. Run from the paper repository:
## OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript checks/order_ctrl.R
## The bundled package source is used directly. No installed undeRfactor is used.
if (sys.nframe() > 0L) {
  script <- normalizePath(sys.frame(1)$ofile, mustWork = TRUE)
} else {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(arg)) stop("Run this file with Rscript.")
  script <- normalizePath(sub("^--file=", "", arg[1]), mustWork = TRUE)
}
paper_dir <- dirname(dirname(script))
package_dir <- normalizePath(file.path(paper_dir, "undeRfactor"), mustWork = TRUE)
thread_vars <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
if (any(Sys.getenv(thread_vars) != "1")) stop("Use the single-thread command printed at the top of this script.")
suppressPackageStartupMessages(library(lavaan))
if (!requireNamespace("numDeriv", quietly = TRUE)) stop("numDeriv is required for independent derivative checks.")
pkg <- new.env(parent = globalenv())
src <- list.files(file.path(package_dir, "R"), "\\.R$", full.names = TRUE)
for (f in src) sys.source(f, envir = pkg)

## Coordinates are (l_1,...,l_p,psi_1,...,psi_p), with latent variance fixed
## to one and observed scales unchanged. Marker fits are mapped to this same
## identified parameterization. Only one common loading-sign flip is applied.
full_ml <- function(l, psi, R) {
  p <- length(l); Sigma <- tcrossprod(l) + diag(psi, p)
  U <- chol(Sigma); Ki <- chol2inv(U); B <- Ki %*% R %*% Ki
  A <- Ki - B
  F <- sum(Ki * R) + 2 * sum(log(diag(U))) - as.numeric(determinant(R, logarithm = TRUE)$modulus) - p
  gradient <- c(2 * as.vector(A %*% l), diag(A))
  hessian <- matrix(0, 2*p, 2*p)
  for (j in seq_len(2*p)) {
    ej <- numeric(p); ej[if (j <= p) j else j-p] <- 1
    dS <- if (j <= p) tcrossprod(ej, l) + tcrossprod(l, ej) else tcrossprod(ej)
    dA <- -Ki %*% dS %*% Ki + Ki %*% dS %*% B + B %*% dS %*% Ki
    hessian[, j] <- c(2 * as.vector(dA %*% l + if (j <= p) A %*% ej else numeric(p)), diag(dA))
  }
  if (max(abs(hessian-t(hessian))) > 1e-9) stop("Analytic Hessian lost symmetry.")
  list(F = F, gradient = gradient, hessian = (hessian+t(hessian))/2, Sigma = Sigma)
}

grad_tol <- 1e-6; eig_tol <- 1e-5; psi_tol <- 1e-8
cov_tol <- 1e-12; objective_tol <- 1e-8; derivative_tol <- 1e-6
pack <- function(x) paste(sprintf("%.12g", x), collapse = ";")
signed_share <- function(a, b, lambda_A, lambda_B, phi) {
  a <- a/lambda_A; b <- b/lambda_B
  if (!is.finite(a+b) || abs(a+b) < 1e-12) return(NA_real_)
  .5 + (1+phi)/(2*(1-phi)) * (a-b)/(a+b)
}

diagnose <- function(l, psi, R, cfg) {
  orientation <- if (sum(l) < 0) -1 else 1
  l <- orientation*l
  f <- full_ml(l, psi, R)
  max_gradient <- max(abs(f$gradient))
  minimum_eigenvalue <- min(eigen(f$hessian, symmetric = TRUE, only.values = TRUE)$values)
  state <- if (min(psi) < -psi_tol) "inadmissible_negative_uniqueness" else
    if (min(psi) <= psi_tol) "uniqueness_boundary_unclassified" else
    if (max_gradient > grad_tol) "nonstationary" else
    if (minimum_eigenvalue < -eig_tol) "saddle_evidence" else
    if (minimum_eigenvalue > eig_tol) "strict_minimum_evidence" else "flat_or_uncertain_stationary"
  A <- seq_len(4); BB <- 4+seq_len(4)
  a <- mean(l[A]); b <- mean(l[BB])
  list(F = f$F, max_full_gradient = max_gradient, min_full_hessian_eigenvalue = minimum_eigenvalue,
       min_uniqueness = min(psi), max_uniqueness = max(psi),
       min_covariance_eigenvalue = min(eigen(f$Sigma, symmetric = TRUE, only.values = TRUE)$values),
       common_sign_flip = orientation, loading_A = a, loading_B = b,
       within_block_spread = max(abs(l[A]-a), abs(l[BB]-b)),
       observed_variance_error = max(abs(diag(f$Sigma)-1)),
       signed_share_A = signed_share(a, b, cfg$lambda_A, cfg$lambda_B, cfg$phi),
       stationary = max_gradient <= grad_tol, numerical_classification = state,
       loadings = pack(l), uniquenesses = pack(psi))
}

make_row <- function(cfg, source, ident, start, order, converged, l, psi, R,
                     latent_variance = 1, lavaan_F = NA_real_, start_error = NA_real_,
                     start_l = NULL, start_psi = NULL, warnings = "") {
  ans <- diagnose(l, psi, R, cfg)
  data.frame(case = cfg$case, phi = cfg$phi, lambda_A = cfg$lambda_A, lambda_B = cfg$lambda_B,
             source = source, identification = ident, start = start, item_order = order,
             lavaan_or_reference_converged = converged, fitted_latent_variance = latent_variance,
             lavaan_F = lavaan_F, lavaan_objective_error = abs(ans$F-lavaan_F),
             starting_covariance_target_error = start_error,
             starting_loadings_unit_latent_variance = pack(start_l), starting_uniquenesses = pack(start_psi),
             as.data.frame(ans, stringsAsFactors = FALSE), warnings = warnings,
             reference_F = NA_real_, objective_gap_to_reference = NA_real_, stringsAsFactors = FALSE)
}

check_derivatives <- function(l, psi, R, label) {
  p <- length(l); z <- c(l, psi)
  fn <- function(x) full_ml(x[seq_len(p)], x[p+seq_len(p)], R)$F
  gr <- function(x) full_ml(x[seq_len(p)], x[p+seq_len(p)], R)$gradient
  exact <- full_ml(l, psi, R)
  ge <- max(abs(numDeriv::grad(fn, z)-exact$gradient))
  he <- max(abs(numDeriv::jacobian(gr, z)-exact$hessian))
  if (ge > derivative_tol || he > derivative_tol) stop("Derivative check failed: ", label)
  cat(sprintf("PASS derivative check %s: max gradient error %.3g; max Hessian error %.3g\n", label, ge, he))
}

## Independent evaluation of the equal-block analytic branches, including
## the symmetric stationary point below threshold, which is a saddle.
equal_reference <- function(lambda, phi, m = 4) {
  L <- lambda^2; n <- m-1; b <- 1+n*L; q <- m*L*phi/b
  phi_star <- b/(m*(2-L))
  stopifnot(phi < phi_star)
  total <- L+(2-L)*q^2
  alpha <- sqrt((total+sqrt(total^2-4*q^2))/2); beta <- q/alpha
  F <- -n*log(1-L)-log(b)+n*log(1-q^2)
  ls <- sqrt(L*(n+m*phi)/(2*m-1))
  list(list(label = "analytic_A", l = rep(c(alpha, beta), each = m), F = F),
       list(label = "analytic_B", l = rep(c(beta, alpha), each = m), F = F),
       list(label = "analytic_symmetric_saddle", l = rep(ls, 2*m), F = NA_real_))
}

fit_control <- function(R, ident, start, order) {
  canonical <- rownames(R); nm <- canonical[order]; S <- R[nm, nm]
  if (start == "default") model <- paste("f =~", paste(nm, collapse = " + ")) else {
    loadings <- if (ident == "std.lv") sprintf("start(.6)*%s", nm) else
      c(nm[1], sprintf("start(1)*%s", nm[-1]))
    model <- paste0("f =~ ", paste(loadings, collapse = " + "), "\n",
                    if (ident == "marker") "f ~~ start(.36)*f\n" else "",
                    paste(sprintf("%s ~~ start(.64)*%s", nm, nm), collapse = "\n"))
  }
  warnings <- character()
  fit <- withCallingHandlers(cfa(model, sample.cov = S, sample.nobs = 1000,
                                 sample.cov.rescale = FALSE, estimator = "ML",
                                 std.lv = ident == "std.lv", se = "none"),
                             warning = function(w) {
                               warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
                             })
  starts <- lavInspect(fit, "start"); est <- lavInspect(fit, "est")
  start_Sigma <- (starts$lambda %*% starts$psi %*% t(starts$lambda) + starts$theta)[canonical, canonical]
  target <- matrix(.36, length(canonical), length(canonical)); diag(target) <- 1
  start_error <- max(abs(start_Sigma-target))
  if (start == "common" && start_error > cov_tol) stop("Common implied starting covariance does not match target.")
  if (est$psi[1,1] <= 0 || starts$psi[1,1] <= 0)
    stop("Nonpositive latent variance requires a separate admissibility analysis.")
  converged <- isTRUE(lavInspect(fit, "converged"))
  lavaan_F <- if (converged) 2*unname(fitMeasures(fit, "fmin")) else NA_real_
  list(converged = converged, l = (est$lambda[,1]*sqrt(est$psi[1,1]))[canonical],
       psi = diag(est$theta)[canonical], factor_variance = est$psi[1,1], lavaan_F = lavaan_F,
       start_error = start_error, start_Sigma = start_Sigma,
       start_l = (starts$lambda[,1]*sqrt(starts$psi[1,1]))[canonical],
       start_psi = diag(starts$theta)[canonical], warnings = unique(warnings))
}

main <- function() {
  git_info <- function(path) {
    if (!file.exists(file.path(path, ".git"))) return("snapshot; see provenance/SOURCES.md")
    rev <- system2("git", c("-C", shQuote(path), "rev-parse", "HEAD"), stdout = TRUE)
    dirty <- system2("git", c("-C", shQuote(path), "status", "--porcelain"), stdout = TRUE)
    paste(rev, if (length(dirty)) "dirty" else "clean")
  }
  seed <- 20261009L; random_starts <- 8L
  cat("Command: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript checks/order_ctrl.R\n")
  cat("UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
  cat("R:", R.version.string, "| lavaan:", as.character(packageVersion("lavaan")),
      "| numDeriv:", as.character(packageVersion("numDeriv")), "\n")
  cat("Paper:", git_info(paper_dir), "| package:", git_info(package_dir), "\n")
  cat("Source directory:", package_dir, "\n")
  cat("Sibling package source version:", read.dcf(file.path(package_dir, "DESCRIPTION"))[1,"Version"], "\n")
  hashed <- c(script, src)
  cat("MD5:", paste(paste(basename(hashed), tools::md5sum(hashed), sep = "="), collapse = " "), "\n")
  cat("Threads:", paste(paste(thread_vars, Sys.getenv(thread_vars), sep = "="), collapse = " "),
      "| seed=", seed, " | random starts per unequal case=", random_starts, "\n", sep = "")
  cat("Population sample.cov.rescale=FALSE; sample.nobs=1000; no generated data.\n")
  cat("Coordinates: unit latent variance, raw observed scales; one common sign flip.\n")
  cat("Thresholds: gradient=", grad_tol, "; Hessian=", eig_tol, "; uniqueness boundary=", psi_tol,
      "; starting covariance=", cov_tol, "; objective=", objective_tol, "; derivative=", derivative_tol, "\n", sep = "")
  cat("Negative uniqueness is inadmissible; boundary fits are unclassified. No constrained KKT conclusion is drawn.\n")
  cat("Unequal-case reference = best strict minimum found numerically, not a global optimality certificate.\n")
  cases <- data.frame(case = 1:6, phi = c(.1, .01, .2, .05, .3, .3),
                      lambda_A = c(.68, .60, .68, .60, .7, .6), lambda_B = c(.72, .80, .72, .80, .7, .8))
  output <- list(); common_errors <- numeric()
  for (ci in seq_len(nrow(cases))) {
    cfg <- cases[ci, ]; tagged <- pkg$two_trait_cor(cfg$phi, cfg$lambda_A, cfg$lambda_B, m_A = 4)
    R <- tagged; attr(R, "two_trait") <- NULL
    cat(sprintf("\nCase %d: phi=%.2f lambda_A=%.2f lambda_B=%.2f\n", ci, cfg$phi, cfg$lambda_A, cfg$lambda_B))
    check_derivatives(rep(.6,8), rep(.64,8), R, paste("case", ci, "common start"))
    refs <- list()
    if (cfg$lambda_A == cfg$lambda_B) {
      sols <- equal_reference(cfg$lambda_A, cfg$phi)
      for (sol in sols) {
        row <- make_row(cfg, "analytic_reference", "unit_latent_variance", sol$label, "canonical", NA,
                        sol$l, 1-sol$l^2, R)
        if (is.finite(sol$F) && abs(row$F-sol$F) > objective_tol) stop("Analytic objective mismatch.")
        refs[[length(refs)+1]] <- row
      }
      stopifnot(refs[[1]]$numerical_classification == "strict_minimum_evidence",
                refs[[2]]$numerical_classification == "strict_minimum_evidence",
                refs[[3]]$numerical_classification == "saddle_evidence")
      check_derivatives(sols[[1]]$l, 1-sols[[1]]$l^2, R, "analytic asymmetric minimum")
      check_derivatives(sols[[3]]$l, 1-sols[[3]]$l^2, R, "analytic symmetric saddle")
    } else {
      cat("Numerical reference seed:", seed+ci, "| five deterministic starts plus", random_starts,
          "random starts; package saddle-escape perturbation retained.\n")
      sols <- pkg$with_seed(seed+ci, pkg$ml_numeric(R, attr(tagged,"two_trait"), nstart = random_starts))
      labels <- c("smc", "A_only", "B_only", "tilt_A", "tilt_B", sprintf("random_%02d", seq_len(random_starts)))
      stopifnot(length(sols) == length(labels))
      for (j in seq_along(sols)) {
        sol <- sols[[j]]
        row <- make_row(cfg, "numerical_reference", "unit_latent_variance", labels[j], "canonical",
                        sol$converged, sol$l, sol$psi, R)
        if (abs(row$F-sol$F) > objective_tol) stop("Concentrated versus full objective mismatch.")
        refs[[length(refs)+1]] <- row
      }
    }
    refs <- do.call(rbind, refs)
    good <- refs$numerical_classification == "strict_minimum_evidence"
    if (!any(good)) stop("No numerically verified strict reference minimum.")
    reference_F <- min(refs$F[good])
    distinct <- list()
    for (j in which(good)) {
      l <- as.numeric(strsplit(refs$loadings[j], ";", fixed = TRUE)[[1]])
      psi <- as.numeric(strsplit(refs$uniquenesses[j], ";", fixed = TRUE)[[1]])
      S <- tcrossprod(l)+diag(psi)
      if (!length(distinct) || all(vapply(distinct, function(x) max(abs(x-S)) > 1e-4, logical(1))))
        distinct[[length(distinct)+1]] <- S
    }
    best <- which.min(ifelse(good, refs$F, Inf))
    cat(sprintf("Reference: F=%.10f signed_share_A=%.6f; %d distinct verified strict minima found; %s\n",
                reference_F, refs$signed_share_A[best], length(distinct),
                if (cfg$lambda_A == cfg$lambda_B) "analytic global value" else "numerical best found"))
    refs$reference_F <- reference_F; refs$objective_gap_to_reference <- refs$F-reference_F
    output[[length(output)+1]] <- refs
    starts_common <- list()
    for (ident in c("std.lv", "marker")) for (start in c("default", "common")) {
      for (ord in c("A_first", "B_first")) {
        order <- if (ord == "A_first") 1:8 else c(5:8,1:4)
        ans <- fit_control(R, ident, start, order)
        row <- make_row(cfg, "lavaan_control", ident, start, ord, ans$converged,
                        ans$l, ans$psi, R, ans$factor_variance, ans$lavaan_F, ans$start_error,
                        ans$start_l, ans$start_psi, paste(ans$warnings, collapse = " | "))
        if (is.finite(row$lavaan_objective_error) && row$lavaan_objective_error > objective_tol)
          stop("lavaan versus direct objective mismatch.")
        row$reference_F <- reference_F; row$objective_gap_to_reference <- row$F-reference_F
        output[[length(output)+1]] <- row
        if (start == "common") { starts_common[[length(starts_common)+1]] <- ans$start_Sigma; common_errors <- c(common_errors,ans$start_error) }
        cat(sprintf("  %-6s %-7s %-7s converged=%-5s F=%.10f gap=% .3g grad=%.2g Hmin=% .3g psi=[%.4f,%.4f] l=[%.5f,%.5f] shareA=%.6f %s\n",
                    ident, start, ord, ans$converged, row$F, row$objective_gap_to_reference,
                    row$max_full_gradient, row$min_full_hessian_eigenvalue, row$min_uniqueness,
                    row$max_uniqueness, row$loading_A, row$loading_B, row$signed_share_A, row$numerical_classification))
        if (length(ans$warnings)) cat("    Warnings:", paste(ans$warnings, collapse = " | "), "\n")
      }
    }
    error <- max(vapply(starts_common, function(x) max(abs(x-starts_common[[1]])), numeric(1)))
    if (error > cov_tol) stop("Common starts differ across identification or item order.")
    cat(sprintf("PASS common starting covariance across four identification/order combinations: max error %.3g\n", error))
  }
  rows <- do.call(rbind, output)
  csv <- file.path(paper_dir,"build","checks","order_ctrl_rerun.csv")
  write.csv(rows, csv, row.names = FALSE, na = "NA")
  cat("\nCSV:", csv, "\n")
  cat("Controlled fits:", sum(rows$source == "lavaan_control"), "| reference candidates:", sum(rows$source != "lavaan_control"), "\n")
  print(with(rows[rows$source == "lavaan_control", ], table(lavaan_or_reference_converged, numerical_classification)))
  cat("Maximum common starting covariance target error:", max(common_errors), "\n")
  cat("PASS: bounded controlled comparisons and all asserted independent checks.\n")
  cat("Convergence, stationarity, curvature, and admissibility are separate diagnostics. Multiple starts do not certify numerical global optimality.\n")
  cat("\nSession info:\n"); print(sessionInfo())
}
if (sys.nframe() == 0L) {
  log <- file.path(paper_dir,"build","checks","order_ctrl_rerun.log")
  sink(log, split = TRUE)
  tryCatch(main(), error = function(e) { cat("\nFAIL:", conditionMessage(e), "\n"); stop(e) }, finally = sink())
}
