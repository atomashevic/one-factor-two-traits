## Bounded population checks. Run with one numerical thread, for example:
## OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript checks/score_share.R
## This sources the bundled package checkout; no installed undeRfactor is used.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(arg)) stop("Run this file with Rscript.")
script <- normalizePath(sub("^--file=", "", arg[1]), mustWork = TRUE)
paper_dir <- dirname(dirname(script))
package_dir <- normalizePath(file.path(paper_dir, "undeRfactor"), mustWork = TRUE)
src <- list.files(file.path(package_dir, "R"), "\\.R$", full.names = TRUE)
if (!length(src)) stop("No bundled package source files found.")
pkg <- new.env(parent = globalenv())
for (f in src) sys.source(f, envir = pkg)

git_info <- function(path) {
    if (!file.exists(file.path(path, ".git"))) return("snapshot; see provenance/SOURCES.md")
  rev <- system2("git", c("-C", shQuote(path), "rev-parse", "--short", "HEAD"), stdout = TRUE)
  dirty <- system2("git", c("-C", shQuote(path), "status", "--porcelain"), stdout = TRUE)
  paste(rev, if (length(dirty)) "dirty" else "clean")
}
cat("R:", R.version.string, "\n")
cat("Paper:", git_info(paper_dir), "| package:", git_info(package_dir), "\n")
cat("Source directory:", package_dir, "\n")
hashed <- c(script, file.path(package_dir, "R", c("pop_fit.R", "factor_share.R", "two_trait_cor.R", "utils.R")))
cat("MD5:", paste(paste(basename(hashed), tools::md5sum(hashed), sep = "="), collapse = " "), "\n")
cat("Threads:", paste(paste(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"),
                         Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")),
                         sep = "="), collapse = " "), "| numerical seed=1\n")

## Algebraic identities use a scaled 1e-10 tolerance in these well-conditioned cases.
assert_close <- function(actual, expected, label, tol = 1e-10) {
  if (length(actual) != length(expected) || any(!is.finite(c(actual, expected))))
    stop(label, ": nonfinite values or different lengths")
  error <- max(abs(actual - expected))
  if (error > tol * max(1, abs(expected)))
    stop(sprintf("%s: maximum absolute error %.3g exceeds tolerance", label, error))
  invisible(error)
}
ratio <- function(v) {
  den <- sum(v)
  if (any(!is.finite(v)) || !is.finite(den) || den == 0) return(NA_real_)
  v[1] / den
}
signed_share <- function(a, b, phi) {
  if (!all(is.finite(c(a, b, phi))) || abs(phi) >= 1 || a + b == 0) return(NA_real_)
  0.5 + ((1 + phi) / (2 * (1 - phi))) * ((a - b) / (a + b))
}

## Independent implementation of the manuscript's equal-block reference.
## Theorem 1 gives alpha*beta=q and alpha^2+beta^2=L+(2-L)*q^2.
## This checks implementation consistency; it is not an independent theorem proof.
equal_reference <- function(m, lambda, phi) {
  L <- lambda^2
  phi_star <- (1 + (m - 1) * L) / (m * (2 - L))
  if (phi < phi_star) {
    q <- m * L * phi / (1 + (m - 1) * L)
    total <- L + (2 - L) * q^2
    alpha <- sqrt((total + sqrt(total^2 - 4 * q^2)) / 2)
    beta <- q / alpha
  } else alpha <- beta <- sqrt(L * (m - 1 + m * phi) / (2 * m - 1))
  c(alpha = alpha, beta = beta,
    diagnostic = signed_share(alpha / lambda, beta / lambda, phi),
    score = ratio(c(alpha / (1 - alpha^2), beta / (1 - beta^2))))
}

score_share <- function(R, fit, label) {
  d <- attr(R, "two_trait"); p <- nrow(R)
  ia <- d$blocks$A; ib <- d$blocks$B
  Lambda <- matrix(0, p, 2)
  Lambda[ia, 1] <- d$lambda_A; Lambda[ib, 2] <- d$lambda_B
  Phi <- matrix(c(1, d$phi, d$phi, 1), 2)
  assert_close(Lambda %*% Phi %*% t(Lambda) + diag(1 - rowSums(Lambda^2)), R,
               paste(label, "generating covariance"))
  l <- fit$loadings; psi <- fit$uniquenesses
  if (any(psi <= 0)) stop(label, ": these score identities require positive uniquenesses")
  Sigma <- tcrossprod(l) + diag(psi)
  invpsi_l <- l / psi
  information <- sum(l * invpsi_l)
  reg <- solve(Sigma, l)
  bart <- invpsi_l / information
  assert_close(reg, invpsi_l / (1 + information), paste(label, "regression identity"))
  assert_close(reg, bart * information / (1 + information), paste(label, "score proportionality"))
  assert_close(sum(bart * l), 1, paste(label, "Bartlett normalization"))
  coefficients <- cbind(regression = as.vector(crossprod(Lambda, reg)),
                        bartlett = as.vector(crossprod(Lambda, bart)),
                        sum = as.vector(crossprod(Lambda, rep(1, p))))
  assert_close(coefficients[, "sum"], c(sum(d$lambda_A), sum(d$lambda_B)),
               paste(label, "sum coefficients"))
  direct <- c(sum(d$lambda_A * l[ia] / psi[ia]), sum(d$lambda_B * l[ib] / psi[ib]))
  assert_close(coefficients[, "regression"], direct / (1 + information),
               paste(label, "regression trait coefficients"))
  shares <- apply(coefficients, 2, ratio)
  assert_close(shares["regression"], shares["bartlett"], paste(label, "score shares"))
  block_l <- c(mean(l[ia]), mean(l[ib]))
  target <- block_l / c(mean(d$lambda_A), mean(d$lambda_B))
  w <- solve(Phi, target)
  assert_close(Phi %*% w, target, paste(label, "covariance matching"))
  diagnostic <- signed_share(target[1], target[2], d$phi)
  assert_close(ratio(w), diagnostic, paste(label, "diagnostic ratio"))
  assert_close(fit$share$share_A, diagnostic, paste(label, "package positive-loading diagnostic"))
  diagnostic_var <- as.numeric(crossprod(w, Phi %*% w))
  block_target <- numeric(p); block_target[ia] <- block_l[1]; block_target[ib] <- block_l[2]
  item_covariance <- as.vector(Lambda %*% Phi %*% w)
  assert_close(c(mean(item_covariance[ia]), mean(item_covariance[ib])), block_l,
               paste(label, "block-mean covariance matching"))
  if (length(unique(d$lambda_A)) == 1 && length(unique(d$lambda_B)) == 1)
    assert_close(item_covariance, block_target, paste(label, "homogeneous item covariance"))
  spread <- max(abs(l - block_target))
  if (length(ia) == length(ib) && length(unique(c(d$lambda_A, d$lambda_B))) == 1) {
    ref <- equal_reference(length(ia), d$lambda_A[1], d$phi)
    assert_close(block_l, ref[c("alpha", "beta")], paste(label, "equal-block loadings"))
    assert_close(c(diagnostic, shares["regression"]), ref[c("diagnostic", "score")],
                 paste(label, "equal-block shares"))
  }
  c(diagnostic = diagnostic, regression = unname(shares["regression"]),
    bartlett = unname(shares["bartlett"]), sum = unname(shares["sum"]),
    diagnostic_variance = diagnostic_var, alpha = block_l[1], beta = block_l[2],
    within_block_spread = spread, item_covariance_mismatch = max(abs(item_covariance - l)))
}

cases <- list(list(.653, .70, 9, 7, "INC .70"), list(.653, .80, 9, 7, "INC .80"),
              list(.912, .70, 9, 7, "SCOPE .70"), list(.912, .80, 9, 7, "SCOPE .80"),
              list(.865, .70, 9, 7, "ESC .70"), list(.865, .80, 9, 7, "ESC .80"),
              list(.3, .7, 4, 4, "exch .3"), list(.2, .7, 4, 4, "exch .2"),
              list(0, .7, 4, 4, "exch 0"), list(.6, .7, 4, 4, "exch .6"))
cat("\nScientific identity checks (diagnostic is covariance matching; scores are model-implied):\n")
for (cfg in cases) {
  R <- pkg$two_trait_cor(cfg[[1]], cfg[[2]], m_A = cfg[[3]], m_B = cfg[[4]])
  fit <- pkg$pop_fit(R, seed = 1)
  ans <- score_share(R, fit, cfg[[5]])
  cat(sprintf("PASS %-9s diagnostic=%.8f regression=%.8f Bartlett=%.8f sum=%.8f variance=%.8f spread=%.2g [%s]\n",
              cfg[[5]], ans["diagnostic"], ans["regression"], ans["bartlett"], ans["sum"],
              ans["diagnostic_variance"], ans["within_block_spread"], fit$method))
  if (cfg[[5]] == "exch .2") {
    assert_close(ans[c("alpha", "beta")], c(.689208, .230271), "rounded assessment loadings", tol = 5e-7)
    assert_close(ans[c("diagnostic", "regression", "diagnostic_variance")],
                 c(.87435, .84372, .98757), "rounded assessment shares and variance", tol = 5e-6)
    cat(sprintf("PASS rounded assessment: alpha=%.9f beta=%.9f\n", ans["alpha"], ans["beta"]))
  }
  if (cfg[[5]] == "exch 0")
    assert_close(ans[c("diagnostic", "regression", "sum", "diagnostic_variance")],
                 c(1, 1, .5, 1), "endpoint shares and diagnostic variance")
  if (cfg[[5]] == "exch .6") {
    ## In the symmetric branch, Var(eta_w) = 1 - (1-phi)/((2m-1)(1+phi)).
    variance_reference <- 1 - (1 - cfg[[1]]) / ((2 * cfg[[3]] - 1) * (1 + cfg[[1]]))
    assert_close(ans["diagnostic_variance"], variance_reference, "symmetric diagnostic variance")
    if (ans["diagnostic_variance"] >= 1) stop("Symmetric diagnostic variance must be below one at phi=.6.")
    cat(sprintf("PASS symmetric variance reference: %.9f < 1\n", variance_reference))
  }
}
cat("Numerical unequal-block fits are best found; these identities do not certify their global optimality.\n")

## A synthetic heterogeneous example, with no optimization. Signed block means
## define w; individual covariances match only when l_i/lambda_i is constant per block.
R_heterogeneous <- pkg$two_trait_cor(.3, c(.6, .8), c(.5, .7))
l_heterogeneous <- c(.55, .65, .35, .55)
fit_heterogeneous <- list(loadings = l_heterogeneous, uniquenesses = 1 - l_heterogeneous^2,
                          share = pkg$factor_share(.6, .45, .7, .6, .3))
heterogeneous <- score_share(R_heterogeneous, fit_heterogeneous, "synthetic heterogeneous")
assert_close(heterogeneous["item_covariance_mismatch"], 1 / 28, "heterogeneous item mismatch")
cat(sprintf("PASS synthetic heterogeneous: block-mean covariances match (.600, .450); individual item mismatch=%.8f\n",
            heterogeneous["item_covariance_mismatch"]))

cat("\nSigned coefficient and magnitude-share domain checks:\n")
phi <- .3; lambda <- .7; loadings <- c(.6, -.2)
Phi <- matrix(c(1, phi, phi, 1), 2)
w_signed <- solve(Phi, loadings / lambda)
signed <- pkg$factor_share(loadings[1], loadings[2], lambda, lambda, phi, mode = "signed")
assert_close(c(signed$w_A, signed$w_B), w_signed, "signed coefficients")
assert_close(lambda * Phi %*% c(signed$w_A, signed$w_B), loadings, "signed covariance")
assert_close(signed$share_A, ratio(w_signed), "signed ratio")
assert_close(signed$magnitude_share_A, abs(w_signed[1]) / sum(abs(w_signed)), "magnitude ratio")
reflected <- pkg$factor_share(-loadings[1], -loadings[2], lambda, lambda, phi, mode = "signed")
assert_close(reflected$share_A, signed$share_A, "global sign invariance")
abs_result <- pkg$factor_share(loadings[1], loadings[2], lambda, lambda, phi, mode = "absolute")
assert_close(lambda * Phi %*% c(abs_result$w_A, abs_result$w_B), abs(loadings), "legacy absolute target")
cat(sprintf("PASS opposite signs: signed ratio=%.8f; magnitude share=%.8f; legacy absolute-target ratio=%.8f\n",
            signed$share_A, signed$magnitude_share_A, abs_result$share_A))
opposite <- pkg$factor_share(.5, -.5, lambda, lambda, phi, mode = "signed")
zero <- pkg$factor_share(0, 0, lambda, lambda, phi, mode = "signed")
stopifnot(is.na(opposite$share_A), opposite$status == "zero_sum", opposite$magnitude_share_A == .5,
          is.na(zero$share_A), is.na(zero$magnitude_share_A), zero$status == "zero_composite")
for (args in list(list(.5, .5, 0, lambda, phi), list(0, .5, 0, lambda, phi),
                  list(.5, .5, lambda, lambda, 1), list(.5, .5, lambda, lambda, -1))) {
  rejected <- tryCatch({ do.call(pkg$factor_share, args); FALSE }, error = function(e) TRUE)
  stopifnot(rejected)
}
cat("PASS zero-sum and zero-composite statuses; nonidentifiable/incompatible reference and singular correlation rejected.\n")

cat("\nNear-unit correlation: stable arithmetic retains intrinsic contrast amplification.\n")
cat("phi condition(Phi) equal_share perturbed_share amplification\n")
for (phi in c(.9, .99, .9999, 1 - 1e-8, 1 - .Machine$double.eps / 2, -.99)) {
  a <- .6 + 1e-6; b <- .6
  stable <- signed_share(a, b, phi)
  result <- pkg$factor_share(a, b, 1, 1, phi, mode = "signed")
  equal <- pkg$factor_share(.6, .6, 1, 1, phi, mode = "signed")
  assert_close(equal$share_A, .5, "near-unit equal loading ratio")
  assert_close(result$share_A, stable, "near-unit contrast identity", tol = 1e-9)
  assert_close(result$condition_number, (1 + abs(phi)) / (1 - abs(phi)), "spectral condition")
  cat(sprintf("%.16g %.3g %.8f %.8f %.3g\n", phi, result$condition_number,
              equal$share_A, result$share_A, (1 + phi) / (2 * (1 - phi))))
}
cat("PASS: all asserted scientific identities, signed-domain checks and rounded reference values.\n")
