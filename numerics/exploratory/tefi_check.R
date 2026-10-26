# Source EGAnet 2.5.0 (GitHub HEAD) internal functions without installing the package
env <- new.env()
for (f in list.files("/tmp/EGAnet/R", full.names = TRUE, pattern = "\\.R$")) {
  try(sys.source(f, envir = env), silent = TRUE)
}
stopifnot(exists("tefi_standard", envir = env), exists("matrix_entropy", envir = env))

vn_eigen <- function(rho) { ev <- eigen(rho, symmetric = TRUE, only.values = TRUE)$values; ev <- ev[ev > 1e-15]; -sum(ev * log(ev)) }
tefi_formula <- function(R, s, H) {           # same formula as EGAnet's tefi_standard, pluggable entropy H
  R <- abs(R); K <- length(unique(s)); S <- H(R / ncol(R))
  Sk <- sapply(sort(unique(s)), function(k) { B <- R[s == k, s == k, drop = FALSE]; H(B / ncol(B)) })
  (mean(Sk) - S) + (S - sum(Sk)) * sqrt(K)
}
pop2 <- function(m, lam, phi) {               # two correlated traits, m items each, equal loadings
  L <- rbind(cbind(rep(lam, m), 0), cbind(0, rep(lam, m))); Phi <- matrix(c(1, phi, phi, 1), 2)
  S <- L %*% Phi %*% t(L); diag(S) <- 1; S
}
m <- 5; lam <- .7
parts <- list(
  "1 community"                 = rep(1, 2*m),
  "true 2 (5+5)"                = rep(1:2, each = m),
  "2, one swap each"            = c(1,1,1,1,2, 2,2,2,2,1),
  "2, two swaps each"           = c(1,1,1,2,2, 2,2,2,1,1),
  "3 (5+3+2, over-split)"       = c(rep(1,5), 2,2,2,3,3),
  "4 (3+2+3+2)"                 = c(1,1,1,2,2, 3,3,3,4,4),
  "10 singletons"               = 1:(2*m))
cat("Check: EGAnet tefi_standard == formula with elementwise entropy:\n")
R <- pop2(m, lam, .5); s <- parts[["true 2 (5+5)"]]
ega <- env$tefi_standard(abs(R), s, FALSE)$VN.Entropy.Fit
mine <- tefi_formula(R, s, function(r) -sum(r * log(r)))
cat(sprintf("  EGAnet %.6f | -sum(rho*log(rho)) %.6f | eigen-VN %.6f\n\n", ega, mine, tefi_formula(R, s, vn_eigen)))
for (phi in c(0, .3, .5, .7, .9, .99)) {
  R <- pop2(m, lam, phi)
  cat(sprintf("phi = %.2f\n", phi))
  for (nm in names(parts)) {
    s <- parts[[nm]]
    e <- env$tefi_standard(abs(R), s, FALSE)$VN.Entropy.Fit
    v <- tefi_formula(R, s, vn_eigen)
    cat(sprintf("  %-24s TEFI(EGAnet) %8.3f   TEFI(true VN) %8.3f\n", nm, e, v))
  }
}
