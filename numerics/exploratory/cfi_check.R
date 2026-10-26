suppressMessages(library(lavaan))
pop <- function(L, Phi) { S <- L %*% Phi %*% t(L); diag(S) <- 1; S }
fit1 <- function(S) {
  p <- ncol(S); v <- paste0("x", 1:p); dimnames(S) <- list(v, v)
  mod <- paste("f =~", paste(v, collapse = " + "))
  cfa(mod, sample.cov = S, sample.nobs = 1e6, std.lv = TRUE, sample.cov.rescale = FALSE)
}
Fml <- function(S, Sig) { A <- solve(Sig) %*% S; sum(diag(A)) - log(det(A)) - ncol(S) }
cat("=== Exchangeable two-trait populations: lavaan 1-factor ML vs closed forms ===\n")
cat(sprintf("%3s %4s %4s | %9s %9s | %8s %8s | %7s %7s | %9s\n","m","lam","phi","F_lavaan","F_closed","tr-p","maxdiag","CFI_lav","CFI_cf","h-gap*2-F"))
for (m in c(2, 4, 8)) for (lam in c(.5, .7)) for (phi in c(.5, .8, .95)) {
  p <- 2*m; L <- rbind(cbind(rep(lam,m),0), cbind(0,rep(lam,m)))
  S <- pop(L, matrix(c(1,phi,phi,1),2)); f <- fit1(S); Sig <- unclass(fitted(f)$cov)
  F <- Fml(S, Sig)
  kap <- m*lam^2*(1-phi)/(1-lam^2); Fc <- (p-1)*log(1+kap/(p-1)) - log(1+kap)
  trc <- sum(diag(solve(Sig) %*% S)) - p
  hgap <- log(det(Sig)) - log(det(S))          # = 2[h(Sigma_hat) - h(S)]
  F0 <- -log(det(S)); CFIc <- 1 - Fc/F0
  cat(sprintf("%3d %4.2f %4.2f | %9.6f %9.6f | %8.1e %8.1e | %7.4f %7.4f | %9.1e\n", m, lam, phi, F, Fc, trc,
      max(abs(diag(Sig)-1)), fitMeasures(f,"cfi"), CFIc, hgap - F))
}
cat("\n=== Non-exchangeable: random loadings, unequal blocks; bound KL_1F <= (kappa2 - log(1+kappa2))/2 ===\n")
set.seed(11); out <- NULL
for (r in 1:60) {
  m1 <- sample(3:8,1); m2 <- sample(3:8,1); phi <- runif(1,.2,.95)
  l1 <- runif(m1,.35,.85); l2 <- runif(m2,.35,.85)
  L <- rbind(cbind(l1,0), cbind(0,l2)); Phi <- matrix(c(1,phi,phi,1),2)
  S <- pop(L, Phi); f <- fit1(S); KL <- Fml(S, unclass(fitted(f)$cov))/2
  s1 <- sum(l1^2/(1-l1^2)); s2 <- sum(l2^2/(1-l2^2))
  k2 <- ((s1+s2) - sqrt((s1-s2)^2 + 4*s1*s2*phi^2))/2
  ub <- (k2 - log(1+k2))/2
  pp <- m1+m2; lb_exch <- ((pp-1)*log(1+k2/(pp-1)) - log(1+k2))/2
  out <- rbind(out, c(m1,m2,phi,k2,KL,ub,KL/ub, lb_exch))
}
colnames(out) <- c("m1","m2","phi","kappa2","KL_1F","upper","KL/upper","exch_formula")
cat(sprintf("violations of upper bound: %d of %d; KL/upper range %.3f to %.3f\n", sum(out[,5] > out[,6] + 1e-9), nrow(out), min(out[,7]), max(out[,7])))
cat(sprintf("KL vs exchangeable formula with same kappa2 and p: ratio range %.3f to %.3f\n", min(out[,5]/out[,8]), max(out[,5]/out[,8])))
print(round(head(out[order(out[,4]),], 8), 4))
