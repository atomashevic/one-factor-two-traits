suppressMessages(library(lavaan))
pop <- function(m, lam, phi){L<-rbind(cbind(rep(lam,m),0),cbind(0,rep(lam,m))); S<-L%*%matrix(c(1,phi,phi,1),2)%*%t(L); diag(S)<-1; v<-paste0("x",1:(2*m)); dimnames(S)<-list(v,v); S}
srmr_cf <- function(m,lam,phi) lam^2*(1-phi)*sqrt(m*(m-1)/(4*m^2-1))
cat("SRMR closed form vs lavaan (ML, phi >= phi*)\n")
for (m in c(2,4,8)) for (lam in c(.5,.7)) for (phi in c(.6,.8,.95)) {
  S <- pop(m,lam,phi); v<-colnames(S)
  f <- cfa(paste("f =~", paste(v,collapse="+")), sample.cov=S, sample.nobs=1e6, std.lv=TRUE, sample.cov.rescale=FALSE)
  fm <- fitMeasures(f, c("srmr","tli","cfi","rmsea"))
  cat(sprintf(" m=%d lam=%.1f phi=%.2f  SRMR lavaan %.5f closed %.5f | TLI %.4f\n", m, lam, phi, fm["srmr"], srmr_cf(m,lam,phi), fm["tli"]))
}
cat("\nULS (MINRES) one-factor fits: block loadings\n")
for (m in c(2,4)) for (phi in c(0,.1,.2,.3,.4)) {
  S <- pop(m,.7,phi); v<-colnames(S)
  best <- NULL
  for (st in list(c(.7,.05), c(.5,.5), c(.05,.7), c(.6,.3))) {
    mod <- paste0("f =~ ", paste0("start(",rep(st,each=m),")*",v,collapse="+"))
    f <- try(cfa(mod, sample.cov=S, sample.nobs=1e5, std.lv=TRUE, estimator="ULS", sample.cov.rescale=FALSE), silent=TRUE)
    if (inherits(f,"try-error")) next
    fv <- f@optim$fx; if (is.null(best) || fv < best$fx) best <- list(fx=fv, l=coef(f)[1:(2*m)])
  }
  cat(sprintf(" m=%d phi=%.1f  loadings block A %.3f, block B %.3f\n", m, phi, mean(abs(best$l[1:m])), mean(abs(best$l[(m+1):(2*m)]))))
}
