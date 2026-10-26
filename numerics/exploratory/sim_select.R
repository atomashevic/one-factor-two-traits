suppressMessages(library(lavaan)); suppressMessages(library(MASS))
env <- new.env(); for (f in list.files("/tmp/EGAnet/R", full.names=TRUE, pattern="\\.R$")) try(sys.source(f, envir=env), silent=TRUE)
vn <- function(r){ev<-eigen(r,TRUE,TRUE)$values; ev<-ev[ev>1e-15]; -sum(ev*log(ev))}
tefi_vn <- function(R,s){R<-abs(R);K<-length(unique(s));S<-vn(R/ncol(R));Sk<-sapply(unique(s),function(k){B<-R[s==k,s==k,drop=FALSE];vn(B/ncol(B))});(mean(Sk)-S)+(S-sum(Sk))*sqrt(K)}
icm_syntax <- function(s){ v<-paste0("x",seq_along(s)); out<-c()
  for(k in sort(unique(s))){ it<-v[s==k]; out<-c(out, paste0("f",k," =~ ",paste(it,collapse="+")))
    if(length(it)==1) out<-c(out, paste0(it," ~~ 0*",it)) }
  paste(out, collapse="\n") }
cands <- list(K1=rep(1,8), true2=rep(1:2,each=4), swap1=c(1,1,1,2,2,2,2,1), split3=c(1,1,1,1,2,2,3,3),
              peel=c(1,2,2,2,2,2,2,2), split4=c(1,1,2,2,3,3,4,4))
syn <- lapply(cands, icm_syntax)
stir <- c(1,127,966,1701); set.seed(2026); lam <- .55; res <- list()
for (phi in c(.3,.5,.7,.9,1)) for (N in c(250,1000)) {
  L <- rbind(cbind(rep(lam,4),0), cbind(0,rep(lam,4))); S <- L%*%matrix(c(1,phi,phi,1),2)%*%t(L); diag(S) <- 1
  truth <- if (phi==1) "K1" else "true2"; pick <- matrix(NA,100,4, dimnames=list(NULL,c("TEFI_EGAnet","TEFI_VN","BIC_ICM","DL")))
  for (r in 1:100) {
    X <- mvrnorm(N, rep(0,8), S); R <- cor(X); v <- paste0("x",1:8); dimnames(R) <- list(v,v)
    te <- sapply(cands, function(s) env$tefi_standard(abs(R), s, FALSE)$VN.Entropy.Fit)
    tv <- sapply(cands, function(s) tefi_vn(R, s))
    bic <- sapply(names(cands), function(nm) { f <- suppressWarnings(try(cfa(syn[[nm]], sample.cov=R, sample.nobs=N, std.lv=TRUE), silent=TRUE))
                   if (inherits(f,"try-error")) Inf else unname(fitMeasures(f,"bic")) })
    names(bic) <- names(cands); dl <- bic + 2*log(stir[sapply(cands,function(z) length(unique(z)))]); pick[r,] <- c(names(which.min(te)), names(which.min(tv)), names(which.min(bic)), names(which.min(dl)))
  }
  acc <- colMeans(pick==truth); modal <- apply(pick,2,function(z) names(sort(table(z),decreasing=TRUE))[1])
  cat(sprintf("phi=%.1f N=%4d truth=%-5s | TEFI(EGAnet) %.2f [%s]  TEFI(VN) %.2f [%s]  BIC %.2f [%s]  DL %.2f [%s]\n",
      phi, N, truth, acc[1], modal[1], acc[2], modal[2], acc[3], modal[3], acc[4], modal[4]))
}
