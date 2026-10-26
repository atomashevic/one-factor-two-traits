suppressMessages(library(lavaan))
pop <- function(m, lam, phi){L<-rbind(cbind(rep(lam,m),0),cbind(0,rep(lam,m))); S<-L%*%matrix(c(1,phi,phi,1),2)%*%t(L); diag(S)<-1; v<-paste0("x",1:(2*m)); dimnames(S)<-list(v,v); S}
best_srmr <- function(m,lam,phi){ S<-pop(m,lam,phi); v<-colnames(S); best<-NULL
  for (st in list(c(lam,.05),c(.5,.5),c(.6,.3))) { mod<-paste0("f =~ ",paste0("start(",rep(st,each=m),")*",v,collapse="+"))
    f<-suppressWarnings(try(cfa(mod,sample.cov=S,sample.nobs=1e6,std.lv=TRUE,sample.cov.rescale=FALSE),silent=TRUE)); if(inherits(f,"try-error")) next
    if(is.null(best)||f@optim$fx<best@optim$fx) best<-f }
  unname(fitMeasures(best,"srmr")) }
for (m in c(2,3,4)) { lam<-.5; g<-seq(0,.5,by=.01); s<-sapply(g,function(f) best_srmr(m,lam,f)); k<-which(s<=.08)[1]
  cat(sprintf("lam=.5 m=%d: SRMR at phi=0: %.3f; passes .08 from phi=%.2f\n", m, s[1], g[k])) }
