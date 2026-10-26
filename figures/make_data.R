## Archived sample exports: no fitting and no installed package required.
## Rscript figures/make_data.R --samples-only
args <- commandArgs(TRUE)
script_arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
script_path <- normalizePath(sub('^--file=', '', script_arg[1]))
repo <- dirname(dirname(script_path))
export_samples <- function() {
  started <- proc.time()[['elapsed']]
  read_summary <- function(name) read.csv(file.path(repo,'results/summary',name),stringsAsFactors=FALSE)
  q1 <- read_summary('q1_fit.csv'); q2 <- read_summary('q2_diagnostics_stage5.csv')
  q3 <- read_summary('q3_minima.csv'); pop <- readRDS(file.path(repo,'results/population.rds'))
  inputs <- c('results/summary/q1_fit.csv','results/summary/q2_diagnostics_stage5.csv',
              'results/summary/q3_minima.csv','results/population.rds')
  stopifnot(!anyDuplicated(q1[c('cell_id','est','model','index')]),
            !anyDuplicated(q2[c('cell_id','est')]), !anyDuplicated(q3[c('cell_id','est')]))
  write_export <- function(x,name) {
    destination <- file.path(repo,'figures',name)
    write.csv(x,destination,row.names=FALSE,na='NA')
    cat('Export',name,nrow(x),'rows; MD5',unname(tools::md5sum(destination)),'\n')
  }
  ## Figure 3: preserve studied range, use exact archived coordinates/rates.
  f3 <- subset(q1,design=='C' & N==250 & model=='1F' & phi>=.5 & phi<=.912 &
                 ((format=='cont' & est=='ML' & index %in% c('cfi','rmsea')) |
                  (format=='ord_skew' & est=='WLSMV' &
                    index %in% c('cfi','cfi_sc','cfi_rob','rmsea','rmsea_sc','rmsea_rob'))))
  f3$metric <- sub('_(sc|rob)$','',f3$index)
  f3$variant <- ifelse(f3$est=='ML','ml',ifelse(grepl('_rob$',f3$index),'rob',
                         ifelse(grepl('_sc$',f3$index),'sc','naive')))
  f3$n_pass <- round(f3$n*f3$p_accept)
  stopifnot(nrow(f3)==64,all(is.finite(f3$p_accept)),all(f3$n>0),
            all(abs(f3$n_pass/f3$n-f3$p_accept)<1e-12),
            !anyDuplicated(f3[c('metric','variant','phi')]))
  f3 <- f3[order(f3$metric,f3$variant,f3$phi),c('cell_id','pop_id','cond','phi','N','format','est',
             'index','metric','variant','n','n_pass','p_accept','cutoff')]
  write_export(f3,'fig3.csv')
  ## Figure 4: mean absolute normalized magnitude contrast, not a detection rate.
  f4 <- subset(q2,cond=='A_m4_l70' & N==2000 & phi<=.8 &
                 ((format=='cont' & est %in% c('ML','ULS')) | (format=='ord_sym' & est=='WLSMV')))
  stopifnot(nrow(f4)==36,all(f4$n_C_magnitude_valid>0),
            all(is.finite(f4$mean_abs_C_magnitude)),all(f4$n_nonpositive_min_psi==0),
            all(f4$n_C_magnitude_valid==f4$n_finite),!anyDuplicated(f4[c('est','phi')]))
  for(e in unique(f4$est)) stopifnot(length(unique(f4$phi[f4$est==e]))==12)
  f4$ml_threshold <- f4$phistar; f4$uls_threshold <- 1/(f4$mA+f4$mB)
  f4 <- f4[order(f4$est,f4$phi),c('cell_id','pop_id','cond','phi','N','format','est',
             'n_attempted','n_converged','n_finite','n_C_magnitude_valid','mean_abs_C_magnitude',
             'mcse_mean_abs_C_magnitude','mean_abs_magnitude_D','n_opposite_sign','p_opposite_sign',
             'ml_threshold','uls_threshold')]
  write_export(f4,'fig4.csv')
  ## Figure 5 left: raw flag counts, not certified sample minimum counts.
  f5l <- subset(q3,N==2000 & format=='cont' & est=='ML' &
             ((cond %in% c('B_l68_l72','B_l60_l80') & phi<=.5) |
              (cond=='C_9_7_l75' & phi<=.5)))
  stopifnot(nrow(f5l)==24)
  f5l$n_flag <- NA_integer_; f5l$n_higher_objective_weak <- NA_integer_
  for(i in seq_len(nrow(f5l))) {
    file <- sprintf('results/raw/cell_%04d.rds',f5l$cell_id[i])
    raw <- readRDS(file.path(repo,file)); inputs <- c(inputs,file)
    ms <- raw$ms[raw$ms$est=='ML',]
    flag <- if(raw$cell$weak=='A') ms$def_on_weak else ms$def_rev_on_weak
    den <- sum(!is.na(flag)); num <- sum(flag,na.rm=TRUE)
    stopifnot(den==1000,abs(num/den-f5l$p_default_on_weak_listed_first[i])<1e-12)
    f5l$n_flag[i] <- den; f5l$n_higher_objective_weak[i] <- num
  }
  ## Fold values are theoretical annotations; C has no recorded fold provenance.
  f5l <- f5l[order(f5l$cond,f5l$phi),c('cell_id','pop_id','cond','phi','N','format','est','weak',
           'n_flag','n_higher_objective_weak','p_default_on_weak_listed_first','phi_fold')]
  write_export(f5l,'fig5_left.csv')
  ## Figure 5 right: pooled selected replications, symmetric ordinal arm only.
  f5r <- do.call(rbind,lapply(c('ML','MLs','ULS','WLSMV'),function(e) {
    x <- subset(q2,design=='A' & N==2000 & est==e &
                  if(e %in% c('ML','ULS')) format=='cont' else format=='ord_sym')
    stopifnot(nrow(x)>0,all(x$n_nonpositive_min_psi==0),
              all(x$n_selected_M_gt_10==x$n_selected_first+x$n_selected_second))
    first <- sum(x$n_selected_first); second <- sum(x$n_selected_second)
    selected <- first+second; attempted <- sum(x$n_attempted)
    stopifnot(selected>0)
    data.frame(est=e,format=if(e %in% c('ML','ULS')) 'cont' else 'ord_sym',N=2000,
       n_cells=nrow(x),n_attempted=attempted,n_selected=selected,n_first=first,n_second=second,
       p_first=first/selected,selection_abs_M=.1,
       source_cell_ids=paste(x$cell_id,collapse=';'))
  }))
  stopifnot(identical(as.integer(f5r$n_first),c(24767L,11508L,12385L,2649L)),
            identical(as.integer(f5r$n_selected),c(28719L,13501L,18344L,4746L)))
  write_export(f5r,'fig5_right.csv')
  ## Table 4: baseline inversion must be valid in every paired replication.
  requests <- data.frame(cond=c(rep('B_l68_l72',5),rep('B_l60_l80',2),rep('C_9_7_l75',2)),
                         phi=c(0,.1,.2,.24,.3,0,.1,0,.1))
  table_rows <- lapply(seq_len(nrow(requests)),function(k) {
    request <- requests[k,]
    population <- subset(pop,cond==request$cond & abs(phi-request$phi)<1e-12 &
                               format=='cont' & est=='ML' & model=='1F')
    stopifnot(nrow(population)==1,is.finite(population$cfi),
              abs(population$cfi-(1-population$F/population$FB))<1e-12)
    out <- data.frame(row=k,pop_id=population$pop_id,cond=request$cond,phi=request$phi,
                      population_CFI=population$cfi,population_F=population$F,population_FB=population$FB)
    for(N in c(500,2000)) {
      lookup <- q3[q3$cond==request$cond & abs(q3$phi-request$phi)<1e-12 &
                   q3$format=='cont' & q3$est=='ML' & q3$N==N,]
      stopifnot(nrow(lookup)==1)
      file <- sprintf('results/raw/cell_%04d.rds',lookup$cell_id)
      raw <- readRDS(file.path(repo,file)); inputs <- c(inputs,file)
      default <- subset(raw$fits,est=='ML' & model=='1F'); ms <- subset(raw$ms,est=='ML')
      stopifnot(!anyDuplicated(default$rep),!anyDuplicated(ms$rep))
      match_row <- match(ms$rep,default$rep); stopifnot(!anyNA(match_row),nrow(ms)==1000)
      default <- default[match_row,]
      model_nc <- default$chisq-default$df
      baseline_nc <- model_nc/(1-default$cfi)
      stopifnot(all(default$converged),all(ms$n_distinct==ms$n_admissible),all(!ms$any_heywood_min),
                all(default$cfi>0 & default$cfi<1),all(model_nc>0),all(baseline_nc>model_nc),
                all(abs(default$chisq-2*N*default$fx)<1e-10),all(is.finite(ms$best_fx)))
      best_nc <- 2*N*ms$best_fx-default$df
      best_cfi <- 1-pmax(best_nc,0)/pmax(best_nc,baseline_nc,0)
      stopifnot(all(((default$fx-ms$best_fx)>1e-7)==ms$def_nonglobal))
      flag <- ms$def_on_weak # Table 4 uses design order, including stronger-first Design C.
      den <- sum(!is.na(flag)); stopifnot(den==1000)
      suffix <- as.character(N)
      out[[paste0('cell_id_',suffix)]] <- lookup$cell_id
      out[[paste0('n_',suffix)]] <- nrow(ms)
      out[[paste0('default_CFI_',suffix)]] <- mean(default$cfi)
      out[[paste0('best_found_CFI_',suffix)]] <- mean(best_cfi)
      out[[paste0('n_higher_objective_weak_',suffix)]] <- sum(flag)
      out[[paste0('p_higher_objective_weak_',suffix)]] <- mean(flag)
      out[[paste0('min_model_nc_',suffix)]] <- min(model_nc)
      out[[paste0('min_baseline_nc_',suffix)]] <- min(baseline_nc)
      out[[paste0('min_baseline_margin_',suffix)]] <- min(baseline_nc-model_nc)
    }
    out
  })
  tab <- do.call(rbind,table_rows); write_export(tab,'table4.csv')
  cat('Table4 paired replication rows:',sum(tab$n_500+tab$n_2000),'\n')
  cat('No fitting performed. Retained clusters have positive recorded uniqueness; stationarity and global optimality are unverified.\n')
  cat('Inputs MD5:\n'); print(tools::md5sum(file.path(repo,unique(inputs))))
  cat('Script MD5:',unname(tools::md5sum(script_path)),'\n')
  cat('Elapsed seconds:',proc.time()[['elapsed']]-started,'\n')
}
if('--samples-only' %in% args) {
  if(!identical(args,'--samples-only')) stop('Only --samples-only is supported for archived exports')
  export_samples(); quit(save='no')
}
if(length(args)) stop('Unknown arguments')

## Run from the repository root: Rscript figures/make_data.R
library(undeRfactor)   # remotes::install_github("atomashevic/undeRfactor")
## Figure 1: block loadings, ML (closed form) and ULS, m = 4, lambda = .7
phis <- sort(unique(c(seq(0, 0.95, by = 0.005), 0.125, phi_star(4, 0.7))))
f1 <- do.call(rbind, lapply(phis, function(p) {
  R <- two_trait_cor(p, 0.7, m_A = 4)
  ml <- pop_fit(R); uls <- pop_fit(R, "ULS", nstart = 5)
  data.frame(phi = p, ml_A = max(ml$minima$l_A[1], ml$minima$l_B[1]), ml_B = min(ml$minima$l_A[1], ml$minima$l_B[1]),
             uls_A = max(uls$minima$l_A[1], uls$minima$l_B[1]), uls_B = min(uls$minima$l_A[1], uls$minima$l_B[1]))
}))
write.csv(f1, "figures/fig1.csv", row.names = FALSE)
## Figure 2: share of the followed trait
grid <- sort(unique(c(seq(0, 0.9, by = 0.01), 0.4089, 0.4411, 0.2415)))
designs <- list(eq = list(0.7, 0.7), mir = list(c(.6, .65, .75, .8), c(.6, .65, .75, .8)),
                mild = list(0.68, 0.72), strong = list(0.60, 0.80))
f2 <- do.call(rbind, lapply(names(designs), function(d) do.call(rbind, lapply(grid, function(p) {
  R <- two_trait_cor(p, designs[[d]][[1]], designs[[d]][[2]], m_A = 4)
  f <- pop_fit(R, nstart = 10)
  loc <- if (nrow(f$minima) > 1 && abs(diff(f$minima$F[1:2])) > 1e-9) f$minima$share_A[2] else NA
  data.frame(design = d, phi = p, share = f$share$share, local = if (is.na(loc)) NA else max(loc, 1 - loc))
}))))
write.csv(f2, "figures/fig2.csv", row.names = FALSE)
cat("check vs widget: mir .3", f2$share[f2$design=="mir" & abs(f2$phi-.3)<1e-9], "(.821) mild .3", f2$share[f2$design=="mild" & abs(f2$phi-.3)<1e-9],
    "(.8501) strong .3", f2$share[f2$design=="strong" & abs(f2$phi-.3)<1e-9], "(.946) mild local .2", f2$local[f2$design=="mild" & abs(f2$phi-.2)<1e-9], "(.8176)\n")
cat("mild local last phi:", max(f2$phi[f2$design=="mild" & !is.na(f2$local)]), "\n")
