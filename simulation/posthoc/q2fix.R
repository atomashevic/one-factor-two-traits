## Signed, descriptive archive diagnostics. Never changes canonical q2_tilt.csv.
## Rscript simulation/posthoc/q2fix.R [--input FILE] [--output FILE]
argv <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), '../..'))
opt <- list(input = file.path(root, 'results/fits1f.rds'),
            output = file.path(root, 'results/summary/q2_diagnostics_stage5.csv'))
if (length(argv) %% 2L) stop('Options require a value')
if (length(argv)) for (i in seq(1L, length(argv), by = 2L)) {
  key <- sub('^--', '', argv[i])
  if (!key %in% names(opt)) stop('Unknown option: ', argv[i])
  opt[[key]] <- argv[i + 1L]
}
if (file.exists(opt$output)) stop('Output exists; preserve it and choose a new --output')
cat('Command:', paste(commandArgs(), collapse = ' '), '\n')
cat('Thread limits:', paste(Sys.getenv(c('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS',
                                        'MKL_NUM_THREADS')), collapse = '/'), '\n')

## Selection concerns magnitudes, irrespective of block polarity. Its nominal
## binomial SE is conditional on the selected count within one cell only.
selection_counts <- function(M, cutoff = .1) {
  first <- sum(M > cutoff); second <- sum(M < -cutoff)
  selected <- first + second
  p <- if (selected) first / selected else NA_real_
  list(n_selected_M_gt_10 = selected, n_selected_first = first,
       n_selected_second = second, p_selected_first = p,
       mcse_p_selected_first = if (selected) sqrt(p * (1 - p) / selected) else NA_real_)
}

contrasts <- function(a, b, mA, mB, rA, rB, phi, tol = 1e-8) {
  finite <- is.finite(a) & is.finite(b)
  orient_sum <- mA * a + mB * b
  flip <- finite & orient_sum < 0
  a[flip] <- -a[flip]; b[flip] <- -b[flip]
  ref_valid <- is.finite(rA) & is.finite(rB) & rA > 0 & rB > 0
  u <- a / rA; v <- b / rB
  margin <- abs(u + v) / (abs(u) + abs(v))
  stable <- finite & ref_valid & is.finite(margin) & margin > tol
  C <- rep(NA_real_, length(a)); C[stable] <- (u[stable] - v[stable]) / (u[stable] + v[stable])
  magnitude_valid <- finite & ref_valid & is.finite(abs(u)+abs(v)) & (abs(u)+abs(v))>0
  C_magnitude <- rep(NA_real_,length(a))
  C_magnitude[magnitude_valid] <- (abs(u[magnitude_valid])-abs(v[magnitude_valid])) /
    (abs(u[magnitude_valid])+abs(v[magnitude_valid]))
  bounded <- stable & a >= 0 & b >= 0 & (a > 0 | b > 0)
  opposite <- finite & ((a < 0 & b > 0) | (a > 0 & b < 0))
  share_valid <- stable & is.finite(phi) & abs(phi) < 1
  share <- rep(NA_real_, length(a))
  share[share_valid] <- .5 + (1 + phi[share_valid]) / (2 * (1 - phi[share_valid])) * C[share_valid]
  data.frame(a, b, D = a - b, C, C_magnitude, magnitude_valid, margin, share, finite, flip,
             orientation_violation = finite & orient_sum < -tol,
             orientation_tie = finite & abs(orient_sum) <= tol,
             ref_valid, stable, bounded, opposite, share_valid)
}
## Mathematical behavior and invalid-domain checks, independent of archived outcomes.
a <- c(.6,.6,.4,.4,0,.4,Inf,.6)
b <- c(.2,-.2,.4,-.4,0,.2,.2,-.6+1e-12)
rA <- rep(.7,8); rA[6] <- 0
z <- contrasts(a,b,rep(4,8),rep(4,8),rA,rep(.7,8),rep(0,8))
n <- contrasts(-a,-b,rep(4,8),rep(4,8),rA,rep(.7,8),rep(0,8))
s <- contrasts(b,a,rep(4,8),rep(4,8),rep(.7,8),rA,rep(0,8))
stopifnot(isTRUE(all.equal(z$C,n$C)), isTRUE(all.equal(z$C,-s$C)),
          z$D[1] != z$D[2], z$C_magnitude[1] == z$C_magnitude[2],
          abs(z$C_magnitude[1]-.5)<1e-12,
          isTRUE(all.equal(z$C_magnitude,n$C_magnitude)),
          isTRUE(all.equal(z$C_magnitude,-s$C_magnitude)), z$bounded[1], !z$bounded[2], z$opposite[2],
          z$C[3] == 0, all(is.na(z$C[4:8])),
          identical(c(9.5 / .5), 19))
M <- abs(z$a[1:5]) - abs(z$b[1:5])
Mn <- abs(n$a[1:5]) - abs(n$b[1:5])
Ms <- abs(s$a[1:5]) - abs(s$b[1:5])
empty_selection <- selection_counts(c(0, .05, -.05))
sel <- selection_counts(M); exchanged <- selection_counts(Ms)
stopifnot(identical(sel, selection_counts(Mn)),
          sel$n_selected_M_gt_10 == exchanged$n_selected_M_gt_10,
          sel$n_selected_first == exchanged$n_selected_second,
          sel$n_selected_second == exchanged$n_selected_first,
          empty_selection$n_selected_M_gt_10 == 0,
          is.na(empty_selection$p_selected_first), is.na(empty_selection$mcse_p_selected_first),
          selection_counts(c(.1, -.1))$n_selected_M_gt_10 == 0)
cat('Contrast invariance, block exchange, sign, reference and cancellation checks passed.\n')
cat('Magnitude-selection sign invariance, block exchange, strict cutoff and zero denominator checks passed.\n')
started <- proc.time()[['elapsed']]
f <- readRDS(opt$input)
stopifnot(!anyDuplicated(f[c('cell_id','est','rep')]))
provenance <- attr(f, 'provenance')
if (!is.null(provenance)) {
  cat('Input provenance:', provenance$note, '\n')
  cat('Input source files:', length(provenance$input_md5), '\n')
  print(unique(provenance$source_runtime[c('R','lavaan')]))
}
rA <- ifelse(f$est == 'MLs', f$refA, f$lamA)
rB <- ifelse(f$est == 'MLs', f$refB, f$lamB)
rp <- ifelse(f$est == 'MLs', f$refphi, f$phi)
z <- contrasts(f$lA,f$lB,f$mA,f$mB,rA,rB,rp)
## Reproduce the former rule solely to quantify its failure and verify extraction.
u <- abs(f$lA) / rA; v <- abs(f$lB) / rB
oldA <- .5 + (1 + rp) / (2 * (1 - rp)) * (u - v) / (u + v)
oldshare <- pmax(oldA, 1 - oldA)
oldD <- abs(f$lA) - abs(f$lB)
old_ok <- f$converged %in% TRUE & is.finite(oldshare)
calkey <- paste(f$cond,f$N,f$format,f$est,sep='|')
cal <- old_ok & abs(f$phi - .9) < 1e-9
q95 <- tapply(oldshare[cal], calkey[cal], quantile, probs=.95)
oldq <- unname(q95[calkey])
key <- paste(f$cell_id,f$est,sep='|')
groups <- split(seq_len(nrow(f)), key)
safe_mean <- function(x) if(length(x)) mean(x) else NA_real_
mcse <- function(x) if(length(x)>1L) sd(x)/sqrt(length(x)) else NA_real_
qq <- function(x,p) if(length(x)) unname(quantile(x,p)) else NA_real_
rows <- lapply(groups, function(i) {
  out <- f[i[1],c('cell_id','pop_id','design','cond','phi','N','format','est',
                  'mA','mB','lamA','lamB','phistar')]
  cv <- i[f$converged[i] %in% TRUE]
  j <- cv[z$finite[cv]]; k <- j[z$bounded[j]]; stable <- j[z$stable[j]]
  opposite <- j[z$opposite[j]]; old <- i[old_ok[i]]
  out$n_attempted <- length(i); out$n_ok <- sum(f$ok[i] %in% TRUE)
  out$n_converged <- length(cv); out$n_finite <- length(j)
  out$n_positive_min_psi <- sum(is.finite(f$min_psi[cv]) & f$min_psi[cv]>0)
  out$n_nonpositive_min_psi <- sum(is.finite(f$min_psi[cv]) & f$min_psi[cv]<=0)
  out$n_unknown_min_psi <- sum(!is.finite(f$min_psi[cv]))
  out$n_common_flips <- sum(z$flip[j]); out$n_orientation_violation <- sum(z$orientation_violation[j])
  out$n_orientation_tie <- sum(z$orientation_tie[j]); out$n_opposite_sign <- length(opposite)
  out$n_invalid_reference <- sum(!z$ref_valid[j]); out$n_unstable_denominator <- sum(z$ref_valid[j] & !z$stable[j])
  out$n_C_stable <- length(stable); out$n_C_bounded <- length(k)
  out$p_opposite_sign <- safe_mean(z$opposite[j]); out$mcse_p_opposite_sign <- mcse(as.numeric(z$opposite[j]))
  out$mean_D <- safe_mean(z$D[j]); out$mcse_mean_D <- mcse(z$D[j])
  out$mean_abs_D <- safe_mean(abs(z$D[j])); out$mcse_mean_abs_D <- mcse(abs(z$D[j]))
  out$q05_abs_D <- qq(abs(z$D[j]),.05); out$median_abs_D <- qq(abs(z$D[j]),.5)
  out$q95_abs_D <- qq(abs(z$D[j]),.95)
  out$mean_abs_D_opposite <- safe_mean(abs(z$D[opposite]))
  magnitude_D <- abs(z$a[j]) - abs(z$b[j])
  out$mean_magnitude_D <- safe_mean(magnitude_D)
  out$mean_abs_magnitude_D <- safe_mean(abs(magnitude_D))
  out$mcse_mean_abs_magnitude_D <- mcse(abs(magnitude_D))
  out$mean_abs_magnitude_D_opposite <- safe_mean(abs(abs(z$a[opposite])-abs(z$b[opposite])))
  selected <- selection_counts(magnitude_D)
  for (field in names(selected)) out[[field]] <- selected[[field]]
  km <- j[z$magnitude_valid[j]]
  out$n_C_magnitude_valid <- length(km)
  out$mean_C_magnitude <- safe_mean(z$C_magnitude[km])
  out$mean_abs_C_magnitude <- safe_mean(abs(z$C_magnitude[km]))
  out$mcse_mean_abs_C_magnitude <- mcse(abs(z$C_magnitude[km]))
  out$median_abs_C_magnitude <- qq(abs(z$C_magnitude[km]),.5)
  out$q95_abs_C_magnitude <- qq(abs(z$C_magnitude[km]),.95)
  out$mean_C_bounded <- safe_mean(z$C[k]); out$mean_abs_C_bounded <- safe_mean(abs(z$C[k]))
  out$mcse_mean_abs_C_bounded <- mcse(abs(z$C[k])); out$median_abs_C_bounded <- qq(abs(z$C[k]),.5)
  out$q95_abs_C_bounded <- qq(abs(z$C[k]),.95)
  out$mean_abs_C_opposite_stable <- safe_mean(abs(z$C[opposite[z$stable[opposite]]]))
  for (d in c(.05,.10,.20)) {
    suffix <- sprintf('%02d',round(100*d)); event <- as.numeric(abs(z$D[j]) > d)
    out[[paste0('p_abs_D_gt_',suffix)]] <- safe_mean(event)
    out[[paste0('mcse_p_abs_D_gt_',suffix)]] <- mcse(event)
    magnitude_event <- as.numeric(abs(magnitude_D) > d)
    out[[paste0('p_abs_M_gt_',suffix)]] <- safe_mean(magnitude_event)
    out[[paste0('mcse_p_abs_M_gt_',suffix)]] <- mcse(magnitude_event)
  }
  ## Preserve the all-finite diagnostics, including the seven improper fits.
  ## Positive recorded uniqueness is a domain filter, not minimum certification.
  positive <- j[is.finite(f$min_psi[j]) & f$min_psi[j] > 0]
  nonpositive <- j[is.finite(f$min_psi[j]) & f$min_psi[j] <= 0]
  unknown <- j[!is.finite(f$min_psi[j])]
  positive_M <- abs(z$a[positive]) - abs(z$b[positive])
  positive_km <- positive[z$magnitude_valid[positive]]
  out$n_positive_min_psi_finite <- length(positive)
  out$n_nonpositive_min_psi_finite <- length(nonpositive)
  out$n_unknown_min_psi_finite <- length(unknown)
  out$n_opposite_sign_positive_min_psi <- sum(z$opposite[positive])
  out$p_opposite_sign_positive_min_psi <- safe_mean(z$opposite[positive])
  out$mcse_p_opposite_sign_positive_min_psi <- mcse(as.numeric(z$opposite[positive]))
  out$mean_abs_M_positive_min_psi <- safe_mean(abs(positive_M))
  out$mcse_mean_abs_M_positive_min_psi <- mcse(abs(positive_M))
  out$n_C_magnitude_valid_positive_min_psi <- length(positive_km)
  out$mean_abs_C_magnitude_positive_min_psi <- safe_mean(abs(z$C_magnitude[positive_km]))
  out$mcse_mean_abs_C_magnitude_positive_min_psi <- mcse(abs(z$C_magnitude[positive_km]))
  positive_selection <- selection_counts(positive_M)
  out$n_selected_admissible_first <- positive_selection$n_selected_first
  out$n_selected_admissible_second <- positive_selection$n_selected_second
  out$n_selected_M_gt_10_positive_min_psi <- positive_selection$n_selected_M_gt_10
  out$p_selected_first_positive_min_psi <- positive_selection$p_selected_first
  out$mcse_p_selected_first_positive_min_psi <- positive_selection$mcse_p_selected_first
  out$n_selected_M_gt_10_nonpositive_min_psi <- sum(abs(abs(z$a[nonpositive])-abs(z$b[nonpositive])) > .1)
  out$n_selected_M_gt_10_unknown_min_psi <- sum(abs(abs(z$a[unknown])-abs(z$b[unknown])) > .1)
  sv <- j[z$share_valid[j]]
  out$n_covariance_share_valid <- length(sv)
  out$n_covariance_share_outside_01 <- sum(z$share[sv]<0 | z$share[sv]>1)
  out$legacy_n <- length(old); out$legacy_q95 <- oldq[i[1]]
  out$legacy_p_tilted <- safe_mean(oldshare[old] > oldq[old])
  out$legacy_mean_share <- safe_mean(oldshare[old]); out$legacy_mean_abs_d <- safe_mean(abs(oldD[old]))
  out$legacy_zero_phi_impossible <- out$phi == 0 && is.finite(out$legacy_q95) && out$legacy_q95 >= 1
  out
})
Q <- do.call(rbind,rows); rownames(Q) <- NULL; Q <- Q[order(Q$cell_id,Q$est),]
legacy <- read.csv(file.path(root,'results/summary/q2_tilt.csv'))
m <- match(paste(Q$cell_id,Q$est),paste(legacy$cell_id,legacy$est))
stopifnot(!anyNA(m),nrow(Q)==nrow(legacy))
legacy_fields <- c(n='legacy_n',q95='legacy_q95',p_tilted='legacy_p_tilted',
                   mean_share='legacy_mean_share',mean_abs_d='legacy_mean_abs_d')
for (original in names(legacy_fields)) {
  field <- legacy_fields[[original]]
  delta <- abs(Q[[field]]-legacy[[original]][m])
  stopifnot(identical(is.na(Q[[field]]),is.na(legacy[[original]][m])), all(delta[is.finite(delta)]<1e-9))
  cat('Legacy reproduction',original,'max absolute error',max(delta,na.rm=TRUE),'\n')
}
stopifnot(all(Q$n_C_bounded+Q$n_opposite_sign <= Q$n_finite),
          all(Q$n_positive_min_psi+Q$n_nonpositive_min_psi+Q$n_unknown_min_psi==Q$n_converged),
          all(Q$n_selected_M_gt_10 == Q$n_selected_first + Q$n_selected_second),
          all(Q$n_selected_M_gt_10 <= Q$n_finite),
          all(Q$n_positive_min_psi_finite + Q$n_nonpositive_min_psi_finite + Q$n_unknown_min_psi_finite == Q$n_finite),
          all(Q$n_selected_M_gt_10_positive_min_psi == Q$n_selected_admissible_first + Q$n_selected_admissible_second),
          all(Q$n_selected_M_gt_10 == Q$n_selected_M_gt_10_positive_min_psi + Q$n_selected_M_gt_10_nonpositive_min_psi + Q$n_selected_M_gt_10_unknown_min_psi),
          identical(is.na(Q$p_selected_first), Q$n_selected_M_gt_10 == 0),
          identical(is.na(Q$mcse_p_selected_first), Q$n_selected_M_gt_10 == 0))
stage4 <- read.csv(file.path(root,'results/summary/q2_diagnostics_stage4.csv'))
m4 <- match(paste(Q$cell_id,Q$est),paste(stage4$cell_id,stage4$est))
stopifnot(nrow(Q) == 1946L, nrow(Q) == nrow(stage4), !anyNA(m4), all(names(stage4) %in% names(Q)))
maximum_stage4_error <- 0
for (field in names(stage4)) {
  previous <- stage4[[field]][m4]; current <- Q[[field]]
  if (is.numeric(previous)) {
    stopifnot(identical(is.na(current),is.na(previous)),
              isTRUE(all.equal(current,previous,tolerance=1e-12,check.attributes=FALSE)))
    error <- abs(current-previous)
    if (any(is.finite(error))) maximum_stage4_error <- max(maximum_stage4_error,error[is.finite(error)])
  } else stopifnot(identical(current,previous))
}
cat('Stage 4 preservation:',length(names(stage4)),'columns; maximum numerical difference',maximum_stage4_error,'\n')
fig5 <- Q[Q$design == 'A' & Q$N == 2000 &
            ((Q$format == 'cont' & Q$est %in% c('ML','ULS')) |
             (Q$format == 'ord_sym' & Q$est %in% c('MLs','WLSMV'))),]
fig5_counts <- aggregate(cbind(n_selected_first,n_selected_second,n_selected_M_gt_10)~est+format,fig5,sum)
expected_first <- c(ML=24767,MLs=11508,ULS=12385,WLSMV=2649)
expected_selected <- c(ML=28719,MLs=13501,ULS=18344,WLSMV=4746)
stopifnot(nrow(fig5_counts)==4L,
          all(fig5_counts$n_selected_first == expected_first[fig5_counts$est]),
          all(fig5_counts$n_selected_M_gt_10 == expected_selected[fig5_counts$est]),
          sum(fig5$n_nonpositive_min_psi_finite) == 0,
          sum(fig5$n_unknown_min_psi_finite) == 0)
fig5_counts$p_selected_first <- fig5_counts$n_selected_first / fig5_counts$n_selected_M_gt_10
cat('Figure 5 right pooled magnitude-selection counts:\n'); print(fig5_counts)
write.csv(Q,opt$output,row.names=FALSE)
cat(nrow(Q),'cell-estimator summaries; elapsed',proc.time()[['elapsed']]-started,'seconds\n')
cat('Input MD5:',unname(tools::md5sum(opt$input)),'\n')
cat('Script MD5:',unname(tools::md5sum(script)),'\n')
cat('Output MD5:',unname(tools::md5sum(opt$output)),'\n')
cat('R:',R.version.string,'\n')
cat('SHA256 input, script and output:\n')
system2('sha256sum',shQuote(c(opt$input,script,opt$output)))
cat('All quantities condition on convergence and their stated domain. No minimum or score certification.\n')
cat('Magnitude imbalance M=abs(a)-abs(b) is distinct from signed D; neither is a score share.\n')
cat('C_magnitude normalizes magnitude imbalance on all valid-reference nonzero pairs, regardless of polarity.\n')
cat('C_bounded requires nonnegative oriented block means; opposite-sign cases remain in D.\n')
cat('MCSE uses independent replications within each cell. Formats are not pooled.\n')
cat('Magnitude selection is abs(M)>.1; first means M>.1, second M<-.1. Conditional proportions use selected counts, not n_finite.\n')
cat('Conditional selection MCSE is the nominal per-cell binomial SE. Pooled heterogeneous-cell bars require their own uncertainty calculation.\n')
cat('Fields named admissible_first/second mean positive recorded min_psi only; they do not certify a minimum.\n')
print(aggregate(cbind(n_attempted,n_converged,n_finite,n_opposite_sign,n_common_flips,
                     n_C_bounded,n_covariance_share_outside_01)~design,Q,sum))
cat('Zero-correlation groups forced to zero by legacy q95 >= 1:',sum(Q$legacy_zero_phi_impossible),
    'of',sum(Q$phi==0),'\n')
