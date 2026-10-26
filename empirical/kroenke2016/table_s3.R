## Section 7 with the published heterogeneous loadings of Kroenke et al. (2016).
##
## Table S3 gives standardized one-factor, two-factor and bifactor loadings for
## each trial, estimated in Mplus with the items treated as ordinal (logistic
## link). Using the two-factor loadings and the reported trait correlations,
## this computes the population one-factor approximation of each trial: from
## one total-correlation summary loading per block (T2, exact quartic) and from
## the full implied matrix (T3, best found). It compares T3's loadings with the published one-factor
## loadings, which come from the trial data. The pooled correlation matrix of
## Table S2 is analysed as a separate check.
##
## Run from the repository root:
##   OPENBLAS_NUM_THREADS=1 Rscript empirical/kroenke2016/table_s3.R
## Outputs: build/empirical/kroenke2016/*.csv

suppressMessages(pkgload::load_all("undeRfactor", quiet = TRUE))
here <- file.path("empirical", "kroenke2016")
out_dir <- file.path("build", "empirical", "kroenke2016")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

s3 <- utils::read.csv(file.path(here, "table_s3_loadings.csv"))
## Kroenke et al. (2016), Table 5: trait correlations and sample sizes
trials <- data.frame(trial = c("SCOPE", "ESCAPE", "INCPAD"), phi = c(0.912, 0.865, 0.653),
                     n = c(250, 241, 405))

rows <- list(); items <- list()
for (k in seq_len(nrow(trials))) {
  tr <- trials$trial[k]; phi <- trials$phi[k]
  d <- s3[s3$trial == tr, ]
  dep <- d$block == "depression"
  Lambda <- cbind(ifelse(dep, d$two_factor, 0), ifelse(dep, 0, d$two_factor))
  rownames(Lambda) <- d$item
  pr <- predict_one_factor(Lambda, phi = phi, lambda_summary = "tc", nstart = 30)
  t3 <- pop_fit(two_trait_cor(phi, d$two_factor[dep], d$two_factor[!dep]), nstart = 30)
  sc <- t3$scores
  for (t in c("T2", "T3")) {
    r <- pr$tiers[pr$tiers$tier == t, ]
    rows[[length(rows) + 1]] <- data.frame(
      trial = tr, phi = phi, tier = t, lambda_dep = pr$blocks$lambda_bar[1],
      lambda_anx = pr$blocks$lambda_bar[2], F = r$F, CFI = r$CFI, TLI = r$TLI,
      RMSEA = r$RMSEA, SRMR = r$SRMR, dep_coefficient_share = r$share_A,
      dep_score_share = r$score_share_A, n_minima = r$n_minima,
      score_cor_dep = if (t == "T3") sc$cor_A[1] else NA,
      score_cor_anx = if (t == "T3") sc$cor_B[1] else NA,
      sum_share_dep = if (t == "T3") sc$share_A[2] else NA,
      sum_cor_dep = if (t == "T3") sc$cor_A[2] else NA,
      sum_cor_anx = if (t == "T3") sc$cor_B[2] else NA)
  }
  items[[k]] <- data.frame(trial = tr, item = d$item, published_one_factor = d$one_factor,
                           predicted_one_factor = unname(t3$loadings))
}
res <- do.call(rbind, rows)
it <- do.call(rbind, items)
agree <- do.call(rbind, lapply(split(it, it$trial), function(x)
  data.frame(trial = x$trial[1], r = stats::cor(x$published_one_factor, x$predicted_one_factor),
             mean_abs_diff = mean(abs(x$published_one_factor - x$predicted_one_factor)),
             max_abs_diff = max(abs(x$published_one_factor - x$predicted_one_factor)),
             published_dep_mean = mean(x$published_one_factor[1:9]),
             published_anx_mean = mean(x$published_one_factor[10:16]),
             predicted_dep_mean = mean(x$predicted_one_factor[1:9]),
             predicted_anx_mean = mean(x$predicted_one_factor[10:16]))))
agree <- agree[match(trials$trial, agree$trial), ]

## Published one-factor loadings alone: their block-mean coefficient share
## against the published two-factor loadings (as the paper's share is defined)
pub <- do.call(rbind, lapply(seq_len(nrow(trials)), function(k) {
  d <- s3[s3$trial == trials$trial[k], ]; dep <- d$block == "depression"
  s <- factor_share(mean(d$one_factor[dep]), mean(d$one_factor[!dep]),
                    mean(d$two_factor[dep]), mean(d$two_factor[!dep]), trials$phi[k],
                    mode = "signed")
  data.frame(trial = trials$trial[k], published_dep_coefficient_share = s$share_A)
}))

## Section 7, Table 5: common-loading population values at loadings .70-.80
common <- do.call(rbind, lapply(seq_len(nrow(trials)), function(k) {
  f <- lapply(c(0.70, 0.80), function(l)
    pop_fit(two_trait_cor(trials$phi[k], l, m_A = 9, m_B = 7)))
  data.frame(trial = trials$trial[k],
             common_CFI = paste(sprintf("%.3f", sapply(f, `[[`, "CFI")), collapse = "-"),
             common_RMSEA = paste(sprintf("%.3f", sapply(f, `[[`, "RMSEA")), collapse = "-"),
             common_dep_share = paste(sprintf("%.2f", sapply(f, function(x) x$share$share_A)),
                                      collapse = "-"))
}))

## Table S2: all trials pooled (n = 896), normal-theory ML on the printed
## correlations. The trials differ in population, so this is a check of the
## workflow on a published matrix, not an estimate for any one trial.
S2 <- as.matrix(utils::read.csv(file.path(here, "table_s2_correlations.csv"), row.names = 1))
nm <- rownames(S2)
f2 <- lavaan::cfa(paste0("dep =~ ", paste(nm[1:9], collapse = " + "), "\nanx =~ ",
                         paste(nm[10:16], collapse = " + ")),
                  sample.cov = S2, sample.nobs = 896, std.lv = TRUE)
pr2 <- predict_one_factor(f2, lambda_summary = "tc", nstart = 30)
oc <- order_check(S2, list(dep = nm[1:9], anx = nm[10:16]), nobs = 896, nstart = 20)
fm2 <- lavaan::fitMeasures(f2, c("cfi", "rmsea", "srmr"))
pooled <- data.frame(
  phi_hat = pr2$phi, two_factor_CFI = fm2[["cfi"]], two_factor_RMSEA = fm2[["rmsea"]],
  T2_F = pr2$tiers$F[pr2$tiers$tier == "T2"], T3_F = pr2$tiers$F[pr2$tiers$tier == "T3"],
  observed_F = 2 * oc$solutions$fmin[1],
  observed_minus_two_factor_F = 2 * oc$solutions$fmin[1] - 2 * lavaan::fitMeasures(f2, "fmin"),
  observed_CFI = oc$solutions$cfi[1], observed_RMSEA = oc$solutions$rmsea[1],
  T3_dep_share = pr2$tiers$share_A[pr2$tiers$tier == "T3"],
  observed_dep_share = oc$solutions$share_A[1], solutions = nrow(oc$solutions),
  order_dependent = oc$order_dependent)

utils::write.csv(res, file.path(out_dir, "table_s3_predictions.csv"), row.names = FALSE)
utils::write.csv(it, file.path(out_dir, "table_s3_item_loadings.csv"), row.names = FALSE)
utils::write.csv(agree, file.path(out_dir, "table_s3_loading_agreement.csv"), row.names = FALSE)
utils::write.csv(pooled, file.path(out_dir, "table_s2_pooled.csv"), row.names = FALSE)
options(width = 140)
print(res, digits = 3, row.names = FALSE)
print(merge(merge(agree, pub), common), digits = 3, row.names = FALSE)
print(pooled, digits = 3, row.names = FALSE)
