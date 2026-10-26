## Table 5 and the Section 7 numbers, recomputed from Kroenke et al. (2016).
## Published two-factor loadings: empirical/kroenke2016/table_s3_loadings.csv,
## extracted from the supplement by empirical/kroenke2016/extract_tables.py.
## Displayed cells are parsed from manuscript/sections/07_phq_illustration.tex.
## Population solutions for the published loadings are best found (multistart);
## with a common loading they come from the exact quartic.
## Run from the repository root: Rscript checks/phq_table5.R
suppressMessages(pkgload::load_all("undeRfactor", quiet = TRUE))
cat("R:", R.version.string, "| lavaan", as.character(utils::packageVersion("lavaan")), "\n")
s3 <- utils::read.csv(file.path("empirical", "kroenke2016", "table_s3_loadings.csv"))
tex <- readLines(file.path("manuscript", "sections", "07_phq_illustration.tex"))
text <- paste(tex, collapse = " ")
trials <- c(SCOPE = 0.912, ESCAPE = 0.865, INCPAD = 0.653)
fails <- 0
check <- function(label, shown, value) {
  ok <- identical(shown, value)
  cat(if (ok) "PASS" else "FAIL", label, "shown", shown, "computed", value, "\n")
  fails <<- fails + !ok
}
r3 <- function(x) sub("^0", "", sprintf("%.3f", x))
r2 <- function(x) sub("^0", "", sprintf("%.2f", x))
range_of <- function(x, f) paste(f(min(x)), "to", f(max(x)))

## Table 5 rows
body <- tex[(grep("label\\{tab:5\\}", tex) + 1):grep("bottomrule", tex)[1]]
rows <- strsplit(grep("&", body, value = TRUE)[-1], "&")
rows <- lapply(rows, function(r) trimws(sub("\\\\\\\\.*$", "", r)))
for (k in seq_along(trials)) {
  tr <- names(trials)[k]; phi <- trials[[k]]
  d <- s3[s3$trial == tr, ]; dep <- d$block == "depression"
  pub <- rows[[2 * k - 1]]; com <- rows[[2 * k]]
  f <- pop_fit(two_trait_cor(phi, d$two_factor[dep], d$two_factor[!dep]), nstart = 30)
  check(paste(tr, "published RMSEA"), pub[6], r3(f$RMSEA))
  check(paste(tr, "published CFI"), pub[7], r3(f$CFI))
  check(paste(tr, "published depression share"), pub[8], r2(f$share$share_A))
  grid <- lapply(seq(0.70, 0.80, by = 0.01), function(l)
    pop_fit(two_trait_cor(phi, l, m_A = 9, m_B = 7)))
  check(paste(tr, "common RMSEA"), com[6], range_of(sapply(grid, `[[`, "RMSEA"), r3))
  check(paste(tr, "common CFI"), com[7], range_of(sapply(grid, `[[`, "CFI"), r3))
  check(paste(tr, "common depression share"), com[8],
        range_of(sapply(grid, function(x) x$share$share_A), r2))
}

## Section 7 prose
found <- function(s) grepl(s, text, fixed = TRUE)
pf <- lapply(names(trials), function(tr) {
  d <- s3[s3$trial == tr, ]; dep <- d$block == "depression"
  list(d = d, dep = dep, fit = pop_fit(two_trait_cor(trials[[tr]], d$two_factor[dep],
                                                     d$two_factor[!dep]), nstart = 30))
})
names(pf) <- names(trials)
means <- sapply(pf, function(x) c(r2(mean(x$d$two_factor[x$dep])), r2(mean(x$d$two_factor[!x$dep]))))
check("block means", found(sprintf("block means of %s and %s in SCOPE, %s and %s in ESCAPE, and %s and %s in INCPAD",
                                    means[1, 1], means[2, 1], means[1, 2], means[2, 2], means[1, 3], means[2, 3])), TRUE)
sh <- sapply(pf, function(x) r2(x$fit$share$share_A))
check("population shares", found(sprintf("coefficient share of %s, %s and %s", sh[1], sh[2], sh[3])), TRUE)
## Section 5.2 criteria for blocks of unequal length: the total correlation of
## each block (decides at phi = 0) and m(m - 1){lambda^2/(1 - lambda^2)}^2 with
## one total-correlation summary loading per block (decides as phi -> 1)
tc_of <- function(l) { R <- tcrossprod(l); diag(R) <- 1
  -0.5 * as.numeric(determinant(R, logarithm = TRUE)$modulus) }
crit <- lapply(names(trials), function(tr) {
  x <- pf[[tr]]; dep <- x$dep; l <- x$d$two_factor
  pr <- predict_one_factor(cbind(ifelse(dep, l, 0), ifelse(dep, 0, l)), phi = trials[[tr]],
                           lambda_summary = "tc")
  lb <- pr$blocks$lambda_bar; m <- c(sum(dep), sum(!dep))
  side <- sapply(seq(0.01, 0.99, by = 0.01), function(ph)
    block_solutions(m[1], m[2], lb[1], lb[2], ph)$share_A[1] > 0.5)
  list(tc = c(tc_of(l[dep]), tc_of(l[!dep])), near1 = m * (m - 1) * (lb^2 / (1 - lb^2))^2,
       side = side, share = x$fit$share$share_A)
})
names(crit) <- names(trials)
tcs <- lapply(crit, function(x) sprintf("%.2f", x$tc))
check("total correlations", found(sprintf("larger total correlation in SCOPE (%s against %s) and INCPAD (%s against %s), whereas in ESCAPE the smaller loading difference leaves the larger total correlation with the longer depression block (%s against %s)",
                                          tcs$SCOPE[2], tcs$SCOPE[1], tcs$INCPAD[2], tcs$INCPAD[1],
                                          tcs$ESCAPE[1], tcs$ESCAPE[2])), TRUE)
for (tr in names(trials)) {
  x <- crit[[tr]]; dep_tc <- x$tc[1] > x$tc[2]
  check(paste(tr, "near-1 index orders the blocks as the total correlation"),
        x$near1[1] > x$near1[2], dep_tc)
  check(paste(tr, "favored block constant over phi and matching the share"),
        all(x$side == dep_tc) && (x$share > 0.5) == dep_tc, TRUE)
}
check("depression larger share only in ESCAPE",
      sapply(crit, function(x) x$share > 0.5), c(SCOPE = FALSE, ESCAPE = TRUE, INCPAD = FALSE))
cors <- sapply(pf, function(x) sprintf("%.3f", stats::cor(x$d$one_factor, x$fit$loadings)))
cors <- sub("^0", "", cors)
check("loading correlations", found(sprintf("correlations are %s, %s and %s", cors[1], cors[2], cors[3])), TRUE)
pub_sh <- sapply(pf, function(x) r2(factor_share(mean(x$d$one_factor[x$dep]), mean(x$d$one_factor[!x$dep]),
                                                mean(x$d$two_factor[x$dep]), mean(x$d$two_factor[!x$dep]),
                                                trials[[x$d$trial[1]]], mode = "signed")$share_A))
check("published one-factor shares", found(sprintf("depression coefficient shares of %s, %s and %s", pub_sh[1], pub_sh[2], pub_sh[3])), TRUE)
sc <- pf$INCPAD$fit$scores
check("INCPAD score correlations", found(sprintf("correlates %s with depression and %s with anxiety, whereas the unit-weighted sum correlates %s and %s",
                                                 r2(sc$cor_A[1]), r2(sc$cor_B[1]), r2(sc$cor_A[2]), r2(sc$cor_B[2]))), TRUE)
check("ESCAPE independence discrepancy", found(sprintf("which is %.1f at ESCAPE's correlation", pf$ESCAPE$fit$F_B)), TRUE)
d <- pf$ESCAPE$d; dep <- pf$ESCAPE$dep
bs <- blind_spot(d$two_factor[dep], d$two_factor[!dep], index = c("CFI", "RMSEA"), step = 0.02, nstart = 10)
check("ESCAPE cutoff crossings", found(sprintf("approximately %s and RMSEA reaches .06 from approximately %s",
                                               r3(bs[["CFI"]]), r3(bs[["RMSEA"]]))), TRUE)
e75 <- pop_fit(two_trait_cor(0.865, 0.75, m_A = 9, m_B = 7))
check("ESCAPE common loading .75", found(sprintf("CFI $=%s$ and RMSEA $=%s$", r3(e75$CFI), r3(e75$RMSEA))), TRUE)
check("SCOPE RMSEA at loading .80", found(sprintf("(%s at loading .80",
  sub("^0", "", sprintf("%.4f", pop_fit(two_trait_cor(0.912, 0.80, m_A = 9, m_B = 7))$RMSEA)))), TRUE)
if (fails) stop(fails, " Table 5 or Section 7 value(s) differ from the recomputation")
cat("PASS: all Table 5 cells and Section 7 values recomputed from Kroenke et al. (2016)\n")
