## Three archived-seed ML refits under the current R/lavaan versions.
## OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript checks/sample_ctrl.R
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(arg)) stop("Run this file with Rscript.")
script <- normalizePath(sub("^--file=", "", arg[1]), mustWork = TRUE)
paper_dir <- dirname(dirname(script))
pop <- new.env(parent = globalenv())
source(file.path(paper_dir, "checks", "order_ctrl.R"), local = pop)
sim <- new.env(parent = globalenv())
sim_source <- file.path(paper_dir, "simulation", "R", "sim_functions.R")
sys.source(sim_source, envir = sim)

main <- function() {
  cat("Command: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript checks/sample_ctrl.R\n")
  cat("UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
  cat("Current R:", R.version.string, "| lavaan:", as.character(packageVersion("lavaan")), "\n")
  cat("Paper HEAD:", system2("git",c("-C",shQuote(paper_dir),"rev-parse","HEAD"),stdout=TRUE), "\n")
  cat("Current RNGkind:", paste(RNGkind(), collapse = ", "), "\n")
  cat("Threads:", paste(paste(pop$thread_vars, Sys.getenv(pop$thread_vars), sep = "="), collapse = " "), "\n")
  cat("MD5:", paste(paste(basename(c(script,pop$script,sim_source)),
                               tools::md5sum(c(script,pop$script,sim_source)),sep="="),collapse=" "), "\n")
  cat("Scope: three continuous-data, first-replication, default-start, one-factor ML fits only.\n")
  cat("Original options retained: std.lv=TRUE, se=none, test=standard, h1/baseline=TRUE, sample.cov.rescale=TRUE.\n")
  cat("The original one-factor fit is the first fit after generation; no preceding estimator or random-start draws are replayed.\n")
  cat("Generator: set.seed(cell$seed + rep); matrix(rnorm(N*p),N,p) %*% chol(pop_sigma(cell)).\n")
  cat("Archive files do not retain generated data, starting matrices, or RNGkind. Exact data identity is unverified.\n")
  cat("Nonzero archival differences are reported explicitly. Absolute ML discrepancy agreement tolerance=1e-8; exact reproduction is separate.\n")
  cat("Derivative coordinates: raw observed scales, unit latent variance; common global sign orientation.\n")
  cat("Stationarity threshold=", pop$grad_tol, "; Hessian threshold=",pop$eig_tol,
      "; uniqueness boundary threshold=",pop$psi_tol,"\n",sep="")
  current_cells <- sim$make_cells()
  rows <- list()
  for (cell_id in c(232L,664L,808L)) {
    archive <- file.path(paper_dir,"results","raw",sprintf("cell_%04d.rds",cell_id))
    old <- readRDS(archive); cell <- old$cell; rep_id <- 1L
    current <- current_cells[current_cells$cell_id == cell_id, ]
    keys <- c("cell_id","pop_id","N","format","seed","mA","mB","lamA","lamB","phi")
    for (key in keys) if (!identical(as.character(cell[[key]]),as.character(current[[key]])))
      stop("Current design differs from archived cell in ",key)
    expected_base <- cell$pop_id*100000 + match(cell$N,sim$NS)*10000
    stopifnot(cell$seed == expected_base,cell$format == "cont",cell$mA == 4,cell$mB == 4)
    seed <- cell$seed+rep_id
    target <- old$fits[old$fits$rep == rep_id & old$fits$est == "ML" & old$fits$model == "1F",]
    stopifnot(nrow(target)==1L,isTRUE(target$converged),is.finite(target$fx))
    set.seed(seed)
    Y <- sim$gen_latent(cell$N,sim$pop_sigma(cell))
    S <- cov(Y); after_generator <- .Random.seed
    warnings <- character()
    fit <- withCallingHandlers(sim$fit_ml_cov(sim$model_1f(ncol(Y)),S,cell$N),
                               warning = function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning") })
    rng_unchanged <- identical(after_generator,.Random.seed)
    effective <- lavInspect(fit,"sampstat")$cov
    rescale_error <- max(abs(effective-S*(cell$N-1)/cell$N))
    if (rescale_error > 1e-12) stop("Effective sample covariance differs from original ML rescaling.")
    est <- lavInspect(fit,"est"); starts <- lavInspect(fit,"start")
    if (est$psi[1,1] <= 0) stop("Nonpositive latent variance in selected empirical fit.")
    l <- est$lambda[,1]*sqrt(est$psi[1,1]); psi <- diag(est$theta)
    cfg <- data.frame(lambda_A=cell$lamA,lambda_B=cell$lamB,phi=cell$phi)
    ans <- pop$diagnose(l,psi,effective,cfg)
    oriented_l <- ans$common_sign_flip*l
    full <- pop$full_ml(oriented_l,psi,effective)
    sd_hat <- sqrt(diag(full$Sigma))
    l_std <- oriented_l/sd_hat; psi_std <- psi/diag(full$Sigma)
    block_std <- c(mean(l_std[1:4]),mean(l_std[5:8]))
    old_block <- c(target$lA,target$lB)
    old_orientation <- if (sum(old_block)<0) -1 else 1
    old_block <- old_orientation*old_block
    new_fx <- as.numeric(lavInspect(fit,"optim")$fx)
    discrepancy_error <- abs(full$F-2*new_fx)
    if (discrepancy_error > 1e-8) stop("Direct full ML discrepancy differs from lavaan objective.")
    archival_difference <- abs(full$F-2*target$fx)
    archived_fx_difference <- abs(new_fx-target$fx)
    exact_objective <- identical(new_fx,unname(target$fx))
    block_difference <- max(abs(block_std-old_block))
    uniqueness_difference <- abs(min(psi_std)-target$min_psi)
    exact_fields <- exact_objective && archival_difference == 0 && block_difference == 0 && uniqueness_difference == 0
    pop$check_derivatives(oriented_l,psi,effective,paste("empirical cell",cell_id))
    starts_Sigma <- starts$lambda %*% starts$psi %*% t(starts$lambda)+starts$theta
    row <- data.frame(cell_id=cell_id,rep=rep_id,condition=cell$cond,N=cell$N,
                      phi=cell$phi,lambda_A=cell$lamA,lambda_B=cell$lamB,
                      seed_base=cell$seed,replication_seed=seed,
                      archive_R=old$R,archive_lavaan=old$lavaan,current_R=R.version.string,
                      current_lavaan=as.character(packageVersion("lavaan")),
                      archived_converged=target$converged,current_lavaan_converged=lavInspect(fit,"converged"),
                      archived_F=2*target$fx,current_F=full$F,absolute_archived_F_difference=archival_difference,
                      absolute_archived_fx_difference=archived_fx_difference,
                      exact_archived_objective_match=exact_objective,
                      exact_all_compared_fields_match=exact_fields,
                      objective_comparison_tolerance=1e-8,
                      archived_objective_within_tolerance=archival_difference<=1e-8,
                      direct_lavaan_objective_error=discrepancy_error,
                      max_full_gradient=ans$max_full_gradient,min_full_hessian_eigenvalue=ans$min_full_hessian_eigenvalue,
                      stationary=ans$stationary,numerical_classification=ans$numerical_classification,
                      min_uniqueness_raw=min(psi),max_uniqueness_raw=max(psi),
                      min_uniqueness_standardized=min(psi_std),max_uniqueness_standardized=max(psi_std),
                      archived_min_uniqueness_standardized=target$min_psi,
                      absolute_archived_min_uniqueness_difference=uniqueness_difference,
                      common_sign_flip=ans$common_sign_flip,
                      loading_A_raw=ans$loading_A,loading_B_raw=ans$loading_B,
                      loading_A_standardized=block_std[1],loading_B_standardized=block_std[2],
                      archived_loading_A_oriented=old_block[1],archived_loading_B_oriented=old_block[2],
                      max_archived_block_loading_difference=block_difference,
                      signed_share_A=pop$signed_share(block_std[1],block_std[2],cell$lamA,cell$lamB,cell$phi),
                      fitted_sample_variance_error=max(abs(diag(full$Sigma)-diag(effective))),
                      rescaling_error=rescale_error,fit_preserved_rng_state=rng_unchanged,
                      generated_sum=sum(Y),generated_sum_squares=sum(Y^2),
                      effective_sample_covariance=pop$pack(as.vector(effective)),
                      current_default_start_covariance=pop$pack(as.vector(starts_Sigma)),
                      loadings_unit_latent_variance=pop$pack(oriented_l),uniquenesses_raw=pop$pack(psi),
                      warnings=paste(unique(warnings),collapse=" | "),stringsAsFactors=FALSE)
    rows[[length(rows)+1]] <- row
    cat(sprintf("\nCell %d rep %d seed %d: %s phi=%.2f N=%d\n",cell_id,rep_id,seed,cell$cond,cell$phi,cell$N))
    cat("Archive:",old$R,"| lavaan",old$lavaan,"| MD5",unname(tools::md5sum(archive)),"\n")
    cat(sprintf("Archived F=%.17g; current F=%.17g; absolute difference=%.3g; exact fx match=%s; within1e-8=%s\n",
                2*target$fx,full$F,archival_difference,exact_objective,archival_difference<=1e-8))
    if (!exact_objective) cat("Exact archival objective reproduction failed; the absolute difference is retained above.\n")
    if (!exact_fields) cat(sprintf("Exact reproduction of all compared floating-point fields failed: block-loading difference %.3g; minimum-uniqueness difference %.3g.\n",
                                  block_difference,uniqueness_difference))
    cat(sprintf("Current lavaan converged=%s; max full gradient=%.3g; min full Hessian eigenvalue=%.6g; %s\n",
                row$current_lavaan_converged,row$max_full_gradient,row$min_full_hessian_eigenvalue,row$numerical_classification))
    cat(sprintf("Raw uniqueness=[%.8f,%.8f]; standardized uniqueness=[%.8f,%.8f]; oriented standardized block means=[%.8f,%.8f]\n",
                min(psi),max(psi),min(psi_std),max(psi_std),block_std[1],block_std[2]))
    cat(sprintf("Direct/lavaan objective error=%.3g; sample rescaling error=%.3g; estimator preserved generator RNG state=%s\n",
                discrepancy_error,rescale_error,rng_unchanged))
    if (length(warnings)) cat("Warnings:",paste(unique(warnings),collapse=" | "),"\n")
  }
  result <- do.call(rbind,rows)
  output <- file.path(paper_dir,"results","summary","sample_ctrl_stage4.csv")
  write.csv(result,output,row.names=FALSE,na="NA")
  cat("\nCSV:",output,"\n")
  cat("Completed exactly three cross-version empirical refits; no simulation source or archived result was changed.\n")
  cat("These selected converged archive rows do not estimate failure rates, establish empirical global minima, or validate all archived classifications.\n")
  cat("PASS: seed/design checks, covariance rescaling, direct objective checks, and independent derivative checks.\n")
  cat("\nSession info:\n");print(sessionInfo())
}
log <- file.path(paper_dir,"results","logs","sample_ctrl_stage4.log")
sink(log,split=TRUE)
tryCatch(main(),error=function(e){cat("\nFAIL:",conditionMessage(e),"\n");stop(e)},finally=sink())
