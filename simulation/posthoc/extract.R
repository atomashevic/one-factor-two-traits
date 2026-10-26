## Extract archived one-factor rows without changing signs or convergence flags.
## Rscript simulation/posthoc/extract.R [--input DIR] [--output FILE] [--limit N]
argv <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), '../..'))
opt <- list(input = Sys.getenv('SIM_RESULTS', file.path(root, 'results/raw')),
            output = file.path(root, 'results/fits1f.rds'), limit = Inf)
if (length(argv) %% 2L) stop('Options require a value')
if (length(argv)) for (i in seq(1L, length(argv), by = 2L)) {
  key <- sub('^--', '', argv[i])
  if (!key %in% names(opt)) stop('Unknown option: ', argv[i])
  opt[[key]] <- argv[i + 1L]
}
limit <- as.numeric(opt$limit)
if (is.na(limit) || limit < 1 || (is.finite(limit) && limit != floor(limit)))
  stop('--limit must be a positive integer')
if (file.exists(opt$output)) stop('Output exists; preserve it and choose a new --output')
fs <- sort(list.files(opt$input, '^cell_[0-9]{4}\\.rds$', full.names = TRUE))
if (!length(fs)) stop('No cell files in ', opt$input)
fs <- head(fs, min(length(fs), limit))
input_md5 <- tools::md5sum(fs)
runtime <- vector('list', length(fs))
required <- c('cell_id', 'pop_id', 'design', 'cond', 'phi', 'N', 'format',
              'est', 'model', 'rep', 'ok', 'converged', 'fx', 'lA', 'lB',
              'shareA', 'min_psi')
started <- proc.time()[['elapsed']]
parts <- lapply(seq_along(fs), function(i) {
  r <- readRDS(fs[i]); cell <- r$cell
  stopifnot(is.data.frame(cell), nrow(cell) == 1L, is.data.frame(r$fits),
            all(required %in% names(r$fits)))
  id <- as.integer(sub('cell_([0-9]+)\\.rds', '\\1', basename(fs[i])))
  stopifnot(identical(as.integer(cell$cell_id), id))
  x <- r$fits[r$fits$model == '1F', , drop = FALSE]
  if (!nrow(x)) stop('Missing one-factor rows: ', fs[i])
  for (k in c('cell_id', 'pop_id', 'design', 'cond', 'phi', 'N', 'format'))
    if (!all(x[[k]] == cell[[k]])) stop('Metadata mismatch: ', fs[i], ' ', k)
  if (anyDuplicated(x[c('est', 'rep')])) stop('Duplicate fit key: ', fs[i])
  for (e in unique(x$est))
    if (!identical(sort(as.integer(x$rep[x$est == e])), seq_len(cell$reps)))
      stop('Incomplete replication sequence: ', fs[i], ' ', e)
  for (k in c('mA', 'mB', 'lamA', 'lamB', 'phistar', 'seed', 'reps')) x[[k]] <- cell[[k]]
  for (k in c('lamA', 'lamB', 'phi')) {
    target <- c(lamA = 'refA', lamB = 'refB', phi = 'refphi')[[k]]
    x[[target]] <- if (is.null(r$ref_scores)) NA_real_ else r$ref_scores[[k]]
  }
  x$source_file <- basename(fs[i])
  runtime[[i]] <<- data.frame(source_file = basename(fs[i]),
                             R = r$R, lavaan = r$lavaan, stringsAsFactors = FALSE)
  if (i %% 200L == 0L) cat('Read', i, 'cells\n')
  x
})
fits <- do.call(rbind, parts); rownames(fits) <- NULL
stopifnot(!anyDuplicated(fits[c('cell_id', 'est', 'rep')]))
attr(fits, 'provenance') <- list(
  command = commandArgs(), R = R.version.string,
  script_md5 = unname(tools::md5sum(script)),
  input_md5 = setNames(unname(input_md5), basename(fs)),
  source_runtime = do.call(rbind, runtime),
  note = 'Signed archived lA/lB, shareA, ok and converged are unchanged; no fit certification inferred.')
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
temporary <- paste0(opt$output, '.tmp-', Sys.getpid())
saveRDS(fits, temporary)
if (!file.rename(temporary, opt$output)) stop('Could not finalize ', temporary)
cat(nrow(fits), 'one-factor fits from', length(fs), 'cells; elapsed',
    proc.time()[['elapsed']] - started, 'seconds\n')
cat('Signed block means:', sum(fits$lA < 0 | fits$lB < 0, na.rm = TRUE),
    'rows with a negative block mean\n')
cat('Output MD5:', unname(tools::md5sum(opt$output)), '\n')
print(unique(attr(fits, 'provenance')$source_runtime[c('R', 'lavaan')]))
