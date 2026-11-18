## Regenerate replications of archived cells with the committed simulation code
## and compare them with the archive. Two modes, so that generation can run on
## a machine without the archive.
##   Rscript simulation/posthoc/reproduce_archive.R --mode generate \
##     --cells 1,34 --reps 5 --out build/report-revision/reproduce_gen.rds
##   Rscript simulation/posthoc/reproduce_archive.R --mode compare \
##     --gen build/report-revision/reproduce_gen.rds [--raw results/raw]
## The committed extract_fit() stores absolute block means for one-factor fits;
## the archive holds signed block means. The comparison therefore uses |lA|,
## |lB| and reports signed agreement separately.
argv <- commandArgs(TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), '../..'))
opt <- list(mode = 'generate', cells = '1', reps = '5', out = '', gen = '',
            raw = file.path(root, 'results/raw'))
if (length(argv) %% 2L) stop('Options require a value')
if (length(argv)) for (i in seq(1L, length(argv), by = 2L)) {
  key <- sub('^--', '', argv[i])
  if (!key %in% names(opt)) stop('Unknown option: ', argv[i])
  opt[[key]] <- argv[i + 1L]
}
if (opt$mode == 'generate') {
  source(file.path(root, 'simulation/R/sim_functions.R'))
  stopifnot(nzchar(opt$out), !file.exists(opt$out))
  cells <- make_cells()
  ids <- as.integer(strsplit(opt$cells, ',')[[1]])
  out <- lapply(ids, function(id) {
    r <- run_cell(cells[cells$cell_id == id, ], reps = as.integer(opt$reps), cores = 1)
    list(cell = r$cell, fits = r$fits, ms = r$ms)
  })
  attr(out, 'runtime') <- list(R = R.version.string,
                               lavaan = as.character(packageVersion('lavaan')),
                               code_md5 = unname(tools::md5sum(file.path(root, 'simulation/R/sim_functions.R'))))
  saveRDS(out, opt$out)
  cat('Generated', length(out), 'cells;', R.version.string, '; lavaan',
      as.character(packageVersion('lavaan')), '\n')
} else if (opt$mode == 'compare') {
  gen <- readRDS(opt$gen)
  print(attr(gen, 'runtime'))
  for (g in gen) {
    a <- readRDS(file.path(opt$raw, sprintf('cell_%04d.rds', g$cell$cell_id)))
    key <- c('est', 'model', 'rep')
    m <- merge(g$fits, a$fits, by = key, suffixes = c('.new', '.arc'))
    stopifnot(nrow(m) == nrow(g$fits))
    d <- function(x, y) max(abs(x - y), na.rm = TRUE)
    f1 <- m[m$model == '1F', ]
    cat(sprintf('cell %d (%s, phi=%s, N=%s, %s; archive R "%s", lavaan %s): %d fit rows\n',
                g$cell$cell_id, g$cell$cond, g$cell$phi, g$cell$N, g$cell$format,
                a$R, a$lavaan, nrow(m)))
    cat(sprintf('  max |diff|: fx %.3g, cfi %.3g, rmsea %.3g, srmr %.3g, min_psi %.3g\n',
                d(m$fx.new, m$fx.arc), d(m$cfi.new, m$cfi.arc), d(m$rmsea.new, m$rmsea.arc),
                d(m$srmr.new, m$srmr.arc), d(m$min_psi.new, m$min_psi.arc)))
    cat(sprintf('  1F rows: max ||lA|-|lA_arc|| %.3g, ||lB|-|lB_arc|| %.3g; archived negative lA %d, lB %d; regenerated negative lA %d, lB %d\n',
                d(abs(f1$lA.new), abs(f1$lA.arc)), d(abs(f1$lB.new), abs(f1$lB.arc)),
                sum(f1$lA.arc < 0), sum(f1$lB.arc < 0), sum(f1$lA.new < 0), sum(f1$lB.new < 0)))
    cat(sprintf('  1F rows: max |shareA - shareA_arc| %.3g (differs where the archive has a negative block mean)\n',
                d(f1$shareA.new, f1$shareA.arc)))
    if (!is.null(g$ms)) {
      mm <- merge(g$ms, a$ms, by = c('est', 'rep'), suffixes = c('.new', '.arc'))
      cat(sprintf('  multistart rows %d: max |best_fx diff| %.3g; n_distinct equal %d/%d; spurious_weak equal %d/%d; def_on_weak equal %d/%d\n',
                  nrow(mm), d(mm$best_fx.new, mm$best_fx.arc),
                  sum(mm$n_distinct.new == mm$n_distinct.arc), nrow(mm),
                  sum(mm$spurious_weak.new == mm$spurious_weak.arc, na.rm = TRUE), nrow(mm),
                  sum(mm$def_on_weak.new == mm$def_on_weak.arc, na.rm = TRUE), nrow(mm)))
    }
  }
} else stop('Unknown mode')
