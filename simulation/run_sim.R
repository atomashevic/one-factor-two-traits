## run_sim.R
## Runs cells of the Section 7.3 design and writes one file per cell,
## <out>/cell_<id>.rds. Finished cells are skipped, so a run can be resumed.
##
## Usage:
##   Rscript run_sim.R --list                       # print the design and exit
##   Rscript run_sim.R --cells 1-50 --cores 8       # run some cells
##   Rscript run_sim.R --task 3 --ntasks 100        # cluster array: task 3 of 100
##   Rscript run_sim.R --pilot --cores 2            # 72 small cells + cost estimate
##
## Options:
##   --out DIR        output directory (default results; pilot: pilot)
##   --cores K        cores for replications within a cell (mclapply)
##   --cells IDS      cell ids, e.g. "1-10,25"
##   --task i --ntasks n   run every cell with (cell_id - 1) %% n == i - 1
##   --reps_scale x   multiply the replication counts (e.g. 0.1 for a dry run)
##   --nstart k       random starts for the multi-start arm (default 20)
##   --force          rerun cells whose output exists

.here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
source(file.path(.here, "R", "sim_functions.R"))
source(file.path(.here, "R", "cli.R"))

args <- parse_args(commandArgs(TRUE), list(out = "", cores = "1", cells = "", task = "",
                                            ntasks = "", reps_scale = "1",
                                            nstart = as.character(NSTART)))
cells <- make_cells()

if (flag(args, "list")) {
  cat(sprintf("%d populations, %d cells\n", length(unique(cells$pop_id)), nrow(cells)))
  print(aggregate(cell_id ~ design + format, data = cells, FUN = length))
  cat(sprintf("replications: %d continuous, %d ordinal; total %d\n",
              sum(cells$reps[cells$format == "cont"]), sum(cells$reps[cells$format != "cont"]),
              sum(cells$reps)))
  cat(sprintf("cells with ULS arm: %d; cells with multi-start arm: %d\n",
              sum(cells$uls), sum(cells$ms)))
  pops <- make_populations()
  for (cd in unique(pops$cond)) cat(sprintf("  %-10s phi = %s\n", cd,
                                            paste(pops$phi[pops$cond == cd], collapse = ", ")))
  quit(save = "no")
}

pilot <- flag(args, "pilot")
out <- if (nzchar(args$out)) args$out else if (pilot) "pilot" else "results"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
cores <- as.integer(args$cores)
nstart <- as.integer(args$nstart)

if (pilot) {
  ## One population per design, every N and format, few replications.
  ## phi = .9 is included so that the tilt criterion of Question 2 can be formed.
  pick <- (cells$cond == "A_m8_l70" & cells$phi %in% c(0.3, 0.9)) |
          (cells$cond == "B_l60_l80" & cells$phi %in% c(0.01, 0.9)) |
          (cells$cond == "C_9_7_l75" & cells$phi %in% c(0.653, 0.9))
  run <- cells[pick, ]
  run$reps <- ifelse(run$format == "cont", 20, 10)
} else {
  run <- cells
  if (nzchar(args$cells)) run <- run[run$cell_id %in% parse_ids(args$cells), ]
  if (nzchar(args$task)) {
    i <- as.integer(args$task); n <- as.integer(args$ntasks)
    run <- run[(run$cell_id - 1) %% n == i - 1, ]
  }
  run$reps <- pmax(1, round(run$reps * as.numeric(args$reps_scale)))
}

for (k in seq_len(nrow(run))) {
  cell <- run[k, ]
  f <- file.path(out, sprintf("cell_%04d.rds", cell$cell_id))
  if (file.exists(f) && !flag(args, "force")) next
  res <- run_cell(cell, reps = cell$reps, cores = cores, nstart = nstart)
  saveRDS(res, f)
  message(sprintf("cell %4d  %-10s phi=%.3f N=%4d %-8s reps=%4d  %.1f s wall, %.2f s/rep",
                  cell$cell_id, cell$cond, cell$phi, cell$N, cell$format, cell$reps,
                  res$wall, mean(res$rep_time)))
}

if (pilot) {
  ## Cost estimate: mean seconds per replication in each pilot cell, applied to
  ## every full-design cell with the same design, format and N (Design A uses
  ## the 16-item condition, so its estimate is conservative).
  done <- lapply(list.files(out, "^cell_.*\\.rds$", full.names = TRUE), readRDS)
  pt <- do.call(rbind, lapply(done, function(r) data.frame(
    design = r$cell$design, format = r$cell$format, N = r$cell$N,
    sec_per_rep = mean(r$rep_time), failed = r$n_failed)))
  pt <- aggregate(cbind(sec_per_rep, failed) ~ design + format + N, data = pt, FUN = mean)
  full <- merge(cells, pt, by = c("design", "format", "N"))
  full$cpu_sec <- full$reps * full$sec_per_rep
  cat("\nSeconds per replication in the pilot:\n")
  print(pt[order(pt$design, pt$format, pt$N), ], row.names = FALSE)
  cat("\nEstimated CPU hours for the full design:\n")
  print(round(tapply(full$cpu_sec, list(full$design, full$format), sum) / 3600, 1))
  cat(sprintf("Total: %.0f CPU hours\n", sum(full$cpu_sec) / 3600))
}
