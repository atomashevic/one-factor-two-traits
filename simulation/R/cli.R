## cli.R: minimal command-line parsing (no extra packages).
## "--key value" sets key; a bare "--flag" sets flag to "TRUE".
parse_args <- function(argv, defaults = list()) {
  out <- defaults
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a)
    key <- gsub("-", "_", sub("^--", "", a))
    if (i < length(argv) && !startsWith(argv[i + 1], "--")) {
      out[[key]] <- argv[i + 1]; i <- i + 2
    } else {
      out[[key]] <- "TRUE"; i <- i + 1
    }
  }
  out
}

## "1-10,15,20-22" -> c(1:10, 15, 20:22)
parse_ids <- function(x) {
  parts <- strsplit(x, ",")[[1]]
  sort(unique(unlist(lapply(parts, function(p) {
    r <- as.integer(strsplit(p, "-")[[1]])
    if (length(r) == 2) r[1]:r[2] else r
  }))))
}

flag <- function(args, key) isTRUE(as.logical(args[[key]]))
