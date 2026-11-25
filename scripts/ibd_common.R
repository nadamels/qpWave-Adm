# Shared helpers for the IBD-based clustering / ancestry-modelling scripts (ibd_0*.R). Sourced, not run.
suppressPackageStartupMessages(library(data.table))

getenv <- function(k, d = NULL) {
  v <- Sys.getenv(k, unset = NA)
  if (is.na(v) || v == "") {
    if (is.null(d)) stop("Environment variable ", k, " is required")
    d
  } else v
}
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# Chromosome labels -> integer (1..22, X = 23)
chr_to_int <- function(x) {
  x <- sub("^chr", "", as.character(x), ignore.case = TRUE)
  x[toupper(x) == "X"] <- "23"
  suppressWarnings(as.integer(x))
}

# Approximate sex-averaged genetic length (cM) of chr1-22; ONLY a fallback for when the IBD table has no cM positions
# (used to turn cM totals into a "fraction of genome" and as jackknife block weights).
FALLBACK_CHR_CM <- c(`1` = 286, `2` = 269, `3` = 223, `4` = 214, `5` = 204, `6` = 192, `7` = 187, `8` = 168, `9` = 166,
                     `10` = 181, `11` = 158, `12` = 175, `13` = 126, `14` = 117, `15` = 121, `16` = 127, `17` = 129,
                     `18` = 121, `19` = 108, `20` = 108, `21` = 60, `22` = 62)

# Leiden community detection on a weighted igraph graph: leidenAlg (as in the methods) if installed, else igraph's own.
leiden_membership <- function(g, resolution, seed = 1) {
  set.seed(seed)
  if (requireNamespace("leidenAlg", quietly = TRUE)) {
    res <- leidenAlg::leiden.community(g, resolution = resolution)
    as.integer(res$membership)
  } else {
    as.integer(igraph::membership(igraph::cluster_leiden(g, objective_function = "modularity",
                                                         weights = igraph::E(g)$weight,
                                                         resolution_parameter = resolution)))
  }
}

# Non-negative least squares with sum-to-one: min ||y - X b||^2  s.t.  b >= 0, sum(b) = 1
# method "lsei": limSolve::lsei (the package named in the methods); falls back to quadprog if limSolve is missing.
# method "nnls": plain limSolve::nnls, weights then rescaled to sum to one.
solve_nnls <- function(X, y, method = "lsei") {
  k <- ncol(X)
  if (k == 1) return(1)
  if (method == "nnls" && requireNamespace("limSolve", quietly = TRUE)) {
    b <- limSolve::nnls(X, y)$X
    return(if (sum(b) > 0) b / sum(b) else rep(NA_real_, k))
  }
  if (requireNamespace("limSolve", quietly = TRUE)) {
    r <- tryCatch(limSolve::lsei(A = X, B = y, E = matrix(1, 1, k), F = 1, G = diag(k), H = rep(0, k), type = 2, verbose = FALSE),
                  error = function(e) NULL)
    if (is.null(r) || r$IsError) return(rep(NA_real_, k))
    b <- as.numeric(r$X); names(b) <- NULL
    return(pmax(b, 0))
  }
  if (!requireNamespace("quadprog", quietly = TRUE)) stop("Neither limSolve nor quadprog is installed")
  D <- 2 * crossprod(X); diag(D) <- diag(D) + 1e-10 * max(diag(D), 1e-12)
  r <- tryCatch(quadprog::solve.QP(D, 2 * crossprod(X, y), cbind(rep(1, k), diag(k)), c(1, rep(0, k)), meq = 1),
                error = function(e) NULL)
  if (is.null(r)) rep(NA_real_, k) else pmax(r$solution, 0)
}

# Read a TSV whose comment lines start with '#' (the header line is the first non-comment line)
read_tsv_nocomment <- function(path) {
  l <- readLines(path); l <- l[nzchar(trimws(l)) & !grepl("^\\s*#", l)]
  fread(text = paste0(paste(l, collapse = "\n"), "\n"), sep = "\t", header = TRUE, colClasses = "character")
}
