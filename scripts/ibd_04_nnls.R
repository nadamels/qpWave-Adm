#!/usr/bin/env Rscript
# =============================================================================
# Step 4: supervised ancestry modelling from IBD painting profiles.
#
#   target (group or single individual)  =  sum_k  w_k * source_k          w_k >= 0, sum w_k = 1
#   solved on the painting profiles over the donor groups by non-negative least squares (limSolve::lsei),
#   SEs = weighted block jackknife leaving out each chromosome in turn (block weight = chromosome genetic length).
#
# Models come from
#   MODELS     config/ibd_models.tsv      columns: target, sources (comma separated) [, model_id]      (explicit models)
#   ENUMERATE  1 -> every combination of 1..MAX_SOURCES groups of SOURCE_POOL (config/ibd_source_pool.tsv) for every
#              target in TARGETS (default: all modern_target groups), ranked by residual sum of squares.
#   A source / target name is a group from ibd_groups.tsv (ancient cluster or modern group) or, for a target, an IID.
#
# Donors of a model = all donor groups EXCEPT the sources and the target's own group (EXCLUDE_SOURCES=1, default), the
# usual SOURCEFIND-style choice: a source would otherwise "copy" mostly from itself, which the target cannot do.
# Profiles are re-normalised to sum to 1 over each model's donor set.
#
# Env: PAINT_RDS OUTDIR [MODELS ENUMERATE=0 SOURCE_POOL TARGETS MAX_SOURCES=3 MIN_WEIGHT=0.02 TOP_N=5 SIZE_NORM=1
#      PROFILE_NORM=sum|none EXCLUDE_SOURCES=1 METHOD=lsei NCORES=1 JACKKNIFE=1]
# =============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "ibd_common.R"))

PAINT_RDS <- getenv("PAINT_RDS"); OUTDIR <- getenv("OUTDIR")
MODELS <- getenv("MODELS", ""); ENUMERATE <- getenv("ENUMERATE", "0") == "1"
SOURCE_POOL <- getenv("SOURCE_POOL", ""); TARGETS <- getenv("TARGETS", "")
MAX_SOURCES <- as.integer(getenv("MAX_SOURCES", 3)); MIN_WEIGHT <- as.numeric(getenv("MIN_WEIGHT", 0.02))
TOP_N <- as.integer(getenv("TOP_N", 5)); SIZE_NORM <- getenv("SIZE_NORM", "1") == "1"
PROFILE_NORM <- getenv("PROFILE_NORM", "sum")                  # sum = each individual profile sums to 1 (default); none = raw mean-IBD vectors
stopifnot(PROFILE_NORM %in% c("sum", "none"))
EXCLUDE_SOURCES <- getenv("EXCLUDE_SOURCES", "1") == "1"; METHOD <- getenv("METHOD", "lsei")
NCORES <- as.integer(getenv("NCORES", 1)); JACKKNIFE <- getenv("JACKKNIFE", "1") == "1"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

PT <- readRDS(PAINT_RDS)
S <- PT$S; Nn <- PT$N; ids <- PT$ids; grp_of <- PT$grp_of; groups <- PT$groups; donors <- PT$donors
chr_w <- as.numeric(PT$chr_cM); nC <- dim(S)[3]
Stot <- apply(S, c(1, 2), sum)
msg("limSolve: ", requireNamespace("limSolve", quietly = TRUE), "  donor groups: ", length(donors), "  chromosomes: ", nC)

members_of <- function(x) if (x %in% groups) which(grp_of == x) else if (x %in% ids) match(x, ids) else stop("unknown group / IID: ", x)

# mean painting profile of a set of individuals over `don`; `left_out` = chromosome index to leave out (NULL = none)
profile <- function(rows, don, left_out = NULL) {
  M <- if (is.null(left_out)) Stot[rows, don, drop = FALSE] else Stot[rows, don, drop = FALSE] - S[rows, don, left_out, drop = FALSE][, , 1]
  if (!is.matrix(M)) M <- matrix(M, nrow = length(rows))
  if (SIZE_NORM) M <- M / Nn[rows, don, drop = FALSE]
  s <- rowSums(M, na.rm = TRUE); keep <- s > 0 & is.finite(s)
  if (!any(keep)) return(rep(NA_real_, length(don)))
  if (PROFILE_NORM == "none") return(colMeans(M[keep, , drop = FALSE]))
  colMeans(M[keep, , drop = FALSE] / s[keep])
}

fit_once <- function(tmem, smem, don, left_out = NULL) {
  y <- profile(tmem, don, left_out)
  X <- sapply(smem, function(m) profile(m, don, left_out))
  if (anyNA(y) || anyNA(X)) return(NULL)
  w <- solve_nnls(X, y, METHOD)
  if (anyNA(w)) return(NULL)
  list(w = w, rss = sum((y - X %*% w)^2), rel_rss = sum((y - X %*% w)^2) / sum(y^2), res = y - X %*% w)
}

# weighted block jackknife (Busing et al. 1999; same estimator as ADMIXTOOLS)
jackknife <- function(tmem, smem, don, est) {
  th <- vapply(seq_len(nC), function(c) { f <- fit_once(tmem, smem, don, c); if (is.null(f)) rep(NA_real_, length(smem)) else f$w }, numeric(length(smem)))
  th <- matrix(th, nrow = length(smem))
  ok <- !is.na(th[1, ]); m <- chr_w[ok]; n <- sum(m); g <- length(m); h <- n / m
  th <- th[, ok, drop = FALSE]
  tj <- g * est - as.numeric(th %*% ((n - m) / n))              # jackknife mean estimate
  tau <- sapply(seq_len(g), function(j) h[j] * est - (h[j] - 1) * th[, j])
  tau <- matrix(tau, nrow = length(smem))
  se <- sqrt(rowSums((tau - tj)^2 / rep(h - 1, each = length(smem))) / g)
  se
}

fit_model <- function(target, sources, jk = JACKKNIFE) {
  tmem <- members_of(target); smem <- lapply(sources, members_of)
  tgrp <- if (target %in% groups) target else grp_of[match(target, ids)]
  don <- if (EXCLUDE_SOURCES) setdiff(donors, c(sources, tgrp)) else setdiff(donors, tgrp)
  don <- don[!is.na(don)]
  f <- fit_once(tmem, smem, don)
  if (is.null(f)) return(NULL)
  se <- if (jk) tryCatch(jackknife(tmem, smem, don, f$w), error = function(e) rep(NA_real_, length(sources))) else rep(NA_real_, length(sources))
  worst <- order(-abs(f$res))[1:min(3, length(don))]
  list(target = target, sources = sources, w = f$w, se = se, rss = f$rss, rel_rss = f$rel_rss, n_donors = length(don), n_target = length(tmem),
       worst = paste(sprintf("%s(%+.3f)", don[worst], f$res[worst]), collapse = ";"))
}

fmt <- function(r, id) data.table(model_id = id, target = r$target, n_target = r$n_target, n_sources = length(r$sources),
  sources = paste(r$sources, collapse = ","), weights = paste(sprintf("%.3f", r$w), collapse = ","), se = paste(sprintf("%.3f", r$se), collapse = ","),
  min_weight = min(r$w), rss = r$rss, rel_rss = r$rel_rss, n_donors = r$n_donors, worst_residual_donors = r$worst)
long <- function(r, id) data.table(model_id = id, target = r$target, source = r$sources, weight = r$w, se = r$se,
  z = ifelse(is.na(r$se) | r$se == 0, NA_real_, r$w / r$se))

res <- list(); lres <- list()

# ---- explicit models -----------------------------------------------------------------------------------------
if (nzchar(MODELS) && file.exists(MODELS)) {
  M <- read_tsv_nocomment(MODELS)
  if (!"model_id" %in% names(M)) M[, model_id := paste0("M", seq_len(.N))]
  msg("explicit models: ", nrow(M))
  for (i in seq_len(nrow(M))) {
    src <- trimws(strsplit(M$sources[i], ",")[[1]]); r <- tryCatch(fit_model(trimws(M$target[i]), src), error = function(e) { msg("model ", M$model_id[i], " failed: ", conditionMessage(e)); NULL })
    if (!is.null(r)) { res[[length(res) + 1]] <- cbind(fmt(r, M$model_id[i]), set = "explicit"); lres[[length(lres) + 1]] <- long(r, M$model_id[i]) }
  }
}

# ---- enumeration over a source pool --------------------------------------------------------------------------
if (ENUMERATE) {
  stopifnot(nzchar(SOURCE_POOL), file.exists(SOURCE_POOL))
  pl <- trimws(sub("#.*", "", readLines(SOURCE_POOL))); pl <- pl[nzchar(pl)]
  pool <- unique(unlist(lapply(pl, function(x) if (startsWith(x, "~")) grep(substring(x, 2), groups, value = TRUE) else x)))
  miss <- setdiff(pool, groups); if (length(miss)) { msg("WARNING: pool entries not found among groups (dropped): ", paste(miss, collapse = ", ")); pool <- intersect(pool, groups) }
  pool <- pool[PT$gsize[pool] >= 1]
  tg <- if (nzchar(TARGETS)) trimws(strsplit(TARGETS, ",")[[1]]) else PT$types[type == "modern_target", group]
  msg("enumeration: ", length(tg), " target(s), pool of ", length(pool), ", up to ", MAX_SOURCES, " sources")
  for (t in tg) {
    cand <- setdiff(pool, t); allm <- list()
    for (k in seq_len(min(MAX_SOURCES, length(cand)))) {
      cmb <- combn(cand, k, simplify = FALSE)
      out <- parallel::mclapply(cmb, function(src) tryCatch(fit_model(t, src, jk = FALSE), error = function(e) NULL), mc.cores = NCORES)
      out <- Filter(Negate(is.null), out)
      d <- rbindlist(lapply(out, fmt, id = NA_character_)); d[, set := "enumerated"]; d[, rank_in_size := frank(rss, ties.method = "first")]
      d[, feasible := min_weight >= MIN_WEIGHT]
      allm[[k]] <- d
      # SEs for the best TOP_N models of this size (and best TOP_N feasible)
      sel <- unique(c(order(d$rss)[seq_len(min(TOP_N, nrow(d)))], which(d$feasible)[order(d$rss[d$feasible])][seq_len(min(TOP_N, sum(d$feasible)))]))
      for (j in sel) { r <- fit_model(t, strsplit(d$sources[j], ",")[[1]], jk = JACKKNIFE); d[j, `:=`(se = fmt(r, NA)$se)]
        lres[[length(lres) + 1]] <- long(r, sprintf("%s|%s", t, d$sources[j])) }
      msg(t, ": ", k, "-source models ", nrow(d), ", feasible (all weights >= ", MIN_WEIGHT, ") ", sum(d$feasible), ", best RSS ", signif(min(d$rss), 3))
    }
    all_t <- rbindlist(allm, fill = TRUE)
    all_t[, model_id := sprintf("%s|%s", target, sources)]
    res[[length(res) + 1]] <- all_t
  }
}

if (!length(res)) stop("no models fitted: provide MODELS and/or ENUMERATE=1")
R <- rbindlist(res, fill = TRUE)
R[, rss := signif(rss, 5)]
setorder(R, target, n_sources, rss)
fwrite(R, file.path(OUTDIR, "ibd_nnls_models.tsv"), sep = "\t")
L <- unique(rbindlist(lres), by = c("model_id", "source"))
fwrite(L, file.path(OUTDIR, "ibd_nnls_weights_long.tsv"), sep = "\t")
if (ENUMERATE) {
  best <- R[set == "enumerated" & feasible == TRUE][, head(.SD, TOP_N), by = .(target, n_sources)]
  fwrite(best, file.path(OUTDIR, "ibd_nnls_best_feasible_by_size.tsv"), sep = "\t")
}
msg("wrote ", file.path(OUTDIR, "ibd_nnls_models.tsv"), " (", nrow(R), " models) and ibd_nnls_weights_long.tsv")

# ---- plot: weights of explicit models / best feasible models ----------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE) && nrow(L)) {
  suppressPackageStartupMessages(library(ggplot2))
  ids_plot <- if (ENUMERATE) { b <- R[set == "enumerated" & feasible == TRUE][order(target, n_sources, rss)][, .SD[1], by = target]; b$model_id } else R[set == "explicit"]$model_id
  D <- L[model_id %in% ids_plot]
  if (nrow(D)) {
    p <- ggplot(D, aes(x = model_id, y = weight, fill = source)) + geom_col(width = 0.7) + coord_flip() + labs(x = NULL, y = "ancestry proportion (IBD-painting NNLS)", fill = "source") +
      theme_minimal(base_size = 10) + theme(legend.position = "bottom")
    ggsave(file.path(OUTDIR, "ibd_nnls_models.png"), p, width = 9, height = max(3, 0.45 * length(unique(D$model_id)) + 2), dpi = 200)
    ggsave(file.path(OUTDIR, "ibd_nnls_models.pdf"), p, width = 9, height = max(3, 0.45 * length(unique(D$model_id)) + 2))
  }
}
