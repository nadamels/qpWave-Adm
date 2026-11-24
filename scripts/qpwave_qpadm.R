#!/usr/bin/env Rscript
# =============================================================================
# qpWave + qpAdm for modern Iranian sub-populations (ADMIXTOOLS 2)
#
#   0. f2 blocks from PLINK (cached)
#   1. qpWave among the target groups (are they distinguishable / how many streams?)
#   2. qpWave screen of candidate SOURCE sets (drop source sets that are not
#      distinguishable = rank-deficient, their qpAdm weights would be meaningless)
#   3. qpAdm "rotating" over every surviving source set, for every target
#      (sources not used on the LEFT go to the RIGHT, on top of a fixed outgroup set)
#   4. model selection: p >= PTHRESH and all weights in [0,1]; simplest model wins
#   5. weights +- SE (block jackknife) for the best models, leave-one-out-RIGHT
#      robustness for the best model per target
#
# All settings come from environment variables (see run_qpwave_qpadm.sh).
# =============================================================================
suppressPackageStartupMessages({
  library(admixtools)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
})

getenv <- function(k, d = NULL) {
  v <- Sys.getenv(k, unset = NA)
  if (is.na(v) || v == "") {
    if (is.null(d)) stop("Environment variable ", k, " is required")
    d
  } else v
}
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

PLINK_PREFIX  <- getenv("PLINK_PREFIX")                       # prefix of bed/bim/fam (fam: FID = population)
OUTDIR        <- getenv("OUTDIR")
GROUP_SUMMARY <- getenv("GROUP_SUMMARY", file.path(OUTDIR, "group_summary.tsv"))
F2DIR         <- getenv("F2DIR", file.path(OUTDIR, "f2_blocks"))
NCORES        <- as.integer(getenv("NCORES", 4))
PTHRESH       <- as.numeric(getenv("PTHRESH", 0.05))          # qpAdm / qpWave significance level
WEIGHT_TOL    <- as.numeric(getenv("WEIGHT_TOL", 0))          # weights must lie in [-tol, 1+tol]
MAX_SOURCES   <- as.integer(getenv("MAX_SOURCES", 3))         # largest source set tried
MIN_SOURCES   <- as.integer(getenv("MIN_SOURCES", 2))         # smallest model that may be SELECTED (1-way fits are still reported separately)
ROTATE        <- getenv("ROTATE", "1") == "1"                 # 1: pool groups not on the LEFT move to the RIGHT (rotating qpAdm)
                                                              # 0: RIGHT = fixed right_base set only (Narasimhan 2019 style)
CI_Z          <- as.numeric(getenv("CI_Z", 0))                # >0: every weight +- CI_Z*SE must lie inside [0,1] (paper: 2)
COMPETE       <- getenv("COMPETE", "0") == "1"                # model competition (needs ROTATE=0): add a competitor's extra sources to the RIGHT
COMPETE_MAX   <- as.integer(getenv("COMPETE_MAX", 300))       # at most this many feasible models per target enter the competition
MIN_WEIGHT    <- as.numeric(getenv("MIN_WEIGHT", 0.05))       # every source must contribute at least this much, else the model is
                                                              # really an (n-1)-source model plus noise
TOP_N         <- as.integer(getenv("TOP_N", 3))               # models kept per target for plots/SEs
MAXMEM_MB     <- as.integer(getenv("MAXMEM_MB", 16000))
BLGSIZE       <- getenv("BLGSIZE", "auto")                    # "auto" | Morgans (0.05) | bp if > 100 (2e6)
FORCE_F2      <- getenv("FORCE_F2", "0") == "1"
SEED          <- as.integer(getenv("SEED", 1))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
options(mc.cores = NCORES)
set.seed(SEED)

pmap_list <- function(x, f) {
  if (NCORES > 1) parallel::mclapply(x, f, mc.cores = NCORES, mc.preschedule = TRUE) else lapply(x, f)
}

# ---------------------------------------------------------------- groups ----
gs <- read_tsv(GROUP_SUMMARY, show_col_types = FALSE)
targets    <- gs$group[gs$role == "target"]
pool       <- gs$group[gs$role == "pool"]
right_base <- gs$group[gs$role == "right_base"]
all_pops   <- c(targets, pool, right_base)
msg("targets (", length(targets), "): ", paste(targets, collapse = ", "))
msg("source pool (", length(pool), "): ", paste(pool, collapse = ", "))
msg("fixed RIGHT (", length(right_base), "): ", paste(right_base, collapse = ", "))
if (length(right_base) < 4) warning("Fewer than 4 fixed right populations; qpAdm will have little power.")
write_tsv(tibble(order = seq_along(pool), source = pool), file.path(OUTDIR, "source_order.tsv"))

# ------------------------------------------------------------ 0. f2 blocks ----
# The cache is trusted ONLY if a completion marker written after a successful extract_f2() matches the current
# population labelling (md5 of the labelled .fam + the population list). A killed/partial run, or a changed
# pool/target definition, therefore triggers a clean recomputation instead of a "file not found" error later.
f2_key <- c(paste(sort(all_pops), collapse = ","), unname(tools::md5sum(paste0(PLINK_PREFIX, ".fam"))))
f2_marker <- file.path(F2DIR, ".f2_complete")
have_f2 <- file.exists(f2_marker) && identical(readLines(f2_marker, warn = FALSE), f2_key)
if (dir.exists(F2DIR) && !have_f2 && !FORCE_F2)
  msg("existing f2 blocks are incomplete or were built for a different population set -> recomputing")
if (!have_f2 || FORCE_F2) {
  if (dir.exists(F2DIR)) unlink(F2DIR, recursive = TRUE)
  bimfile <- paste0(PLINK_PREFIX, ".bim")
  bim_head <- read_table(bimfile, col_names = FALSE, col_types = "cccccc", n_max = 200000, progress = FALSE)
  no_genetic_map <- all(as.numeric(bim_head$X3) == 0)
  blg <- if (BLGSIZE == "auto") { if (no_genetic_map) 2e6 else 0.05 } else as.numeric(BLGSIZE)
  msg("extract_f2: blgsize = ", blg, if (blg > 100) " bp (no genetic map in .bim)" else " Morgan")
  # All genotypes are diploid (ancient samples are imputed), so do NOT auto-detect
  # pseudohaploid samples. maxmiss = 0 keeps only SNPs present in every population.
  extract_f2(PLINK_PREFIX, F2DIR, pops = all_pops, blgsize = blg, maxmiss = 0,
             adjust_pseudohaploid = FALSE, auto_only = TRUE, n_cores = NCORES,
             maxmem = MAXMEM_MB, overwrite = TRUE, verbose = TRUE)
  writeLines(f2_key, f2_marker)
} else msg("re-using f2 blocks in ", F2DIR, " (set FORCE_F2=1 to recompute)")
f2 <- f2_from_precomp(F2DIR, pops = all_pops, afprod = TRUE, verbose = FALSE)
lost <- setdiff(all_pops, dimnames(f2)[[1]])
if (length(lost)) stop("Populations missing from f2 blocks: ", paste(lost, collapse = ", "))
msg("f2 array: ", paste(dim(f2), collapse = " x "), "  (pops x pops x SNP blocks)")

# ------------------------------------------------ 1. qpWave among targets ----
msg("qpWave: pairwise clade test between target groups")
right_all <- c(right_base, pool)
tp <- qpwave_pairs(f2, targets, right_all)
write_tsv(tp, file.path(OUTDIR, "qpwave_targets_pairwise.tsv"))
joint <- qpwave(f2, left = targets, right = right_all, verbose = FALSE)$rankdrop
write_tsv(joint, file.path(OUTDIR, "qpwave_targets_joint_rank.tsv"))
msg("qpWave (all targets jointly): rank ", paste0(joint$f4rank, " p=", signif(joint$p, 3), collapse = " | "))

# ---------------------------------------- 2. qpWave screen of source sets ----
combos <- unlist(lapply(seq_len(min(MAX_SOURCES, length(pool))),
                        function(k) combn(pool, k, simplify = FALSE)), recursive = FALSE)
right_for <- function(left) if (ROTATE) c(right_base, setdiff(pool, left)) else right_base
msg("candidate source sets: ", length(combos), " (1..", MAX_SOURCES, " sources)")

screen_one <- function(left) {
  k <- length(left)
  if (k == 1) return(tibble(left = left, nsrc = 1L, rank_tested = NA_integer_, chisq = NA_real_,
                            dof = NA_integer_, p_rankdrop = NA_real_, distinguishable = TRUE))
  o <- tryCatch(qpadm_p(f2, left, right_for(left), target = NULL, rnk = k - 2), error = function(e) NULL)
  if (is.null(o)) return(NULL)
  tibble(left = paste(left, collapse = "+"), nsrc = k, rank_tested = k - 2L,
         chisq = o$chisq[1], dof = o$dof[1], p_rankdrop = o$p[1],
         distinguishable = o$p[1] < PTHRESH)
}
msg("qpWave screen of source sets (rank k-2 must be REJECTED => sources are distinguishable)")
scr <- pmap_list(combos, screen_one) %>% bind_rows()
scr$left[scr$nsrc == 1] <- unlist(combos[lengths(combos) == 1])
write_tsv(scr, file.path(OUTDIR, "qpwave_source_screen.tsv"))
ok_sets <- scr$left[scr$distinguishable %in% TRUE]
combos <- combos[map_chr(combos, paste, collapse = "+") %in% ok_sets]
msg("source sets passing the qpWave screen: ", length(combos), " / ", nrow(scr))

# -------------------------------------------------------- 3. qpAdm rotating ----
fit_model <- function(target, left, right = right_for(left)) {
  k <- length(left)
  o <- tryCatch(qpadm_p(f2, left, right, target = target, weights = TRUE), error = function(e) NULL)
  if (is.null(o)) return(NULL)
  wcols <- setdiff(names(o), c("f4rank", "dof", "chisq", "p", "feasible"))
  w <- if (k == 1) 1 else as.numeric(unlist(o[1, wcols]))
  tibble(target = target, nsrc = k, left = paste(left, collapse = "+"),
         n_right = length(right), chisq = o$chisq[1], dof = o$dof[1], p = o$p[1],
         weights = paste(signif(w, 4), collapse = ";"), min_weight = min(w),
         feasible = all(w >= -WEIGHT_TOL & w <= 1 + WEIGHT_TOL))
}
jobs <- expand_grid(target = targets, i = seq_along(combos))
msg("qpAdm: ", nrow(jobs), " models (", length(targets), " targets x ", length(combos), " source sets)")
res <- pmap_list(seq_len(nrow(jobs)), function(j) fit_model(jobs$target[j], combos[[jobs$i[j]]])) %>% bind_rows()
# --- optional paper-style criterion: every weight +- CI_Z SE inside [0,1] (needs jackknife SEs, so only for models that already pass p and weights)
res$w_lo <- NA_real_; res$w_hi <- NA_real_
if (CI_Z > 0) {
  cand <- which(res$p >= PTHRESH & res$feasible)
  msg("jackknife SEs for ", length(cand), " models with p >= ", PTHRESH, " and weights in [0,1] (criterion: weight +- ", CI_Z, " SE inside [0,1])")
  if (length(cand)) {
    ci <- pmap_list(cand, function(i) {
      left <- strsplit(res$left[i], "\\+")[[1]]
      if (length(left) == 1) return(c(1, 1))
      o <- tryCatch(qpadm(f2, left = left, right = right_for(left), target = res$target[i], verbose = FALSE), error = function(e) NULL)
      if (is.null(o) || is.null(o$weights)) return(c(NA_real_, NA_real_))
      c(min(o$weights$weight - CI_Z * o$weights$se), max(o$weights$weight + CI_Z * o$weights$se))
    })
    ci <- do.call(rbind, ci)
    res$w_lo[cand] <- ci[, 1]; res$w_hi[cand] <- ci[, 2]
  }
}
res$ci_ok <- if (CI_Z > 0) (!is.na(res$w_lo) & res$w_lo >= 0 & res$w_hi <= 1) else TRUE
res <- res %>%
  mutate(status = case_when(p >= PTHRESH & feasible & min_weight >= MIN_WEIGHT & ci_ok ~ "FEASIBLE",
                            p >= PTHRESH & feasible & min_weight >= MIN_WEIGHT          ~ "WEIGHT_CI_FAIL",
                            p >= PTHRESH & feasible                                     ~ "NEGLIGIBLE_SOURCE",
                            feasible                                                    ~ "REJECTED_p",
                            TRUE                                                        ~ "INFEASIBLE_weights"))

# --- optional model competition (Narasimhan 2019, S4.3): for two feasible models A and B, the sources that only B has are added to the
# RIGHT set of A and A is re-tested. A is "eliminated" if it then fails while B survives the symmetric test; "contested" if both fail.
res$competition <- NA_character_
if (COMPETE && ROTATE) msg("COMPETE ignored: with ROTATE=1 every unused candidate is already on the right")
if (COMPETE && !ROTATE) {
  for (t in targets) {
    fi <- which(res$target == t & res$status == "FEASIBLE" & res$nsrc >= MIN_SOURCES)
    fi <- head(fi[order(res$nsrc[fi], -res$p[fi])], COMPETE_MAX)
    n <- length(fi)
    if (n < 2) { res$competition[fi] <- "winner"; next }
    srcs <- lapply(res$left[fi], function(x) strsplit(x, "\\+")[[1]])
    pr <- expand.grid(a = seq_len(n), b = seq_len(n)); pr <- pr[pr$a != pr$b, ]
    pr <- pr[vapply(seq_len(nrow(pr)), function(r) length(setdiff(srcs[[pr$b[r]]], srcs[[pr$a[r]]])) > 0, NA), ]
    msg("model competition for ", t, ": ", n, " feasible models, ", nrow(pr), " pairwise re-tests")
    fl <- unlist(pmap_list(seq_len(nrow(pr)), function(r) {
      A <- srcs[[pr$a[r]]]
      f <- fit_model(t, A, c(right_for(A), setdiff(srcs[[pr$b[r]]], A)))
      is.null(f) || !(f$p >= PTHRESH && f$feasible)
    }))
    Fm <- matrix(NA, n, n); Fm[cbind(pr$a, pr$b)] <- fl
    others <- function(a) setdiff(seq_len(n), a)
    elim   <- vapply(seq_len(n), function(a) any(vapply(others(a), function(b) isTRUE(Fm[a, b]) && !isTRUE(Fm[b, a]), NA)), NA)
    mutual <- vapply(seq_len(n), function(a) any(vapply(others(a), function(b) isTRUE(Fm[a, b]) &&  isTRUE(Fm[b, a]), NA)), NA)
    res$competition[fi] <- ifelse(elim, "eliminated", ifelse(mutual, "contested", "winner"))
  }
  msg("competition outcome: ", paste(names(table(res$competition)), table(res$competition), collapse = ", "))
}
write_tsv(res, file.path(OUTDIR, "qpadm_all_models.tsv"))
msg("qpAdm done: ", sum(res$status == "FEASIBLE" & res$nsrc >= MIN_SOURCES), " feasible models with >= ", MIN_SOURCES,
    " sources (p >= ", PTHRESH, ", every weight in [", MIN_WEIGHT, ",1])")

# 1-source fits: the target is statistically indistinguishable from ONE candidate w.r.t. the right set.
# Reported for transparency; not selected (MIN_SOURCES) because such a fit says "clade with X", not "admixed".
one <- res %>% filter(nsrc == 1, status == "FEASIBLE") %>% arrange(target, desc(p)) %>% select(target, source = left, p, chisq, dof)
write_tsv(one, file.path(OUTDIR, "qpadm_one_source_fits.tsv"))
if (nrow(one)) for (t in unique(one$target))
  msg("NOTE 1-source fit(s) for ", t, ": ", paste0(one$source[one$target == t], " (p=", signif(one$p[one$target == t], 2), ")", collapse = ", "))

# ------------------------------------------------------------ 4. selection ----
pick_models <- function(d) {
  d <- d %>% filter(nsrc >= MIN_SOURCES)
  feas <- d %>% filter(status == "FEASIBLE") %>% arrange(competition %in% "eliminated", nsrc, desc(p)) %>% mutate(selection = "feasible")
  if (nrow(feas)) return(feas)
  fb <- d %>% filter(feasible) %>% arrange(desc(p)) %>% mutate(selection = "NO_FIT_best_p_with_valid_weights")
  if (nrow(fb)) return(fb)
  d %>% arrange(desc(p)) %>% mutate(selection = "NO_FIT_best_p_invalid_weights")
}
ranked <- res %>% group_by(target) %>% group_modify(~ pick_models(.x)) %>% ungroup()
top <- ranked %>% group_by(target) %>% slice_head(n = TOP_N) %>% mutate(rank = row_number()) %>% ungroup()

# how often is each source part of a feasible model? (source-selection evidence)
feas_sel <- res %>% filter(status == "FEASIBLE", nsrc >= MIN_SOURCES)
nfeas <- feas_sel %>% count(target, name = "n_feasible_total")
src_freq <- feas_sel %>%
  mutate(src = strsplit(left, "\\+")) %>% unnest(src) %>%
  count(target, src, name = "n_feasible_models") %>%
  left_join(nfeas, by = "target") %>%
  mutate(share_of_feasible = n_feasible_models / n_feasible_total) %>%
  arrange(target, desc(n_feasible_models))
write_tsv(src_freq, file.path(OUTDIR, "qpadm_source_frequency_in_feasible_models.tsv"))
n_by_size <- res %>% filter(status == "FEASIBLE") %>% count(target, nsrc, name = "n_feasible")
write_tsv(n_by_size, file.path(OUTDIR, "qpadm_feasible_models_by_size.tsv"))

# best model at EACH size (FEASIBLE first, then models whose weights are all in [0,1], then by p). A high p with
# weights_valid = FALSE (negative / huge weights) is NOT a fit: the target lies outside the span of those sources.
# shows what 2-, 3-, 4- and 5-way models achieve,
# even when a simpler size already fits or nothing fits at all
best_by_size <- res %>% group_by(target, nsrc) %>%
  summarise(n_models = n(), n_feasible = sum(status == "FEASIBLE"),
            .groups = "drop") %>%
  left_join(res %>% group_by(target, nsrc) %>% arrange(desc(status == "FEASIBLE"), desc(feasible), desc(p)) %>% slice_head(n = 1) %>%
              ungroup() %>% select(target, nsrc, best_model = left, p, status, weights_valid = feasible, weights),
            by = c("target", "nsrc")) %>%
  arrange(target, nsrc)
write_tsv(best_by_size, file.path(OUTDIR, "qpadm_best_by_size.tsv"))
cat("\n---- best model at each size (p, status) ----\n")
print(best_by_size %>% mutate(p = signif(p, 3)) %>% select(target, nsrc, n_feasible, p, status, weights_valid, best_model), n = Inf, width = 200)

# --------------------------------- 5. weights +- SE and LOO-right robustness ----
msg("re-fitting the top ", TOP_N, " models per target with block-jackknife SEs")
refit <- function(i) {
  r <- top[i, ]
  left <- strsplit(r$left, "\\+")[[1]]
  right <- right_for(left)
  o <- tryCatch(qpadm(f2, left = left, right = right, target = r$target, verbose = FALSE),
                error = function(e) NULL)
  w <- if (!is.null(o) && !is.null(o$weights)) o$weights else
    tibble(target = r$target, left = left, weight = if (length(left) == 1) 1 else NA_real_, se = NA_real_, z = NA_real_)
  tibble(target = r$target, rank = r$rank, nsrc = r$nsrc, p = r$p, chisq = r$chisq, dof = r$dof,
         selection = r$selection, competition = r$competition, source = w$left, weight = w$weight, se = if (r$nsrc == 1) NA_real_ else w$se, z = w$z,
         right = paste(right, collapse = ","))
}
sel <- pmap_list(seq_len(nrow(top)), refit) %>% bind_rows()
write_tsv(sel, file.path(OUTDIR, "qpadm_selected_models.tsv"))

msg("leave-one-RIGHT-out robustness of the best model per target")
loo_one <- function(t) {
  r <- top %>% filter(target == t, rank == 1)
  left <- strsplit(r$left, "\\+")[[1]]
  right <- right_for(left)
  map_dfr(right, function(drop) {
    rr <- setdiff(right, drop)
    f <- fit_model(t, left, rr)
    if (is.null(f)) return(NULL)
    tibble(target = t, model = r$left, dropped_right = drop, p = f$p, feasible = f$feasible,
           still_fits = f$p >= PTHRESH & f$feasible)
  })
}
loo <- pmap_list(targets, loo_one) %>% bind_rows()
write_tsv(loo, file.path(OUTDIR, "qpadm_best_loo_right.tsv"))

# ------------------------------------------------------------ console summary ----
cat("\n================ BEST MODEL PER TARGET ================\n")
for (t in targets) {
  s <- sel %>% filter(target == t, rank == 1)
  l <- loo %>% filter(target == t)
  cat(sprintf("\n%s   [%s]   p = %.4g   (%d-way)   LOO-right: %d/%d still fit\n", t, s$selection[1], s$p[1],
              s$nsrc[1], sum(l$still_fits), nrow(l)))
  cat(sprintf("   %-18s %6.1f%%  +- %s\n", s$source, 100 * s$weight, ifelse(is.na(s$se), "NA", sprintf("%.1f", 100 * s$se))), sep = "")
}
cat("\nResults in: ", OUTDIR, "\n", sep = "")
