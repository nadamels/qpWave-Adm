#!/usr/bin/env Rscript
# =============================================================================
# Homogeneity QC of the pooled groups (targets and candidate sources)
#
# Every target/pool group is split into sub-labels (published label x sampling region, see make_pops.py --sublabels). For each group,
# all pairs of sub-labels are tested with qpWave rank 0 ("do these two sub-labels form a clade?") against the design's fixed right set,
# the same logic as Narasimhan et al. 2019 SM S3.1 (individuals from a site tested against each other with the RightAll set).
# A group whose sub-labels are mutually rejected is not a homogeneous source/target and should be split or have outliers removed.
#
# Limits: this is group-level (sub-label pairs), not individual-level outlier detection; sub-labels with 1 individual have low power;
# groups with a single sub-label cannot be tested. The expected rejection rate for a homogeneous group is about HOMOG_P.
# =============================================================================
suppressPackageStartupMessages({ library(admixtools); library(dplyr); library(tidyr); library(readr); library(purrr); library(tibble); library(ggplot2) })
getenv <- function(k, d = NULL) { v <- Sys.getenv(k, unset = NA); if (is.na(v) || v == "") { if (is.null(d)) stop("Environment variable ", k, " is required"); d } else v }
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
PLINK_PREFIX_SUB <- getenv("PLINK_PREFIX_SUB"); OUTDIR <- getenv("OUTDIR")
SUBGROUPS     <- getenv("SUBGROUPS", file.path(OUTDIR, "sample_subgroups.tsv"))
GROUP_SUMMARY <- getenv("GROUP_SUMMARY", file.path(OUTDIR, "group_summary.tsv"))
F2DIR  <- getenv("F2DIR_SUB", file.path(OUTDIR, "f2_blocks_sub"))
NCORES <- as.integer(getenv("NCORES", 4)); MAXMEM_MB <- as.integer(getenv("MAXMEM_MB", 16000))
HOMOG_P <- as.numeric(getenv("HOMOG_P", 0.01)); BLGSIZE <- getenv("BLGSIZE", "auto"); FORCE_F2 <- getenv("FORCE_F2", "0") == "1"
OUT <- file.path(OUTDIR, "homogeneity"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

sub <- read_tsv(SUBGROUPS, show_col_types = FALSE)
gs  <- read_tsv(GROUP_SUMMARY, show_col_types = FALSE)
right <- gs$group[gs$role == "right_base"]
subsum <- sub %>% count(group, role, sublabel, name = "n")
pops <- c(subsum$sublabel, right)
msg("sub-labels: ", nrow(subsum), " in ", n_distinct(subsum$group), " groups; fixed right set (", length(right), "): ", paste(right, collapse = ", "))

key <- c(paste(sort(pops), collapse = ","), unname(tools::md5sum(paste0(PLINK_PREFIX_SUB, ".fam"))))
marker <- file.path(F2DIR, ".f2_complete")
ok <- file.exists(marker) && identical(readLines(marker, warn = FALSE), key)
if (!ok || FORCE_F2) {
  if (dir.exists(F2DIR)) unlink(F2DIR, recursive = TRUE)
  bim <- read_table(paste0(PLINK_PREFIX_SUB, ".bim"), col_names = FALSE, col_types = "cccccc", n_max = 200000, progress = FALSE)
  blg <- if (BLGSIZE == "auto") { if (all(as.numeric(bim$X3) == 0)) 2e6 else 0.05 } else as.numeric(BLGSIZE)
  extract_f2(PLINK_PREFIX_SUB, F2DIR, pops = pops, blgsize = blg, maxmiss = 0, adjust_pseudohaploid = FALSE,
             auto_only = TRUE, n_cores = NCORES, maxmem = MAXMEM_MB, overwrite = TRUE, verbose = FALSE)
  writeLines(key, marker)
} else msg("re-using f2 blocks in ", F2DIR)
f2 <- f2_from_precomp(F2DIR, pops = pops, afprod = TRUE, verbose = FALSE)

testable <- subsum %>% group_by(group) %>% filter(n() >= 2) %>% ungroup()
untestable <- setdiff(unique(subsum$group), unique(testable$group))
msg("groups tested: ", n_distinct(testable$group), "; groups with a single sub-label (cannot be tested): ", paste(untestable, collapse = ", "))
pairs <- map_dfr(unique(testable$group), function(g) {
  l <- testable$sublabel[testable$group == g]
  qpwave_pairs(f2, l, right) %>% filter(pop1 < pop2) %>% mutate(group = g)
}) %>%
  left_join(subsum %>% select(pop1 = sublabel, n1 = n), by = "pop1") %>% left_join(subsum %>% select(pop2 = sublabel, n2 = n), by = "pop2") %>%
  mutate(reject = p < HOMOG_P)
write_tsv(pairs, file.path(OUT, "homogeneity_pairs.tsv"))

by_group <- pairs %>% group_by(group) %>%
  summarise(n_sublabels = n_distinct(c(pop1, pop2)), n_pairs = n(), n_rejected = sum(reject), frac_rejected = mean(reject), min_p = min(p), .groups = "drop") %>%
  mutate(verdict = case_when(n_rejected == 0 ~ "no pair rejected",
                             frac_rejected <= 0.10 ~ "few pairs rejected (check which)",
                             TRUE ~ "HETEROGENEOUS: many pairs rejected")) %>% arrange(desc(frac_rejected))
by_sub <- bind_rows(pairs %>% transmute(group, sublabel = pop1, reject), pairs %>% transmute(group, sublabel = pop2, reject)) %>%
  group_by(group, sublabel) %>% summarise(n_pairs = n(), n_rejected = sum(reject), frac = mean(reject), .groups = "drop") %>%
  left_join(subsum %>% select(sublabel, n), by = "sublabel") %>%
  mutate(flag = case_when(n_pairs >= 3 & frac >= 0.5 & n >= 2 ~ "outlier candidate",
                          frac >= 0.5 ~ "check (few partners or n=1)", TRUE ~ "")) %>% arrange(group, desc(frac))
write_tsv(by_group, file.path(OUT, "homogeneity_by_group.tsv")); write_tsv(by_sub, file.path(OUT, "homogeneity_by_sublabel.tsv"))

cat("\n================ HOMOGENEITY (qpWave rank-0 between sub-labels, p <", HOMOG_P, "= rejected; ~", 100 * HOMOG_P, "% expected by chance) ================\n")
print(by_group %>% mutate(min_p = signif(min_p, 3), frac_rejected = round(frac_rejected, 2)), n = Inf, width = 200)
fl <- by_sub %>% filter(flag != "")
if (nrow(fl)) { cat("\nSub-labels rejected against >= 50% of their group:\n"); print(fl %>% select(group, sublabel, n, n_pairs, n_rejected, flag), n = Inf, width = 200) }
if (length(untestable)) cat("\nNot testable (one sub-label): ", paste(untestable, collapse = ", "), "\n")

# heat map
shorten <- function(x, g) substring(x, nchar(g) + 3)
pl <- bind_rows(pairs %>% transmute(group, a = shorten(pop1, group), b = shorten(pop2, group), p),
                pairs %>% transmute(group, a = shorten(pop2, group), b = shorten(pop1, group), p)) %>% mutate(lp = pmin(-log10(pmax(p, 1e-300)), 10))
g <- ggplot(pl, aes(a, b, fill = lp)) + geom_tile(colour = "white") +
  facet_wrap(~group, scales = "free", ncol = 3) +
  scale_fill_gradient(low = "#e8f0fb", high = "#1f5fae", limits = c(0, 10), name = "-log10(p)") +
  labs(x = NULL, y = NULL, title = "Homogeneity of pooled groups: qpWave rank-0 between sub-labels",
       subtitle = sprintf("dark = the two sub-labels are NOT a clade against the fixed right set (p < %g is rejected)", HOMOG_P)) +
  theme_minimal(base_size = 8) + theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom", plot.title.position = "plot")
ng <- n_distinct(pl$group); ggsave(file.path(OUT, "homogeneity_heatmaps.pdf"), g, width = 14, height = 3.5 * ceiling(ng / 3) + 1)
tryCatch(ggsave(file.path(OUT, "homogeneity_heatmaps.png"), g, width = 14, height = 3.5 * ceiling(ng / 3) + 1, dpi = 200), error = function(e) NULL)
msg("done: ", OUT)
