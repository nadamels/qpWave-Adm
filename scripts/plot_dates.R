#!/usr/bin/env Rscript
# =============================================================================
# Archaeological dates of the ancient Iranian target groups: one bar per individual, mean +- SD per group
#
#   Rscript scripts/plot_dates.R [out_dir=figures]
#
# Dates come from fam_meta/sample_meta.tsv ("Archeological Date"). Two kinds are present and drawn differently:
#   * radiocarbon: "4780-3518 calBCE (4780+-30 BP, lab)" -> the 95% calibrated interval (solid bar)
#   * archaeological context: "2600-2500 BCE" -> the stated range (dashed, lighter)
# Per individual the midpoint of the interval is used. Per group: mean of the midpoints and the SD ACROSS individuals (n-1; undefined for n = 1).
# The SD describes how spread out in time the pooled individuals are (a homogeneity check); it is not a measurement error.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2); library(stringr); library(tibble) })
a <- commandArgs(TRUE); OUT <- if (length(a) >= 1) a[1] else "figures"
REPO <- Sys.getenv("REPO_DIR", "."); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

tmp <- file.path(tempdir(), "grp")
st <- system2("python3", c("-I", file.path(REPO, "scripts/make_pops.py"), "--meta", file.path(REPO, "fam_meta/sample_meta.tsv"),
  "--fam", file.path(REPO, "fam_meta/all_anc_mod.run3.LDpruned.pca.fam"), "--pools", file.path(REPO, "config/pop_pools_ancient_iran_distal.tsv"),
  "--targets", "config", "--outdir", tmp), stdout = FALSE)
if (st != 0) stop("make_pops.py failed")
grp  <- read_tsv(file.path(tmp, "sample_groups.tsv"), show_col_types = FALSE) %>% filter(role == "target") %>% select(IID, group, population_label)
meta <- read_tsv(file.path(REPO, "fam_meta/sample_meta.tsv"), show_col_types = FALSE) %>% select(IID, date_raw = `Archeological Date`)
d <- grp %>% left_join(meta, by = "IID")

rngm <- str_match(d$date_raw, "^\\s*([0-9]{3,5})\\s*-\\s*([0-9]{3,5})\\s*(cal)?BCE")
bpm <- str_match(d$date_raw, "([0-9]{3,5})\\s*±\\s*([0-9]{1,3})\\s*BP")
d <- d %>% mutate(older = -as.numeric(rngm[, 2]), younger = -as.numeric(rngm[, 3]),           # negative = BCE, so older is further left
                  mid = (older + younger) / 2, halfwidth = (younger - older) / 2,
                  kind = ifelse(!is.na(bpm[, 2]), "radiocarbon (95% calibrated interval)", "archaeological context range"),
                  bp = as.numeric(bpm[, 2]), bp_sd = as.numeric(bpm[, 3]))
if (any(is.na(d$mid))) stop("unparsed dates: ", paste(d$IID[is.na(d$mid)], collapse = ", "))

# groups asked for (the two Shahr-i-Sokhta groups are kept separate, as in the qpAdm config)
keep <- c("Iran_ChL_HajjiFiruz", "Iran_ChL_TepeHissar", "Iran_BA_ShahrISokhta_BA1", "Iran_BA_ShahrISokhta_IRN", "Iran_BA_Hasanlu", "Iran_IA_Hasanlu")
d <- d %>% filter(group %in% keep)

summ <- d %>% group_by(group) %>% summarise(
  n = n(), n_radiocarbon = sum(grepl("^radiocarbon", kind)),
  mean_mid = mean(mid), sd_mid = ifelse(n() > 1, sd(mid), NA_real_),
  mean_mid_radiocarbon = ifelse(any(grepl("^radiocarbon", kind)), mean(mid[grepl("^radiocarbon", kind)]), NA_real_),
  sd_mid_radiocarbon = ifelse(sum(grepl("^radiocarbon", kind)) > 1, sd(mid[grepl("^radiocarbon", kind)]), NA_real_),
  earliest = min(older), latest = max(younger), .groups = "drop") %>% arrange(mean_mid)
write_tsv(d %>% arrange(group, mid) %>% select(group, IID, population_label, kind, older_BCE = older, younger_BCE = younger, midpoint = mid, bp, bp_sd, date_raw),
          file.path(OUT, "ancient_iran_dates_by_individual.tsv"))
write_tsv(summ, file.path(OUT, "ancient_iran_dates_by_group.tsv"))
fy <- function(x) sprintf("%.0f", abs(x))
cat("\nGroup summary (midpoints; negative = BCE):\n"); print(summ %>% mutate(across(where(is.numeric), ~ round(.x, 0))), width = 200)

lev <- summ$group
d$group <- factor(d$group, levels = lev)
d <- d %>% arrange(group, mid) %>% group_by(group) %>% mutate(row = row_number()) %>% ungroup()
strip <- setNames(sprintf("%s  (n=%d)\nmean %s BCE%s", summ$group, summ$n, fy(summ$mean_mid),
                          ifelse(is.na(summ$sd_mid), "  (SD undefined, n=1)", sprintf("  +/- %s yr SD", fy(summ$sd_mid)))), summ$group)
sg <- summ %>% mutate(group = factor(group, levels = lev), lo = mean_mid - ifelse(is.na(sd_mid), 0, sd_mid), hi = mean_mid + ifelse(is.na(sd_mid), 0, sd_mid))
hues <- c("#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#6250d6", "#eda100"); cols <- setNames(hues, keep)
ink <- "#0b0b0b"; ink2 <- "#52514e"

p <- ggplot(d) +
  geom_rect(data = sg %>% filter(!is.na(sd_mid)), aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf, fill = group), alpha = 0.13, inherit.aes = FALSE) +
  geom_vline(data = sg, aes(xintercept = mean_mid, colour = group), linewidth = 0.7, linetype = "solid", alpha = 0.9) +
  geom_segment(aes(x = older, xend = younger, y = row, yend = row, colour = group, linetype = kind), linewidth = 1.6, lineend = "butt", alpha = 0.85) +
  geom_point(aes(x = mid, y = row, colour = group), size = 1.8, shape = 21, fill = "#fcfcfb", stroke = 0.9) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", labeller = labeller(group = strip), switch = "y") +
  scale_colour_manual(values = cols, guide = "none") + scale_fill_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = c("radiocarbon (95% calibrated interval)" = "solid", "archaeological context range" = "22"), name = NULL) +
  scale_x_continuous(breaks = seq(-6000, -500, 500), labels = function(x) abs(x), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(add = 0.7)) +
  labs(x = "Years BCE (calibrated)", y = NULL, title = "Archaeological dates of the ancient Iranian target groups",
       subtitle = "one bar per individual (interval) with its midpoint; vertical line = group mean of midpoints; shaded band = mean +/- 1 SD across individuals") +
  theme_minimal(base_size = 10.5) +
  theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA), panel.background = element_rect(fill = "#fcfcfb", colour = NA),
        text = element_text(colour = ink), axis.text = element_text(colour = ink2), axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_line(colour = "#e6e5e1", linewidth = 0.4),
        strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold", size = 9), strip.placement = "outside",
        panel.spacing.y = unit(0.5, "lines"), legend.position = "bottom", plot.title.position = "plot", plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = ink2), plot.margin = margin(8, 14, 8, 8))
ht <- 2.2 + 0.19 * nrow(d) + 0.5 * length(lev)
ggsave(file.path(OUT, "ancient_iran_dates.pdf"), p, width = 12.5, height = ht)
tryCatch(ggsave(file.path(OUT, "ancient_iran_dates.png"), p, width = 12.5, height = ht, dpi = 250), error = function(e) message("PNG not written: ", conditionMessage(e)))
message("written to ", OUT)
