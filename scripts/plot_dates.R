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

lev <- rev(summ$group)                                   # oldest group at the top of the plot
d$group <- factor(d$group, levels = lev)
d <- d %>% arrange(group, mid) %>% group_by(group) %>%
  mutate(off = if (n() > 1) seq(-0.30, 0.30, length.out = n()) else 0) %>% ungroup() %>% mutate(y = as.numeric(group) + off)
sg <- summ %>% mutate(group = factor(group, levels = lev), y = as.numeric(group), sd0 = ifelse(is.na(sd_mid), 0, sd_mid))
ylab_txt <- setNames(sprintf("%s\nn=%d, %s BCE%s", summ$group, summ$n, fy(summ$mean_mid),
                             ifelse(is.na(summ$sd_mid), "", sprintf(" +/- %s", fy(summ$sd_mid)))), summ$group)
# colours match the user's legend (Office theme colours); both Shahr-i-Sokhta groups share the purple, both Hasanlu groups share the green
cols <- c(Iran_ChL_HajjiFiruz = "#156082", Iran_ChL_TepeHissar = "#e97132", Iran_BA_ShahrISokhta_BA1 = "#a02b93",
          Iran_BA_ShahrISokhta_IRN = "#a02b93", Iran_BA_Hasanlu = "#4ea72e", Iran_IA_Hasanlu = "#4ea72e")
ink <- "#0b0b0b"; ink2 <- "#52514e"

p <- ggplot() +
  geom_rect(data = sg %>% filter(!is.na(sd_mid)), aes(xmin = mean_mid - sd0, xmax = mean_mid + sd0, ymin = y - 0.42, ymax = y + 0.42, fill = group), alpha = 0.16) +
  geom_segment(data = sg, aes(x = mean_mid, xend = mean_mid, y = y - 0.42, yend = y + 0.42, colour = group), linewidth = 0.9) +
  geom_segment(data = d, aes(x = older, xend = younger, y = y, yend = y, colour = group, linetype = kind), linewidth = 1.1, lineend = "butt", alpha = 0.85) +
  geom_point(data = d, aes(x = mid, y = y, colour = group), size = 1.1, shape = 21, fill = "#fcfcfb", stroke = 0.7) +
  scale_colour_manual(values = cols, guide = "none") + scale_fill_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = c("radiocarbon (95% calibrated interval)" = "solid", "archaeological context range" = "22"), name = NULL) +
  scale_x_continuous(breaks = seq(-6000, -500, 500), labels = function(x) abs(x), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(breaks = seq_along(lev), labels = ylab_txt[lev], expand = expansion(add = 0.5)) +
  labs(x = "Years BCE (calibrated)", y = NULL, title = "Archaeological dates of the ancient Iranian target groups",
       subtitle = "bars = individual intervals (dots = midpoints); line = group mean; band = mean +/- 1 SD across individuals") +
  theme_minimal(base_size = 9) +
  theme(plot.background = element_rect(fill = "#fcfcfb", colour = NA), panel.background = element_rect(fill = "#fcfcfb", colour = NA),
        text = element_text(colour = ink), axis.text = element_text(colour = ink2), axis.text.y = element_text(colour = ink, hjust = 0, size = 7.5, lineheight = 0.95),
        panel.grid.major.y = element_line(colour = "#efeeea", linewidth = 0.3), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "#e6e5e1", linewidth = 0.4), legend.position = "bottom", legend.key.width = unit(1.4, "lines"),
        legend.margin = margin(0, 0, 0, 0), plot.title.position = "plot", plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(colour = ink2, size = 8), plot.margin = margin(6, 10, 4, 6))
ggsave(file.path(OUT, "ancient_iran_dates.pdf"), p, width = 8.2, height = 3.9)
tryCatch(ggsave(file.path(OUT, "ancient_iran_dates.png"), p, width = 8.2, height = 3.9, dpi = 300), error = function(e) message("PNG not written: ", conditionMessage(e)))
message("written to ", OUT)
