#!/usr/bin/env Rscript
# =============================================================================
# PSMC curves of the modern Iranians, coloured by sub-population (the five target groups)
#
#   Rscript scripts/plot_psmc.R [psmc_dir=psmc] [out_dir=figures] [g=25] [mu=1.25e-8] [s=100] [iter=20] [xmin=1e3] [xmax=1e5]
#
#   g   generation time in years
#   mu  mutation rate per site per generation
#   s   bin size in bp used when the .psmcfa files were made (fq2psmcfa -s, PSMC default 100)
#   iter  iteration to plot (psmc_plot.pl default 20; the psmc runs here used -N25)
#   xmin, xmax  time window in years (default 1e3 to 1e5, i.e. 10^3-10^5 years)
#
# Scaling follows psmc_plot.pl:  N0 = theta0 / (4 * mu * s);  time(years) = 2 * N0 * t_k * g;  Ne = N0 * lambda_k,
# as psmc_plot.pl does: iteration `iter`, and the final group of tied time intervals (the unreliable oldest bins) is not drawn. Files are matched to groups by sample ID (file name) through
# scripts/make_pops.py; files without a group (not in the sample metadata, or not one of the 68 clustered Iranians) are left out and listed.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(purrr); library(ggplot2); library(tibble) })
a <- commandArgs(TRUE)
PSMC_DIR <- if (length(a) >= 1) a[1] else "psmc"; OUT <- if (length(a) >= 2) a[2] else "figures"
G  <- if (length(a) >= 3) as.numeric(a[3]) else 25
MU <- if (length(a) >= 4) as.numeric(a[4]) else 1.25e-8
S  <- if (length(a) >= 5) as.numeric(a[5]) else 100
ITER <- if (length(a) >= 6) as.integer(a[6]) else 20
XMIN <- if (length(a) >= 7) as.numeric(a[7]) else 1e3
XMAX <- if (length(a) >= 8) as.numeric(a[8]) else 1e5
REPO <- Sys.getenv("REPO_DIR", "."); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- group of every sample (same definitions as the qpAdm pipeline)
tmp <- file.path(tempdir(), "grp")
status <- system2("python3", c("-I", file.path(REPO, "scripts/make_pops.py"), "--meta", file.path(REPO, "fam_meta/sample_meta.tsv"),
  "--fam", file.path(REPO, "fam_meta/all_anc_mod.run3.LDpruned.pca.fam"), "--pools", file.path(REPO, "config/pop_pools.tsv"), "--outdir", tmp),
  stdout = FALSE)
if (status != 0) stop("make_pops.py failed")
grp <- read_tsv(file.path(tmp, "sample_groups.tsv"), show_col_types = FALSE) %>% filter(role == "target") %>% select(IID, group)

# ---- parse the last iteration of each .psmc file
parse_psmc <- function(f) {
  x <- readLines(f)
  # last group of tied intervals is not plotted by psmc_plot.pl (unless -L): derive it from the -p pattern
  pat <- sub(".*pattern:([^,]+),.*", "\\1", grep("^MM\tpattern:", x, value = TRUE)[1])
  parts <- strsplit(gsub("\\s", "", pat), "\\+")[[1]]
  sizes <- unlist(lapply(parts, function(p) if (grepl("\\*", p)) rep(as.integer(sub(".*\\*", "", p)), as.integer(sub("\\*.*", "", p))) else as.integer(p)))
  last_group_n <- sizes[length(sizes)]; n_intv <- sum(sizes)
  max_k <- n_intv - last_group_n - 1                       # keep k <= max_k
  rd <- grep("^RD\t", x); rd_i <- grep(paste0("^RD\t", ITER, "$"), x)
  start <- if (length(rd_i)) rd_i[1] else rd[length(rd)]  # fall back to the last iteration
  nxt <- rd[rd > start]; blk <- x[start:(if (length(nxt)) nxt[1] - 1 else length(x))]
  tr <- strsplit(grep("^TR\t", blk, value = TRUE)[1], "\t")[[1]]; theta0 <- as.numeric(tr[2])
  rs <- do.call(rbind, strsplit(grep("^RS\t", blk, value = TRUE), "\t"))
  n0 <- theta0 / (4 * MU * S)
  tibble(IID = sub("\\.psmc$", "", basename(f)), k = as.integer(rs[, 2]), t_years = 2 * n0 * as.numeric(rs[, 3]) * G,
         Ne = n0 * as.numeric(rs[, 4]), N0 = n0, theta0 = theta0) %>% filter(k <= max_k)
}
files <- list.files(PSMC_DIR, pattern = "\\.psmc$", full.names = TRUE)
cur <- map_dfr(files, parse_psmc) %>% left_join(grp, by = "IID")
left_out <- unique(cur$IID[is.na(cur$group)])
missing_psmc <- setdiff(grp$IID, unique(cur$IID))
cur <- cur %>% filter(!is.na(group))
message(length(files), " PSMC files; ", n_distinct(cur$IID), " assigned to a group.")
if (length(left_out)) message("No group (left out): ", paste(left_out, collapse = ", "))
if (length(missing_psmc)) message("Clustered Iranians without a PSMC file: ", paste(missing_psmc, collapse = ", "))
message(sprintf("N0 (present-day scale) across samples: median %.0f, range %.0f-%.0f  [theta0 %.4f median; mu=%g, s=%g]",
                median(cur$N0), min(cur$N0), max(cur$N0), median(cur$theta0), MU, S))

# ---- step functions on a common log-time grid; group median across individuals
grid <- 10^seq(log10(XMIN), log10(XMAX), length.out = 400)
step_at <- function(d, t) { i <- findInterval(t, d$t_years); d$Ne[pmax(i, 1)] }
per_ind <- cur %>% group_by(IID, group) %>% group_modify(~ tibble(t = grid, Ne = step_at(arrange(.x, k), grid))) %>% ungroup()
med <- per_ind %>% group_by(group, t) %>% summarise(Ne = median(Ne), .groups = "drop")
write_tsv(cur, file.path(OUT, "psmc_curves_by_individual.tsv")); write_tsv(med, file.path(OUT, "psmc_group_median.tsv"))

# ---- plot
order <- c("Mod_Central_Iranian", "Mod_Zoroastrian", "Mod_Jewish_branch", "Mod_Armenian_Assyrian", "Mod_Turk")
nice  <- c(Mod_Central_Iranian = "Central Iranian", Mod_Zoroastrian = "Zoroastrian", Mod_Jewish_branch = "Jewish branch",
           Mod_Armenian_Assyrian = "Armenian/Assyrian", Mod_Turk = "Turk")
cols  <- c("#2a78d6", "#eb6834", "#1baf7a", "#e87ba4", "#6250d6"); names(cols) <- order
ns <- cur %>% distinct(IID, group) %>% count(group)
lab <- setNames(sprintf("%s (n=%d)", nice[order], ns$n[match(order, ns$group)]), order)
per_ind$group <- factor(per_ind$group, order); med$group <- factor(med$group, order)
inwin <- med %>% filter(t >= XMIN, t <= XMAX); unit <- if (max(inwin$Ne) < 1e4) 1e3 else 1e4
ymax <- ceiling(max(inwin$Ne) * 1.5 / unit * 2) / 2      # headroom for the legend and labels; individual spikes may exceed it
ink <- "#0b0b0b"; ink2 <- "#52514e"
p <- ggplot() +
  annotate("rect", xmin = 1e4, xmax = min(XMAX, 1.2e5), ymin = -Inf, ymax = Inf, fill = "#ececea") +
  annotate("text", x = sqrt(XMIN * 1e4), y = ymax * 0.99, label = "Holocene\n0-10 kya", vjust = 1, size = 3.3, colour = ink, lineheight = 0.95) +
  annotate("text", x = sqrt(1e4 * XMAX), y = ymax * 0.99, label = "Last glaciation\n10-120 kya", vjust = 1, size = 3.3, colour = ink, lineheight = 0.95) +
  geom_step(data = per_ind, aes(t, Ne / unit, group = IID, colour = group), linewidth = 0.3, alpha = 0.28, show.legend = FALSE) +
  geom_step(data = med, aes(t, Ne / unit, colour = group), linewidth = 1.2) +
  scale_colour_manual(values = cols, labels = lab, name = NULL, breaks = order) +
  scale_x_log10(limits = c(XMIN, XMAX), breaks = 10^(3:5), labels = c(expression(10^3), expression(10^4), expression(10^5)), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, ymax), clip = "on") +
  scale_y_continuous(expand = c(0, 0)) +
  annotation_logticks(sides = "b", size = 0.3, colour = ink2) +
  labs(x = bquote(paste("Years (g=", .(G), ", ", mu, "=", .(MU), ", s=", .(S), ")")),
       y = if (unit == 1e4) expression("Effective population size (" %*% 10^4 * ")") else expression("Effective population size (" %*% 10^3 * ")"),
       title = "PSMC of modern Iranians by sub-population", subtitle = "thin lines: individuals; thick lines: group median") +
  theme_classic(base_size = 12) +
  theme(legend.position = c(0.80, 0.70), legend.background = element_blank(), legend.key = element_blank(), legend.text = element_text(size = 10.5),
        axis.line = element_line(colour = ink, linewidth = 0.5), axis.text = element_text(colour = ink), plot.title.position = "plot", plot.margin = margin(8, 22, 8, 8),
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = ink2), panel.border = element_blank()) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.6, alpha = 1)))
ggsave(file.path(OUT, "psmc_by_subpopulation.pdf"), p, width = 9, height = 5.6)
tryCatch(ggsave(file.path(OUT, "psmc_by_subpopulation.png"), p, width = 9, height = 5.6, dpi = 300), error = function(e) message("PNG not written: ", conditionMessage(e)))
message("written to ", OUT)
