#!/usr/bin/env Rscript
# Plots for the qpWave / qpAdm results written by qpwave_qpadm.R
#   qpadm_best_models        stacked ancestry proportions of the best model per target
#   qpadm_top_models         the top-N competing models per target (model uncertainty)
#   qpwave_targets_pairwise  heat map of pairwise qpWave (clade) p-values between target groups
# Usage: OUTDIR=<results dir> Rscript plot_qpadm.R
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2) })

OUTDIR  <- Sys.getenv("OUTDIR"); if (OUTDIR == "") stop("set OUTDIR")
PTHRESH <- as.numeric(ifelse(Sys.getenv("PTHRESH") == "", 0.05, Sys.getenv("PTHRESH")))
rd <- function(f) read_tsv(file.path(OUTDIR, f), show_col_types = FALSE)
sel   <- rd("qpadm_selected_models.tsv")
order <- rd("source_order.tsv")
gs    <- rd("group_summary.tsv")
pairs <- rd("qpwave_targets_pairwise.tsv")

# ---- design tokens (validated categorical palette, light mode; fixed order, never cycled) ----
surface <- "#fcfcfb"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e6e5e1"; other <- c("#8f8e88", "#c4c3bd", "#62615c")   # neutral shades for sources beyond the 8 hues
hues <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")

# colour follows the SOURCE (config order), so a source keeps its colour across every figure.
best_src <- order$source[order$source %in% unique(sel$source[sel$rank == 1])]          # best-model sources get the 8 hues first
used <- c(best_src, order$source[order$source %in% setdiff(unique(sel$source), best_src)])
extra <- rep(other, length.out = max(0, length(used) - length(hues)))
cols <- setNames(c(hues, extra)[seq_along(used)], used)
if (length(used) > length(hues))
  message("NOTE: ", length(used) - length(hues), " source(s) beyond 8 are drawn in neutral greys (named in-bar).")

nlab <- setNames(sprintf("%s  (n=%d)", gs$group, gs$n), gs$group)
tg_order <- gs$group[gs$role == "target"]
tg_lab   <- function(x) factor(nlab[x], levels = rev(nlab[tg_order]))

theme_q <- function() {
  theme_minimal(base_size = 10) +
    theme(plot.background = element_rect(fill = surface, colour = NA),
          panel.background = element_rect(fill = surface, colour = NA),
          text = element_text(colour = ink), axis.text = element_text(colour = ink2),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = grid, linewidth = 0.3),
          legend.position = "bottom", legend.title = element_blank(),
          plot.title.position = "plot", plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(colour = ink2), strip.text = element_text(face = "bold", hjust = 0))
}
save2 <- function(p, stem, w, h) {
  ggsave(file.path(OUTDIR, paste0(stem, ".pdf")), p, width = w, height = h)
  tryCatch(ggsave(file.path(OUTDIR, paste0(stem, ".png")), p, width = w, height = h, dpi = 300),
           error = function(e) message("PNG not written for ", stem, ": ", conditionMessage(e)))
}
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
seg_label <- function(d) with(d, ifelse(weight >= 0.10,
      ifelse(is.na(se), sprintf("%s\n%.0f%%", source, 100 * weight),
             sprintf("%s\n%.0f%% +/- %.0f", source, 100 * weight, 100 * se)), ""))

# =============== 1. best model per target =================
best <- sel %>% filter(rank == 1) %>%
  mutate(target_f = tg_lab(target), source = factor(source, levels = names(cols)))
ann <- best %>% distinct(target_f, p, nsrc, selection) %>%
  mutate(txt = paste0(fmt_p(p), "  |  ", nsrc, "-way", ifelse(selection == "feasible", "", "\nNO FITTING MODEL")))
p1 <- ggplot(best, aes(x = weight, y = target_f, fill = source)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.62, colour = surface, linewidth = 0.8) +
  geom_text(aes(label = seg_label(best)), position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 2.7, colour = ink, lineheight = 0.9) +
  geom_text(data = ann, aes(x = 1.02, y = target_f, label = txt), inherit.aes = FALSE,
            hjust = 0, size = 2.9, colour = ink2, lineheight = 0.95) +
  scale_fill_manual(values = cols, drop = TRUE) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, .25), limits = c(0, 1.32),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Ancestry proportion", y = NULL, title = "qpAdm: best model per modern Iranian group",
       subtitle = sprintf("simplest model with p >= %g and all weights in [0,1]; +/- = block-jackknife SE", PTHRESH)) +
  guides(fill = guide_legend(nrow = 2)) + theme_q()
save2(p1, "qpadm_best_models", 10, 1.6 + 0.95 * length(unique(best$target)))

# =============== 2. top-N models per target =================
topn <- sel %>% mutate(source = factor(source, levels = names(cols)),
                       model = sprintf("#%d  %s", rank, fmt_p(p)),
                       target_t = factor(target, levels = tg_order)) %>%
  group_by(target_t) %>% mutate(model = factor(model, levels = rev(unique(model[order(rank)])))) %>% ungroup()
p2 <- ggplot(topn, aes(x = weight, y = model, fill = source)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.62, colour = surface, linewidth = 0.8) +
  geom_text(aes(label = seg_label(topn)), position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 2.5, colour = ink, lineheight = 0.9) +
  facet_wrap(~target_t, ncol = 1, scales = "free_y", labeller = labeller(target_t = function(x) nlab[x])) +
  scale_fill_manual(values = cols, drop = TRUE) +
  scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, .25), expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Ancestry proportion", y = NULL, title = "qpAdm: competing models (top per target)",
       subtitle = "ranked by parsimony (fewest sources), then p-value; models flagged in the table if none fit") +
  guides(fill = guide_legend(nrow = 2)) + theme_q()
save2(p2, "qpadm_top_models", 10, 1.2 + 1.45 * length(unique(topn$target)))

# =============== 3. qpWave pairwise target heat map =================
pw <- pairs %>% mutate(a = factor(pop1, levels = tg_order), b = factor(pop2, levels = rev(tg_order)),
                       lp = pmin(-log10(pmax(p, 1e-300)), 12))
p3 <- ggplot(pw, aes(a, b, fill = lp)) +
  geom_tile(colour = surface, linewidth = 1.2) +
  geom_text(aes(label = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)), colour = ifelse(lp > 6, "#ffffff", ink)), size = 3) +
  scale_colour_identity() +
  scale_fill_gradient(low = "#e8f0fb", high = "#1f5fae", name = expression(-log[10](p)), limits = c(0, 12)) +
  labs(x = NULL, y = NULL, title = "qpWave: are the target groups distinguishable?",
       subtitle = sprintf("rank-0 (clade) test, pairwise; p < %g (darker) = the two groups are NOT a clade w.r.t. the right set", PTHRESH)) +
  theme_q() + theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.major.x = element_blank(),
                    legend.position = "right")
save2(p3, "qpwave_targets_pairwise", 7.5, 5)

message("plots written to ", OUTDIR)
