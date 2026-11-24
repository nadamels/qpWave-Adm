#!/usr/bin/env Rscript
# Plots for explore_admixture_pca.py:  OUTDIR=<explore dir> Rscript plot_admixture_pca.R
#   admixture_group_profiles   mean ADMIXTURE ancestry of targets, candidate sources, outgroups
#   pca_targets_vs_candidates  PC1-PC2 and PC3-PC4: target individuals vs. candidate/outgroup centroids
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2) })
OUTDIR <- Sys.getenv("OUTDIR"); if (OUTDIR == "") stop("set OUTDIR")
rd <- function(f) read_tsv(file.path(OUTDIR, f), show_col_types = FALSE)
gsum <- rd("admixture_pca_group_summary.tsv"); smp <- rd("admixture_pca_samples.tsv"); lab <- rd("admixture_component_labels.tsv")

surface <- "#fcfcfb"; ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e6e5e1"
hues <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
theme_q <- function() theme_minimal(base_size = 10) + theme(
  plot.background = element_rect(fill = surface, colour = NA), panel.background = element_rect(fill = surface, colour = NA),
  text = element_text(colour = ink), axis.text = element_text(colour = ink2), panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = grid, linewidth = 0.3), legend.position = "bottom", legend.title = element_blank(),
  plot.title.position = "plot", plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(colour = ink2),
  strip.text = element_text(face = "bold", hjust = 0))
save2 <- function(p, stem, w, h) {
  ggsave(file.path(OUTDIR, paste0(stem, ".pdf")), p, width = w, height = h)
  tryCatch(ggsave(file.path(OUTDIR, paste0(stem, ".png")), p, width = w, height = h, dpi = 300),
           error = function(e) message("PNG not written: ", conditionMessage(e)))
}
K <- sum(grepl("^Q[0-9]+$", names(gsum)))
short <- vapply(seq_len(K), function(k) paste0("K", k, "  ", paste(head(strsplit(lab$top_groups_by_mean_fraction[k], ", ")[[1]], 2), collapse = " / ")), "")
role_lab <- c(target = "Modern Iranian targets", pool = "Candidate sources (pool)", ignore = "Ignored (post-Neolithic Iran)", right_base = "Fixed RIGHT outgroups")

# ---- 1. ADMIXTURE profiles
d <- gsum %>% select(group, role, n, starts_with("Q")) %>% select(-starts_with("sdQ")) %>%
  pivot_longer(starts_with("Q"), names_to = "comp", values_to = "frac") %>%
  mutate(comp = factor(comp, levels = paste0("Q", seq_len(K)), labels = short),
         role = factor(role, levels = names(role_lab), labels = role_lab),
         glab = factor(sprintf("%s (n=%d)", group, n), levels = rev(unique(sprintf("%s (n=%d)", gsum$group, gsum$n)))))
p1 <- ggplot(d, aes(frac, glab, fill = comp)) +
  geom_col(width = 0.72, colour = surface, linewidth = 0.4) +
  facet_grid(role ~ ., scales = "free_y", space = "free_y", labeller = label_wrap_gen(16)) +
  scale_fill_manual(values = setNames(hues[seq_len(K)], levels(d$comp))) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Mean ADMIXTURE ancestry", y = NULL, title = sprintf("Unsupervised ADMIXTURE (K=%d): who carries which component", K),
       subtitle = "hypothesis generation only: components are drifted clusters, not ancestral populations") +
  guides(fill = guide_legend(nrow = 4)) + theme_q() + theme(panel.grid.major.y = element_blank())
save2(p1, "admixture_group_profiles", 10, 2.5 + 0.2 * nrow(gsum))

# ---- 2. PCA
tg <- unique(gsum$group[gsum$role == "target"])
tcol <- setNames(hues[seq_along(tg)], tg)
pairs <- list(c("PC1", "PC2"), c("PC3", "PC4"))
long <- function(df, id) bind_rows(lapply(pairs, function(p) df %>% transmute(!!id := .data[[id]], x = .data[[p[1]]], y = .data[[p[2]]],
                                                                                role, panel = paste(p, collapse = " vs "))))
ind <- long(smp %>% filter(role == "target") %>% mutate(group = group) %>% select(group, role, PC1:PC4), "group")
cen <- long(gsum %>% filter(role != "target") %>% select(group, role, PC1:PC4), "group")
# zoom on the targets: keep only centroids within +-1x the targets' range of the target cloud (far outgroups such as Mbuti
# or Karitiana would otherwise squeeze all Iranians into one corner); dropped centroids are listed in the subtitle
win <- ind %>% group_by(panel) %>% summarise(xl = min(x) - diff(range(x)), xh = max(x) + diff(range(x)),
                                             yl = min(y) - diff(range(y)), yh = max(y) + diff(range(y)), .groups = "drop")
cen <- cen %>% left_join(win, by = "panel") %>% mutate(inside = x >= xl & x <= xh & y >= yl & y <= yh)
outside <- sort(unique(cen$group[!cen$inside])); outside <- setdiff(outside, unique(cen$group[cen$inside]))
cen <- cen %>% filter(inside)
p2 <- ggplot() +
  geom_point(data = cen, aes(x, y, shape = role), colour = ink2, size = 2.2) +
  geom_text(data = cen, aes(x, y, label = group), colour = ink2, size = 2.4, vjust = -0.9, check_overlap = TRUE) +
  geom_point(data = ind, aes(x, y, colour = group), size = 2, alpha = 0.85) +
  facet_wrap(~panel, scales = "free") +
  scale_colour_manual(values = tcol) + scale_shape_manual(values = c(pool = 17, ignore = 4, right_base = 15),
                                                           labels = c(pool = "candidate (centroid)", ignore = "ignored (centroid)", right_base = "right outgroup (centroid)")) +
  labs(x = NULL, y = NULL, title = "PCA: modern Iranian individuals vs ancient / outgroup centroids",
       subtitle = paste0("zoomed on the targets; not shown (far away in both panels): ", paste(outside, collapse = ", "))) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.001), expand = expansion(mult = 0.1)) + scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  theme_q()
save2(p2, "pca_targets_vs_candidates", 12, 5.5)
message("plots written to ", OUTDIR)
