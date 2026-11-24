#!/usr/bin/env Rscript
# Colour legend for the K=7 / K=8 ADMIXTURE bar plots of the modern Iranians.
# Colours are ColorBrewer "Set1" (what the bar plot uses). Labels are checked against the K=8 .Q file:
# for every component the populations with the highest mean fraction are printed next to the label.
# usage: Rscript scripts/make_admixture_legend.R [repo_dir] [out_dir]
suppressPackageStartupMessages({ library(ggplot2); library(readr); library(dplyr) })
args <- commandArgs(TRUE); repo <- if (length(args) >= 1) args[1] else "."; out <- if (length(args) >= 2) args[2] else file.path(repo, "figures")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# column of the .Q file  ->  colour in the figure  ->  label
leg <- tribble(
  ~comp, ~colour,   ~colour_name, ~label,
  5,     "#f781bf", "pink",       "Neolithic Iran / Turkmenistan",
  4,     "#4daf4a", "green",      "Anatolian Neolithic farmer",
  3,     "#a65628", "brown",      "Hunter-gatherer: EEHG / WHG (also Steppe)",
  1,     "#984ea3", "purple",     "ASI (South Asian)",
  6,     "#e41a1c", "red",        "East Asian",
  2,     "#377eb8", "blue",       "African",
  8,     "#ff7f00", "orange",     "Indigenous American (IAM)",
  7,     "#ffff33", "yellow",     "Oceanian (Papuan)")

meta <- read_tsv(file.path(repo, "fam_meta/sample_meta.tsv"), show_col_types = FALSE)
fam  <- read.table(file.path(repo, "fam_meta/all_anc_mod.run3.LDpruned.pca.fam"), stringsAsFactors = FALSE)$V2
Q    <- as.matrix(read.table(file.path(repo, "fam_meta/all_anc_mod.run3.LDpruned.8.Q")))
stopifnot(nrow(Q) == length(fam))
pop  <- meta$Population[match(fam, meta$IID)]
support <- function(c, k = 4, minn = 3) {
  d <- tibble(pop = pop, q = Q[, c]) %>% group_by(pop) %>% summarise(m = mean(q), n = n(), .groups = "drop") %>%
    filter(n >= minn) %>% arrange(desc(m)) %>% slice_head(n = k)
  paste0(d$pop, " (", sprintf("%.2f", d$m), ")", collapse = ", ")
}
leg$highest_in <- vapply(leg$comp, support, "")
write_tsv(leg %>% transmute(colour_name, hex = colour, label, K8_component = paste0("K", comp), highest_in_K8_data = highest_in),
          file.path(out, "admixture_legend.tsv"))

ink <- "#0b0b0b"; ink2 <- "#52514e"; surface <- "#ffffff"
leg$y <- rev(seq_len(nrow(leg)))
note <- paste0(
  "Check against the K=8 .Q file (the 70 modern Iranians in the figure): pink 41%, green 39%, brown 8%, ASI (purple) only 4%.\n",
  "K=7 has no pink. The K=7 purple band (about a quarter of these genomes in the figure) is much larger than true ASI (4% at K=8),\n",
  "so it probably also carries Iranian-Neolithic-related ancestry (the component that becomes pink at K=8). Treat the K=7 'ASI' label with care.\n",
  "Green: the data support 'Anatolian Neolithic farmer' (max in Barcin_N and Anatolia_N, also early European farmers), not specifically PPN:\n",
  "the only Pre-Pottery Anatolian sample (Mardin PPN) carries 45% green + 29% pink.   The K=7 .Q file is not in the repo, so K=7 was not re-checked.")
note <- gsub("\\\\n", "\n", note)
p <- ggplot(leg) +
  geom_tile(aes(x = 0, y = y), width = 0.55, height = 0.7, fill = leg$colour, colour = ink2, linewidth = 0.4) +
  geom_text(aes(x = 0.5, y = y + 0.06, label = label), hjust = 0, vjust = 0, size = 4, fontface = "bold", colour = ink) +
  geom_text(aes(x = 0.5, y = y - 0.04, label = paste0("highest in: ", chartr("\u0131\u0130\u015f\u015e\u00e7\u00c7\u011f\u011e\u00f6\u00d6\u00fc\u00dc", "iIsScCgGoOuU", highest_in))), hjust = 0, vjust = 1, size = 2.7, colour = ink2) +
  annotate("text", x = -0.3, y = -0.55, label = note, hjust = 0, vjust = 1, size = 2.8, colour = ink2, lineheight = 1.05) +
  scale_x_continuous(limits = c(-0.35, 10.6)) + scale_y_continuous(limits = c(-3.3, 8.7)) +
  labs(title = "ADMIXTURE colour legend: modern Iranians", subtitle = "K=8 components (Set1 colours as in the bar plot)") +
  theme_void(base_size = 11) +
  theme(plot.background = element_rect(fill = surface, colour = NA), plot.margin = margin(10, 14, 8, 10),
        plot.title = element_text(face = "bold", size = 13, colour = ink), plot.subtitle = element_text(colour = ink2))
ggsave(file.path(out, "admixture_legend.pdf"), p, width = 10.5, height = 5.8)
ggsave(file.path(out, "admixture_legend.png"), p, width = 10.5, height = 5.8, dpi = 300)
message("written to ", out)
