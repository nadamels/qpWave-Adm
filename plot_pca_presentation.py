#!/usr/bin/env python3
"""Presentation-quality PCA plot (PC1 vs PC2) of ancient + modern samples.

Reads  fam_meta/all_anc_mod.run3.LDpruned.pca.{eigenvec,eigenval}
       fam_meta/mega_sample_table.xlsx
Writes figures/pca_presentation.{png,svg,pdf}

Encoding
  colour  = region (super-population); Iran is red
  small translucent dots = modern reference individuals
  large outlined circles = ancient individuals (bigger = older)
  red markers (circle/square/diamond/triangle/star) = modern Iranians by cluster
"""
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
META = ROOT / "fam_meta"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
STEM = "all_anc_mod.run3.LDpruned.pca"

# ----------------------------------------------------------------- data
ev = pd.read_csv(META / f"{STEM}.eigenvec", sep=r"\s+", header=None, dtype={1: str})
ev.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, ev.shape[1] - 1)]
val = np.loadtxt(META / f"{STEM}.eigenval")
# Only the top-20 eigenvalues are stored, so normalise by the total variance
# (= number of samples for standardised genotypes): 313.36/5348 = 5.86 %, 241.46/5348 = 4.51 %.
var = val / len(ev) * 100
tab = pd.read_excel(META / "mega_sample_table.xlsx", sheet_name="mega_table", dtype={"IID": str})
d = ev.merge(tab, on="IID", how="left")
d["anc"] = d["Undergone Imputation"].eq("Yes")
d["sp"] = d["Super-population"].replace({"Caucuses": "Caucasus"})
d["cluster"] = d["Empirical cluster (Study Pop only)"]


def age_bp(s):
    """Midpoint of a date range such as '5898-5531 calBCE (..)' -> years BP."""
    m = re.match(r"\s*(\d+)\s*(BCE|CE)?\s*-\s*(\d+)\s*(?:cal)?(BCE|CE)", str(s))
    if not m:
        return np.nan
    a, ea, b, eb = int(m[1]), m[2] or m[4], int(m[3]), m[4]
    ya = -a if ea == "BCE" else a
    yb = -b if eb == "BCE" else b
    return 1950 - (ya + yb) / 2


d["age"] = d["Archeological Date"].map(age_bp)

# ----------------------------------------------------------------- style
INK, MUTED = "#1f2430", "#6b7280"
plt.rcParams.update({
    "font.family": "Inter", "font.size": 13, "text.color": INK,
    "axes.edgecolor": "#c9ccd3", "axes.labelcolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

COL = {  # fixed, entity-based colours, matched to the original interactive PCA viewer
    "Western-Europe": "#c80fd4", "Northern-Europe": "#eaa8f0",
    "Eastern-Europe": "#e040e6", "Southern-Europe": "#1f0b3d",
    "Caucasus": "#a9d4ec", "Anatolia": "#8d5b4c",
    "Levant": "#6d7dc4", "Arabia": "#4b5b9e",
    "Central-Asia": "#d9a900", "South-Asia": "#f58a1f",
    "North-Africa": "#a8283a", "Africa": "#8e1b2a", "Iran": "#ff1a1a",
}
OTHER = "#d4d7dd"  # Africa, Americas, East/SE Asia, Oceania (context only)
# region -> (text, label x, label y, ha, colour, arrow target (x, y) or None = median of points in view)
LABEL = {
    "Western-Europe": ("Western Europe", -0.0040, 0.0200, "right", "#a800b8", (-0.0004, 0.0182)),
    "Southern-Europe": ("Southern Europe", 0.0040, 0.0200, "left", "#1f0b3d", (0.0006, 0.0176)),
    "Anatolia": ("Anatolia", 0.0048, 0.0186, "left", "#8d5b4c", (0.0012, 0.0168)),
    "Caucasus": ("Caucasus", 0.0048, 0.0172, "left", "#4a95bd", (0.0002, 0.0155)),
    "Northern-Europe": ("Northern\nEurope", -0.0084, 0.0166, "left", "#c070d0", (-0.0017, 0.0150)),
    "Eastern-Europe": ("Eastern\nEurope", -0.0084, 0.0130, "left", "#b010c0", (-0.0032, 0.0121)),
    "Central-Asia": ("Central Asia", -0.0092, 0.0089, "left", "#a67c00", (-0.0040, 0.0089)),
    "South-Asia": ("South Asia", -0.0020, 0.0042, "left", "#c2610c", (-0.0031, 0.0058)),
    "Levant": ("Levant", 0.0060, 0.0154, "left", "#4d5fb0", (0.0032, 0.0140)),
    "Arabia": ("Arabian\nPeninsula", 0.0088, 0.0128, "left", "#3b4a8c", (0.0060, 0.0114)),
    "North-Africa": ("North Africa", 0.0098, 0.0074, "left", "#a8283a", (0.0085, 0.0095)),
}
# Modern Iranians use the marker shapes of the original PCA viewer (v, X, D, P).
# The viewer's own "Iranian group" assignment is not in the repo, so shapes are mapped onto the
# empirical clusters of mega_sample_table.xlsx: X = cluster 0 (upper right, next to the Caucasus),
# v = the large central group (clusters 1, 3, 4), D = cluster 2, P = no cluster assigned.
IRAN_GROUPS = [  # (matplotlib marker, legend text, clusters or None = unassigned, marker area)
    ("v", "Clusters 1, 3, 4", ["Iran_cluster_1", "Iran_cluster_3", "Iran_cluster_4"], 120),
    ("X", "Cluster 0", ["Iran_cluster_0"], 130),
    ("D", "Cluster 2", ["Iran_cluster_2"], 100),
    ("P", "Unassigned", None, 120),
]
XLIM, YLIM = (-0.0098, 0.0128), (0.0033, 0.0208)


def size_from_age(a):  # marker area (pt^2): older = bigger
    return 34 + 120 * np.clip(a.fillna(6000), 0, 12000) / 12000


PCX, PCY = "PC1", "PC2"


def draw(ax, pcx=PCX, pcy=PCY):
    mod, anc = d[~d.anc], d[d.anc]
    # context: groups we don't label
    ctx = mod[~mod.sp.isin(COL)]
    ax.scatter(ctx[pcx], ctx[pcy], s=8, c=OTHER, alpha=.6, lw=0, zorder=1, rasterized=True)
    # modern reference, small translucent dots
    for sp, g in mod[mod.sp.isin(COL) & (mod.sp != "Iran")].groupby("sp"):
        ax.scatter(g[pcx], g[pcy], s=18, c=COL[sp], alpha=.7, lw=0, zorder=2, rasterized=True)
    # ancient: larger, outlined, size = age (draw oldest first so young sit on top)
    a = anc[anc.sp.isin(COL) & (anc.sp != "Iran")].sort_values("age", ascending=False, na_position="first")
    ax.scatter(a[pcx], a[pcy], s=size_from_age(a.age), c=a.sp.map(COL), alpha=.85,
               edgecolor="white", lw=.9, zorder=3)
    # ancient Iran: darker red, white outline (same outline rule as all ancients)
    ai = anc[anc.sp == "Iran"].sort_values("age", ascending=False, na_position="first")
    ax.scatter(ai[pcx], ai[pcy], s=size_from_age(ai.age), c="#d01010", alpha=.95,
               edgecolor="white", lw=.9, zorder=4)
    # modern Iranians
    mi = mod[mod.sp == "Iran"]
    for mk, _, cls, sz in IRAN_GROUPS:
        g = mi[mi.cluster.isin(cls)] if cls else mi[mi.cluster.isna()]
        ax.scatter(g[pcx], g[pcy], marker=mk, s=sz, c="#ff1a1a", edgecolor="#000000", lw=1,
                   zorder=6)


def finish(ax, labels=True):
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_xlabel(f"PC1  ({var[0]:.1f}% of variance)", fontsize=13, labelpad=8)
    ax.set_ylabel(f"PC2  ({var[1]:.1f}% of variance)", fontsize=13, labelpad=8)
    ax.set_xticks([]); ax.set_yticks([])  # PC loadings are unitless; keep the plot calm
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if labels:
        vis = d[(d[PCX].between(*XLIM)) & (d[PCY].between(*YLIM))]
        for sp, (txt, x, y, ha, c, tgt) in LABEL.items():
            if tgt is None:
                g = vis[(vis.sp == sp) & ~vis.anc]
                tgt = (g[PCX].median(), g[PCY].median())
            ax.annotate(txt, xy=tgt, xytext=(x, y), ha=ha, va="center", fontsize=15,
                        fontweight="semibold", color=c, zorder=10, linespacing=.95,
                        arrowprops=dict(arrowstyle="-", color=c, lw=1, alpha=.7,
                                        shrinkA=4, shrinkB=3, relpos=(1 if ha == "right" else 0, .5)))
        ax.text(0.0018, 0.0112, "Iranians", ha="left", va="center", fontsize=22,
                fontweight="bold", color="#d01010", zorder=10)


def legends(fig, ax):
    # modern Iranian clusters
    mi = d[(~d.anc) & (d.sp == "Iran")]
    h = []
    for mk, txt, cls, _ in IRAN_GROUPS:
        n = int((mi.cluster.isin(cls) if cls else mi.cluster.isna()).sum())
        h.append(Line2D([], [], ls="", marker=mk, ms=10, mfc="#ff1a1a", mec="#000000", mew=1,
                        label=f"{txt}  (n={n})"))
    l1 = ax.legend(handles=h, title="Modern Iranians", loc="upper right", frameon=False,
                   ncol=1, fontsize=12, title_fontsize=12.5, labelcolor=INK,
                   bbox_to_anchor=(1.0, 1.0), handletextpad=.4, labelspacing=.5)
    l1._legend_box.align = "left"; l1.get_title().set_color(MUTED)
    ax.add_artist(l1)
    # ancient vs modern + size key
    h2 = [Line2D([], [], ls="", marker="o", ms=4, mfc="#9aa0aa", mec="none", label="Modern reference"),
          Line2D([], [], ls="", marker="o", ms=10, mfc="#9aa0aa", mec="white", mew=1,
                 label="Ancient (bigger = older)"),
          Line2D([], [], ls="", marker="o", ms=10, mfc="#d01010", mec="white", mew=1,
                 label="Ancient Iranian")]
    l2 = ax.legend(handles=h2, loc="upper right", frameon=False, fontsize=12, labelcolor=INK,
                   bbox_to_anchor=(1.0, 0.76), handletextpad=.4, labelspacing=.6)
    ax.add_artist(l2)


def main():
    fig, ax = plt.subplots(figsize=(16, 9), dpi=200)
    fig.patch.set_facecolor("white")
    draw(ax); finish(ax); legends(fig, ax)
    fig.subplots_adjust(left=.06, right=.97, top=.97, bottom=.09)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"pca_presentation.{ext}", facecolor="white")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
