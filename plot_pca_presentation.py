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

COL = {  # exact colours of the interactive PCA viewer
    "Western-Europe": "#c000db", "Northern-Europe": "#f0b5ff",
    "Eastern-Europe": "#d200e5", "Southern-Europe": "#1b002f",
    "Caucasus": "#d4effc", "Anatolia": "#8c564b",
    "Levant": "#4a6795", "Arabia": "#4a6795",
    "Central-Asia": "#E1AD01", "South-Asia": "#fb8400",
    "North-Africa": "#960016", "Africa": "#960016", "Iran": "#ff0000",
}
OTHER = "#d4d7dd"  # Africa, Americas, East/SE Asia, Oceania (context only)
# region -> (text, label x, label y, ha, colour, arrow target (x, y) or None = median of points in view)
LABEL = {
    "Western-Europe": ("Western Europe", -0.0040, 0.0200, "right", "#c000db", (-0.0004, 0.0182)),
    "Southern-Europe": ("Southern Europe", 0.0040, 0.0200, "left", "#1b002f", (0.0006, 0.0176)),
    "Anatolia": ("Anatolia", 0.0048, 0.0186, "left", "#8c564b", (0.0012, 0.0168)),
    "Caucasus": ("Caucasus", 0.0048, 0.0172, "left", "#5b9bbd", (0.0002, 0.0155)),
    "Northern-Europe": ("Northern\nEurope", -0.0084, 0.0166, "left", "#c77fe0", (-0.0017, 0.0150)),
    "Eastern-Europe": ("Eastern\nEurope", -0.0084, 0.0130, "left", "#d200e5", (-0.0032, 0.0121)),
    "Central-Asia": ("Central Asia", -0.0092, 0.0089, "left", "#a67c00", (-0.0040, 0.0089)),
    "South-Asia": ("South Asia", -0.0020, 0.0042, "left", "#d96f00", (-0.0031, 0.0058)),
    "Levant": ("Levant", 0.0060, 0.0154, "left", "#4a6795", (0.0032, 0.0140)),
    "Arabia": ("Arabian\nPeninsula", 0.0046, 0.0090, "right", "#4a6795", (0.0060, 0.0114)),
    "North-Africa": ("North Africa", 0.0098, 0.0074, "left", "#960016", (0.0085, 0.0095)),
}
# Modern Iranian groups (marker shapes follow the original PCA viewer). Grouping rules:
#   Zoroastrian            all samples labelled Zoroastrian
#   Assyrian / Armenian    Assyrian, Armenian and Assyrian_Armenian
#   Turk                   the single Turk sample
#   Jewish-related branch  samples labelled Jewish + IGP_19, IGP_25, IGP_57 (they sit in the same branch)
#   Central Iranian        everyone else (incl. the 2 SGDP "Iranian" reference individuals)
JEWISH_EXTRA = {"IGP_19", "IGP_25", "IGP_57"}


def iran_group(row):
    pop = row["Population"]
    if pop == "Zoroastrian":
        return "Zoroastrian"
    if pop in ("Armenian", "Assyrian", "Assyrian_Armenian"):
        return "Assyrian / Armenian"
    if pop == "Turk":
        return "Turk"
    if pop == "Jewish" or row["IID"] in JEWISH_EXTRA:
        return "Jewish-related branch"
    return "Central Iranian"


IRAN_GROUPS = [  # (group, matplotlib marker, marker area); drawn in this order
    ("Central Iranian", "v", 120),
    ("Jewish-related branch", "X", 130),
    ("Zoroastrian", "^", 125),
    ("Assyrian / Armenian", "s", 105),
    ("Turk", "D", 120),
]
d.loc[~d.anc & (d.sp == "Iran"), "igroup"] = d[~d.anc & (d.sp == "Iran")].apply(iran_group, axis=1)
XLIM, YLIM = (-0.0098, 0.0128), (0.0033, 0.0208)


# Ancient samples: marker area (pt^2) by age tier -- older = much bigger.
# Cut-offs (years BP) follow the archaeological period labels in the sample table
# (median age: Neolithic ~7,500, Chalcolithic ~5,950, Early/Middle Bronze ~3,800-4,850, Iron Age ~2,900).
AGE_TIERS = [  # (min age BP, legend text, marker area)
    (7000, "Neolithic & older  (> 7,000 BP)", 640),
    (5200, "Chalcolithic  (5,200\u20137,000)", 300),
    (3500, "Bronze Age  (3,500\u20135,200)", 150),
    (2500, "Iron Age  (2,500\u20133,500)", 80),
    (0, "Historical  (< 2,500 BP)", 42),
]


def size_from_age(a):
    a = a.fillna(3000)  # 4 undated samples -> Iron Age size
    out = pd.Series(AGE_TIERS[-1][2], index=a.index, dtype=float)
    for lo, _, sz in reversed(AGE_TIERS):
        out[a >= lo] = sz
    return out


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
    ax.scatter(a[pcx], a[pcy], s=size_from_age(a.age), c=a.sp.map(COL), alpha=.9,
               edgecolor="black", lw=.7, zorder=3)
    # ancient Iran: red circles, sized by age like all ancients
    ai = anc[anc.sp == "Iran"].sort_values("age", ascending=False, na_position="first")
    ax.scatter(ai[pcx], ai[pcy], s=size_from_age(ai.age), c="#ff0000", alpha=.95,
               edgecolor="black", lw=.8, zorder=4)
    # modern Iranians
    mi = mod[mod.sp == "Iran"]
    for grp, mk, sz in IRAN_GROUPS:
        g = mi[mi.igroup == grp]
        ax.scatter(g[pcx], g[pcy], marker=mk, s=sz, c="#ff0000", edgecolor="#000000", lw=1,
                   zorder=6)


def finish(ax, labels=True):
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_xlabel("PC1", fontsize=13, labelpad=8)
    ax.set_ylabel("PC2", fontsize=13, labelpad=8)
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
                fontweight="bold", color="#e00000", zorder=10)


def legends(fig, ax):
    # modern Iranian clusters
    mi = d[(~d.anc) & (d.sp == "Iran")]
    h = [Line2D([], [], ls="", marker=mk, ms=10, mfc="#ff0000", mec="#000000", mew=1,
                label=f"{grp}  (n={int((mi.igroup == grp).sum())})") for grp, mk, _ in IRAN_GROUPS]
    l1 = ax.legend(handles=h, title="Modern Iranians", loc="upper right", frameon=False,
                   ncol=1, fontsize=12, title_fontsize=12.5, labelcolor=INK,
                   bbox_to_anchor=(1.0, 1.0), handletextpad=.4, labelspacing=.5)
    l1._legend_box.align = "left"; l1.get_title().set_color(MUTED)
    ax.add_artist(l1)
    # ancient vs modern + size key
    h2 = [Line2D([], [], ls="", marker="o", ms=np.sqrt(sz), mfc="#d9d9de", mec="black", mew=.7,
                  label=txt) for _, txt, sz in AGE_TIERS]
    h2 += [Line2D([], [], ls="", marker="o", ms=10, mfc="#ff0000", mec="black", mew=.8,
                  label="Ancient Iranian (same size scale)")]
    l2 = ax.legend(handles=h2, title="Ancient samples (size = age)", loc="upper right",
                   frameon=False, fontsize=11.5, title_fontsize=12.5, labelcolor=INK,
                   bbox_to_anchor=(1.0, 0.74), handletextpad=.6, labelspacing=.45,
                   borderaxespad=.4)
    l2._legend_box.align = "left"; l2.get_title().set_color(MUTED)
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
