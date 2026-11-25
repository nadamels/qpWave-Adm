#!/usr/bin/env python3
"""Build the population-labelled PLINK .fam used by ADMIXTOOLS 2 (FID = population).

Inputs
  --meta   fam_meta/sample_meta.tsv        (IID, Population, Super-population, ... )
  --fam    <prefix>.fam                    (original fam; only IIDs/order are used)
  --pools  config/pop_pools.tsv            (group, role, population_label)
  --outdir <dir>                           (writes plink/data.fam, sample_groups.tsv, group_summary.tsv)

With --targets config, groups whose role in the pools file is "target" (population labels) are the targets, e.g. ancient Iranians.
Target (modern Iranian) groups follow the hierarchical-clustering result (n = 68 modern Iranians):
  Mod_Zoroastrian        Population == Zoroastrian
  Mod_Turk               Population == Turk               (the single Turk sample; NOT Turk_Azeri)
  Mod_Jewish_branch      Population == Jewish + IGP_57 + IGP_19 + IGP_25
  Mod_Armenian_Assyrian  Population in Armenian / Assyrian / Assyrian_Armenian
  Mod_Central_Iranian    every other modern Iranian
"""
import argparse, csv, os, sys
from collections import Counter, defaultdict

JEWISH_BRANCH_EXTRA = {"IGP_57", "IGP_19", "IGP_25"}
ARM_ASSYR = {"Armenian", "Assyrian", "Assyrian_Armenian"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--fam", required=True)
    ap.add_argument("--pools", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--targets", choices=["modern", "config", "both"], default="modern",
                    help="modern: the 68 modern Iranians in five groups (default); config: groups with role 'target' in the pools file "
                         "(e.g. ancient Iranians); both: union")
    ap.add_argument("--sublabels", action="store_true",
                    help="also write plink_sub/data.fam + sample_subgroups.tsv: target/pool groups split into sub-labels "
                         "(published label x sampling region) for the homogeneity test")
    ap.add_argument("--keep-ignored", action="store_true",
                    help="also label groups whose role is 'ignore' (used only by explore_admixture_pca.py)")
    a = ap.parse_args()

    with open(a.meta, newline="", encoding="utf-8") as f:
        meta = {r["IID"]: r for r in csv.DictReader(f, delimiter="\t")}

    fam = [l.split() for l in open(a.fam, encoding="utf-8") if l.strip()]
    missing = [r[1] for r in fam if r[1] not in meta]
    if missing:
        sys.exit(f"ERROR: {len(missing)} fam IIDs not found in {a.meta}, e.g. {missing[:5]}")

    # ---- targets: modern Iranians with an Iran_cluster_* label (= the 68 in the dendrogram)
    target = {}
    for iid, r in (meta.items() if a.targets in ("modern", "both") else []):
        if (r["Super-population"] == "Iran" and r["Undergone Imputation"] == "No"
                and r["Empirical cluster (Study Pop only)"].startswith("Iran_cluster")):
            pop = r["Population"]
            if pop == "Zoroastrian":
                g = "Mod_Zoroastrian"
            elif pop == "Turk":
                g = "Mod_Turk"
            elif pop in ARM_ASSYR:
                g = "Mod_Armenian_Assyrian"
            elif pop == "Jewish" or iid in JEWISH_BRANCH_EXTRA:
                g = "Mod_Jewish_branch"
            else:
                g = "Mod_Central_Iranian"
            target[iid] = g
    if a.targets in ("modern", "both"):
        if len(target) != 68:
            print(f"WARNING: expected 68 modern Iranians, found {len(target)}", file=sys.stderr)
        for need in JEWISH_BRANCH_EXTRA:
            if need not in target:
                sys.exit(f"ERROR: {need} (Jewish-branch sample) is not among the modern Iranian targets")

    # ---- pools / right_base from the pools file (exact Population labels)
    label2group, group_role, group_order = {}, {}, []
    with open(a.pools, encoding="utf-8") as f:
        rows = [l.rstrip("\n") for l in f if l.strip() and not l.startswith("#")]
    rd = csv.DictReader(rows, delimiter="\t")
    for r in rd:
        g, role, lab = r["group"].strip(), r["role"].strip(), r["population_label"].strip()
        if role not in ("target", "pool", "right_base", "ignore"):
            sys.exit(f"ERROR: bad role '{role}' for group {g}")
        if g not in group_order:
            group_order.append(g)
        if group_role.setdefault(g, role) != role:
            sys.exit(f"ERROR: group {g} has more than one role")
        if lab in label2group and label2group[lab] != g:
            sys.exit(f"ERROR: label {lab} assigned to two groups")
        label2group[lab] = g
    known = {r["Population"] for r in meta.values()}
    bad = sorted(l for l in label2group if l not in known)
    if bad:
        sys.exit("ERROR: population labels not found in metadata:\n  " + "\n  ".join(bad))

    # ---- assign every fam sample
    assign, role_of = {}, {}
    for r in fam:
        iid = r[1]
        if iid in target:
            assign[iid] = target[iid]; role_of[assign[iid]] = "target"
        elif meta[iid]["Super-population"] != "Iran" and meta[iid]["Population"] in label2group:
            g = label2group[meta[iid]["Population"]]
            if a.keep_ignored or group_role[g] != "ignore":
                assign[iid] = g; role_of[g] = group_role[g]
        elif meta[iid]["Super-population"] == "Iran" and meta[iid]["Undergone Imputation"] == "Yes" \
                and meta[iid]["Population"] in label2group:   # ancient Iranians live in super-pop 'Iran'
            g = label2group[meta[iid]["Population"]]
            if a.keep_ignored or group_role[g] != "ignore":
                assign[iid] = g; role_of[g] = group_role[g]

    os.makedirs(os.path.join(a.outdir, "plink"), exist_ok=True)
    with open(os.path.join(a.outdir, "plink", "data.fam"), "w") as out:
        for r in fam:
            out.write(" ".join([assign.get(r[1], "unused")] + r[1:]) + "\n")

    with open(os.path.join(a.outdir, "sample_groups.tsv"), "w") as out:
        out.write("IID\tgroup\trole\tpopulation_label\n")
        for r in fam:
            g = assign.get(r[1])
            if g:
                out.write(f"{r[1]}\t{g}\t{role_of[g]}\t{meta[r[1]]['Population']}\n")

    if a.sublabels:
        import re
        san = lambda x: re.sub(r"[^A-Za-z0-9]+", "_", x).strip("_")
        os.makedirs(os.path.join(a.outdir, "plink_sub"), exist_ok=True)
        with open(os.path.join(a.outdir, "plink_sub", "data.fam"), "w") as out, \
             open(os.path.join(a.outdir, "sample_subgroups.tsv"), "w") as tab:
            tab.write("IID\tgroup\trole\tsublabel\tpopulation_label\tsuper_population\n")
            for r in fam:
                g = assign.get(r[1])
                if g is None or role_of[g] == "ignore":
                    fid = "unused"
                elif role_of[g] in ("target", "pool"):
                    fid = f"{g}__{san(meta[r[1]]['Population'])}__{san(meta[r[1]]['Super-population'])}"
                    tab.write(f"{r[1]}\t{g}\t{role_of[g]}\t{fid}\t{meta[r[1]]['Population']}\t{meta[r[1]]['Super-population']}\n")
                else:
                    fid = g
                out.write(" ".join([fid] + r[1:]) + "\n")
    n = Counter(assign.values())
    TGT_ORDER = ["Mod_Zoroastrian", "Mod_Turk", "Mod_Jewish_branch", "Mod_Armenian_Assyrian", "Mod_Central_Iranian"]

    def sort_key(g):
        r = ["target", "pool", "right_base", "ignore"].index(role_of[g])
        pos = TGT_ORDER.index(g) if g in TGT_ORDER else group_order.index(g)
        return (r, pos)

    configured = {g for g, ro in group_role.items() if a.keep_ignored or ro != "ignore"}
    with open(os.path.join(a.outdir, "group_summary.tsv"), "w") as out:
        out.write("group\trole\tn\n")
        for g in sorted(n, key=sort_key):
            out.write(f"{g}\t{role_of[g]}\t{n[g]}\n")
    empty = sorted(configured - set(n))
    print("group\trole\tn")
    for g in sorted(n, key=sort_key):
        print(f"{g}\t{role_of[g]}\t{n[g]}")
    if empty:
        print("WARNING: configured groups with 0 samples in this fam (dropped): " + ", ".join(empty), file=sys.stderr)


if __name__ == "__main__":
    main()
