#!/usr/bin/env python3
"""Sample table for the IBD pipeline, built from fam_meta/sample_meta.tsv.

Output columns (tab separated)
  IID, is_ancient, population, super_population, date_mid_BP, lat, lon, modern_group

  is_ancient    TRUE = imputed ancient genome ("Undergone Imputation" == Yes)
  date_mid_BP   midpoint of "Archeological Date" in years BP (1950 = 0); NA if not parseable
  modern_group  target / donor label of a MODERN individual (ancient individuals get NA; their group is the IBD cluster):
                  modern Iranians with an Iran_cluster_* label -> Mod_Zoroastrian / Mod_Turk / Mod_Jewish_branch /
                  Mod_Armenian_Assyrian / Mod_Central_Iranian   (same rules as make_pops.py, 68 samples)
                  the 2 modern Iranians without a cluster label -> NA (ignored, as in the qpAdm pipeline)
                  every other modern sample -> its Population label (modern donor group)
"""
import argparse, csv, re, sys

JEWISH_BRANCH_EXTRA = {"IGP_57", "IGP_19", "IGP_25"}
ARM_ASSYR = {"Armenian", "Assyrian", "Assyrian_Armenian"}


def parse_date_bp(s):
    """'5898-5531 calBCE', '1100-1300 CE', '776-545 calBCE (2510±20 BP, ...)', '800 BCE - 200 CE' -> midpoint in years BP."""
    if not s:
        return None
    s = re.sub(r"\([^)]*\)|\[[^\]]*\]", " ", s)          # drop lab codes / remarks
    toks = re.findall(r"(\d+)\s*(?:cal)?\s*(BCE|CE)?", s)
    if not toks:
        return None
    years, era = [], None
    for num, e in reversed(toks):                        # an era applies to the numbers before it ("5898-5531 calBCE")
        era = e or era
        if era is None:
            return None
        n = int(num)
        years.append(-(n - 1) if era == "BCE" else n)    # no year 0
    return 1950 - sum(years) / len(years)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    n_anc = n_mod = 0
    with open(a.meta, newline="", encoding="utf-8") as f, open(a.out, "w") as out:
        out.write("IID\tis_ancient\tpopulation\tsuper_population\tdate_mid_BP\tlat\tlon\tmodern_group\n")
        for r in csv.DictReader(f, delimiter="\t"):
            iid, pop, sup = r["IID"], r["Population"], r["Super-population"]
            anc = r["Undergone Imputation"] == "Yes"
            grp = "NA"
            if not anc:
                if sup == "Iran":
                    if r["Empirical cluster (Study Pop only)"].startswith("Iran_cluster"):
                        if pop == "Zoroastrian":
                            grp = "Mod_Zoroastrian"
                        elif pop == "Turk":
                            grp = "Mod_Turk"
                        elif pop in ARM_ASSYR:
                            grp = "Mod_Armenian_Assyrian"
                        elif pop == "Jewish" or iid in JEWISH_BRANCH_EXTRA:
                            grp = "Mod_Jewish_branch"
                        else:
                            grp = "Mod_Central_Iranian"
                else:
                    grp = pop
            d = parse_date_bp(r["Archeological Date"]) if anc else None
            out.write("\t".join([iid, "TRUE" if anc else "FALSE", pop, sup, "NA" if d is None else f"{d:.0f}",
                                 r["Latitude"] or "NA", r["Longitude"] or "NA", grp]) + "\n")
            n_anc += anc; n_mod += (not anc)
    print(f"wrote {a.out}: {n_anc} ancient, {n_mod} modern samples", file=sys.stderr)


if __name__ == "__main__":
    main()
