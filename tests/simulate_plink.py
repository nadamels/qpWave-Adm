#!/usr/bin/env python3
"""Simulate a small diploid PLINK dataset with the same group labels as the real analysis
(for smoke-testing the pipeline only - the numbers mean nothing biologically).

usage: simulate_plink.py <sample_groups.tsv from make_pops.py> <out_prefix> [n_snps]
Writes <out_prefix>.bed/.bim/.fam  (FID = group, as the pipeline expects).
"""
import sys, numpy as np

groups_tsv, out, nsnp = sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 120000
rng = np.random.default_rng(7)
rows = [l.rstrip("\n").split("\t") for l in open(groups_tsv)][1:]
iid = [r[0] for r in rows]; grp = [r[1] for r in rows]

def drift(p, F):
    F = max(F, 1e-4)
    a = np.clip(p, 1e-3, 1 - 1e-3) * (1 - F) / F
    b = (1 - np.clip(p, 1e-3, 1 - 1e-3)) * (1 - F) / F
    return rng.beta(a, b)

base = rng.beta(0.4, 0.4, nsnp)
nonafr = drift(base, .05); afr = drift(base, .10)
ea = drift(nonafr, .08); weur = drift(nonafr, .04)
whg = drift(weur, .05); eehg = drift(weur, .05); irann = drift(weur, .05)
anat = 0.8 * drift(weur, .03) + 0.2 * whg
lev = drift(weur, .03)
steppe = 0.6 * eehg + 0.4 * irann
mix = lambda *a: sum(w * c for w, c in a)
comp = {
    "Mbuti": drift(afr, .01), "Papuan": drift(ea, .08), "East_Asia_Han": drift(ea, .02), "Karitiana": drift(ea, .10),
    "Onge": drift(ea, .12), "Chukchi": drift(ea, .05), "WHG": whg, "EEHG": eehg,
    "Iran_N": drift(irann, .01), "Iran_ChL": drift(mix((.7, irann), (.3, anat)), .01),
    "Iran_BA_SE": drift(mix((.6, irann), (.2, anat), (.2, ea)), .01),
    "Iran_BIA_NW": drift(mix((.55, irann), (.25, anat), (.2, steppe)), .01),
    "Anatolia_N": drift(anat, .01), "Anatolia_ChL_BA": drift(mix((.8, anat), (.2, irann)), .01),
    "Levant_N": drift(lev, .01), "Natufian": drift(lev, .02), "Armenia_N": drift(mix((.5, irann), (.5, anat)), .01),
    "Armenia_EBA": drift(mix((.5, irann), (.4, anat), (.1, steppe)), .01),
    "Armenia_LBA_IA": drift(mix((.45, irann), (.3, anat), (.25, steppe)), .01),
    "Steppe_EMBA": drift(steppe, .01), "Steppe_MLBA": drift(mix((.8, steppe), (.2, anat)), .01),
    "Turan_N": drift(mix((.8, irann), (.2, ea)), .01), "BMAC": drift(mix((.6, irann), (.2, anat), (.2, ea)), .01),
    "South_Asia_BA": drift(mix((.6, irann), (.4, ea)), .01),
    "Mod_Zoroastrian": mix((.80, irann), (.10, steppe), (.10, anat)),
    "Mod_Central_Iranian": mix((.55, irann), (.25, anat), (.20, steppe)),
    "Mod_Jewish_branch": mix((.35, irann), (.35, lev), (.30, anat)),
    "Mod_Armenian_Assyrian": mix((.30, irann), (.45, anat), (.25, steppe)),
    "Mod_Turk": mix((.40, irann), (.30, anat), (.20, steppe), (.10, ea)),
}
freq = {g: np.clip(drift(v, .005) if g.startswith("Mod_") else v, 0, 1) for g, v in comp.items()}
missing = sorted(set(grp) - set(freq))
if missing: sys.exit(f"no simulation model for groups: {missing}")

G = np.empty((nsnp, len(iid)), dtype=np.uint8)
for j, g in enumerate(grp):
    G[:, j] = rng.binomial(2, freq[g])
code = np.array([0b00, 0b10, 0b11], dtype=np.uint8)[G]      # A2 count 0/1/2 -> plink 2-bit code
pad = (-code.shape[1]) % 4
if pad: code = np.hstack([code, np.zeros((nsnp, pad), np.uint8)])
b = code[:, 0::4] | (code[:, 1::4] << 2) | (code[:, 2::4] << 4) | (code[:, 3::4] << 6)
with open(out + ".bed", "wb") as f:
    f.write(bytes([0x6C, 0x1B, 0x01])); f.write(b.astype(np.uint8).tobytes())
chrom = np.sort(rng.integers(1, 23, nsnp)); pos = np.zeros(nsnp, int)
for c in range(1, 23):
    m = chrom == c; pos[m] = np.sort(rng.choice(150_000_000, m.sum(), replace=False))
with open(out + ".bim", "w") as f:
    for i in range(nsnp): f.write(f"{chrom[i]}\trs{i}\t0\t{pos[i]}\tA\tG\n")
with open(out + ".fam", "w") as f:
    for i, g in zip(iid, grp): f.write(f"{g} {i} 0 0 0 -9\n")
print(f"wrote {out}.bed/.bim/.fam : {len(iid)} samples x {nsnp} SNPs")
