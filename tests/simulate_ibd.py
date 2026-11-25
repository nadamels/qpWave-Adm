#!/usr/bin/env python3
"""Fake IBD dataset to check that the IBD pipeline (scripts/ibd_0*.R) runs and recovers known proportions.

  python3 tests/simulate_ibd.py OUTDIR
writes  OUTDIR/sample_meta.tsv   (same columns as fam_meta/sample_meta.tsv)
        OUTDIR/fake.final.with_cM.tsv.gz   (IBDseq-like segments + cM columns, with header)

Ancient: three "populations" A, B, C (A = two nearby sub-groups A1/A2 that share a bit less IBD with each other than within),
  D = an unrelated group.  Modern targets: Mod_Central_Iranian = 70% A + 30% B,  Mod_Zoroastrian = 100% C (+ nothing else),
  Mod_Turk = 50% B + 50% D.  IBD rate between two individuals = sum_ab q_i[a] q_j[b] L[a,b].
"""
import gzip, os, sys, random
import numpy as np

out = sys.argv[1]; os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(1); random.seed(1)
anc_groups = {"A1": 40, "A2": 40, "B": 60, "C": 50, "D": 40}
base = {"A1": "A", "A2": "A", "B": "B", "C": "C", "D": "D"}
ancestry_comp = ["A1", "A2", "B", "C", "D"]
L = np.full((5, 5), 0.25)
for i in range(5): L[i, i] = 3.0
L[0, 1] = L[1, 0] = 1.6
geo = {"A1": ("Armenia_BA", 40.0, 44.0, 4000), "A2": ("Armenia_BA", 40.5, 44.5, 3900), "B": ("Anatolia_N", 38.0, 32.0, 8000),
       "C": ("Steppe_BA", 50.0, 60.0, 4800), "D": ("Han_like", 35.0, 110.0, 3000)}
mod = {"Mod_Central_Iranian": ("Fars", 40, {"A1": .35, "A2": .35, "B": .3}), "Mod_Zoroastrian": ("Zoroastrian", 12, {"C": 1.0}),
       "Mod_Turk": ("Turk", 8, {"B": .5, "D": .5})}

samples = []     # (iid, q vector)
rows = []
def q_of(d):
    q = np.zeros(5)
    for k, v in d.items(): q[ancestry_comp.index(k)] = v
    return q
for g, n in anc_groups.items():
    pop, lat, lon, bp = geo[g]
    for k in range(n):
        iid = f"{g}_{k}"; samples.append((iid, q_of({g: 1.0})))
        rows.append(dict(IID=iid, pop=pop, sup="Caucuses", imp="Yes", lat=lat + rng.normal(0, .2), lon=lon + rng.normal(0, .2),
                         date=f"{bp+100+k}-{bp-100} calBCE", cl="NA"))
for g, (pop, n, d) in mod.items():
    for k in range(n):
        iid = f"{g}_{k}"; samples.append((iid, q_of(d)))
        rows.append(dict(IID=iid, pop=pop, sup="Iran", imp="No", lat="", lon="", date="", cl="Iran_cluster_1"))
hdr = ["IID", "Source Publication(s)", "Empirical cluster (Study Pop only)", "Population", "Super-population", "Sex", "Archeological Period",
       "Archeological Date", "Original Genome Coverage (ancient samples only)", "Undergone Imputation", "Latitude", "Longitude"]
with open(f"{out}/sample_meta.tsv", "w") as f:
    f.write("\t".join(hdr) + "\n")
    for r in rows:
        f.write("\t".join([r["IID"], "sim", r["cl"], r["pop"], r["sup"], "U", "", r["date"], "", r["imp"], str(r["lat"]), str(r["lon"])]) + "\n")

chr_cM = [280, 260, 220, 210, 200, 190, 180, 165, 160, 180, 155, 170, 125, 115, 120, 125, 125, 120, 105, 105, 60, 62]
cum = np.cumsum([0] + chr_cM)
n = len(samples); Q = np.array([q for _, q in samples])
rate = Q @ L @ Q.T
ids = [s for s, _ in samples]
with gzip.open(f"{out}/fake.final.with_cM.tsv.gz", "wt") as f:
    f.write("sample1\thap1\tsample2\thap2\tchromosome\tstart\tend\tLOD\tstart_cM\tend_cM\tlength_cM\n")
    for i in range(n):
        lam = rate[i, i + 1:] * 2.0
        for j, l in zip(range(i + 1, n), lam):
            k = rng.poisson(l)
            for _ in range(k):
                c = rng.choice(22, p=np.array(chr_cM) / sum(chr_cM)); ln = 2 + rng.exponential(2.5)
                s = rng.uniform(0, chr_cM[c] - ln - 0.01) if chr_cM[c] > ln + 0.1 else 0
                h1, h2 = rng.integers(1, 3, 2)
                f.write(f"{ids[i]}\t{h1}\t{ids[j]}\t{h2}\t{c+1}\t{int(s*1e6)}\t{int((s+ln)*1e6)}\t{rng.uniform(3,40):.2f}\t{s:.3f}\t{s+ln:.3f}\t{ln:.3f}\n")
    # a few within-individual (ROH) segments
    for i in range(0, n, 7):
        f.write(f"{ids[i]}\t1\t{ids[i]}\t2\t3\t1000000\t6000000\t10.0\t10.000\t15.000\t5.000\n")
print("simulated", n, "individuals")
