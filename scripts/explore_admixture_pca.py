#!/usr/bin/env python3
"""Hypothesis generation from unsupervised ADMIXTURE + PCA, BEFORE qpWave / qpAdm.

What this is for
  * label the ADMIXTURE components (which groups carry each one),
  * show which components each modern Iranian target group carries,
  * check whether the candidate source pool (config/pop_pools.tsv) *covers* those components,
  * check whether the fixed RIGHT set spans the components (qpAdm needs differential relatedness),
  * rank candidate sources by distance to each target in ADMIXTURE and PCA space.

What this is NOT
  * It is not evidence for any admixture model. ADMIXTURE components are drifted clusters, not
    ancestral populations: similar Q profiles do not prove a shared source, and K is arbitrary.
  * Its numbers never enter qpWave/qpAdm. The pool is FROZEN (config/pop_pools.tsv) before qpAdm
    is run, and only qpAdm p-values / weights / SEs are reported as results.

usage:
  explore_admixture_pca.py --meta fam_meta/sample_meta.tsv --fam <fam> --pools config/pop_pools.tsv \
      --q <file>.Q --eigenvec <file>.eigenvec --outdir <dir> [--npcs 4]
Row order of the .Q file must equal the .fam order (ADMIXTURE keeps the input order).
"""
import argparse, csv, os, subprocess, sys
from collections import defaultdict
import numpy as np

COMP_FRAC_MIN = 0.05      # a target "carries" a component if its mean Q is at least this
CARRIER_MIN = 0.25        # a group "carries" a component if the component is >= this share of its mean profile
REDUNDANT_Q = 0.08        # candidate pairs closer than this (Euclid., Q space) are likely to be merged by the qpWave screen
SD_FLAG = 0.10            # within-group SD of a component above this = heterogeneous group


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True); ap.add_argument("--fam", required=True)
    ap.add_argument("--pools", required=True); ap.add_argument("--q", required=True)
    ap.add_argument("--eigenvec", required=True); ap.add_argument("--outdir", required=True)
    ap.add_argument("--npcs", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    here = os.path.dirname(os.path.abspath(__file__))
    tmp = os.path.join(a.outdir, "_groups")
    subprocess.run([sys.executable, "-I", os.path.join(here, "make_pops.py"), "--meta", a.meta, "--fam", a.fam,
                    "--pools", a.pools, "--outdir", tmp, "--keep-ignored"], check=True, stdout=subprocess.DEVNULL)

    fam = [l.split()[1] for l in open(a.fam, encoding="utf-8") if l.strip()]
    Q = np.loadtxt(a.q)
    ev = [l.split() for l in open(a.eigenvec)]
    if Q.shape[0] != len(fam) or len(ev) != len(fam):
        sys.exit(f"ERROR: rows differ: fam {len(fam)}, Q {Q.shape[0]}, eigenvec {len(ev)}")
    if [e[1] for e in ev] != fam:
        sys.exit("ERROR: eigenvec IIDs are not in fam order")
    K = Q.shape[1]
    PC = np.array([[float(x) for x in e[2:2 + a.npcs]] for e in ev])
    row = {iid: i for i, iid in enumerate(fam)}

    meta = {r["IID"]: r for r in csv.DictReader(open(a.meta, encoding="utf-8"), delimiter="\t")}
    grp, role = {}, {}
    for r in csv.DictReader(open(os.path.join(tmp, "sample_groups.tsv")), delimiter="\t"):
        grp[r["IID"]] = r["group"]; role[r["group"]] = r["role"]
    order = {}
    for g in role:
        order.setdefault(g, len(order))

    members = defaultdict(list)
    for iid, g in grp.items():
        members[g].append(row[iid])
    # reference groups (everyone else) by Population label, only to NAME the components
    refs = defaultdict(list)
    for iid in fam:
        if iid not in grp:
            refs["ref:" + meta[iid]["Population"]].append(row[iid])
    refs = {k: v for k, v in refs.items() if len(v) >= 5}

    def mean_q(idx): return Q[idx].mean(0)

    # ---- 1. component names: the 4 reference/ancient groups with the highest mean of each component
    names = []
    allg = {**{g: v for g, v in members.items()}, **refs}
    for c in range(K):
        top = sorted(((Q[v, c].mean(), g) for g, v in allg.items() if len(v) >= 3), reverse=True)[:4]
        names.append(f"K{c+1}: " + ", ".join(g.replace("ref:", "") for _, g in top))
    with open(os.path.join(a.outdir, "admixture_component_labels.tsv"), "w") as f:
        f.write("component\ttop_groups_by_mean_fraction\n")
        for c in range(K):
            f.write(f"K{c+1}\t{names[c].split(': ', 1)[1]}\n")

    # ---- 2. group summary (mean Q, within-group SD, PC centroids)
    gl = sorted(members, key=lambda g: (["target", "pool", "right_base", "ignore"].index(role[g]), order[g]))
    with open(os.path.join(a.outdir, "admixture_pca_group_summary.tsv"), "w") as f:
        f.write("group\trole\tn\t" + "\t".join(f"Q{c+1}" for c in range(K)) + "\t"
                + "\t".join(f"sdQ{c+1}" for c in range(K)) + "\t" + "\t".join(f"PC{p+1}" for p in range(a.npcs)) + "\n")
        for g in gl:
            v = members[g]
            f.write("\t".join([g, role[g], str(len(v))] + [f"{x:.4f}" for x in Q[v].mean(0)]
                              + [f"{x:.4f}" for x in Q[v].std(0)] + [f"{x:.5f}" for x in PC[v].mean(0)]) + "\n")
    with open(os.path.join(a.outdir, "admixture_pca_samples.tsv"), "w") as f:
        f.write("IID\tgroup\trole\t" + "\t".join(f"Q{c+1}" for c in range(K)) + "\t"
                + "\t".join(f"PC{p+1}" for p in range(a.npcs)) + "\n")
        for g in gl:
            for i in members[g]:
                f.write("\t".join([fam[i], g, role[g]] + [f"{x:.4f}" for x in Q[i]] + [f"{x:.5f}" for x in PC[i]]) + "\n")

    mq = {g: mean_q(members[g]) for g in gl}
    mpc = {g: PC[members[g]].mean(0) for g in gl}
    tg = [g for g in gl if role[g] == "target"]
    pool = [g for g in gl if role[g] == "pool"]
    right = [g for g in gl if role[g] == "right_base"]
    cand = [g for g in gl if role[g] in ("pool", "ignore")]

    # ---- 3. component coverage per target
    carriers = lambda g, c: mq[g][c] / mq[g].sum() >= CARRIER_MIN
    rows = []
    for t in tg:
        sd = Q[members[t]].std(0) if len(members[t]) > 1 else np.full(K, np.nan)
        for c in range(K):
            if mq[t][c] < COMP_FRAC_MIN:
                continue
            pc = [g for g in pool if carriers(g, c)]
            ig = [g for g in gl if role[g] == "ignore" and carriers(g, c)]
            rc = [g for g in right if carriers(g, c)]
            status = "covered by pool" if pc else ("only in IGNORED groups" if ig else "NOT COVERED by any pool group")
            rows.append([t, f"K{c+1}", f"{mq[t][c]:.3f}", "" if np.isnan(sd[c]) else f"{sd[c]:.3f}", names[c].split(": ", 1)[1],
                         ",".join(pc), ",".join(ig), ",".join(rc), status])
    with open(os.path.join(a.outdir, "admixture_component_coverage.tsv"), "w") as f:
        f.write("target\tcomponent\ttarget_mean_fraction\ttarget_within_group_sd\tcomponent_top_groups\tpool_carriers\t"
                "ignored_carriers\tright_base_carriers\tstatus\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    # ---- 4. candidate ranking by distance (ADMIXTURE and PCA space)
    with open(os.path.join(a.outdir, "admixture_pca_candidate_screen.tsv"), "w") as f:
        f.write("target\tcandidate\trole\tn\tq_distance\tq_rank\tpc_distance\tpc_rank\n")
        for t in tg:
            qd = {g: float(np.linalg.norm(mq[t] - mq[g])) for g in cand}
            pd_ = {g: float(np.linalg.norm(mpc[t] - mpc[g])) for g in cand}
            qr = {g: i + 1 for i, g in enumerate(sorted(qd, key=qd.get))}
            pr = {g: i + 1 for i, g in enumerate(sorted(pd_, key=pd_.get))}
            for g in sorted(cand, key=lambda g: qr[g]):
                f.write(f"{t}\t{g}\t{role[g]}\t{len(members[g])}\t{qd[g]:.4f}\t{qr[g]}\t{pd_[g]:.5f}\t{pr[g]}\n")

    # ---- 5. redundancy among pool groups (hint only; the qpWave screen is what decides)
    red = []
    for i, g in enumerate(pool):
        for h in pool[i + 1:]:
            d = float(np.linalg.norm(mq[g] - mq[h]))
            if d < REDUNDANT_Q:
                red.append((g, h, d))
    with open(os.path.join(a.outdir, "admixture_pool_redundancy.tsv"), "w") as f:
        f.write("group_a\tgroup_b\tq_distance\n")
        for g, h, d in sorted(red, key=lambda x: x[2]):
            f.write(f"{g}\t{h}\t{d:.4f}\n")

    # ---- console summary
    print("\nADMIXTURE components (K=%d):" % K)
    for n in names: print("  " + n)
    print("\nMean ancestry per target (components >= %.0f%%):" % (100 * COMP_FRAC_MIN))
    for t in tg:
        parts = ", ".join(f"K{c+1}={mq[t][c]:.2f}" for c in np.argsort(-mq[t]) if mq[t][c] >= COMP_FRAC_MIN)
        print(f"  {t:<24} n={len(members[t]):<3} {parts}")
    print("\nCoverage problems (target carries a component that no pool group carries):")
    bad = [r for r in rows if r[-1] != "covered by pool"]
    for r in bad:
        print(f"  {r[0]:<24} {r[1]} ({float(r[2]):.2f})  {r[-1]}   ignored carriers: {r[6] or '-'}   right carriers: {r[7] or '-'}")
    if not bad: print("  none")
    print("\nComponents with no carrier in the fixed RIGHT set (covered only if a pool group rotates onto the right):")
    for c in range(K):
        if not any(carriers(g, c) for g in right):
            print(f"  K{c+1}  {names[c].split(': ', 1)[1]}")
    print("\nNearest candidates (ADMIXTURE | PCA), top 3 each:")
    for t in tg:
        qd = sorted(cand, key=lambda g: np.linalg.norm(mq[t] - mq[g]))[:3]
        pd_ = sorted(cand, key=lambda g: np.linalg.norm(mpc[t] - mpc[g]))[:3]
        print(f"  {t:<24} Q: {', '.join(qd)} | PC: {', '.join(pd_)}")
    print("\nHeterogeneous targets (within-group SD of a component > %.2f):" % SD_FLAG)
    any_h = False
    for t in tg:
        if len(members[t]) < 3: continue
        sd = Q[members[t]].std(0)
        hh = [f"K{c+1} sd={sd[c]:.2f}" for c in range(K) if sd[c] > SD_FLAG]
        if hh: any_h = True; print(f"  {t:<24} " + ", ".join(hh))
    if not any_h: print("  none")
    print("\nLikely-redundant pool pairs (Q distance < %.2f): %s" % (REDUNDANT_Q, "; ".join(f"{g}~{h} ({d:.3f})" for g, h, d in red) or "none"))
    print("\nOutputs in", a.outdir)


if __name__ == "__main__":
    main()
