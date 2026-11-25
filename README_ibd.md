# IBD-based clustering and ancestry modelling (IBDseq → Leiden clusters → painting profiles → NNLS)

Continues from the files you already made with IBDseq + the paper filters
(`*.paper_filtered.final.with_cM.tsv.gz`: LOD ≥ 3, ≥ 2 cM, excess-long-IBD regions removed). It does **not** redo IBDseq or the
GenomicRanges excess-IBD filter.

```bash
cd /scratch/nme3923/IranGen/modern_anc_unfilt_PLUS54/merged_datasets/ancient_modern_all/qpWave-Adm
git pull
bash setup_ibd_packages.sh                        # once, LOGIN node: data.table, igraph, limSolve, leidenAlg, ggplot2
IBD_DIR=/path/to/dir/with/the/ibd/files STEPS="samples pairs cluster" sbatch --export=ALL run_ibd_ancestry.sh
#   -> look at ../ibd_ancestry_out/ibd_cluster_summary.tsv and ibd_cluster_hierarchy_nodes.tsv, then fill in
#      config/ibd_source_pool.tsv and/or config/ibd_models.tsv with the cluster names, and:
STEPS="paint nnls" sbatch --export=ALL run_ibd_ancestry.sh
```
`IBD_DIR` (or `IBD_FILE`) must point at the `.with_cM.tsv.gz`; `DATA_DIR`/`OUTDIR` default as in `run_qpwave_qpadm.sh`.
Check the first lines of `ibd_pair_totals.log`: it prints the header it found and the columns it chose. If the header is not
recognised set `COL_ID1 COL_ID2 COL_CHR COL_LEN_CM` (name or 1-based index; or `COL_START_CM`+`COL_END_CM`). If the IDs in the IBD file differ
from `fam_meta/sample_meta.tsv` (e.g. a haplotype suffix) use `ID_SUB='_[12]$'`; the script stops if > 50 % of segments have unknown IDs.

| step | script | what it does (methods paragraph) |
|---|---|---|
| samples | `scripts/ibd_make_samples.py` | ancient/modern flag (`Undergone Imputation`), date midpoint (BP), lat/lon, modern target groups (`Mod_*`, same rules as `make_pops.py`) and modern donor groups (= `Population`) |
| pairs | `scripts/ibd_01_pair_totals.R` | streams the 3.8 GB table through `zcat` in 4 M-line chunks → total cM per (pair, chromosome). Also within-individual segments (`roh_by_individual.tsv`, only ≥ 2 cM here) |
| cluster | `scripts/ibd_02_cluster.R` | igraph network of the **ancient** individuals, weight = cM shared / genome cM; iterative Leiden (r = 0.5, +0.5 up to 3; split into communities and recurse, each sub-network restarts at 0.5); then top-down collapse of the hierarchy by date span (`MAX_DATE_SPAN` = 1000 yr) and geographic span (`MAX_KM` = 500 km), with manual overrides in `config/ibd_cluster_overrides.tsv` |
| paint | `scripts/ibd_03_painting.R` | IBD of every individual with every donor group (ancient clusters + modern groups), per chromosome, leave-one-out (a member of a group sees n−1 donors in its own group) |
| nnls | `scripts/ibd_04_nnls.R` | target = Σ w·source, w ≥ 0, Σw = 1, `limSolve::lsei`; SE = weighted block jackknife leaving out each chromosome (block weight = chromosome cM, ADMIXTOOLS estimator) |

Outputs in `ibd_ancestry_out/`: `ibd_clusters.tsv` (per individual), `ibd_cluster_summary.tsv`, `ibd_cluster_hierarchy_nodes.tsv`,
`ibd_groups.tsv`, `group_mean_ibd_cM.tsv`, `painting_profiles_by_group.tsv`, `nnls_{sum,none}/ibd_nnls_models.tsv`
(`weights`, `se`, `rss`, `min_weight`, `worst_residual_donors`), `…_weights_long.tsv`, `…_best_feasible_by_size.tsv`, `ibd_nnls_models.png|pdf`.

## Modelling choices you should know about
* **Models.** Explicit (`config/ibd_models.tsv`) and/or enumerated (`config/ibd_source_pool.tsv`: all 1..`MAX_SOURCES` combinations for every `Mod_*` target,
  ranked by RSS; "feasible" = every weight ≥ `MIN_WEIGHT`). NNLS has no p-value: more sources always lower the RSS, so compare models of the same size and look at SEs.
* **Donors of a model** = all donor groups (≥ `DONOR_MIN_N` members) except the model's sources and the target's own group (`EXCLUDE_SOURCES=1`). Needs enough *other* groups that
  are differentially related to the sources; `EXCLUDE_SOURCES=0` keeps the sources as donors (more information, but copying from the source itself then dominates).
* **Profile of an individual** = mean cM shared per donor individual (`SIZE_NORM=1`, so large groups are not favoured), then a group profile = mean over its members.
* **Normalisation (`PROFILE_NORM`).** The methods say profiles are "normalised to [0,1]". `sum` (each individual's profile sums to 1) does that, but it is not linear in
  ancestry: a mixture's weight is pulled towards the source that has more total IBD with the donors (small / inbred / high-coverage sources). `none` keeps raw mean-IBD vectors
  (linear, unbiased, but sensitive to differences in overall IBD detectability between individuals, e.g. coverage). The driver runs **both** (`nnls_sum/`, `nnls_none/`);
  if they disagree, the difference is that bias, not noise. On the simulated test (`bash tests/run_ibd_test.sh`; truth 70/30 and 50/50) `none` gave 0.68/0.32 and 0.54/0.46,
  `sum` 0.50/0.50 and 0.64/0.36.
* **Spatiotemporal collapse** is rule-based, not manual as in the paper; it is a starting point. Read `ibd_cluster_hierarchy_nodes.tsv` (n, populations, date/distance spans per node) and
  override nodes in `config/ibd_cluster_overrides.tsv`. Cluster names are `<top Population>_c<node path>`. Unclustered = ancient individuals with no IBD with any other ancient individual.
* **Not reproduced exactly:** the clustering in the methods re-runs IBDseq on ancient samples only with INFO > 0.5 and ≥ 1 cM. Your file is the ≥ 2 cM, LOD ≥ 3 set of ancient + modern Iran,
  so the network is built from that (set `MIN_CM` only to be *stricter*). If clusters look too coarse for low-coverage individuals, re-run IBDseq as described and point `IBD_FILE` at it
  (the modern individuals are ignored by the clustering anyway). ROH ≥ 1 cM likewise needs that lower-cutoff file.
* `leidenAlg` is used when installed; otherwise `igraph::cluster_leiden` (modularity, same resolution parameter) is used and the log says so. `limSolve` missing → `quadprog` fallback (testing only).

## Settings (environment variables, defaults)
`STEPS`, `IBD_DIR/IBD_FILE`, `MIN_CM=0`, `ID_SUB`, `R_START=0.5 R_STEP=0.5 R_MAX=3`, `MIN_COMMUNITY=1`, `MAX_DATE_SPAN=1000`, `MAX_KM=500`, `OVERRIDES`, `SEED=1`, `DONOR_MIN_N=2`,
`MAX_SOURCES=3`, `MIN_WEIGHT=0.02`, `TOP_N=5`, `SIZE_NORM=1`, `EXCLUDE_SOURCES=1`, `METHOD=lsei|nnls`, `PROFILE_NORMS="sum none"`, `ENUMERATE`, `TARGETS`, `JACKKNIFE=1`.
Run settings and git commit are appended to `run_info.txt`.

## Caveats
* The test data are simulated; the real file's header/IDs were not available while writing this, hence the auto-detection and the checks above. Read `ibd_pair_totals.log` first.
* Memory: the aggregated pair × chromosome table is small compared with the input, but `painting.rds` holds individuals × groups × 22 numbers; 64 GB is ample.
* Same caveat as the qpAdm pipeline: ancient genomes are imputed, targets are defined from the same data, and single-genome groups (`Mod_Turk`) have no real SE.
