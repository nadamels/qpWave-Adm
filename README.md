# qpWave / qpAdm for modern Iranian sub-populations

> IBD-based clustering + painting/NNLS ancestry modelling (IBDseq output → Leiden clusters → `limSolve`) lives in [`README_ibd.md`](README_ibd.md) (`run_ibd_ancestry.sh`).

```bash
cd /scratch/nme3923/IranGen/modern_anc_unfilt_PLUS54/merged_datasets/ancient_modern_all/qpWave-Adm
git pull
bash setup_admixtools.sh          # once, on a LOGIN node (needs internet): installs ADMIXTOOLS 2 into the glimpse2 env
sbatch run_qpwave_qpadm.sh        # the whole analysis; results -> ../qpwave_qpadm_out/
```

Run settings are environment variables (`MAX_SOURCES=4 PTHRESH=0.01 sbatch --export=ALL run_qpwave_qpadm.sh`):
`BED_PREFIX` (default `all_anc_mod.run3`, all SNPs; `all_anc_mod.run3.LDpruned` also works), `PTHRESH` (0.05),
`MIN_SOURCES` (2), `MAX_SOURCES` (3), `MIN_WEIGHT` (0.05), `TOP_N` (3), `WEIGHT_TOL` (0), `FORCE_F2` (0),
`RUN_EXPLORE` (1), `CONDA_ENV`, `DATA_DIR`, `OUTDIR`.

## Pipeline
| step | script | what |
|---|---|---|
| 0 | `scripts/explore_admixture_pca.py` + `plot_admixture_pca.R` | **hypothesis generation** from unsupervised ADMIXTURE (K=8) + PCA: component labels, what each target carries, is every component covered by the candidate pool, does the right set span the components. Never feeds qpAdm. |
| 1 | `scripts/make_pops.py` | joins `fam_meta/sample_meta.tsv` to the fam; writes a copy of the fam with **FID = population** (bed/bim symlinked, never modified) |
| 2 | `scripts/qpwave_qpadm.R` | `extract_f2` (cached) -> **qpWave** between targets -> **qpWave screen of source sets** -> **qpAdm rotation** -> selection -> weights ± SE -> leave-one-right-out check |
| 3 | `scripts/plot_qpadm.R` | best-model bars, top-N competing models, qpWave heat map |

**Targets** (68 modern Iranians, from the Ward clustering): `Mod_Zoroastrian` (6), `Mod_Turk` (1; not Turk_Azeri),
`Mod_Jewish_branch` (10 = 7 Jewish + IGP_57, IGP_19, IGP_25), `Mod_Armenian_Assyrian` (5 = Armenian, Assyrian, Assyrian_Armenian),
`Mod_Central_Iranian` (46 = everyone else, incl. Fars_Jewish IGP_16). The 2 `Iranian:LP6005…` samples have no cluster and are ignored.

**Sources / right set** live in `config/pop_pools.tsv` (edit freely):
* `right_base` = always RIGHT, never a source: Mbuti, Papuan, Karitiana, Onge, Chukchi (deep outgroups, o9a-style) +
  **WHG and EEHG** (regional discriminators that help separate CHG-like/Iran_N ancestry from Anatolia_N ancestry).
* `pool` = 14 rotating candidates: Iran_N (= Ganj_Dareh_N, the only Iranian source), Anatolia_N, Anatolia_ChL_BA, Levant_N, Natufian,
  Armenia_N, Armenia_EBA, Armenia_LBA_IA, Steppe_EMBA, Steppe_MLBA, Turan_N, BMAC, South_Asia_BA, **East_Asia_Han**.
  For every model, pool populations not used on the LEFT are moved to the RIGHT ("rotating" qpAdm, Harney et al. 2021).
* `ignore` = Iran_ChL, Iran_BA_SE, Iran_BIA_NW (post-Neolithic Iranians): removed completely, on purpose, because they are so close to modern
  Iranians that they give 1-source fits. Set a group's role back to `pool` to use it.
* East_Asia_Han is a candidate (not a fixed outgroup) because step 0 showed the Turk sample carries ~13% of the East Asian component (K6)
  and nothing else in the pool carries it; a fixed right population that is the source would absorb that ancestry.

**Source selection** = (i) qpWave on every candidate source set: the rank k-2 hypothesis must be *rejected*, otherwise the sources are not
distinguishable and the set is dropped; (ii) qpAdm on the survivors; (iii) a model is *feasible* if p >= `PTHRESH` **and every weight is in
[`MIN_WEIGHT`, 1]**; (iv) only models with >= `MIN_SOURCES` sources can be selected; the best is the one with the **fewest sources, then highest p**.
`MIN_WEIGHT` stops a near-clade fit from being "improved" into a 2-way model by a negligible second source.
1-source fits are not selected but are listed in `qpadm_one_source_fits.tsv` (they mean "indistinguishable from X", which is itself a result).
If nothing fits, the best-p model is reported and flagged `NO_FIT_...`.

## Using ADMIXTURE / PCA together with qpWave / qpAdm: the rules
ADMIXTURE and PCA are **discovery** tools; qpWave/qpAdm are **hypothesis tests**. To keep the tests valid:
1. **Hypotheses only.** ADMIXTURE/PCA are used to propose candidate sources and to sanity-check the design. Q-fractions or PC positions are never
   inputs to qpAdm and never reported as ancestry proportions; only qpAdm weights ± SE and p-values are results.
2. **Components are not populations.** A component is a drifted cluster; K is arbitrary; ADMIXTURE cannot tell admixture from shared drift, and
   similar Q profiles do not prove a shared source (at K=8, Natufian and Levant_N are 0.077 apart in Q space; the qpWave screen decides, not Q).
3. **Freeze the design before looking at qpAdm output.** The pool and right set are in git (`config/pop_pools.tsv`); every run saves the
   exact file used plus the git commit (`pop_pools.used.tsv`, `run_info.txt`). Changing the pool after seeing results is a new analysis: report
   how many models were tried (`qpadm_all_models.tsv`) and why the pool changed. Do not tune the right set until a model fits.
4. **The right set must be differentially related to the sources.** Step 0 lists components with no carrier in the fixed right set
   (K5, the Iran_N/Turan component, is represented only through rotation, so Iran_N, Turan_N and BMAC must not all be on the left).
5. **Targets defined from the same data.** The target groups come from clustering the same samples; qpAdm p-values do not account for
   that, so treat groupings as descriptive, and note `Mod_Turk` is a single genome.
6. **Multiple testing.** Hundreds of models per target are tested at p >= 0.05 without correction; read the top-N alternatives,
   the leave-one-right-out table and source frequencies, not just the best model.

## Outputs (`qpwave_qpadm_out/`)
`explore/` (ADMIXTURE/PCA: `admixture_group_profiles.png`, `pca_targets_vs_candidates.png`, `admixture_component_coverage.tsv`,
`admixture_pca_candidate_screen.tsv`, `admixture_pool_redundancy.tsv`, …) · `qpwave_targets_pairwise.tsv|png` · `qpwave_targets_joint_rank.tsv` ·
`qpwave_source_screen.tsv` · `qpadm_all_models.tsv` · `qpadm_selected_models.tsv` (top-N, weights ± SE) · `qpadm_one_source_fits.tsv` ·
`qpadm_best_loo_right.tsv` · `qpadm_source_frequency_in_feasible_models.tsv` · `qpadm_feasible_models_by_size.tsv` ·
`qpadm_best_models.pdf|png` · `qpadm_top_models.pdf|png` · `pop_pools.used.tsv` · `run_info.txt`.

## Caveats (read before interpreting)
* **Ancient genomes here are imputed diploid calls.** `extract_f2` is run with `adjust_pseudohaploid=FALSE`. Imputation pulls genotypes toward the
  reference panel and can shrink differences between ancient groups; treat results as exploratory and confirm key models on non-imputed data if possible.
* `Mod_Turk` is one genome; `Armenia_N` (2), `Natufian` (2), `South_Asia_BA` (2), `WHG` (2) are tiny.
* Even without the post-Neolithic Iranian groups, `Armenia_EBA` / `Armenia_LBA_IA` / `Armenia_N` are very close to several targets in ADMIXTURE
  space (Q-distance 0.06-0.10), which is why `MIN_SOURCES=2` and `MIN_WEIGHT` exist.
* Pooled ancient groups are heterogeneous (`Armenia_LBA_IA`, `Anatolia_ChL_BA`); change the pooling in `pop_pools.tsv`.
* No CHG, Ust_Ishim, Kostenki14 or MA1 genomes are in this dataset. Han is a southern East Asian proxy; northern alternatives (Buryat, Mongolian,
  Yakut) are listed in the config.
* `tests/simulate_plink.py` makes a fake diploid dataset with the same labels; it only checks that the plumbing runs.

References: Harney et al. 2021 *Genetics* (qpAdm user guide); Maier et al. 2023 *eLife* (ADMIXTOOLS 2); admixr rotation vignette;
Agranat-Tamir et al. 2020 *Cell* supplement (o9a right set).
