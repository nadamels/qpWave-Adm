# qpWave / qpAdm for modern Iranian sub-populations

```bash
cd /scratch/nme3923/IranGen/modern_anc_unfilt_PLUS54/merged_datasets/ancient_modern_all/qpWave-Adm
bash setup_admixtools.sh          # once, on a LOGIN node (needs internet): installs ADMIXTOOLS 2 into the glimpse2 env
sbatch run_qpwave_qpadm.sh        # the whole analysis; results -> ../qpwave_qpadm_out/
```

Run settings are environment variables (`MAX_SOURCES=4 PTHRESH=0.01 sbatch --export=ALL run_qpwave_qpadm.sh`):
`BED_PREFIX` (default `all_anc_mod.run3`, all SNPs; `all_anc_mod.run3.LDpruned` also works), `MAX_SOURCES` (3), `PTHRESH` (0.05),
`TOP_N` (3), `WEIGHT_TOL` (0), `FORCE_F2` (0), `CONDA_ENV`, `DATA_DIR`, `OUTDIR`.

## What the pipeline does
| step | script | what |
|---|---|---|
| 1 | `scripts/make_pops.py` | joins `fam_meta/sample_meta.tsv` to the fam; writes a copy of the fam with **FID = population** (bed/bim are symlinked, never modified) |
| 2 | `scripts/qpwave_qpadm.R` | `extract_f2` (cached) -> **qpWave** between targets -> **qpWave screen of source sets** -> **qpAdm rotation** -> model selection -> weights ± SE -> leave-one-right-out check |
| 3 | `scripts/plot_qpadm.R` | best-model bars, top-N competing models, qpWave heat map |

**Targets** (68 modern Iranians, from the Ward clustering): `Mod_Zoroastrian` (6), `Mod_Turk` (1; not Turk_Azeri),
`Mod_Jewish_branch` (10 = 7 Jewish + IGP_57, IGP_19, IGP_25), `Mod_Armenian_Assyrian` (5 = Armenian, Assyrian, Assyrian_Armenian),
`Mod_Central_Iranian` (46 = everyone else, incl. Fars_Jewish IGP_16). The 2 `Iranian:LP6005…` samples have no cluster and are ignored.

**Sources / right set** live in `config/pop_pools.tsv` (edit freely):
* `right_base` = always RIGHT, never a source: Mbuti, Papuan, Karitiana, Han, Onge, Chukchi (deep outgroups, o9a-style) +
  **WHG and EEHG** (regional discriminators that help separate CHG-like/Iran_N ancestry from Anatolia_N ancestry).
* `pool` = 13 candidate sources (Iran_N = Ganj_Dareh_N is the only Iranian source; Anatolia_N, Anatolia_ChL_BA, Levant_N, Natufian, Armenia_N,
  Armenia_EBA, Armenia_LBA_IA, Steppe_EMBA, Steppe_MLBA, Turan_N, BMAC, South_Asia_BA). Post-Neolithic Iranian groups (Iran_ChL, Iran_BA_SE, Iran_BIA_NW) are `ignore`d on purpose (they are so close to modern Iranians that they give 1-source models); set their role back to `pool` in the TSV to use them. For every model, pool populations not used on the
  LEFT are moved to the RIGHT ("rotating" qpAdm, Harney et al. 2021), so each model has its own right set.

**Source selection** = (i) qpWave on each candidate source set (right = outgroups + remaining pool): the rank k-2 hypothesis must be
*rejected*, otherwise the sources are not distinguishable and the set is dropped; (ii) qpAdm on the survivors; (iii) a model is
*feasible* if p >= `PTHRESH` and all weights are in [0,1]; the best model is the one with the **fewest sources, then highest p**.
If nothing fits, the best-p model is reported and flagged `NO_FIT_...`.

## Outputs (`qpwave_qpadm_out/`)
`qpwave_targets_pairwise.tsv|png` · `qpwave_targets_joint_rank.tsv` · `qpwave_source_screen.tsv` · `qpadm_all_models.tsv` ·
`qpadm_selected_models.tsv` (top-N, weights ± SE) · `qpadm_best_loo_right.tsv` · `qpadm_source_frequency_in_feasible_models.tsv` ·
`qpadm_feasible_models_by_size.tsv` · `qpadm_best_models.pdf|png` · `qpadm_top_models.pdf|png`.

## Caveats (read before interpreting)
* **Ancient genomes here are imputed diploid calls.** `extract_f2` is run with `adjust_pseudohaploid=FALSE`. Imputation pulls genotypes toward the
  reference panel and can shrink differences between ancient groups; treat results as exploratory and confirm key models on non-imputed data if possible.
* `Mod_Turk` is a single genome; its weights are valid for that individual only. `Armenia_N` (2), `Natufian` (2), `South_Asia_BA` (2), `WHG` (2) are tiny.
* Pooled ancient groups are heterogeneous (`Armenia_LBA_IA`, `Anatolia_ChL_BA`); change the pooling in `pop_pools.tsv`.
* ~700 models per target are tested at p >= 0.05 without multiple-testing correction; look at the top-N alternatives, the
  leave-one-right-out table and source frequencies, not only the single best model. Many p-values near 0.05-0.1 with wide SEs means weak support.
* No CHG, Ust_Ishim, Kostenki14 or MA1 genomes are in this dataset.
* `tests/simulate_plink.py` makes a fake diploid dataset with the same labels; it only checks that the plumbing runs.

References: Harney et al. 2021 *Genetics* (qpAdm user guide); Maier et al. 2023 *eLife* (ADMIXTOOLS 2); admixr rotation vignette;
Agranat-Tamir et al. 2020 *Cell* supplement (o9a right set).
