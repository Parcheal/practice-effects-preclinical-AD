# Early practice effects and clinical progression in preclinical Alzheimer's disease: analysis code

R code for the three-cohort analyses (A4, PREVENT-AD, ADNI) of early practice effects (PE) and clinical
progression in cognitively unimpaired older adults.

**No data are included in this repository.** It contains no participant-level data, identifiers, or
derived result files. All cohort data are governed by data-use agreements and must be obtained directly from
each cohort (see [Data access](#data-access)).

## What the code does

For each cohort the scripts:

1. estimate individual PE as the within-person slope of cognition across repeated early assessments
   (`lme4::lmer`, fixed + random intercept and slope; empirical Bayes slope), residualized on baseline
   cognition and standardized within the PE-estimable sample. Scores used: A4 Cogstate Composite (C3) at
   Visits 1, 3, 9, 15; PREVENT-AD RBANS total scale index at baseline, month 12, month 24; ADNI PHC_MEM at
   baseline, month 6, month 12 (month 24 if month 12 is missing);
2. fit landmark Cox models for post-landmark clinical progression (A4: confirmed CDR-GS 0 to 0.5 or higher
   after Week 48; PREVENT-AD: incident MCI after month 24; ADNI: DXSUM CN to MCI/dementia after the terminal
   PE visit), first with demographics, APOE e4 and baseline cognition, then adding amyloid/tau biomarkers;
3. run the sensitivity, validation (bootstrap and cross-validated delta C-index), proportional-hazards,
   interaction and subgroup analyses reported in the Supplement;
4. across cohorts, run the random-effects meta-analysis, competing-risk checks, decision-curve analysis and
   the A4 two-step plasma p-tau217 + PE enrichment analysis.

Plotting scripts are not included.

## Run order

`R/00_run_all.R` prints this order. By default it is a dry run; set `PE_RUN_ANALYSES=true` to source the
scripts in sequence. Each script can also be run on its own with `Rscript R/<cohort>/<script>.R`.

| Step | Script | Cohort |
|---|---|---|
| 1 | `R/a4/01_run_corrected_analysis.R` | A4 PE estimation, PACC LMM, CDR landmark Cox, bootstrap/CV |
| 2 | `R/a4/02_generate_corrected_tables_and_sensitivity.R` | A4 tables and sensitivity analyses |
| 3 | `R/a4/03_run_c2_pe_covariate_interactions.R` | A4 PE x covariate interactions |
| 4 | `R/preventad/01_run_prevent_ad_replication.R` | PREVENT-AD PE estimation and MCI landmark Cox |
| 5 | `R/preventad/02_run_prevent_ad_supplementary.R` | PREVENT-AD bootstrap/CV supplement |
| 6 | `R/adni/01_run_adni_replication_v3.R` | ADNI PE estimation and landmark Cox |
| 7 | `R/adni/02_run_adni_replication_v3_supplement.R` | ADNI supplement (primary cohort table, PH checks, MCI endpoint) |
| 8 | `R/adni/03_run_adni_replication_v3_phase2.R` | ADNI biomarker timing, IPW, A/T strata, interactions |
| 9 | `R/a4/05_run_incremental_full.R` | A4 incremental-value, enrichment and extended analyses |
| 10 | `R/common/06_run_upgrade_D1_D7.R` | Cross-cohort meta-analysis, competing risk, DCA, two-step enrichment |

`R/verify_key_values.R` writes the headline values expected from a rerun to `verify/expected_key_values.csv`.

### Where outputs go

Scripts write reruns under `verify/` (created on first run, git-ignored). Steps 2, 3, 9 and 10 read upstream
results from `outputs/locked/`, so copy the step 1 and steps 4 to 8 results from `verify/` into
`outputs/locked/` before running them:

| Read by downstream scripts | Written by |
|---|---|
| `outputs/locked/a4/analysis_baseline_w1_15_strict_preW48_lmer.csv`, `learningcurve_PE_lmer_w1_15_strict_preW48.csv`, `surv_primary_complete_w1_15_strict_preW48_lmer.rds`, `surv_group_complete_w1_15_strict_preW48_lmer.rds`, `survival_bootstrap_delta_c_w1_15_strict_preW48_lmer.csv`, `survival_cv_delta_c_w1_15_strict_preW48_lmer.csv` | step 1 (`verify/a4/locked_outputs_corrected/`) |
| `outputs/locked/preventad/survival_dataset_primary.csv` | step 4 (`verify/preventad/`) |
| `outputs/locked/adni/v3/exposure_subject_level.csv` | step 6 (`verify/adni/v3/`) |
| `outputs/locked/adni/supplement/Table1_primary_cohort_dataset.csv` | step 7 (`verify/adni/supplement/`) |
| `outputs/locked/adni/v3/adni_dxsum_primary_452.csv`, `outputs/locked/adni/mediation/analytic_dx23_main_biomarker.csv`, `outputs/locked/adni/mediation/mediation_results.csv` | not produced by scripts in this repository (earlier ADNI primary-replication and mediation scripts); step 9 reads `adni_dxsum_primary_452.csv` unconditionally and stops if it is missing, and runs the mediation block only when the two mediation files are present |

These intermediate files are participant-level. Keep them local; `.gitignore` blocks `outputs/`, `verify/`
and all tabular/R data formats.

## Display-item map

`MANIFEST.csv` lists, for each table and figure of the analysis package, the producing script and the
intermediate files it was built from. Paths under `outputs/` refer to local, participant-derived files that
are not part of this repository. Display labels follow the analysis package (eTable/eFigure = Supplementary
Table/Figure) and may not match the numbering of the published Supplementary Information.

## Data access

| Cohort | Access |
|---|---|
| A4 Study | Alzheimer's Clinical Trials Consortium platform, https://www.a4studydata.org/ |
| ADNI | ADNI database via LONI IDA, https://adni.loni.usc.edu/ |
| PREVENT-AD | StoP-AD Centre Registered PREVENT-AD resource, https://registeredpreventad.loris.ca/ |

## Expected data layout

Set `PE_DATA_ROOT` to your data folder (default: `<repo>/data`). File names must match the cohort exports
used for the paper (ADNI LONI exports dated 10Apr2026; PREVENT-AD registered release 2024-08-01).

```
<PE_DATA_ROOT>/
  a4/
    Derived Data/SUBJINFO.csv
    Derived Data/PACC.csv
    Derived Data/SV.csv
    Derived Data/DS.csv
    External Data/imaging_SUVR_amyloid.csv
    External Data/biomarker_pTau217.csv
    External Data/cogstate_battery.csv
    Raw Data/cdr.csv
  preventad/
    Clinical_data/RBANS_Registered_PREVENTAD_2024-08-01.csv
    Clinical_data/Demographics_Registered_PREVENTAD_2024-08-01.csv
    Clinical_data/Genetics_Registered_PREVENTAD.csv
    Clinical_data/PET_NAV_SUVR_ref-cerebellumCortex_Registered_PREVENTAD_2024-08-01.csv
    Clinical_data/CSF_proteins_Registered_PREVENTAD.csv
    Clinical_data/EL_CDR_MoCA_Registered_PREVENTAD.csv
    Clinical_data/MCI_Registered_PREVENTAD_2024-08-01.csv
  adni/
    DXSUM_10Apr2026.csv
    CDR_10Apr2026.csv
    PTDEMOG_10Apr2026.csv
    APOERES_10Apr2026.csv
    ADSP_PHC_COGN_10Apr2026.csv
    UWNPSYCHSUM_10Apr2026.csv
    EMBICqCP_10Apr2026.csv
    ADNI_EMBICDCB_10Apr2026.csv
    UPENN_PLASMA_FUJIREBIO_QUANTERIX_10Apr2026.csv
    ADSP_PHC_PET_Amyloid_Simple_10Apr2026.csv
    UPENNBIOMK_ROCHE_ELECSYS_10Apr2026.csv
    ADSP_PHC_BIOMARKER_10Apr2026.csv
    ADSP_PHC_T1_FS_10Apr2026.csv
    ADSP_PHC_PET_Tau_Detailed_10Apr2026.csv
    ADVERSE_10Apr2026.csv
    RECADV_10Apr2026.csv
    ADNIMERGE2/ADNIMERGE2/data/UCSFFSX7.rda
```

## Software

Primary analyses: R 4.5.2 (Windows) with survival 3.8-3, lme4 1.1-37, lmerTest 3.1-3, Matrix 1.7-4,
dplyr 1.1.4, readr 2.1.5, tidyr 1.3.1, tibble 3.3.0, stringr 1.6.0, broom 1.0.11, broom.mixed 0.2.9.6,
ggplot2 4.0.0, nlme 3.1-168 (`env/sessionInfo.txt`). The ADNI v3 and supplement runs used R 4.3.3 on Ubuntu
(`env/sessionInfo_adni_v3.txt`, `env/sessionInfo_adni_supplement.txt`).

Also used: metafor, cmprsk, survminer, patchwork (see `env/packages.txt`); steps 9 and 10 additionally call
riskRegression, mice, mediation, EValue, lcmm, JM, JMbayes2, scales, ragg and svglite when installed, and
record a module as not estimable if a package is missing.

## Citation

Pei J, Wang Y, Jiang D, Wu D, Guo T, Lin Z, for the Alzheimer's Disease Neuroimaging Initiative.
Early practice effects and clinical progression in preclinical Alzheimer's disease. *npj Dementia* (under review).

Data used in preparation of this work were obtained from the A4 Study, PREVENT-AD and ADNI; please follow each
cohort's acknowledgement and publication policy.

## License

MIT; see `LICENSE`.
