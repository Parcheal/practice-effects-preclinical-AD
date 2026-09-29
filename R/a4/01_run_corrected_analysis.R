

# PE_analysis_code_package path adapter.
script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg)) script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
pkg_root <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(pkg_root, "MANIFEST.csv"))) {
  cur <- pkg_root
  for (i in seq_len(8)) {
    if (file.exists(file.path(cur, "MANIFEST.csv"))) { pkg_root <- cur; break }
    cur <- dirname(cur)
  }
}
source(file.path(pkg_root, "R", "config.R"))
setwd(pkg_root)


# ============================================================================
# A4 Corrected Analysis — fixes time-axis bugs in CDR Cox and PACC LMM,
# and adopts the strict pre-W48 PE window (V1, V3, V9, V15).
#
# Bugs fixed (relative to practice_effects_pe_study/run_learningcurve_from_raw.R):
#
#   B1. CDR Cox (lines 410-437 in original)
#       Original:  filter VISCODE_num >= 51, time = (first_conversion - 51) / 52
#       Problem:   cdr.csv VISCODE is SoA visit number (V18=W48, V48=W168, V54=W192).
#                  VISCODE 51 does not exist; the filter excluded V18, V36, V48 entirely;
#                  treating VISCODE differences as week differences compressed the
#                  real ~3 yr median event time into a fake 0.12 yr.
#       Fix:       Use CDADTC_DAYS_T0 directly. Anchor at each subject's V18 day
#                  (W48 protocol visit). Keep all CDR records with CDADTC_DAYS_T0
#                  greater than or equal to 0. Compute time as real days / 365.25.
#
#   B2. PACC LMM time variable (lines 121-122, 354, 367, 382, 395 in original)
#       Original:  Week = as.numeric(VISCODE), filter Week > 48,
#                  time_years = (Week - 48) / 52
#       Problem:   PACC.csv VISCODE is the SoA visit number, not the calendar week.
#                  V18 = Week 48 (verified against SV.csv visit labels).
#                  filter `Week > 48` (with Week = VISCODE) keeps V54+ which is
#                  Week 192+, not Week 48+.
#       Fix:       Compute true Week = (VISCODE - 6) * 4 (verified by SV.csv).
#                  Filter true_Week > 48, time_years = (true_Week - 48) / 52.
#
#   D1. PE window (decision)
#       Original:  w1_21 = c(1, 3, 9, 15, 21) — visits span Week -14 to Week +60,
#                  which leaks into the post-W48 outcome window via VISCODE 21.
#       Fix:       Adopt strict pre-W48 window w1_15 = c(1, 3, 9, 15) — visits
#                  span Week -14 to Week +36, no overlap with outcome window.
#
# What is NOT changed:
#   - PE estimation method itself (LMM with TDate_DAYS_T0) — it was already correct
#     because cogstate_battery.csv has TDate_DAYS_T0 and the existing code used it.
#   - Cohort eligibility logic.
#   - Covariate set, residualization, z-score grouping thresholds (>0.2, <-0.2).
#   - PACC LMM formula. Cox PH formula.
#
# Outputs land in: source_data/locked_outputs_corrected/
# (the original locked_outputs/ directory is never overwritten.)
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(survival)
  library(broom)
  library(broom.mixed)
})

# ----------------------------------------------------------------------------
# Resolve paths. This script is invoked from the PE_JPAD_clean_package root,
# but reads raw A4 data from the parent Clinical/ directory.
# ----------------------------------------------------------------------------
clinical_dir <- file.path(data_root, "a4")  # packaged A4 input root
pkg_dir      <- pkg_root

raw_path <- function(rel) file.path(clinical_dir, rel)

corrected_dir <- file.path(pkg_dir, "verify", "a4", "locked_outputs_corrected")
if (!dir.exists(corrected_dir)) dir.create(corrected_dir, recursive = TRUE)

cat("============================================\n")
cat("A4 CORRECTED ANALYSIS\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Reads from :", clinical_dir, "\n")
cat("Writes to  :", corrected_dir, "\n")
cat("============================================\n\n")

# ----------------------------------------------------------------------------
# 1) Load A4 data files (identical to original analysis)
# ----------------------------------------------------------------------------
cfg_amyloid_centiloid_cutpoint <- 20

subjinfo <- read_csv(raw_path("Derived Data/SUBJINFO.csv"),
                    show_col_types = FALSE, progress = FALSE) %>%
  filter(SUBSTUDY == "A4") %>%
  filter(MITTFL == 1) %>%
  filter(!is.na(AMYLCENT) & AMYLCENT >= cfg_amyloid_centiloid_cutpoint) %>%
  mutate(
    APOE4 = if_else(APOEGN %in% c("E4/E4", "E3/E4", "E2/E4"), 1L, 0L),
    Sex   = case_when(SEX == 1 ~ "Male", SEX == 2 ~ "Female", TRUE ~ NA_character_)
  ) %>%
  transmute(BID, SUBSTUDY, MITTFL, SUBJITTTR,
            Treatment = TX, Age = AGEYR, Education = EDCCNTU, Sex,
            APOE4, SUVRCER_subjinfo = SUVRCER, AMYLCENT,
            MMSE_V1 = MMSETSV1)

cat("SUBJINFO (A4 + MITTFL + amyloid+):", nrow(subjinfo), "\n")

amyloid <- read_csv(raw_path("External Data/imaging_SUVR_amyloid.csv"),
                    show_col_types = FALSE, progress = FALSE) %>%
  mutate(VISCODE_num = suppressWarnings(as.numeric(as.character(VISCODE)))) %>%
  filter(VISCODE_num == 2, brain_region == "Composite_Summary") %>%
  group_by(BID) %>% slice(1) %>% ungroup() %>%
  transmute(BID, SUVR = suvr_cer)

ptau <- read_csv(raw_path("External Data/biomarker_pTau217.csv"),
                 show_col_types = FALSE, progress = FALSE) %>%
  mutate(VISCODE_num = suppressWarnings(as.numeric(as.character(VISCODE)))) %>%
  filter(VISCODE_num == 6, TEST == "Phosphorylated-tau-217") %>%
  mutate(pTau217 = suppressWarnings(as.numeric(ORRES))) %>%
  filter(!is.na(pTau217)) %>%
  group_by(BID) %>% slice(1) %>% ungroup() %>%
  transmute(BID, pTau217)

pacc_raw <- read_csv(raw_path("Derived Data/PACC.csv"),
                     show_col_types = FALSE, progress = FALSE)

pacc_baseline <- pacc_raw %>%
  mutate(VISCODE_num = suppressWarnings(as.numeric(as.character(VISCODE)))) %>%
  filter(VISCODE_num == 1) %>%
  transmute(BID, baseline_pacc = PACC)

# CORRECTED: build per-subject (BID, VISCODE) -> day_T0 lookup from SV.csv,
# then join to PACC. This is more robust than the formula (VISCODE-6)*4 because:
#   - It uses each subject's actual visit date (not the median).
#   - Sentinel/special VISCODE values (701, 997, 999) are not present in SV.csv,
#     so they fall out of the join automatically.
#   - The mapping is self-documenting (SV.csv labels e.g. V18 = "wk48 Clinic").
sv <- read_csv(raw_path("Derived Data/SV.csv"),
               show_col_types = FALSE, progress = FALSE) %>%
  mutate(VISCODE_num = suppressWarnings(as.numeric(as.character(VISITCD)))) %>%
  filter(!is.na(VISCODE_num) & !is.na(SVSTDTC_DAYS_T0)) %>%
  group_by(BID, VISCODE_num) %>%
  summarise(day_T0 = as.numeric(min(SVSTDTC_DAYS_T0)), .groups = "drop")

cat("SV.csv per-subject visit-date lookup rows:", nrow(sv), "\n")

# Anchor at protocol-defined Week 48 = 336 days post-T0 (fixed, same for all)
W48_DAYS <- 336

pacc_long <- pacc_raw %>%
  mutate(VISCODE_num = suppressWarnings(as.numeric(as.character(VISCODE)))) %>%
  filter(!is.na(VISCODE_num) & !is.na(PACC)) %>%
  inner_join(sv, by = c("BID", "VISCODE_num")) %>%
  filter(day_T0 > W48_DAYS) %>%   # post-Week 48
  mutate(time_years = (day_T0 - W48_DAYS) / 365.25,
         Week = day_T0 / 7) %>%
  transmute(BID, VISCODE = VISCODE_num, day_T0, Week, time_years, PACC)

cat("PACC long-format (post-W48, SV-joined):", nrow(pacc_long),
    "obs across", n_distinct(pacc_long$BID), "subjects\n")
cat("PACC time_years (real days from W48) range:",
    round(range(pacc_long$time_years), 3), "years\n")

cog <- read_csv(raw_path("External Data/cogstate_battery.csv"),
                show_col_types = FALSE, progress = FALSE)
stopifnot("TDate_DAYS_T0" %in% names(cog))

cog_c3 <- cog %>%
  mutate(VISCODE_num = suppressWarnings(as.numeric(as.character(VISCODE)))) %>%
  filter(!is.na(VISCODE_num)) %>%
  filter(!is.na(C3Comp))

# Pull baseline (V1) and W48 (V51) Cogstate for cohort eligibility - same as original
cog_bl_w48 <- cog_c3 %>%
  filter(VISCODE_num %in% c(1, 51)) %>%
  group_by(BID, VISCODE_num) %>%
  summarise(C3Comp = mean(C3Comp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = VISCODE_num, values_from = C3Comp, names_prefix = "VIS") %>%
  rename(baseline_C3Comp = VIS1, week48_C3Comp = VIS51)

cohort <- subjinfo %>%
  inner_join(cog_bl_w48, by = "BID") %>%
  inner_join(pacc_baseline, by = "BID") %>%
  semi_join(pacc_long, by = "BID") %>%
  left_join(amyloid, by = "BID") %>%
  left_join(ptau, by = "BID") %>%
  mutate(
    Age_z         = as.numeric(scale(Age)),
    Education_z   = as.numeric(scale(Education)),
    SUVR_z        = as.numeric(scale(SUVR)),
    pTau217_z     = as.numeric(scale(pTau217)),
    baseline_C3_z = as.numeric(scale(baseline_C3Comp)),
    Sex_F         = factor(Sex, levels = c("Male", "Female")),
    APOE4_carrier = factor(APOE4, levels = c(0, 1), labels = c("Non-carrier", "Carrier"))
  )

cat("Cohort N:", nrow(cohort), "\n")
write_csv(cohort, file.path(corrected_dir, "cohort_baseline.csv"))

# ----------------------------------------------------------------------------
# 2) PE estimation — strict pre-W48 window (V1, V3, V9, V15) using real days
# ----------------------------------------------------------------------------
EARLY_VISITS <- c(1, 3, 9, 15)
WINDOW_NAME  <- "w1_15_strict_preW48"

cog_agg <- cog_c3 %>%
  filter(BID %in% cohort$BID) %>%
  filter(VISCODE_num %in% EARLY_VISITS) %>%
  group_by(BID, VISCODE_num) %>%
  summarise(
    C3Comp        = mean(C3Comp, na.rm = TRUE),
    TDate_DAYS_T0 = mean(TDate_DAYS_T0, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  group_by(BID) %>%
  mutate(
    t0_days    = suppressWarnings(min(TDate_DAYS_T0, na.rm = TRUE)),
    time_years = (TDate_DAYS_T0 - t0_days) / 365.25
  ) %>%
  ungroup() %>%
  select(-t0_days)

cat("\nPE window:", WINDOW_NAME, "(visits", paste(EARLY_VISITS, collapse = ","), ")\n")
cat("Time spans (real Week, by visit):\n")
print(cog_agg %>%
        group_by(VISCODE_num) %>%
        summarise(median_week = round(median(TDate_DAYS_T0, na.rm = TRUE) / 7, 1),
                  n = n(), .groups = "drop"))

# LMM PE (the manuscript's primary PE estimator)
lmer_pe_dat <- cog_agg %>%
  group_by(BID) %>%
  filter(n() >= 2) %>%
  ungroup()

m_lc <- lmer(C3Comp ~ time_years + (time_years | BID),
             data = lmer_pe_dat, REML = TRUE)
re <- ranef(m_lc)$BID
fe_slope <- unname(fixef(m_lc)["time_years"])
slope_subject <- fe_slope + re[, "time_years"]

learning_lmer <- tibble(
  BID = rownames(re),
  learning_slope_raw_lmer = as.numeric(slope_subject)
)

base_c3 <- cog_agg %>%
  group_by(BID) %>%
  summarise(
    n_early_visits      = n(),
    baseline_C3_from_cog = {x <- C3Comp[VISCODE_num == 1]; if (length(x) == 0) NA_real_ else x[1]},
    .groups = "drop"
  )

learning_lmer <- base_c3 %>%
  inner_join(learning_lmer, by = "BID") %>%
  filter(!is.na(learning_slope_raw_lmer) & !is.na(baseline_C3_from_cog))

res_mod <- lm(learning_slope_raw_lmer ~ baseline_C3_from_cog, data = learning_lmer)
learning_lmer <- learning_lmer %>%
  mutate(
    learning_slope_resid_lmer        = resid(res_mod),
    learning_slope_resid_z_lmer      = as.numeric(scale(learning_slope_resid_lmer)),
    learning_slope_resid_z_decr_lmer = -learning_slope_resid_z_lmer,
    LC_group_lmer = case_when(
      learning_slope_resid_z_lmer > 0.2 ~ "High_PE",
      learning_slope_resid_z_lmer < -0.2 ~ "Low_PE",
      TRUE ~ "Moderate_PE"
    )
  )

write_csv(learning_lmer,
          file.path(corrected_dir,
                    paste0("learningcurve_PE_lmer_", WINDOW_NAME, ".csv")))

cat("\nPE estimable subjects:", nrow(learning_lmer), "\n")
cat("PE group counts:\n")
print(table(learning_lmer$LC_group_lmer))

dat_lmer <- cohort %>%
  left_join(learning_lmer %>%
              select(BID, n_early_visits, learning_slope_raw_lmer,
                     learning_slope_resid_z_lmer, learning_slope_resid_z_decr_lmer,
                     LC_group_lmer),
            by = "BID") %>%
  mutate(LC_group_lmer = factor(LC_group_lmer,
                                levels = c("High_PE", "Moderate_PE", "Low_PE")))

write_csv(dat_lmer,
          file.path(corrected_dir,
                    paste0("analysis_baseline_", WINDOW_NAME, "_lmer.csv")))

# ----------------------------------------------------------------------------
# 3) PACC LMM — corrected time variable
# ----------------------------------------------------------------------------
pacc_long_lc <- pacc_long %>%
  inner_join(dat_lmer %>% filter(!is.na(LC_group_lmer)), by = "BID")
# time_years already computed correctly upstream from day_T0 / 365.25

cat("\nPACC LMM analytic data:", nrow(pacc_long_lc),
    "obs across", n_distinct(pacc_long_lc$BID), "subjects\n")
cat("PACC time_years range (corrected):",
    round(range(pacc_long_lc$time_years), 3), "years\n")

m_pacc_group <- lmer(
  PACC ~ time_years * LC_group_lmer + Age_z + Sex_F + Education_z +
    APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z +
    (1 + time_years | BID),
  data = pacc_long_lc %>% filter(!is.na(SUVR_z) & !is.na(pTau217_z)),
  REML = FALSE
)

g_tidy <- broom.mixed::tidy(m_pacc_group, conf.int = TRUE) %>%
  filter(str_detect(term, "^time_years:LC_group"))
write_csv(g_tidy, file.path(corrected_dir,
          paste0("PACC_LMM_group_model5_tidy_", WINDOW_NAME, "_lmer.csv")))

pacc_long_cont <- pacc_long %>%
  inner_join(dat_lmer %>% filter(!is.na(learning_slope_resid_z_decr_lmer)), by = "BID")
# time_years already computed correctly upstream from day_T0 / 365.25

m_pacc_cont <- lmer(
  PACC ~ time_years * learning_slope_resid_z_decr_lmer + Age_z + Sex_F + Education_z +
    APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z +
    (1 + time_years | BID),
  data = pacc_long_cont %>% filter(!is.na(SUVR_z) & !is.na(pTau217_z)),
  REML = FALSE
)

c_tidy <- broom.mixed::tidy(m_pacc_cont, conf.int = TRUE) %>%
  filter(term == "time_years:learning_slope_resid_z_decr_lmer")
write_csv(c_tidy, file.path(corrected_dir,
          paste0("PACC_LMM_cont_model5_term_", WINDOW_NAME, "_lmer.csv")))

cat("\nPACC LMM continuous PE term (corrected time scale):\n")
print(c_tidy)

# ----------------------------------------------------------------------------
# 4) CDR Cox — corrected time variable using CDADTC_DAYS_T0
# ----------------------------------------------------------------------------
cdr_raw <- read_csv(raw_path("Raw Data/cdr.csv"),
                    show_col_types = FALSE, progress = FALSE) %>%
  filter(!is.na(CDGLOBAL)) %>%
  mutate(VISCODE_num = suppressWarnings(as.numeric(as.character(VISCODE)))) %>%
  # Keep all post-randomization CDR records; exclude screening (negative days)
  # and sentinel codes (701/702/997/999).
  filter(!is.na(VISCODE_num)) %>%
  filter(!is.na(CDADTC_DAYS_T0) & CDADTC_DAYS_T0 >= 0) %>%
  filter(VISCODE_num <= 200) %>%
  transmute(BID, VISCODE = VISCODE_num, day_T0 = CDADTC_DAYS_T0,
            CDR_GS = CDGLOBAL)

cat("\nCDR records (post-T0, non-sentinel):", nrow(cdr_raw),
    "across", n_distinct(cdr_raw$BID), "subjects\n")

# Per-subject V18 anchor (Week 48); fall back to fixed 336 days when V18 missing.
v18_anchor <- cdr_raw %>%
  filter(VISCODE == 18) %>%
  group_by(BID) %>%
  summarise(anchor_day = min(day_T0), .groups = "drop")

cat("Subjects with V18 (W48) CDR:", nrow(v18_anchor), "\n")

# Find first qualifying confirmed conversion per subject
cdr_conv <- cdr_raw %>%
  arrange(BID, day_T0, VISCODE) %>%
  group_by(BID) %>%
  mutate(
    n_visits             = n(),
    CDR_05               = (CDR_GS >= 0.5),
    confirmed_conversion = (CDR_05 & lead(CDR_05)),
    conversion_day       = ifelse(confirmed_conversion, day_T0, NA_real_)
  ) %>%
  summarise(
    n_visits         = first(n_visits),
    converted        = ifelse(n_visits >= 2, any(confirmed_conversion, na.rm = TRUE), FALSE),
    first_conv_day   = ifelse(any(confirmed_conversion, na.rm = TRUE),
                              min(conversion_day, na.rm = TRUE), NA_real_),
    last_visit_day   = max(day_T0),
    first_conv_visit = ifelse(any(confirmed_conversion, na.rm = TRUE),
                              VISCODE[which(day_T0 == min(conversion_day, na.rm = TRUE))[1]],
                              NA_real_),
    .groups = "drop"
  ) %>%
  left_join(v18_anchor, by = "BID") %>%
  mutate(
    anchor_day            = if_else(is.na(anchor_day), 336, anchor_day),
    time_to_event_days    = if_else(converted & !is.na(first_conv_day),
                                    first_conv_day - anchor_day,
                                    last_visit_day - anchor_day)
  )

cat("\nCDR conversion summary (corrected):\n")
cat("  n subjects:", nrow(cdr_conv), "\n")
cat("  converted (events):", sum(cdr_conv$converted), "\n")
cat("  first_conv_visit distribution (events only):\n")
print(table(cdr_conv$first_conv_visit[cdr_conv$converted]))
cat("  time_to_event_years quantiles (events only):\n")
print(quantile(pmax(cdr_conv$time_to_event_days[cdr_conv$converted] / 365.25, 0.01),
               c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE))

# Identify prevalent W48 cases: CDR >= 0.5 at the V18 (W48) visit itself,
# or a confirmed conversion whose first qualifying CDR visit is at/before the
# subject-specific W48 anchor. The second condition guards against a missing
# V18 CDR record masking a pre-landmark confirmed conversion.
v18_status <- cdr_raw %>%
  filter(VISCODE == 18) %>%
  group_by(BID) %>%
  summarise(v18_cdr05 = any(CDR_GS >= 0.5, na.rm = TRUE), .groups = "drop")

cdr_conv <- cdr_conv %>%
  left_join(v18_status, by = "BID") %>%
  mutate(
    v18_cdr05 = ifelse(is.na(v18_cdr05), FALSE, v18_cdr05),
    pre_anchor_confirmed = converted & !is.na(first_conv_day) &
      first_conv_day <= anchor_day,
    prevalent_at_anchor = v18_cdr05 | pre_anchor_confirmed
  )

cat("\nCDR conv breakdown:\n")
cat("  Total events: ", sum(cdr_conv$converted), "\n")
cat("  Events with V18 CDR>=0.5: ", sum(cdr_conv$converted & cdr_conv$v18_cdr05), "\n")
cat("  Events confirmed at/before anchor: ", sum(cdr_conv$pre_anchor_confirmed), "\n")
cat("  Events prevalent at anchor: ", sum(cdr_conv$converted & cdr_conv$prevalent_at_anchor), "\n")
cat("  Events incident after anchor: ", sum(cdr_conv$converted & !cdr_conv$prevalent_at_anchor & cdr_conv$time_to_event_days > 0), "\n")

# Primary analysis = incident progression after W48. Exclude prevalent cases
# and require positive raw follow-up time before applying a 0.01-year floor for
# Cox model compatibility. This prevents subjects with only a W48 CDR and no
# post-landmark CDR follow-up from entering as 0.01-year censored observations.
surv_dat_lmer <- dat_lmer %>%
  left_join(cdr_conv, by = "BID") %>%
  filter(!is.na(time_to_event_days)) %>%
  filter(!prevalent_at_anchor) %>%
  filter(time_to_event_days > 0) %>%
  mutate(time_to_event_years = pmax(time_to_event_days / 365.25, 0.01))

surv_group_lmer <- surv_dat_lmer %>% filter(!is.na(LC_group_lmer))
surv_cont_lmer  <- surv_dat_lmer %>% filter(!is.na(learning_slope_resid_z_decr_lmer))
primary_cox_dat <- surv_group_lmer %>% filter(!is.na(SUVR_z) & !is.na(pTau217_z))

saveRDS(surv_group_lmer,
        file.path(corrected_dir,
                  paste0("surv_group_complete_", WINDOW_NAME, "_lmer.rds")))
saveRDS(primary_cox_dat,
        file.path(corrected_dir,
                  paste0("surv_primary_complete_", WINDOW_NAME, "_lmer.rds")))

cat("\nA4 Cox sample (corrected):", nrow(primary_cox_dat),
    "subjects, events =",
    sum(primary_cox_dat$converted), "\n")

cox_g_lmer <- coxph(
  Surv(time_to_event_years, converted) ~ LC_group_lmer + Age_z + Sex_F +
    Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z,
  data = primary_cox_dat
)

cox_c_lmer <- coxph(
  Surv(time_to_event_years, converted) ~ learning_slope_resid_z_decr_lmer + Age_z +
    Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z,
  data = surv_cont_lmer %>% filter(!is.na(SUVR_z) & !is.na(pTau217_z))
)

cox_g_tidy <- tidy(cox_g_lmer, exponentiate = TRUE, conf.int = TRUE)
cox_c_tidy <- tidy(cox_c_lmer, exponentiate = TRUE, conf.int = TRUE)
write_csv(cox_g_tidy, file.path(corrected_dir,
          paste0("CDR_Cox_group_fullyadjusted_", WINDOW_NAME, "_lmer.csv")))
write_csv(cox_c_tidy, file.path(corrected_dir,
          paste0("CDR_Cox_cont_fullyadjusted_", WINDOW_NAME, "_lmer.csv")))

cat("\nCDR Cox (corrected) — group contrast:\n")
print(cox_g_tidy %>% filter(str_detect(term, "LC_group")))
cat("\nCDR Cox (corrected) — continuous PE:\n")
print(cox_c_tidy %>% filter(term == "learning_slope_resid_z_decr_lmer"))

# ----------------------------------------------------------------------------
# 4b) V18 prevalent cases — descriptive summary + one-line sensitivity Cox
# Rationale: these subjects progressed before the post-W48 risk window starts
# and are excluded from the primary analysis. We describe who they are and
# quantify (in a single supplementary line) what happens if they are included.
# ----------------------------------------------------------------------------
prev_full <- dat_lmer %>%
  left_join(cdr_conv, by = "BID") %>%
  filter(!is.na(v18_cdr05))

prev_subset <- prev_full %>% filter(prevalent_at_anchor)

# Descriptive summary
descr_overall <- tibble(
  metric = c("N total with CDR",
             "N V18-prevalent (CDR>=0.5 at W48)",
             "N V18-prevalent in biomarker-complete cohort"),
  value = c(nrow(prev_full),
            nrow(prev_subset),
            nrow(prev_subset %>% filter(!is.na(SUVR_z) & !is.na(pTau217_z))))
)

descr_pe_dist <- prev_subset %>%
  filter(!is.na(LC_group_lmer)) %>%
  count(LC_group_lmer) %>%
  mutate(pct = round(100 * n / sum(n), 1),
         label = paste0(LC_group_lmer, ": ", n, " (", pct, "%)"))

descr_continuous <- prev_subset %>%
  filter(!is.na(LC_group_lmer)) %>%
  summarise(
    Age_mean        = round(mean(Age, na.rm = TRUE), 2),
    Age_sd          = round(sd(Age, na.rm = TRUE), 2),
    Edu_mean        = round(mean(Education, na.rm = TRUE), 2),
    Edu_sd          = round(sd(Education, na.rm = TRUE), 2),
    SUVR_mean       = round(mean(SUVR, na.rm = TRUE), 3),
    SUVR_sd         = round(sd(SUVR, na.rm = TRUE), 3),
    pTau217_mean    = round(mean(pTau217, na.rm = TRUE), 3),
    pTau217_sd      = round(sd(pTau217, na.rm = TRUE), 3),
    APOE4_carrier_pct = round(100 * mean(APOE4 == 1, na.rm = TRUE), 1),
    Female_pct      = round(100 * mean(Sex == "Female", na.rm = TRUE), 1)
  )

write_csv(descr_overall,
          file.path(corrected_dir, "v18_prevalent_descriptive_counts.csv"))
write_csv(descr_pe_dist,
          file.path(corrected_dir, "v18_prevalent_PE_distribution.csv"))
write_csv(descr_continuous,
          file.path(corrected_dir, "v18_prevalent_baseline_characteristics.csv"))

cat("\nV18 prevalent subgroup (descriptive):\n")
cat("  N V18-prevalent: ", nrow(prev_subset), "\n")
cat("  PE group distribution:\n")
print(descr_pe_dist)
cat("  Baseline characteristics:\n")
print(descr_continuous)

# One-line sensitivity Cox: include prevalent cases as t = 0+ events
surv_dat_sens <- dat_lmer %>%
  left_join(cdr_conv, by = "BID") %>%
  filter(!is.na(time_to_event_days)) %>%
  mutate(time_to_event_years = pmax(time_to_event_days / 365.25, 0.01))

surv_group_sens <- surv_dat_sens %>% filter(!is.na(LC_group_lmer))
surv_cont_sens  <- surv_dat_sens %>% filter(!is.na(learning_slope_resid_z_decr_lmer))

cox_g_sens <- coxph(
  Surv(time_to_event_years, converted) ~ LC_group_lmer + Age_z + Sex_F +
    Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z,
  data = surv_group_sens %>% filter(!is.na(SUVR_z) & !is.na(pTau217_z))
)
cox_c_sens <- coxph(
  Surv(time_to_event_years, converted) ~ learning_slope_resid_z_decr_lmer + Age_z +
    Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z,
  data = surv_cont_sens %>% filter(!is.na(SUVR_z) & !is.na(pTau217_z))
)

sens_g_tidy <- tidy(cox_g_sens, exponentiate = TRUE, conf.int = TRUE)
sens_c_tidy <- tidy(cox_c_sens, exponentiate = TRUE, conf.int = TRUE)

sens_row <- tibble(
  analysis = "V18 landmark, INCLUDE prevalent-at-anchor cases as t=0+ events (sensitivity)",
  n_subjects = nrow(surv_group_sens %>% filter(!is.na(SUVR_z) & !is.na(pTau217_z))),
  n_events = sum((surv_group_sens %>% filter(!is.na(SUVR_z) & !is.na(pTau217_z)))$converted),
  HR_low_vs_high_PE = round((sens_g_tidy %>% filter(str_detect(term, "Low_PE")))$estimate[1], 3),
  HR_low_vs_high_PE_ci_lo = round((sens_g_tidy %>% filter(str_detect(term, "Low_PE")))$conf.low[1], 3),
  HR_low_vs_high_PE_ci_hi = round((sens_g_tidy %>% filter(str_detect(term, "Low_PE")))$conf.high[1], 3),
  HR_low_vs_high_PE_p = signif((sens_g_tidy %>% filter(str_detect(term, "Low_PE")))$p.value[1], 3),
  HR_continuous_per_SD = round((sens_c_tidy %>% filter(term == "learning_slope_resid_z_decr_lmer"))$estimate[1], 3),
  HR_continuous_ci_lo = round((sens_c_tidy %>% filter(term == "learning_slope_resid_z_decr_lmer"))$conf.low[1], 3),
  HR_continuous_ci_hi = round((sens_c_tidy %>% filter(term == "learning_slope_resid_z_decr_lmer"))$conf.high[1], 3),
  HR_continuous_p = signif((sens_c_tidy %>% filter(term == "learning_slope_resid_z_decr_lmer"))$p.value[1], 3)
)
write_csv(sens_row, file.path(corrected_dir, "v18_prevalent_sensitivity_one_line.csv"))

cat("\nOne-line sensitivity (INCLUDE V18-prevalent):\n")
print(sens_row)

# ----------------------------------------------------------------------------
# 5) Incremental discrimination (apparent + bootstrap + nested CV)
# ----------------------------------------------------------------------------
cox_dat <- surv_cont_lmer %>%
  filter(!is.na(SUVR_z) & !is.na(pTau217_z))

base_form <- Surv(time_to_event_years, converted) ~ Age_z + Sex_F + Education_z +
  APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z
plus_form <- update(base_form, . ~ . + learning_slope_resid_z_decr_lmer)

c_index <- function(fit, dat) {
  lp <- predict(fit, newdata = dat, type = "lp")
  s  <- Surv(dat$time_to_event_years, dat$converted)
  unname(survival::concordance(s ~ lp, reverse = TRUE)$concordance)
}

base_fit <- coxph(base_form, data = cox_dat)
plus_fit <- coxph(plus_form, data = cox_dat)
c_base <- c_index(base_fit, cox_dat)
c_plus <- c_index(plus_fit, cox_dat)
delta_c_apparent <- c_plus - c_base

cat("\nIncremental discrimination (corrected, apparent in-sample):\n")
cat("  base  C-index:", round(c_base, 4), "\n")
cat("  +PE   C-index:", round(c_plus, 4), "\n")
cat("  delta C      :", round(delta_c_apparent, 4), "\n")

# Bootstrap delta C (1000 paired resamples)
set.seed(42)
B <- 1000
boot_dc <- numeric(B)
for (b in seq_len(B)) {
  idx <- sample.int(nrow(cox_dat), replace = TRUE)
  d_b <- cox_dat[idx, , drop = FALSE]
  bf  <- tryCatch(coxph(base_form, data = d_b), error = function(e) NULL)
  pf  <- tryCatch(coxph(plus_form, data = d_b), error = function(e) NULL)
  if (is.null(bf) || is.null(pf)) { boot_dc[b] <- NA_real_; next }
  boot_dc[b] <- c_index(pf, d_b) - c_index(bf, d_b)
}
boot_dc <- boot_dc[!is.na(boot_dc)]
boot_summary <- tibble(
  B_kept   = length(boot_dc),
  mean_dC  = mean(boot_dc),
  ci_lo    = quantile(boot_dc, 0.025),
  ci_hi    = quantile(boot_dc, 0.975),
  pct_pos  = mean(boot_dc > 0)
)
write_csv(boot_summary, file.path(corrected_dir,
          paste0("survival_bootstrap_delta_c_", WINDOW_NAME, "_lmer.csv")))
cat("\nBootstrap delta C: mean=", round(boot_summary$mean_dC, 4),
    " 95% CI [", round(boot_summary$ci_lo, 4), ",",
    round(boot_summary$ci_hi, 4), "] pct>0=", round(boot_summary$pct_pos, 3), "\n")

# Nested 5-fold CV (R = 10 repeats), re-deriving PE z within each training fold
nested_cv_dc <- function(K = 5, R = 10) {
  out <- numeric()
  for (rep_i in seq_len(R)) {
    set.seed(100 + rep_i)
    folds <- sample(rep(seq_len(K), length.out = nrow(cox_dat)))
    for (k in seq_len(K)) {
      tr <- which(folds != k); te <- which(folds == k)
      tr_dat <- cox_dat[tr, , drop = FALSE]
      te_dat <- cox_dat[te, , drop = FALSE]
      # Re-z PE within training set, apply same transform to test
      mu <- mean(tr_dat$learning_slope_resid_z_decr_lmer, na.rm = TRUE)
      sd <- sd(tr_dat$learning_slope_resid_z_decr_lmer, na.rm = TRUE)
      tr_dat$pe_cv <- (tr_dat$learning_slope_resid_z_decr_lmer - mu) / sd
      te_dat$pe_cv <- (te_dat$learning_slope_resid_z_decr_lmer - mu) / sd
      bf <- tryCatch(coxph(update(base_form, . ~ .), data = tr_dat),
                     error = function(e) NULL)
      pf <- tryCatch(coxph(update(base_form, . ~ . + pe_cv), data = tr_dat),
                     error = function(e) NULL)
      if (is.null(bf) || is.null(pf)) next
      out <- c(out, c_index(pf, te_dat) - c_index(bf, te_dat))
    }
  }
  out
}

cv_dc <- nested_cv_dc()
cv_summary <- tibble(
  n_folds = length(cv_dc),
  mean_dC = mean(cv_dc, na.rm = TRUE),
  min_dC  = min(cv_dc, na.rm = TRUE),
  max_dC  = max(cv_dc, na.rm = TRUE)
)
write_csv(cv_summary, file.path(corrected_dir,
          paste0("survival_cv_delta_c_", WINDOW_NAME, "_lmer.csv")))
cat("\nNested CV delta C: mean=", round(cv_summary$mean_dC, 4),
    " range [", round(cv_summary$min_dC, 4), ",",
    round(cv_summary$max_dC, 4), "]\n")

# ----------------------------------------------------------------------------
# 6) Old vs new headline comparison
# ----------------------------------------------------------------------------
# Read the old locked outputs from source_data/locked_outputs/ and tabulate.

old_dir <- file.path(pkg_dir, "outputs", "locked", "a4_old_not_packaged")

read_safe <- function(path) {
  if (file.exists(path)) suppressMessages(read_csv(path, show_col_types = FALSE))
  else NULL
}

old_pacc_g <- read_safe(file.path(old_dir, "PACC_LMM_group_model5_tidy_w1_21_lmer.csv"))
old_pacc_c <- read_safe(file.path(old_dir, "PACC_LMM_cont_model5_term_w1_21_lmer.csv"))
old_cox_g  <- read_safe(file.path(old_dir, "CDR_Cox_group_fullyadjusted_w1_21_lmer.csv"))
old_cox_c  <- read_safe(file.path(old_dir, "CDR_Cox_cont_fullyadjusted_w1_21_lmer.csv"))
old_inc    <- read_safe(file.path(old_dir, "incremental_value_w1_21_lmer.csv"))

old_boot_all <- read_safe(file.path(old_dir, "survival_bootstrap_delta_c.csv"))
old_boot <- if (!is.null(old_boot_all)) {
  old_boot_all %>%
    filter(window == "w1_21", grepl("lmer", method), type == "group") %>%
    rename(mean_dC = delta_c_boot_mean,
           ci_lo   = delta_c_boot_p025,
           ci_hi   = delta_c_boot_p975) %>%
    slice(1)
} else NULL

old_cv_all <- read_safe(file.path(old_dir, "survival_cv_delta_c.csv"))
old_cv <- if (!is.null(old_cv_all)) {
  old_cv_all %>%
    filter(window == "w1_21", grepl("lmer", method), type == "group") %>%
    rename(mean_dC = delta_c_cv_mean) %>%
    slice(1)
} else NULL

old_surv_path <- file.path(old_dir, "surv_group_complete_w1_21_lmer.rds")
old_surv <- if (file.exists(old_surv_path)) readRDS(old_surv_path) else NULL

pull_estimate <- function(df, term_pattern) {
  if (is.null(df)) return(rep(NA_real_, 3))
  row <- df %>% filter(str_detect(term, term_pattern)) %>% slice(1)
  if (nrow(row) == 0) return(rep(NA_real_, 3))
  c(estimate = row$estimate,
    conf.low = row$conf.low,
    conf.high = row$conf.high)
}

new_pacc_g <- g_tidy %>% filter(str_detect(term, "Low_PE"))
new_pacc_c <- c_tidy
new_cox_g  <- cox_g_tidy %>% filter(str_detect(term, "LC_group_lmerLow_PE"))
new_cox_c  <- cox_c_tidy %>% filter(term == "learning_slope_resid_z_decr_lmer")

n_events_old <- if (!is.null(old_surv)) sum(old_surv$converted) else NA
n_subj_old   <- if (!is.null(old_surv)) sum(!is.na(old_surv$SUVR_z) & !is.na(old_surv$pTau217_z)) else NA
median_t_old <- if (!is.null(old_surv)) {
                  ev <- old_surv$time_to_event_years[old_surv$converted]
                  median(ev, na.rm = TRUE)
                } else NA

# Reuse the biomarker-complete primary Cox sample saved above as a standalone
# RDS so downstream lock files can cite a file whose row count matches the
# manuscript's 719-subject primary analysis directly.
cox_n_new <- nrow(primary_cox_dat)
cox_e_new <- sum(primary_cox_dat$converted)
# Median event time for the PRIMARY Cox sample (biomarker-complete, V18-prevalent
# excluded); NOT over all raw CDR conversions in cdr_conv.
median_t_new <- median(primary_cox_dat$time_to_event_years[primary_cox_dat$converted], na.rm = TRUE)

comp <- tibble(
  metric = c(
    "Cox sample n",
    "Cox events n",
    "Median event time (years)",
    "PACC LMM beta: Low PE x time (per yr)",
    "PACC LMM beta: continuous PE x time (per yr per SD)",
    "CDR Cox HR: Low vs High PE",
    "CDR Cox HR: continuous PE per SD",
    "Bootstrap mean delta C",
    "Bootstrap delta C 95% CI lo",
    "Bootstrap delta C 95% CI hi",
    "Nested CV mean delta C"
  ),
  old = c(
    n_subj_old,
    n_events_old,
    median_t_old,
    pull_estimate(old_pacc_g, "Low_PE")["estimate"],
    pull_estimate(old_pacc_c, "learning_slope")["estimate"],
    pull_estimate(old_cox_g,  "Low_PE")["estimate"],
    pull_estimate(old_cox_c,  "learning_slope")["estimate"],
    if (!is.null(old_boot)) old_boot$mean_dC[1] else NA,
    if (!is.null(old_boot)) old_boot$ci_lo[1]   else NA,
    if (!is.null(old_boot)) old_boot$ci_hi[1]   else NA,
    if (!is.null(old_cv))   old_cv$mean_dC[1]   else NA
  ),
  new = c(
    cox_n_new,
    cox_e_new,
    median_t_new,
    if (nrow(new_pacc_g)) new_pacc_g$estimate[1] else NA,
    if (nrow(new_pacc_c)) new_pacc_c$estimate[1] else NA,
    if (nrow(new_cox_g))  new_cox_g$estimate[1]  else NA,
    if (nrow(new_cox_c))  new_cox_c$estimate[1]  else NA,
    boot_summary$mean_dC,
    boot_summary$ci_lo,
    boot_summary$ci_hi,
    cv_summary$mean_dC
  )
) %>%
  mutate(
    pct_change = (new - old) / abs(old) * 100,
    abs_diff   = new - old
  )

write_csv(comp, file.path(corrected_dir, "old_vs_new_comparison.csv"))

cat("\n============================================\n")
cat("OLD vs NEW HEADLINE COMPARISON\n")
cat("============================================\n")
print(comp, n = Inf)

cat("\nFinished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("All corrected outputs in:", corrected_dir, "\n")
