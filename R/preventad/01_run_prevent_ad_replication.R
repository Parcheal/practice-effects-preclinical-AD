

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
# PREVENT-AD external validation: primary analyses + sensitivity enhancements
# - Rebuild cohort and PE estimates from raw PREVENT-AD data
# - Clarify sample accounting and censoring proxy
# - Update LMM / Cox analyses and high-priority sensitivity analyses
# ============================================================================

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(readr)
  library(stringr)
  library(survival)
  library(tibble)
  library(tidyr)
})

if (!requireNamespace("broom.mixed", quietly = TRUE)) {
  stop("Package 'broom.mixed' is required but not installed.")
}

out_dir <- file.path(pkg_root, "verify", "preventad")
data_dir <- file.path(data_root, "preventad", "Clinical_data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pet_threshold <- 1.54
csf_abeta42_threshold <- 970

parse_visit_months <- function(x) {
  out <- suppressWarnings(as.numeric(str_extract(x, "\\d+")))
  out[!is.na(x) & x == "BL00"] <- 0
  out
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, formatC(x, digits = digits, format = "f"))
}

write_tidy_csv <- function(model, path, exponentiate = FALSE) {
  if (inherits(model, "cox.zph")) {
    tab <- as.data.frame(model$table) %>%
      rownames_to_column("term") %>%
      transmute(
        term = term,
        chisq = `chisq`,
        df = `df`,
        p.value = `p`
      )
    write_csv(tab, path)
    return(invisible(tab))
  }

  tidy_out <- if (inherits(model, "lmerMod") || inherits(model, "lmerModLmerTest")) {
    broom.mixed::tidy(model, effects = c("fixed", "ran_pars"), conf.int = TRUE)
  } else {
    broom::tidy(model, conf.int = TRUE, exponentiate = exponentiate)
  }
  write_csv(tidy_out, path)
  invisible(tidy_out)
}

harrell_c <- function(model) {
  as.numeric(summary(model)$concordance[1])
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[[1]]
}

build_survival_dataset <- function(dat, time_var, label, strict_landmark = TRUE) {
  out <- dat %>%
    mutate(
      converted = ifelse(is.na(converted), FALSE, converted),
      censor_source = label,
      time_to_event_months = ifelse(converted, mci_visit_months - 24, .data[[time_var]] - 24),
      time_to_event_years = pmax(time_to_event_months / 12, 0.01)
    ) %>%
    filter(time_to_event_months > 0)
  if (isTRUE(strict_landmark)) {
    # Drop anyone with an MCI diagnosis at or before the FU24 landmark.
    # Without this, FU24-prevalent MCI cases are silently censored as
    # event = 0 and remain in the risk-set, which is not a clean landmark.
    out <- out %>%
      filter(is.na(mci_visit_months) | mci_visit_months > 24)
  }
  out
}

fit_primary_cox <- function(data, formula_rhs) {
  coxph(
    as.formula(paste0("Surv(time_to_event_years, converted) ~ ", formula_rhs)),
    data = data
  )
}

tidy_cox_exp <- function(model) {
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
}

make_landmark_dataset <- function(data, landmark_years) {
  data %>%
    filter(time_to_event_years > landmark_years) %>%
    mutate(
      time_to_event_years = time_to_event_years - landmark_years
    )
}

# ----------------------------
# Load source data
# ----------------------------
rbans <- read_csv(file.path(data_dir, "RBANS_Registered_PREVENTAD_2024-08-01.csv"), show_col_types = FALSE) %>%
  mutate(
    visit_months = parse_visit_months(Visit_label),
    visit_years = visit_months / 12
  )

demo <- read_csv(file.path(data_dir, "Demographics_Registered_PREVENTAD_2024-08-01.csv"), show_col_types = FALSE) %>%
  transmute(
    CONP_ID,
    CONP_CandID,
    Sex,
    Education = as.numeric(Education_years),
    INTREPAD_arm = INTREPAD_Tx_assignment
  )

genetics <- read_csv(file.path(data_dir, "Genetics_Registered_PREVENTAD.csv"), show_col_types = FALSE) %>%
  transmute(
    CONP_ID,
    APOE_raw = APOE,
    APOE4 = case_when(
      is.na(APOE) ~ NA_real_,
      str_detect(APOE, "(^|\\s)4($|\\s)") ~ 1,
      TRUE ~ 0
    )
  )

pet <- read_csv(file.path(data_dir, "PET_NAV_SUVR_ref-cerebellumCortex_Registered_PREVENTAD_2024-08-01.csv"), show_col_types = FALSE) %>%
  transmute(
    CONP_ID,
    amyloid_SUVR = as.numeric(amyloid_index_SUVR),
    amyloid_pos_pet = ifelse(is.na(amyloid_SUVR), NA, amyloid_SUVR >= pet_threshold)
  )

csf <- read_csv(file.path(data_dir, "CSF_proteins_Registered_PREVENTAD.csv"), show_col_types = FALSE) %>%
  mutate(
    visit_months = parse_visit_months(Visit_label),
    date_taken = as.Date(paste0(Date_taken, "-01"))
  ) %>%
  arrange(CONP_ID, visit_months, date_taken) %>%
  group_by(CONP_ID) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    CONP_ID,
    CSF_ptau181 = as.numeric(ptau),
    CSF_Abeta42 = as.numeric(Amyloid_beta_1_42),
    CSF_tau = as.numeric(tau),
    CSF_first_visit_months = visit_months,
    amyloid_pos_csf = ifelse(is.na(CSF_Abeta42), NA, CSF_Abeta42 < csf_abeta42_threshold)
  )

cmoca <- read_csv(file.path(data_dir, "EL_CDR_MoCA_Registered_PREVENTAD.csv"), show_col_types = FALSE) %>%
  mutate(
    visit_months = parse_visit_months(Visit_label),
    date_taken = as.Date(paste0(Date_taken, "-01"))
  )

moca_bl <- cmoca %>%
  filter(visit_months == 0) %>%
  group_by(CONP_ID) %>%
  summarise(
    baseline_MOCA = first_non_missing(as.numeric(MOCA_total_score)),
    baseline_CDR = first_non_missing(as.numeric(CDR_total_score)),
    .groups = "drop"
  )

mci_first <- read_csv(file.path(data_dir, "MCI_Registered_PREVENTAD_2024-08-01.csv"), show_col_types = FALSE) %>%
  mutate(
    mci_visit_months = parse_visit_months(initial_MCI_visit_label),
    review_visit_months = parse_visit_months(MCI_review_visit_label)
  ) %>%
  arrange(CONP_ID, mci_visit_months, review_visit_months) %>%
  group_by(CONP_ID) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    CONP_ID,
    mci_visit_months = mci_visit_months,
    converted = !is.na(mci_visit_months) & mci_visit_months > 24
  )

# ----------------------------
# Build analytic cohort and PE estimates
# ----------------------------
rbans_bl <- rbans %>%
  filter(visit_months == 0, !is.na(total_scale_index_score)) %>%
  group_by(CONP_ID) %>%
  summarise(
    CONP_CandID = first_non_missing(CONP_CandID),
    baseline_RBANS_total = first_non_missing(as.numeric(total_scale_index_score)),
    baseline_age_years = first_non_missing(as.numeric(Candidate_Age)),
    .groups = "drop"
  )

rbans_early <- rbans %>%
  filter(visit_months %in% c(0, 12, 24), !is.na(total_scale_index_score)) %>%
  mutate(total_scale_index_score = as.numeric(total_scale_index_score))

early_summary <- rbans_early %>%
  group_by(CONP_ID) %>%
  summarise(
    n_pe_visits = n_distinct(visit_months),
    has_bl = any(visit_months == 0),
    has_fu12 = any(visit_months == 12),
    has_fu24 = any(visit_months == 24),
    has_all_3 = has_bl & has_fu12 & has_fu24,
    .groups = "drop"
  )

cohort_baseline <- early_summary %>%
  filter(n_pe_visits >= 2) %>%
  left_join(rbans_bl, by = "CONP_ID") %>%
  left_join(demo, by = "CONP_ID") %>%
  left_join(genetics, by = "CONP_ID") %>%
  left_join(pet, by = "CONP_ID") %>%
  left_join(csf, by = "CONP_ID") %>%
  left_join(moca_bl, by = "CONP_ID") %>%
  mutate(
    amyloid_positive = case_when(
      !is.na(amyloid_pos_pet) ~ amyloid_pos_pet,
      !is.na(amyloid_pos_csf) ~ amyloid_pos_csf,
      TRUE ~ NA
    ),
    has_amyloid_data = !is.na(amyloid_pos_pet) | !is.na(amyloid_pos_csf),
    Age_z = safe_scale(baseline_age_years),
    Education_z = safe_scale(Education),
    baseline_RBANS_z = safe_scale(baseline_RBANS_total),
    baseline_MOCA_z = safe_scale(baseline_MOCA),
    SUVR_z = safe_scale(amyloid_SUVR),
    CSF_ptau181_z = safe_scale(CSF_ptau181),
    CSF_Abeta42_z = safe_scale(CSF_Abeta42),
    Sex_F = factor(Sex, levels = c("Male", "Female")),
    APOE4_carrier = factor(ifelse(APOE4 == 1, "Carrier", "Non-carrier"), levels = c("Non-carrier", "Carrier"))
  ) %>%
  arrange(CONP_ID)

write_csv(cohort_baseline, file.path(out_dir, "cohort_baseline.csv"))

pe_ids <- cohort_baseline %>%
  filter(has_bl, n_pe_visits >= 2, !is.na(baseline_RBANS_total)) %>%
  pull(CONP_ID)

rbans_pe <- rbans_early %>%
  filter(CONP_ID %in% pe_ids) %>%
  mutate(time_years = visit_months / 12) %>%
  select(CONP_ID, time_years, total_scale_index_score)

pe_input <- cohort_baseline %>%
  filter(CONP_ID %in% pe_ids) %>%
  select(CONP_ID, baseline_RBANS_from_pe = baseline_RBANS_total)

pe_ols <- rbans_pe %>%
  group_by(CONP_ID) %>%
  group_modify(~ {
    fit <- lm(total_scale_index_score ~ time_years, data = .x)
    tibble(pe_slope_raw = coef(fit)[["time_years"]])
  }) %>%
  ungroup() %>%
  left_join(pe_input, by = "CONP_ID")

ols_resid_fit <- lm(pe_slope_raw ~ baseline_RBANS_from_pe, data = pe_ols)
pe_ols <- pe_ols %>%
  mutate(
    pe_resid = resid(ols_resid_fit),
    pe_resid_z = safe_scale(pe_resid),
    pe_resid_z_decr = -pe_resid_z,
    PE_group = case_when(
      pe_resid_z > 0.2 ~ "High_PE",
      pe_resid_z < -0.2 ~ "Low_PE",
      TRUE ~ "Moderate_PE"
    )
  )

write_csv(pe_ols, file.path(out_dir, "PE_OLS.csv"))

lmer_pe_fit <- lmer(total_scale_index_score ~ time_years + (1 + time_years | CONP_ID), data = rbans_pe, REML = FALSE)
lmer_re <- ranef(lmer_pe_fit)$CONP_ID %>%
  rownames_to_column("CONP_ID") %>%
  transmute(CONP_ID, rand_slope = time_years)
fixed_slope <- unname(fixef(lmer_pe_fit)[["time_years"]])

pe_lmer <- pe_input %>%
  left_join(lmer_re, by = "CONP_ID") %>%
  mutate(pe_slope_raw_lmer = fixed_slope + rand_slope)

lmer_resid_fit <- lm(pe_slope_raw_lmer ~ baseline_RBANS_from_pe, data = pe_lmer)
pe_lmer <- pe_lmer %>%
  mutate(
    pe_resid_lmer = resid(lmer_resid_fit),
    pe_resid_z_lmer = safe_scale(pe_resid_lmer),
    pe_resid_z_decr_lmer = -pe_resid_z_lmer,
    PE_group_lmer = case_when(
      pe_resid_z_lmer > 0.2 ~ "High_PE",
      pe_resid_z_lmer < -0.2 ~ "Low_PE",
      TRUE ~ "Moderate_PE"
    )
  )

write_csv(pe_lmer, file.path(out_dir, "PE_LMER.csv"))

dat_lmer <- cohort_baseline %>%
  inner_join(pe_lmer, by = "CONP_ID") %>%
  mutate(
    PE_group_lmer = factor(PE_group_lmer, levels = c("High_PE", "Moderate_PE", "Low_PE")),
    PE_group_tertile = factor(
      c("Low_PE", "Moderate_PE", "High_PE")[ntile(pe_resid_z_lmer, 3)],
      levels = c("High_PE", "Moderate_PE", "Low_PE")
    )
  ) %>%
  arrange(CONP_ID)

write_csv(dat_lmer, file.path(out_dir, "analysis_baseline_lmer.csv"))

# ----------------------------
# Sample accounting and descriptive tables
# ----------------------------
rbans_last <- rbans %>%
  filter(CONP_ID %in% cohort_baseline$CONP_ID) %>%
  group_by(CONP_ID) %>%
  summarise(last_rbans_visit_months = max(visit_months, na.rm = TRUE), .groups = "drop")

cmoca_last <- cmoca %>%
  filter(CONP_ID %in% cohort_baseline$CONP_ID) %>%
  group_by(CONP_ID) %>%
  summarise(last_cmoca_visit_months = max(visit_months, na.rm = TRUE), .groups = "drop")

contact_proxy <- full_join(rbans_last, cmoca_last, by = "CONP_ID") %>%
  mutate(
    across(c(last_rbans_visit_months, last_cmoca_visit_months), ~ ifelse(is.infinite(.x), NA, .x)),
    last_contact_proxy_months = pmax(last_rbans_visit_months, last_cmoca_visit_months, na.rm = TRUE),
    last_contact_proxy_months = ifelse(is.infinite(last_contact_proxy_months), NA, last_contact_proxy_months),
    proxy_equals_rbans = dplyr::near(last_contact_proxy_months, last_rbans_visit_months)
  )

write_csv(contact_proxy, file.path(out_dir, "last_contact_proxy.csv"))

surv_seed <- cohort_baseline %>%
  left_join(pe_lmer, by = "CONP_ID") %>%
  left_join(dat_lmer %>% select(CONP_ID, PE_group_tertile), by = "CONP_ID") %>%
  left_join(mci_first, by = "CONP_ID") %>%
  left_join(contact_proxy, by = "CONP_ID")

surv_primary <- build_survival_dataset(surv_seed, "last_contact_proxy_months", "last_contact_proxy")
surv_rbans <- build_survival_dataset(surv_seed, "last_rbans_visit_months", "last_rbans")

write_csv(surv_primary, file.path(out_dir, "survival_dataset_primary.csv"))
write_csv(surv_rbans, file.path(out_dir, "survival_dataset_lastcontact.csv"))

sample_accounting <- tibble(
  step = c(
    "All registered subjects",
    "Have baseline RBANS",
    "Have >=2 early RBANS visits",
    "PE estimable (LMER + baseline residualization)",
    "Follow-up sample (>FU24)",
    "Primary Cox complete-case sample"
  ),
  n = c(
    n_distinct(rbans$CONP_ID),
    n_distinct(rbans$CONP_ID[rbans$visit_months == 0]),
    nrow(cohort_baseline),
    nrow(dat_lmer),
    nrow(surv_primary),
    sum(
      with(
        surv_primary,
        !is.na(Age_z) & !is.na(Education_z) & !is.na(Sex_F) &
          !is.na(APOE4_carrier) & !is.na(baseline_RBANS_z) & !is.na(PE_group_lmer)
      )
    )
  )
)
write_csv(sample_accounting, file.path(out_dir, "sample_accounting.csv"))

cohort_flow <- tibble(
  step = c(
    "All registered subjects",
    "Have baseline RBANS",
    "Have >= 2 early RBANS visits (BL/FU12/FU24)",
    "PE estimable (requires BL + >=2 early visits)",
    "Survival sample (follow-up > FU24)",
    "MCI events post-FU24",
    "Primary Cox complete-case sample"
  ),
  n = c(
    n_distinct(rbans$CONP_ID),
    n_distinct(rbans$CONP_ID[rbans$visit_months == 0]),
    nrow(cohort_baseline),
    nrow(dat_lmer),
    nrow(surv_primary),
    sum(surv_primary$converted, na.rm = TRUE),
    sum(
      with(
        surv_primary,
        !is.na(Age_z) & !is.na(Education_z) & !is.na(Sex_F) &
          !is.na(APOE4_carrier) & !is.na(baseline_RBANS_z) & !is.na(PE_group_lmer)
      )
    )
  )
)
write_csv(cohort_flow, file.path(out_dir, "cohort_flow.csv"))

censoring_proxy_summary <- tibble(
  metric = c(
    "n_with_rbans_followup",
    "n_with_cmoca_followup",
    "n_proxy_later_than_rbans",
    "pct_proxy_equals_rbans",
    "median_proxy_minus_rbans_months"
  ),
  value = c(
    sum(!is.na(contact_proxy$last_rbans_visit_months)),
    sum(!is.na(contact_proxy$last_cmoca_visit_months)),
    sum(contact_proxy$last_contact_proxy_months > contact_proxy$last_rbans_visit_months, na.rm = TRUE),
    round(100 * mean(contact_proxy$proxy_equals_rbans, na.rm = TRUE), 1),
    median(contact_proxy$last_contact_proxy_months - contact_proxy$last_rbans_visit_months, na.rm = TRUE)
  )
)
write_csv(censoring_proxy_summary, file.path(out_dir, "censoring_proxy_summary.csv"))

table1_complete <- bind_rows(
  dat_lmer %>%
    summarise(
      PE_group_lmer = "Overall",
      N = n(),
      Age_mean = round(mean(baseline_age_years, na.rm = TRUE), 1),
      Age_sd = round(sd(baseline_age_years, na.rm = TRUE), 1),
      Pct_female = round(100 * mean(Sex == "Female", na.rm = TRUE), 1),
      Education_mean = round(mean(Education, na.rm = TRUE), 1),
      Education_sd = round(sd(Education, na.rm = TRUE), 1),
      Pct_APOE4 = round(100 * mean(APOE4 == 1, na.rm = TRUE), 1),
      RBANS_baseline_mean = round(mean(baseline_RBANS_total, na.rm = TRUE), 1),
      RBANS_baseline_sd = round(sd(baseline_RBANS_total, na.rm = TRUE), 1),
      N_amyloid_data = sum(has_amyloid_data, na.rm = TRUE),
      Pct_amyloid_pos = round(100 * mean(amyloid_positive, na.rm = TRUE), 1),
      SUVR_mean = round(mean(amyloid_SUVR, na.rm = TRUE), 3),
      SUVR_sd = round(sd(amyloid_SUVR, na.rm = TRUE), 3)
    ),
  dat_lmer %>%
    group_by(PE_group_lmer) %>%
    summarise(
      N = n(),
      Age_mean = round(mean(baseline_age_years, na.rm = TRUE), 1),
      Age_sd = round(sd(baseline_age_years, na.rm = TRUE), 1),
      Pct_female = round(100 * mean(Sex == "Female", na.rm = TRUE), 1),
      Education_mean = round(mean(Education, na.rm = TRUE), 1),
      Education_sd = round(sd(Education, na.rm = TRUE), 1),
      Pct_APOE4 = round(100 * mean(APOE4 == 1, na.rm = TRUE), 1),
      RBANS_baseline_mean = round(mean(baseline_RBANS_total, na.rm = TRUE), 1),
      RBANS_baseline_sd = round(sd(baseline_RBANS_total, na.rm = TRUE), 1),
      N_amyloid_data = sum(has_amyloid_data, na.rm = TRUE),
      Pct_amyloid_pos = round(100 * mean(amyloid_positive, na.rm = TRUE), 1),
      SUVR_mean = round(mean(amyloid_SUVR, na.rm = TRUE), 3),
      SUVR_sd = round(sd(amyloid_SUVR, na.rm = TRUE), 3),
      .groups = "drop"
    )
)

write_csv(table1_complete, file.path(out_dir, "Table1_complete.csv"))
write_csv(filter(table1_complete, PE_group_lmer != "Overall"), file.path(out_dir, "Table1_by_PE_group.csv"))

# ----------------------------
# Longitudinal RBANS analyses
# ----------------------------
rbans_long <- rbans %>%
  filter(visit_months > 24, !is.na(total_scale_index_score)) %>%
  transmute(
    CONP_ID,
    visit_months,
    time_years = (visit_months - 24) / 12,
    RBANS_version = factor(RBANS_version),
    total_scale_index_score = as.numeric(total_scale_index_score),
    immediate_memory_index_score = as.numeric(immediate_memory_index_score),
    visuospatial_constructional_index_score = as.numeric(visuospatial_constructional_index_score),
    language_index_score = as.numeric(language_index_score),
    attention_index_score = as.numeric(attention_index_score),
    delayed_memory_index_score = as.numeric(delayed_memory_index_score)
  ) %>%
  inner_join(dat_lmer, by = "CONP_ID")

lmm_group_basic <- lmer(
  total_scale_index_score ~ time_years * PE_group_lmer +
    Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z +
    (1 + time_years | CONP_ID),
  data = rbans_long,
  REML = FALSE
)
write_tidy_csv(lmm_group_basic, file.path(out_dir, "RBANS_LMM_group_basic_tidy.csv"))

lmm_cont_basic <- lmer(
  total_scale_index_score ~ time_years * pe_resid_z_decr_lmer +
    Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z +
    (1 + time_years | CONP_ID),
  data = rbans_long,
  REML = FALSE
)
write_tidy_csv(lmm_cont_basic, file.path(out_dir, "RBANS_LMM_cont_basic_tidy.csv"))

lmm_group_basic_versionadj <- lmer(
  total_scale_index_score ~ time_years * PE_group_lmer + RBANS_version +
    Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z +
    (1 + time_years | CONP_ID),
  data = rbans_long,
  REML = FALSE
)
write_tidy_csv(lmm_group_basic_versionadj, file.path(out_dir, "RBANS_LMM_group_basic_versionadj_tidy.csv"))

rbans_long_pet <- rbans_long %>% filter(!is.na(SUVR_z))
lmm_group_amyloid <- lmer(
  total_scale_index_score ~ time_years * PE_group_lmer +
    Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + baseline_RBANS_z +
    (1 + time_years | CONP_ID),
  data = rbans_long_pet,
  REML = FALSE
)
write_tidy_csv(lmm_group_amyloid, file.path(out_dir, "RBANS_LMM_group_amyloid_tidy.csv"))

rbans_long_amyloidpos <- rbans_long %>% filter(amyloid_positive %in% TRUE)
lmm_group_amyloidpos <- lmer(
  total_scale_index_score ~ time_years * PE_group_lmer +
    Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z +
    (1 + time_years | CONP_ID),
  data = rbans_long_amyloidpos,
  REML = FALSE
)
write_tidy_csv(lmm_group_amyloidpos, file.path(out_dir, "RBANS_LMM_group_amyloidpos_tidy.csv"))

rbans_long_amyloiddat <- rbans_long %>% filter(!is.na(amyloid_positive))
lmm_pe_amyloid_interaction <- lmer(
  total_scale_index_score ~ time_years * pe_resid_z_decr_lmer * amyloid_positive +
    Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z +
    (1 + time_years | CONP_ID),
  data = rbans_long_amyloiddat,
  REML = FALSE
)
write_tidy_csv(lmm_pe_amyloid_interaction, file.path(out_dir, "RBANS_LMM_PE_amyloid_interaction.csv"))

subtest_specs <- tribble(
  ~subtest, ~var,
  "Immediate Memory", "immediate_memory_index_score",
  "Visuospatial", "visuospatial_constructional_index_score",
  "Language", "language_index_score",
  "Attention", "attention_index_score",
  "Delayed Memory", "delayed_memory_index_score"
)

subtest_results <- lapply(seq_len(nrow(subtest_specs)), function(i) {
  var_name <- subtest_specs$var[[i]]
  subtest_name <- subtest_specs$subtest[[i]]
  dat_sub <- rbans_long %>% filter(!is.na(.data[[var_name]]))
  fit_sub <- lmer(
    as.formula(paste0(
      var_name,
      " ~ time_years * pe_resid_z_decr_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z + (1 + time_years | CONP_ID)"
    )),
    data = dat_sub,
    REML = FALSE
  )
  broom.mixed::tidy(fit_sub, effects = "fixed", conf.int = TRUE) %>%
    filter(term == "time_years:pe_resid_z_decr_lmer") %>%
    mutate(subtest = subtest_name)
})

write_csv(bind_rows(subtest_results), file.path(out_dir, "RBANS_subtest_PE_effects.csv"))

# ----------------------------
# Survival analyses
# ----------------------------
surv_complete <- surv_primary %>%
  mutate(
    PE_group_lmer = factor(PE_group_lmer, levels = c("High_PE", "Moderate_PE", "Low_PE")),
    PE_group_tertile = factor(PE_group_tertile, levels = c("High_PE", "Moderate_PE", "Low_PE"))
  ) %>%
  filter(
    !is.na(Age_z),
    !is.na(Education_z),
    !is.na(Sex_F),
    !is.na(APOE4_carrier),
    !is.na(baseline_RBANS_z),
    !is.na(PE_group_lmer)
  )

cox_group_basic <- fit_primary_cox(
  surv_complete,
  "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z"
)
write_csv(tidy_cox_exp(cox_group_basic), file.path(out_dir, "MCI_Cox_group_basic.csv"))

cox_cont_basic <- fit_primary_cox(
  surv_complete,
  "pe_resid_z_decr_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z"
)
write_csv(tidy_cox_exp(cox_cont_basic), file.path(out_dir, "MCI_Cox_cont_basic.csv"))

surv_amyloid <- surv_complete %>% filter(!is.na(SUVR_z))
cox_group_amyloid <- fit_primary_cox(
  surv_amyloid,
  "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + baseline_RBANS_z"
)
write_csv(tidy_cox_exp(cox_group_amyloid), file.path(out_dir, "MCI_Cox_group_amyloid.csv"))

cox_cont_amyloid <- fit_primary_cox(
  surv_amyloid,
  "pe_resid_z_decr_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + baseline_RBANS_z"
)
write_csv(tidy_cox_exp(cox_cont_amyloid), file.path(out_dir, "MCI_Cox_cont_amyloid.csv"))

cox_group_basic_ph <- cox.zph(cox_group_basic)
write_tidy_csv(cox_group_basic_ph, file.path(out_dir, "MCI_Cox_group_basic_PH.csv"))

cox_group_amyloid_ph <- cox.zph(cox_group_amyloid)
write_tidy_csv(cox_group_amyloid_ph, file.path(out_dir, "MCI_Cox_group_amyloid_PH.csv"))

cox_cont_amyloid_ph <- cox.zph(cox_cont_amyloid)
write_tidy_csv(cox_cont_amyloid_ph, file.path(out_dir, "MCI_Cox_cont_amyloid_PH.csv"))

surv_moca <- surv_complete %>% filter(!is.na(baseline_MOCA_z))
cox_group_moca <- fit_primary_cox(
  surv_moca,
  "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z + baseline_MOCA_z"
)
write_csv(tidy_cox_exp(cox_group_moca), file.path(out_dir, "MCI_Cox_group_moca.csv"))

surv_all3 <- surv_complete %>% filter(has_all_3)
cox_group_all3 <- fit_primary_cox(
  surv_all3,
  "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z"
)
write_csv(tidy_cox_exp(cox_group_all3), file.path(out_dir, "MCI_Cox_group_all3.csv"))

cox_group_tertile <- fit_primary_cox(
  surv_complete,
  "PE_group_tertile + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z"
)
write_csv(tidy_cox_exp(cox_group_tertile), file.path(out_dir, "MCI_Cox_group_tertile.csv"))

landmark1 <- make_landmark_dataset(surv_complete, 1)
landmark2 <- make_landmark_dataset(surv_complete, 2)

cox_group_landmark1 <- fit_primary_cox(
  landmark1,
  "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z"
)
write_csv(tidy_cox_exp(cox_group_landmark1), file.path(out_dir, "MCI_Cox_group_landmark1y.csv"))

cox_group_landmark2 <- fit_primary_cox(
  landmark2,
  "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z"
)
write_csv(tidy_cox_exp(cox_group_landmark2), file.path(out_dir, "MCI_Cox_group_landmark2y.csv"))

age_cut <- median(surv_complete$baseline_age_years, na.rm = TRUE)
surv_age <- surv_complete %>%
  mutate(age_stratum = factor(ifelse(baseline_age_years >= age_cut, "Older", "Younger")))
cox_group_age_strat <- coxph(
  Surv(time_to_event_years, converted) ~
    PE_group_lmer + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z + strata(age_stratum),
  data = surv_age
)
write_csv(tidy_cox_exp(cox_group_age_strat), file.path(out_dir, "MCI_Cox_group_age_stratified.csv"))
write_csv(
  surv_age %>%
    group_by(age_stratum) %>%
    summarise(n = n(), n_events = sum(converted), age_min = min(baseline_age_years), age_max = max(baseline_age_years), .groups = "drop"),
  file.path(out_dir, "MCI_Cox_group_age_stratified_summary.csv")
)

surv_lastcontact <- surv_rbans %>%
  mutate(PE_group_lmer = factor(PE_group_lmer, levels = c("High_PE", "Moderate_PE", "Low_PE"))) %>%
  filter(
    !is.na(Age_z),
    !is.na(Education_z),
    !is.na(Sex_F),
    !is.na(APOE4_carrier),
    !is.na(baseline_RBANS_z),
    !is.na(PE_group_lmer)
  )
cox_group_lastcontact <- fit_primary_cox(
  surv_lastcontact,
  "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z"
)
write_csv(tidy_cox_exp(cox_group_lastcontact), file.path(out_dir, "MCI_Cox_group_lastcontact.csv"))

surv_csf <- surv_complete %>%
  filter(!is.na(CSF_ptau181_z), !is.na(CSF_Abeta42_z), !is.na(CSF_first_visit_months), CSF_first_visit_months <= 24)
cox_group_csf <- fit_primary_cox(
  surv_csf,
  "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z + CSF_ptau181_z + CSF_Abeta42_z"
)
write_csv(tidy_cox_exp(cox_group_csf), file.path(out_dir, "MCI_Cox_group_csf_biomarkers.csv"))

surv_intrepad <- surv_complete %>% filter(!is.na(INTREPAD_arm))
cox_group_intrepad <- fit_primary_cox(
  surv_intrepad,
  "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z + INTREPAD_arm"
)
write_csv(tidy_cox_exp(cox_group_intrepad), file.path(out_dir, "MCI_Cox_group_intrepad_adjusted.csv"))
write_csv(
  tibble(
    note = "PE x INTREPAD interaction not promoted because arm-known survival subset is small.",
    n = nrow(surv_intrepad),
    n_events = sum(surv_intrepad$converted)
  ),
  file.path(out_dir, "INTREPAD_interaction_note.csv")
)

surv_pet_pos <- surv_complete %>% filter(amyloid_pos_pet %in% TRUE)
if (nrow(surv_pet_pos) >= 25 && sum(surv_pet_pos$converted) >= 10) {
  cox_group_pet_pos <- fit_primary_cox(
    surv_pet_pos,
    "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z"
  )
  write_csv(tidy_cox_exp(cox_group_pet_pos), file.path(out_dir, "MCI_Cox_group_pet_positive.csv"))
} else {
  write_csv(
    tibble(
      note = "PET-positive subgroup survival model not fit because sample size/events were too limited.",
      n = nrow(surv_pet_pos),
      n_events = sum(surv_pet_pos$converted)
    ),
    file.path(out_dir, "MCI_Cox_group_pet_positive.csv")
  )
}

base_model <- fit_primary_cox(
  surv_complete,
  "Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z"
)

incremental_discrimination <- tibble(
  dataset = "PREVENT-AD",
  n = nrow(surv_complete),
  n_events = sum(surv_complete$converted),
  C_base = round(harrell_c(base_model), 4),
  C_pe = round(harrell_c(cox_group_basic), 4),
  delta_C = round(harrell_c(cox_group_basic) - harrell_c(base_model), 4),
  LRT_chisq = as.numeric(2 * (logLik(cox_group_basic) - logLik(base_model))),
  LRT_df = attr(logLik(cox_group_basic), "df") - attr(logLik(base_model), "df"),
  LRT_p = as.numeric(pchisq(LRT_chisq, df = LRT_df, lower.tail = FALSE))
)
write_csv(incremental_discrimination, file.path(out_dir, "incremental_discrimination_MCI.csv"))

landmark_summary <- tibble(
  analysis = c("primary", "landmark_1y", "landmark_2y"),
  n = c(nrow(surv_complete), nrow(landmark1), nrow(landmark2)),
  n_events = c(sum(surv_complete$converted), sum(landmark1$converted), sum(landmark2$converted))
)
write_csv(landmark_summary, file.path(out_dir, "MCI_landmark_summary.csv"))

sensitivity_summary <- bind_rows(
  tidy_cox_exp(cox_group_basic) %>% filter(term == "PE_group_lmerLow_PE") %>% mutate(analysis = "primary"),
  tidy_cox_exp(cox_group_moca) %>% filter(term == "PE_group_lmerLow_PE") %>% mutate(analysis = "baseline_moca_adjusted"),
  tidy_cox_exp(cox_group_all3) %>% filter(term == "PE_group_lmerLow_PE") %>% mutate(analysis = "all_three_visits"),
  tidy_cox_exp(cox_group_landmark1) %>% filter(term == "PE_group_lmerLow_PE") %>% mutate(analysis = "landmark_1y"),
  tidy_cox_exp(cox_group_landmark2) %>% filter(term == "PE_group_lmerLow_PE") %>% mutate(analysis = "landmark_2y"),
  tidy_cox_exp(cox_group_age_strat) %>% filter(term == "PE_group_lmerLow_PE") %>% mutate(analysis = "age_stratified"),
  tidy_cox_exp(cox_group_lastcontact) %>% filter(term == "PE_group_lmerLow_PE") %>% mutate(analysis = "rbans_only_censoring"),
  tidy_cox_exp(cox_group_tertile) %>% filter(term == "PE_group_tertileLow_PE") %>% mutate(analysis = "tertile_low_vs_high")
) %>%
  select(analysis, term, estimate, conf.low, conf.high, p.value)
write_csv(sensitivity_summary, file.path(out_dir, "MCI_sensitivity_summary.csv"))

# ----------------------------
# Sensitivity: include FU24-prevalent MCI (i.e. mci_visit_months <= 24)
# in the risk-set as event = 0. The strict landmark default excludes them
# because keeping them is not a clean landmark; this block is provided
# only so the directional comparison can be reported as a supplement note.
# ----------------------------
surv_primary_inclprevalent <- build_survival_dataset(
  surv_seed,
  "last_contact_proxy_months",
  "last_contact_proxy",
  strict_landmark = FALSE
)

surv_complete_inclprevalent <- surv_primary_inclprevalent %>%
  mutate(
    PE_group_lmer = factor(PE_group_lmer, levels = c("High_PE", "Moderate_PE", "Low_PE"))
  ) %>%
  filter(
    !is.na(Age_z), !is.na(Education_z), !is.na(Sex_F),
    !is.na(APOE4_carrier), !is.na(baseline_RBANS_z), !is.na(PE_group_lmer)
  )

cox_group_basic_inclprevalent <- fit_primary_cox(
  surv_complete_inclprevalent,
  "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z"
)
cox_cont_basic_inclprevalent <- fit_primary_cox(
  surv_complete_inclprevalent,
  "pe_resid_z_decr_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z"
)
write_csv(
  tidy_cox_exp(cox_group_basic_inclprevalent),
  file.path(out_dir, "MCI_Cox_group_basic_inclprevalent.csv")
)
write_csv(
  tidy_cox_exp(cox_cont_basic_inclprevalent),
  file.path(out_dir, "MCI_Cox_cont_basic_inclprevalent.csv")
)

surv_amyloid_inclprevalent <- surv_complete_inclprevalent %>% filter(!is.na(SUVR_z))
cox_group_amyloid_inclprevalent <- fit_primary_cox(
  surv_amyloid_inclprevalent,
  "PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + baseline_RBANS_z"
)
write_csv(
  tidy_cox_exp(cox_group_amyloid_inclprevalent),
  file.path(out_dir, "MCI_Cox_group_amyloid_inclprevalent.csv")
)

write_csv(
  tibble(
    analysis = c("primary_strict_landmark", "include_FU24_prevalent_as_nonevents"),
    n = c(nrow(surv_complete), nrow(surv_complete_inclprevalent)),
    n_events = c(sum(surv_complete$converted), sum(surv_complete_inclprevalent$converted)),
    pet_n = c(
      sum(!is.na(surv_complete$SUVR_z)),
      sum(!is.na(surv_complete_inclprevalent$SUVR_z))
    ),
    pet_n_events = c(
      sum(!is.na(surv_complete$SUVR_z) & surv_complete$converted),
      sum(!is.na(surv_complete_inclprevalent$SUVR_z) & surv_complete_inclprevalent$converted)
    )
  ),
  file.path(out_dir, "MCI_riskset_comparison.csv")
)

# ----------------------------
# Figures and supporting summaries
# ----------------------------
p_pe <- ggplot(dat_lmer, aes(x = pe_resid_z_lmer, fill = PE_group_lmer)) +
  geom_histogram(binwidth = 0.25, color = "white", alpha = 0.8) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
  scale_fill_manual(values = c("High_PE" = "#1b9e77", "Moderate_PE" = "#7570b3", "Low_PE" = "#d95f02")) +
  labs(x = "Residualized PE (z)", y = "Count", fill = "PE group", title = "PREVENT-AD PE distribution") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "Figure1_PE_distribution.pdf"), p_pe, width = 7, height = 5)
ggsave(file.path(out_dir, "Figure1_PE_distribution.png"), p_pe, width = 7, height = 5, dpi = 300)

traj_summary <- rbans %>%
  filter(visit_months %in% c(0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120), !is.na(total_scale_index_score)) %>%
  inner_join(dat_lmer %>% select(CONP_ID, PE_group_lmer), by = "CONP_ID") %>%
  group_by(PE_group_lmer, visit_months) %>%
  summarise(
    mean_rbans = mean(as.numeric(total_scale_index_score), na.rm = TRUE),
    se_rbans = sd(as.numeric(total_scale_index_score), na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_traj <- ggplot(traj_summary, aes(x = visit_months / 12, y = mean_rbans, color = PE_group_lmer, fill = PE_group_lmer)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean_rbans - 1.96 * se_rbans, ymax = mean_rbans + 1.96 * se_rbans), alpha = 0.15, color = NA) +
  geom_vline(xintercept = 2, linetype = "dashed") +
  scale_color_manual(values = c("High_PE" = "#1b9e77", "Moderate_PE" = "#7570b3", "Low_PE" = "#d95f02")) +
  scale_fill_manual(values = c("High_PE" = "#1b9e77", "Moderate_PE" = "#7570b3", "Low_PE" = "#d95f02")) +
  labs(x = "Study year", y = "RBANS total score", color = "PE group", fill = "PE group", title = "RBANS trajectories by PE group") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "Figure2_RBANS_trajectories.pdf"), p_traj, width = 7, height = 5)
ggsave(file.path(out_dir, "Figure2_RBANS_trajectories.png"), p_traj, width = 7, height = 5, dpi = 300)

km_fit <- survfit(Surv(time_to_event_years, converted) ~ PE_group_lmer, data = surv_complete)
km_df <- broom::tidy(km_fit) %>%
  mutate(strata = str_replace(strata, "PE_group_lmer=", ""))

p_km <- ggplot(km_df, aes(x = time, y = estimate, color = strata)) +
  geom_step(linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("High_PE" = "#1b9e77", "Moderate_PE" = "#7570b3", "Low_PE" = "#d95f02")) +
  labs(x = "Years since FU24", y = "MCI-free survival", color = "PE group", title = "MCI conversion by PE group") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "Figure3_KM_MCI_conversion.pdf"), p_km, width = 7, height = 5)
ggsave(file.path(out_dir, "Figure3_KM_MCI_conversion.png"), p_km, width = 7, height = 5, dpi = 300)

a4_vs_preventad <- tibble(
  dimension = c(
    "Discovery cohort", "External validation cohort", "Primary clinical endpoint",
    "Median event time (years)", "PE window", "Cognitive platform",
    "Amyloid coverage", "p-tau217 availability", "Key validation message"
  ),
  A4 = c(
    "A4 placebo arm", "N/A", "Confirmed CDR progression",
    "0.29", "21 weeks", "Cogstate", "All PET-positive", "Available",
    "Discovery of PE-prognosis association"
  ),
  PREVENT_AD = c(
    "N/A", "PREVENT-AD", "Incident MCI after FU24",
    fmt_num(median(surv_primary$time_to_event_years[surv_primary$converted], na.rm = TRUE), 2),
    "24 months", "RBANS", "Partial PET/CSF", "Not available",
    "Partial external validation of clinical progression association"
  )
)
write_csv(a4_vs_preventad, file.path(out_dir, "A4_vs_PREVENTAD_comparison.csv"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
cat("PREVENT-AD replication analyses complete. Outputs written to:", out_dir, "\n")
