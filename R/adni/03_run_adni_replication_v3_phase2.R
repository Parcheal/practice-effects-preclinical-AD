

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
# ADNI replication v3 Phase 2 supplementary analyses
# - Biomarker interactions (PE x CSF_PTAU, PE x CSF_ABETA42, PE x plasma, PE x amyloid)
# - AT stratification (Elecsys cutoffs)
# - Included vs excluded comparison (CSF availability)
# - Biomarker timing summary
# - IPW sensitivity for CSF selection
# - Stricter timing sensitivity (+/- 365 days)
# - Outputs appended to adni_replication_v3_supplement/
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(readr)
  library(stringr)
  library(survival)
  library(tibble)
})

out_dir <- file.path(pkg_root, "verify", "adni", "supplement")
data_dir <- file.path(data_root, "adni")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── utility functions (same as Phase 1) ────────────────────────────────────

pdate <- function(x) as.Date(x)

vis_month <- function(x) {
  x <- tolower(as.character(x))
  out <- suppressWarnings(as.numeric(str_extract(x, "\\d+")))
  out[!is.na(x) & x %in% c("bl", "sc", "f")] <- 0
  out
}

sscale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

first_nm <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[[1]]
}

pe_group <- function(z) {
  factor(
    dplyr::case_when(
      is.na(z) ~ NA_character_,
      z > 0.2 ~ "High_PE",
      z < -0.2 ~ "Low_PE",
      TRUE ~ "Moderate_PE"
    ),
    levels = c("High_PE", "Moderate_PE", "Low_PE")
  )
}

rhs <- function(...) paste(c(...), collapse = " + ")

keep_complete <- function(df, vars) {
  df %>% filter(if_all(all_of(vars), ~ !is.na(.)))
}

fit_cox <- function(df, rhs_text) {
  tryCatch(
    coxph(as.formula(paste0("Surv(time_to_event_years, event) ~ ", rhs_text)), data = df),
    error = function(e) NULL
  )
}

extract_term <- function(model, term, meta) {
  base <- as_tibble(meta)
  if (is.null(model)) {
    return(base %>% mutate(term = term, estimate = NA_real_, conf.low = NA_real_,
                           conf.high = NA_real_, p.value = NA_real_))
  }
  s <- summary(model)
  if (!(term %in% rownames(s$coefficients))) {
    return(base %>% mutate(term = term, estimate = NA_real_, conf.low = NA_real_,
                           conf.high = NA_real_, p.value = NA_real_))
  }
  ci <- tryCatch(suppressMessages(confint(model)), error = function(e) NULL)
  est <- exp(unname(s$coefficients[term, "coef"]))
  p   <- unname(s$coefficients[term, "Pr(>|z|)"])
  lo  <- if (!is.null(ci)) exp(ci[term, 1]) else NA_real_
  hi  <- if (!is.null(ci)) exp(ci[term, 2]) else NA_real_
  base %>% mutate(term = term, estimate = est, conf.low = lo, conf.high = hi, p.value = p)
}

extract_beta <- function(model, term, meta) {
  base <- as_tibble(meta)
  if (is.null(model)) {
    return(base %>% mutate(term = term, beta = NA_real_, se = NA_real_,
                           conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_))
  }
  s <- summary(model)
  if (!(term %in% rownames(s$coefficients))) {
    return(base %>% mutate(term = term, beta = NA_real_, se = NA_real_,
                           conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_))
  }
  ci <- tryCatch(suppressMessages(confint(model)), error = function(e) NULL)
  b  <- unname(s$coefficients[term, "coef"])
  se <- unname(s$coefficients[term, "se(coef)"])
  p  <- unname(s$coefficients[term, "Pr(>|z|)"])
  lo <- if (!is.null(ci)) ci[term, 1] else NA_real_
  hi <- if (!is.null(ci)) ci[term, 2] else NA_real_
  base %>% mutate(term = term, beta = b, se = se, conf.low = lo, conf.high = hi, p.value = p)
}

harrell_c <- function(model) {
  out <- tryCatch(summary(model)$concordance[1], error = function(e) NA_real_)
  as.numeric(out)
}

nearest_bio_full <- function(df, date_col, analysis_df, req) {
  df %>%
    filter(if_all(all_of(req), ~ !is.na(.))) %>%
    inner_join(analysis_df %>% select(RID, anchor_date, analysis_end_date), by = "RID") %>%
    filter(.data[[date_col]] <= analysis_end_date) %>%
    mutate(diff_days = as.numeric(.data[[date_col]] - anchor_date),
           abs_diff_days = abs(diff_days)) %>%
    arrange(RID, abs_diff_days, diff_days) %>%
    group_by(RID) %>%
    slice(1) %>%
    ungroup()
}

fit_blup_slope_model <- function(use_df) {
  lmer(
    score ~ time_years + (1 + time_years | RID),
    data = use_df,
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
}

predict_blup_slopes <- function(fit, target_df) {
  beta <- lme4::fixef(fit)
  random_cov <- as.matrix(lme4::VarCorr(fit)$RID)
  sigma2 <- sigma(fit)^2

  bind_rows(lapply(split(target_df, target_df$RID), function(dat) {
    dat <- dat %>% arrange(time_years)
    x_mat <- cbind(`(Intercept)` = 1, time_years = dat$time_years)
    y_vec <- dat$score
    resid_vec <- as.numeric(y_vec - x_mat %*% beta[c("(Intercept)", "time_years")])
    v_mat <- x_mat %*% random_cov %*% t(x_mat) + diag(sigma2, nrow(x_mat))
    b_hat <- random_cov %*% t(x_mat) %*% solve(v_mat, resid_vec)
    tibble(RID = dat$RID[[1]], raw_slope_blup = unname(beta[["time_years"]] + b_hat[2, 1]))
  }))
}

select_pe_window <- function(long_df, window_mode) {
  long_df <- long_df %>% mutate(visit_months = vis_month(VISCODE2))

  use_df <- if (window_mode == "exact") {
    ids <- long_df %>%
      filter(VISCODE2 %in% c("bl", "m06", "m12"), !is.na(score)) %>%
      distinct(RID, VISCODE2) %>%
      count(RID) %>%
      filter(n == 3) %>%
      select(RID)

    long_df %>%
      semi_join(ids, by = "RID") %>%
      filter(VISCODE2 %in% c("bl", "m06", "m12"))
  } else {
    elig <- long_df %>% filter(VISCODE2 %in% c("bl", "m06", "m12", "m24"), !is.na(score))
    ids <- elig %>%
      group_by(RID) %>%
      summarise(
        has_bl = any(VISCODE2 == "bl"),
        has_m06 = any(VISCODE2 == "m06"),
        has_m12 = any(VISCODE2 == "m12"),
        has_m24 = any(VISCODE2 == "m24"),
        chosen_end = case_when(
          has_bl & has_m06 & has_m12 ~ "m12",
          has_bl & has_m06 & has_m24 ~ "m24",
          TRUE ~ NA_character_
        ),
        .groups = "drop"
      ) %>%
      filter(!is.na(chosen_end))

    elig %>%
      inner_join(ids, by = "RID") %>%
      filter(VISCODE2 %in% c("bl", "m06") | VISCODE2 == chosen_end)
  }

  use_df %>%
    filter(!is.na(score), !is.na(visit_months)) %>%
    mutate(time_years = visit_months / 12)
}

compute_pe <- function(long_df, exposure_id, window_mode) {
  use_df <- select_pe_window(long_df, window_mode)
  blup_fit <- fit_blup_slope_model(use_df)
  blup_slopes <- predict_blup_slopes(blup_fit, use_df)

  subj <- use_df %>%
    arrange(RID, visit_months) %>%
    group_by(RID) %>%
    summarise(
      anchor_date = dplyr::last(source_date),
      baseline_score = score[which.min(visit_months)],
      raw_slope_ols = coef(lm(score ~ I(visit_months / 12)))[["I(visit_months/12)"]],
      .groups = "drop"
    ) %>%
    left_join(blup_slopes, by = "RID")

  resid_fit <- lm(raw_slope_blup ~ baseline_score, data = subj)
  ols_resid_fit <- lm(raw_slope_ols ~ baseline_score, data = subj)

  subj %>%
    mutate(
      exposure_id = exposure_id,
      window_mode = window_mode,
      pe_estimator = "BLUP",
      pe_resid_z = sscale(resid(resid_fit)),
      pe_resid_z_decr = -pe_resid_z,
      PE_group = pe_group(pe_resid_z),
      ols_pe_resid_z = sscale(resid(ols_resid_fit)),
      ols_pe_resid_z_decr = -ols_pe_resid_z,
      PE_group_ols = pe_group(ols_pe_resid_z),
      baseline_score_z = sscale(baseline_score)
    )
}

build_dx23 <- function(subject_df, dxsum_df, exclude_anchor_prevalent = TRUE) {
  if (isTRUE(exclude_anchor_prevalent)) {
    prevalent_ids <- dxsum_df %>%
      inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>%
      filter(
        !is.na(EXAMDATE), !is.na(DIAGNOSIS),
        EXAMDATE <= anchor_date,
        DIAGNOSIS %in% c(2, 3)
      ) %>%
      distinct(RID)
    subject_df <- subject_df %>% anti_join(prevalent_ids, by = "RID")
  }

  post <- dxsum_df %>% inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>%
    filter(!is.na(EXAMDATE), EXAMDATE > anchor_date) %>% arrange(RID, EXAMDATE)
  ev <- post %>% group_by(RID) %>%
    summarise(first_event_date = if (any(DIAGNOSIS %in% c(2, 3), na.rm = TRUE))
      EXAMDATE[which(DIAGNOSIS %in% c(2, 3))[1]] else as.Date(NA),
      last_followup_date = max(EXAMDATE, na.rm = TRUE), .groups = "drop") %>%
    mutate(event = !is.na(first_event_date),
           analysis_end_date = if_else(event, first_event_date, last_followup_date))
  subject_df %>% left_join(ev, by = "RID") %>% filter(!is.na(last_followup_date)) %>%
    mutate(outcome_id = "DXSUM_23",
           time_to_event_years = pmax(as.numeric(analysis_end_date - anchor_date) / 365.25, 0.01))
}

build_mci_only <- function(subject_df, dxsum_df) {
  post <- dxsum_df %>% inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>%
    filter(!is.na(EXAMDATE), EXAMDATE > anchor_date) %>% arrange(RID, EXAMDATE)
  ev <- post %>% group_by(RID) %>%
    summarise(first_mci_date = if (any(DIAGNOSIS == 2, na.rm = TRUE))
      EXAMDATE[which(DIAGNOSIS == 2)[1]] else as.Date(NA),
      first_dem_date = if (any(DIAGNOSIS == 3, na.rm = TRUE))
        EXAMDATE[which(DIAGNOSIS == 3)[1]] else as.Date(NA),
      last_followup_date = max(EXAMDATE, na.rm = TRUE), .groups = "drop") %>%
    mutate(event = !is.na(first_mci_date) & (is.na(first_dem_date) | first_mci_date <= first_dem_date),
           analysis_end_date = case_when(event ~ first_mci_date,
                                         !is.na(first_dem_date) ~ first_dem_date,
                                         TRUE ~ last_followup_date))
  subject_df %>% left_join(ev, by = "RID") %>% filter(!is.na(last_followup_date)) %>%
    mutate(outcome_id = "DXSUM_2_MCI_only",
           time_to_event_years = pmax(as.numeric(analysis_end_date - anchor_date) / 365.25, 0.01))
}

# ── load data ──────────────────────────────────────────────────────────────

cat("Loading data...\n")

dxsum <- read_csv(file.path(data_dir, "DXSUM_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(EXAMDATE = pdate(EXAMDATE))
cdr <- read_csv(file.path(data_dir, "CDR_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(VISDATE = pdate(VISDATE))
ptdemog <- read_csv(file.path(data_dir, "PTDEMOG_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(VISDATE = pdate(VISDATE))
apoe_raw <- read_csv(file.path(data_dir, "APOERES_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(APTESTDT = pdate(APTESTDT))
phc_cogn <- read_csv(file.path(data_dir, "ADSP_PHC_COGN_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(EXAMDATE = pdate(EXAMDATE))
plasma_fuji <- read_csv(file.path(data_dir, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_10Apr2026.csv"),
                        show_col_types = FALSE) %>%
  mutate(EXAMDATE = pdate(EXAMDATE))
amyloid_pet <- read_csv(file.path(data_dir, "ADSP_PHC_PET_Amyloid_Simple_10Apr2026.csv"),
                        show_col_types = FALSE) %>%
  mutate(PHC_SCANDATE = pdate(PHC_SCANDATE))
csf_elecsys <- read_csv(file.path(data_dir, "UPENNBIOMK_ROCHE_ELECSYS_10Apr2026.csv"),
                        show_col_types = FALSE) %>%
  mutate(EXAMDATE = pdate(EXAMDATE))

baseline_cu <- dxsum %>% filter(VISCODE2 == "bl", DIAGNOSIS == 1) %>% distinct(RID)

ptdemog_subject <- ptdemog %>% arrange(RID, VISDATE) %>% group_by(RID) %>%
  summarise(PTGENDER = first_nm(as.numeric(PTGENDER[PTGENDER >= 0])),
            education_years = first_nm(as.numeric(PTEDUCAT[PTEDUCAT >= 0])),
            .groups = "drop")
age_subject <- phc_cogn %>% filter(VISCODE2 == "bl") %>% group_by(RID) %>%
  summarise(baseline_age_years = first_nm(as.numeric(PHC_Age_Cognition)), .groups = "drop")
apoe_subject <- apoe_raw %>% arrange(RID, APTESTDT) %>% group_by(RID) %>%
  summarise(APOE4 = ifelse(str_detect(first_nm(GENOTYPE), "4"), 1, 0), .groups = "drop")

subject_cov <- baseline_cu %>%
  left_join(age_subject, by = "RID") %>%
  left_join(ptdemog_subject, by = "RID") %>%
  left_join(apoe_subject, by = "RID") %>%
  mutate(Age_z = sscale(baseline_age_years),
         Sex_F = factor(ifelse(PTGENDER == 2, "Female", "Male"), levels = c("Male", "Female")),
         Education_z = sscale(education_years))

phc_mem_long <- phc_cogn %>% semi_join(baseline_cu, by = "RID") %>%
  transmute(RID, VISCODE2, source_date = EXAMDATE,
            score = if_else(as.numeric(PHC_MEM_PreciseFilter) == 1, as.numeric(PHC_MEM), NA_real_))

plasma_prepped <- plasma_fuji %>%
  transmute(RID, biomarker_date = EXAMDATE,
            pT217_F = as.numeric(pT217_F), pT217_F_z = sscale(as.numeric(pT217_F)))
amy_prepped <- amyloid_pet %>%
  transmute(RID, biomarker_date = PHC_SCANDATE,
            amyloid_centiloids = as.numeric(PHC_CENTILOIDS),
            amyloid_centiloids_z = sscale(as.numeric(PHC_CENTILOIDS)))
csf_prepped <- csf_elecsys %>%
  transmute(RID, biomarker_date = EXAMDATE,
            CSF_PTAU = as.numeric(PTAU), CSF_PTAU_z = sscale(as.numeric(PTAU)),
            CSF_ABETA42 = as.numeric(ABETA42), CSF_ABETA42_z = sscale(as.numeric(ABETA42)))

base_rhs <- rhs("Age_z", "Sex_F", "Education_z", "APOE4", "baseline_score_z")
csf_rhs  <- rhs(base_rhs, "CSF_PTAU_z", "CSF_ABETA42_z")
mainbio_rhs <- rhs(base_rhs, "pT217_F_z", "amyloid_centiloids_z")

# ── build primary cohort (PHC_MEM flex) ────────────────────────────────────

cat("Building primary cohort...\n")

phc_mem_flex <- compute_pe(phc_mem_long, "PHC_MEM", "flex") %>%
  left_join(subject_cov, by = "RID")

dx23 <- build_dx23(phc_mem_flex, dxsum)
mci  <- build_mci_only(phc_mem_flex, dxsum)

# match biomarkers with timing info preserved
csf_match_full <- nearest_bio_full(
  csf_prepped, "biomarker_date", dx23, c("CSF_PTAU_z", "CSF_ABETA42_z")
) %>% transmute(RID, CSF_PTAU, CSF_PTAU_z, CSF_ABETA42, CSF_ABETA42_z,
                csf_diff_days = diff_days)

plasma_match_full <- nearest_bio_full(
  plasma_prepped, "biomarker_date", dx23, "pT217_F_z"
) %>% transmute(RID, pT217_F, pT217_F_z, plasma_diff_days = diff_days)

amy_match_full <- nearest_bio_full(
  amy_prepped, "biomarker_date", dx23, "amyloid_centiloids_z"
) %>% transmute(RID, amyloid_centiloids, amyloid_centiloids_z, amy_diff_days = diff_days)

dx23_full <- dx23 %>%
  left_join(csf_match_full, by = "RID") %>%
  left_join(plasma_match_full, by = "RID") %>%
  left_join(amy_match_full, by = "RID")

csf_vars <- c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z",
              "APOE4", "baseline_score_z", "pe_resid_z_decr", "PE_group",
              "CSF_PTAU_z", "CSF_ABETA42_z")
bio_vars <- c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z",
              "APOE4", "baseline_score_z", "pe_resid_z_decr", "PE_group",
              "pT217_F_z", "amyloid_centiloids_z")

dx23_csf <- keep_complete(dx23_full, csf_vars)
dx23_bio <- keep_complete(dx23_full, bio_vars)

# MCI-only with CSF
mci_csf_match <- nearest_bio_full(
  csf_prepped, "biomarker_date", mci, c("CSF_PTAU_z", "CSF_ABETA42_z")
) %>% transmute(RID, CSF_PTAU, CSF_PTAU_z, CSF_ABETA42, CSF_ABETA42_z)

mci_csf <- mci %>%
  left_join(mci_csf_match, by = "RID") %>%
  keep_complete(csf_vars)

cat("  DXSUM_23 CSF:  n=", nrow(dx23_csf), " events=", sum(dx23_csf$event), "\n")
cat("  DXSUM_23 bio:  n=", nrow(dx23_bio), " events=", sum(dx23_bio$event), "\n")
cat("  MCI-only CSF:  n=", nrow(mci_csf),  " events=", sum(mci_csf$event), "\n")

# ============================================================================
# PART 1: BIOMARKER INTERACTIONS
# ============================================================================

cat("\n== Part 1: Biomarker interactions ==\n")

interaction_rows <- list()

# 1a. PE x CSF_PTAU (DXSUM_23)
int_rhs_ptau <- rhs("pe_resid_z_decr", csf_rhs, "pe_resid_z_decr:CSF_PTAU_z")
m_int_ptau <- fit_cox(dx23_csf, int_rhs_ptau)
interaction_rows[[length(interaction_rows) + 1]] <- extract_beta(
  m_int_ptau, "pe_resid_z_decr:CSF_PTAU_z",
  list(outcome_id = "DXSUM_23", interaction = "PE_x_CSF_PTAU",
       biomarker_layer = "csf_elecsys", n = nrow(dx23_csf), n_events = sum(dx23_csf$event)))

# 1b. PE x CSF_ABETA42 (DXSUM_23)
int_rhs_ab42 <- rhs("pe_resid_z_decr", csf_rhs, "pe_resid_z_decr:CSF_ABETA42_z")
m_int_ab42 <- fit_cox(dx23_csf, int_rhs_ab42)
interaction_rows[[length(interaction_rows) + 1]] <- extract_beta(
  m_int_ab42, "pe_resid_z_decr:CSF_ABETA42_z",
  list(outcome_id = "DXSUM_23", interaction = "PE_x_CSF_ABETA42",
       biomarker_layer = "csf_elecsys", n = nrow(dx23_csf), n_events = sum(dx23_csf$event)))

# 1c. PE x pT217_F (DXSUM_23, plasma+amyloid layer)
int_rhs_pt217 <- rhs("pe_resid_z_decr", mainbio_rhs, "pe_resid_z_decr:pT217_F_z")
m_int_pt217 <- fit_cox(dx23_bio, int_rhs_pt217)
interaction_rows[[length(interaction_rows) + 1]] <- extract_beta(
  m_int_pt217, "pe_resid_z_decr:pT217_F_z",
  list(outcome_id = "DXSUM_23", interaction = "PE_x_plasma_pT217",
       biomarker_layer = "main_biomarker", n = nrow(dx23_bio), n_events = sum(dx23_bio$event)))

# 1d. PE x amyloid centiloids (DXSUM_23, plasma+amyloid layer)
int_rhs_amy <- rhs("pe_resid_z_decr", mainbio_rhs, "pe_resid_z_decr:amyloid_centiloids_z")
m_int_amy <- fit_cox(dx23_bio, int_rhs_amy)
interaction_rows[[length(interaction_rows) + 1]] <- extract_beta(
  m_int_amy, "pe_resid_z_decr:amyloid_centiloids_z",
  list(outcome_id = "DXSUM_23", interaction = "PE_x_amyloid_centiloids",
       biomarker_layer = "main_biomarker", n = nrow(dx23_bio), n_events = sum(dx23_bio$event)))

# 1e. PE x CSF_PTAU (MCI-only)
m_int_ptau_mci <- fit_cox(mci_csf, int_rhs_ptau)
interaction_rows[[length(interaction_rows) + 1]] <- extract_beta(
  m_int_ptau_mci, "pe_resid_z_decr:CSF_PTAU_z",
  list(outcome_id = "DXSUM_2_MCI_only", interaction = "PE_x_CSF_PTAU",
       biomarker_layer = "csf_elecsys", n = nrow(mci_csf), n_events = sum(mci_csf$event)))

# 1f. PE x CSF_ABETA42 (MCI-only)
m_int_ab42_mci <- fit_cox(mci_csf, int_rhs_ab42)
interaction_rows[[length(interaction_rows) + 1]] <- extract_beta(
  m_int_ab42_mci, "pe_resid_z_decr:CSF_ABETA42_z",
  list(outcome_id = "DXSUM_2_MCI_only", interaction = "PE_x_CSF_ABETA42",
       biomarker_layer = "csf_elecsys", n = nrow(mci_csf), n_events = sum(mci_csf$event)))

interaction_results <- bind_rows(interaction_rows)
write_csv(interaction_results, file.path(out_dir, "phase2_interaction_results.csv"))
cat("  Interactions saved.\n")

# ============================================================================
# PART 2: AT STRATIFICATION
# ============================================================================

cat("\n== Part 2: AT stratification ==\n")

# Elecsys cutoffs (published): ABETA42 < 880 = A+, PTAU > 24 = T+
dx23_csf <- dx23_csf %>%
  mutate(
    A_positive = CSF_ABETA42 < 880,
    T_positive = CSF_PTAU > 24,
    AT_group = factor(
      paste0(ifelse(A_positive, "A+", "A-"), ifelse(T_positive, "T+", "T-")),
      levels = c("A-T-", "A-T+", "A+T-", "A+T+")
    ),
    A_status = factor(ifelse(A_positive, "A+", "A-"), levels = c("A-", "A+")),
    T_status = factor(ifelse(T_positive, "T+", "T-"), levels = c("T-", "T+"))
  )

# AT group descriptive table
at_desc <- dx23_csf %>%
  group_by(AT_group) %>%
  summarise(
    n = n(),
    n_events = sum(event),
    event_rate_pct = 100 * mean(event),
    median_pe_z = median(pe_resid_z, na.rm = TRUE),
    mean_pe_z = mean(pe_resid_z, na.rm = TRUE),
    n_low_pe = sum(PE_group == "Low_PE"),
    n_high_pe = sum(PE_group == "High_PE"),
    .groups = "drop"
  )
write_csv(at_desc, file.path(out_dir, "phase2_AT_descriptive.csv"))

# Stratified Cox: PE HR within each AT group
at_strat_rows <- list()

for (grp in levels(dx23_csf$AT_group)) {
  sub <- dx23_csf %>% filter(AT_group == grp)
  if (nrow(sub) >= 20 && sum(sub$event) >= 10) {
    m <- fit_cox(sub, rhs("pe_resid_z_decr", base_rhs))
    at_strat_rows[[length(at_strat_rows) + 1]] <- extract_term(
      m, "pe_resid_z_decr",
      list(AT_group = grp, n = nrow(sub), n_events = sum(sub$event)))
  } else {
    at_strat_rows[[length(at_strat_rows) + 1]] <- tibble(
      AT_group = grp, n = nrow(sub), n_events = sum(sub$event),
      term = "pe_resid_z_decr", estimate = NA_real_,
      conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_,
      note = "insufficient sample")
  }
}

# Collapsed: A+ vs A-
for (grp_col in c("A_status", "T_status")) {
  for (lev in levels(dx23_csf[[grp_col]])) {
    sub <- dx23_csf %>% filter(.data[[grp_col]] == lev)
    if (nrow(sub) >= 30 && sum(sub$event) >= 15) {
      m <- fit_cox(sub, rhs("pe_resid_z_decr", base_rhs))
      at_strat_rows[[length(at_strat_rows) + 1]] <- extract_term(
        m, "pe_resid_z_decr",
        list(AT_group = lev, n = nrow(sub), n_events = sum(sub$event)))
    }
  }
}

at_strat_results <- bind_rows(at_strat_rows)
write_csv(at_strat_results, file.path(out_dir, "phase2_AT_stratified_cox.csv"))
cat("  AT stratification saved.\n")

# ============================================================================
# PART 3: INCLUDED VS EXCLUDED (CSF AVAILABILITY)
# ============================================================================

cat("\n== Part 3: Included vs excluded ==\n")

dx23_all <- dx23_full %>%
  keep_complete(c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z",
                  "APOE4", "baseline_score_z", "pe_resid_z_decr", "PE_group"))

dx23_all <- dx23_all %>%
  mutate(has_csf = !is.na(CSF_PTAU_z) & !is.na(CSF_ABETA42_z))

incl_excl_rows <- list()

# Continuous variables: t-test
for (var in c("baseline_age_years", "education_years", "baseline_score",
              "pe_resid_z_decr", "time_to_event_years")) {
  w <- dx23_all %>% filter(has_csf) %>% pull(!!sym(var))
  wo <- dx23_all %>% filter(!has_csf) %>% pull(!!sym(var))
  w <- w[!is.na(w)]; wo <- wo[!is.na(wo)]
  tt <- tryCatch(t.test(w, wo), error = function(e) NULL)
  incl_excl_rows[[length(incl_excl_rows) + 1]] <- tibble(
    variable = var,
    with_csf_mean = mean(w), with_csf_sd = sd(w), with_csf_n = length(w),
    without_csf_mean = mean(wo), without_csf_sd = sd(wo), without_csf_n = length(wo),
    test = "t-test",
    p_value = if (!is.null(tt)) tt$p.value else NA_real_
  )
}

# Binary variables: chi-square / fisher
for (var_info in list(
  list(var = "Sex_F", val = "Female", label = "female"),
  list(var = "APOE4", val = 1, label = "APOE4_carrier"),
  list(var = "event", val = TRUE, label = "event")
)) {
  tbl <- table(
    dx23_all$has_csf,
    if (is.character(dx23_all[[var_info$var]])) {
      dx23_all[[var_info$var]] == var_info$val
    } else {
      dx23_all[[var_info$var]] == var_info$val
    }
  )
  test_res <- tryCatch(
    if (any(tbl < 5)) fisher.test(tbl) else chisq.test(tbl),
    error = function(e) NULL
  )
  w_sub <- dx23_all %>% filter(has_csf)
  wo_sub <- dx23_all %>% filter(!has_csf)
  w_pct <- mean(w_sub[[var_info$var]] == var_info$val, na.rm = TRUE)
  wo_pct <- mean(wo_sub[[var_info$var]] == var_info$val, na.rm = TRUE)
  incl_excl_rows[[length(incl_excl_rows) + 1]] <- tibble(
    variable = var_info$label,
    with_csf_mean = w_pct, with_csf_sd = NA_real_, with_csf_n = nrow(w_sub),
    without_csf_mean = wo_pct, without_csf_sd = NA_real_, without_csf_n = nrow(wo_sub),
    test = "chi-square/fisher",
    p_value = if (!is.null(test_res)) test_res$p.value else NA_real_
  )
}

incl_excl_results <- bind_rows(incl_excl_rows)
write_csv(incl_excl_results, file.path(out_dir, "phase2_included_vs_excluded.csv"))
cat("  Included vs excluded saved.\n")

# ============================================================================
# PART 4: BIOMARKER TIMING SUMMARY
# ============================================================================

cat("\n== Part 4: Biomarker timing ==\n")

timing_rows <- list()

# CSF timing
csf_timing <- dx23_full %>% filter(!is.na(csf_diff_days))
timing_rows[[1]] <- tibble(
  biomarker = "CSF_Elecsys",
  n_matched = nrow(csf_timing),
  median_diff_days = median(csf_timing$csf_diff_days),
  q25_diff_days = quantile(csf_timing$csf_diff_days, 0.25),
  q75_diff_days = quantile(csf_timing$csf_diff_days, 0.75),
  mean_abs_diff_days = mean(abs(csf_timing$csf_diff_days)),
  pct_after_anchor = 100 * mean(csf_timing$csf_diff_days > 0),
  pct_within_365 = 100 * mean(abs(csf_timing$csf_diff_days) <= 365),
  pct_within_180 = 100 * mean(abs(csf_timing$csf_diff_days) <= 180)
)

# Plasma timing
plasma_timing <- dx23_full %>% filter(!is.na(plasma_diff_days))
timing_rows[[2]] <- tibble(
  biomarker = "Plasma_pT217",
  n_matched = nrow(plasma_timing),
  median_diff_days = median(plasma_timing$plasma_diff_days),
  q25_diff_days = quantile(plasma_timing$plasma_diff_days, 0.25),
  q75_diff_days = quantile(plasma_timing$plasma_diff_days, 0.75),
  mean_abs_diff_days = mean(abs(plasma_timing$plasma_diff_days)),
  pct_after_anchor = 100 * mean(plasma_timing$plasma_diff_days > 0),
  pct_within_365 = 100 * mean(abs(plasma_timing$plasma_diff_days) <= 365),
  pct_within_180 = 100 * mean(abs(plasma_timing$plasma_diff_days) <= 180)
)

# Amyloid PET timing
amy_timing <- dx23_full %>% filter(!is.na(amy_diff_days))
timing_rows[[3]] <- tibble(
  biomarker = "Amyloid_PET",
  n_matched = nrow(amy_timing),
  median_diff_days = median(amy_timing$amy_diff_days),
  q25_diff_days = quantile(amy_timing$amy_diff_days, 0.25),
  q75_diff_days = quantile(amy_timing$amy_diff_days, 0.75),
  mean_abs_diff_days = mean(abs(amy_timing$amy_diff_days)),
  pct_after_anchor = 100 * mean(amy_timing$amy_diff_days > 0),
  pct_within_365 = 100 * mean(abs(amy_timing$amy_diff_days) <= 365),
  pct_within_180 = 100 * mean(abs(amy_timing$amy_diff_days) <= 180)
)

timing_summary <- bind_rows(timing_rows)
write_csv(timing_summary, file.path(out_dir, "phase2_biomarker_timing.csv"))
cat("  Biomarker timing saved.\n")

# ============================================================================
# PART 5: STRICTER TIMING SENSITIVITY (+/- 365 DAYS)
# ============================================================================

cat("\n== Part 5: Stricter timing sensitivity ==\n")

dx23_csf_strict <- dx23_full %>%
  filter(!is.na(csf_diff_days), abs(csf_diff_days) <= 365) %>%
  keep_complete(csf_vars)

strict_timing_rows <- list()

# Continuous PE
m_strict <- fit_cox(dx23_csf_strict, rhs("pe_resid_z_decr", csf_rhs))
strict_timing_rows[[1]] <- extract_term(
  m_strict, "pe_resid_z_decr",
  list(outcome_id = "DXSUM_23", timing = "strict_365", pe_model = "continuous",
       n = nrow(dx23_csf_strict), n_events = sum(dx23_csf_strict$event)))

# Grouped PE
m_strict_g <- fit_cox(dx23_csf_strict, rhs("PE_group", csf_rhs))
strict_timing_rows[[2]] <- extract_term(
  m_strict_g, "PE_groupLow_PE",
  list(outcome_id = "DXSUM_23", timing = "strict_365", pe_model = "grouped",
       n = nrow(dx23_csf_strict), n_events = sum(dx23_csf_strict$event)))

# Also run for MCI-only
mci_csf_strict <- mci %>%
  left_join(
    nearest_bio_full(csf_prepped, "biomarker_date", mci, c("CSF_PTAU_z", "CSF_ABETA42_z")) %>%
      transmute(RID, CSF_PTAU, CSF_PTAU_z, CSF_ABETA42, CSF_ABETA42_z,
                csf_diff_days = diff_days),
    by = "RID"
  ) %>%
  filter(!is.na(csf_diff_days), abs(csf_diff_days) <= 365) %>%
  keep_complete(csf_vars)

m_mci_strict <- fit_cox(mci_csf_strict, rhs("pe_resid_z_decr", csf_rhs))
strict_timing_rows[[3]] <- extract_term(
  m_mci_strict, "pe_resid_z_decr",
  list(outcome_id = "DXSUM_2_MCI_only", timing = "strict_365", pe_model = "continuous",
       n = nrow(mci_csf_strict), n_events = sum(mci_csf_strict$event)))

strict_timing_results <- bind_rows(strict_timing_rows)
write_csv(strict_timing_results, file.path(out_dir, "phase2_strict_timing_sensitivity.csv"))
cat("  Strict timing: n=", nrow(dx23_csf_strict), " events=", sum(dx23_csf_strict$event),
    " (original n=", nrow(dx23_csf), ")\n")

# ============================================================================
# PART 6: IPW SENSITIVITY FOR CSF SELECTION
# ============================================================================

cat("\n== Part 6: IPW sensitivity ==\n")

# Propensity model: P(has CSF | covariates)
ipw_data <- dx23_all %>%
  keep_complete(c("baseline_age_years", "Sex_F", "education_years", "APOE4",
                  "baseline_score_z", "pe_resid_z_decr"))

prop_model <- glm(has_csf ~ baseline_age_years + Sex_F + education_years + APOE4 +
                    baseline_score_z + pe_resid_z_decr,
                  data = ipw_data, family = binomial)

ipw_data$ps <- predict(prop_model, type = "response")
# Stabilized IPW: weight = P(has_csf) / ps for those with CSF
prev <- mean(ipw_data$has_csf)
ipw_data$ipw <- ifelse(ipw_data$has_csf, prev / ipw_data$ps, (1 - prev) / (1 - ipw_data$ps))
# Truncate weights at 1st and 99th percentile
ipw_q <- quantile(ipw_data$ipw[ipw_data$has_csf], c(0.01, 0.99), na.rm = TRUE)
ipw_data$ipw_trunc <- pmin(pmax(ipw_data$ipw, ipw_q[1]), ipw_q[2])

# IPW-weighted Cox on CSF subset
ipw_csf <- ipw_data %>%
  filter(has_csf) %>%
  keep_complete(c(csf_vars, "ipw_trunc"))

m_ipw <- tryCatch(
  coxph(Surv(time_to_event_years, event) ~ pe_resid_z_decr + Age_z + Sex_F +
          Education_z + APOE4 + baseline_score_z + CSF_PTAU_z + CSF_ABETA42_z,
        data = ipw_csf, weights = ipw_trunc),
  error = function(e) NULL
)

ipw_rows <- list()
ipw_rows[[1]] <- extract_term(
  m_ipw, "pe_resid_z_decr",
  list(outcome_id = "DXSUM_23", model = "IPW_weighted", pe_model = "continuous",
       n = nrow(ipw_csf), n_events = sum(ipw_csf$event)))

# Also unweighted for direct comparison
m_unw <- fit_cox(ipw_csf, rhs("pe_resid_z_decr", csf_rhs))
ipw_rows[[2]] <- extract_term(
  m_unw, "pe_resid_z_decr",
  list(outcome_id = "DXSUM_23", model = "unweighted_reference", pe_model = "continuous",
       n = nrow(ipw_csf), n_events = sum(ipw_csf$event)))

ipw_results <- bind_rows(ipw_rows)
write_csv(ipw_results, file.path(out_dir, "phase2_ipw_sensitivity.csv"))

# Propensity model summary
ps_summary <- tibble(
  variable = c(names(coef(prop_model))),
  coefficient = unname(coef(prop_model)),
  p_value = unname(summary(prop_model)$coefficients[, "Pr(>|z|)"]),
  prevalence_csf = prev,
  ipw_mean = mean(ipw_csf$ipw_trunc),
  ipw_median = median(ipw_csf$ipw_trunc),
  ipw_min = min(ipw_csf$ipw_trunc),
  ipw_max = max(ipw_csf$ipw_trunc)
)
write_csv(ps_summary, file.path(out_dir, "phase2_propensity_model.csv"))
cat("  IPW saved. PS prevalence=", round(prev, 3), "\n")

# ============================================================================
# PART 7: PHASE 2 KEY RESULTS SUMMARY
# ============================================================================

cat("\n== Saving key results ==\n")

phase2_key <- bind_rows(
  interaction_results %>% mutate(analysis = "interaction"),
  at_strat_results %>% mutate(analysis = "AT_stratification"),
  strict_timing_results %>% mutate(analysis = "strict_timing"),
  ipw_results %>% mutate(analysis = "IPW")
)
write_csv(phase2_key, file.path(out_dir, "phase2_key_results.csv"))
capture.output(sessionInfo(), file = file.path(out_dir, "phase2_session_info.txt"))

cat("\nPhase 2 complete. All outputs in:", out_dir, "\n")
