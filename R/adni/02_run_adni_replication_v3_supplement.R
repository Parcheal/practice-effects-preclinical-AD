

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
# ADNI replication v3 supplementary analyses for the PE manuscript
# - Phase 1 summaries and v3-based MCI-only analyses
# - Writes only to adni_replication_v3_supplement
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
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
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
    return(base %>% mutate(term = term, estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_))
  }
  s <- summary(model)
  if (!(term %in% rownames(s$coefficients))) {
    return(base %>% mutate(term = term, estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_))
  }
  ci <- tryCatch(suppressMessages(confint(model)), error = function(e) NULL)
  est <- exp(unname(s$coefficients[term, "coef"]))
  p <- unname(s$coefficients[term, "Pr(>|z|)"])
  lo <- if (!is.null(ci)) exp(ci[term, 1]) else NA_real_
  hi <- if (!is.null(ci)) exp(ci[term, 2]) else NA_real_
  base %>% mutate(term = term, estimate = est, conf.low = lo, conf.high = hi, p.value = p)
}

harrell_c <- function(model) {
  out <- tryCatch(summary(model)$concordance[1], error = function(e) NA_real_)
  as.numeric(out)
}

lrt_row <- function(df, base_rhs, add_term, outcome_id, pe_model) {
  m0 <- fit_cox(df, base_rhs)
  m1 <- fit_cox(df, rhs(add_term, base_rhs))
  if (is.null(m0) || is.null(m1)) {
    return(tibble(outcome_id = outcome_id, pe_model = pe_model, n = nrow(df), n_events = sum(df$event), C_base = NA_real_, C_add = NA_real_, delta_C = NA_real_, LRT_chisq = NA_real_, LRT_df = NA_real_, LRT_p = NA_real_))
  }
  ll0 <- as.numeric(logLik(m0))
  ll1 <- as.numeric(logLik(m1))
  df0 <- attr(logLik(m0), "df")
  df1 <- attr(logLik(m1), "df")
  chisq <- 2 * (ll1 - ll0)
  df_lrt <- df1 - df0
  tibble(
    outcome_id = outcome_id,
    pe_model = pe_model,
    n = nrow(df),
    n_events = sum(df$event),
    C_base = harrell_c(m0),
    C_add = harrell_c(m1),
    delta_C = harrell_c(m1) - harrell_c(m0),
    LRT_chisq = chisq,
    LRT_df = df_lrt,
    LRT_p = pchisq(chisq, df = df_lrt, lower.tail = FALSE)
  )
}

nearest_bio <- function(df, date_col, analysis_df, req) {
  df %>%
    filter(if_all(all_of(req), ~ !is.na(.))) %>%
    inner_join(analysis_df %>% select(RID, anchor_date, analysis_end_date), by = "RID") %>%
    filter(.data[[date_col]] <= analysis_end_date) %>%
    mutate(diff_days = as.numeric(.data[[date_col]] - anchor_date), abs_diff_days = abs(diff_days)) %>%
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

  post <- dxsum_df %>% inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>% filter(!is.na(EXAMDATE), EXAMDATE > anchor_date) %>% arrange(RID, EXAMDATE)
  ev <- post %>% group_by(RID) %>% summarise(first_event_date = if (any(DIAGNOSIS %in% c(2, 3), na.rm = TRUE)) EXAMDATE[which(DIAGNOSIS %in% c(2, 3))[1]] else as.Date(NA), last_followup_date = max(EXAMDATE, na.rm = TRUE), .groups = "drop") %>% mutate(event = !is.na(first_event_date), analysis_end_date = if_else(event, first_event_date, last_followup_date))
  subject_df %>% left_join(ev, by = "RID") %>% filter(!is.na(last_followup_date)) %>% mutate(outcome_id = "DXSUM_23", time_to_event_years = pmax(as.numeric(analysis_end_date - anchor_date) / 365.25, 0.01))
}

build_mci_only <- function(subject_df, dxsum_df) {
  post <- dxsum_df %>% inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>% filter(!is.na(EXAMDATE), EXAMDATE > anchor_date) %>% arrange(RID, EXAMDATE)
  ev <- post %>% group_by(RID) %>% summarise(first_mci_date = if (any(DIAGNOSIS == 2, na.rm = TRUE)) EXAMDATE[which(DIAGNOSIS == 2)[1]] else as.Date(NA), first_dem_date = if (any(DIAGNOSIS == 3, na.rm = TRUE)) EXAMDATE[which(DIAGNOSIS == 3)[1]] else as.Date(NA), last_followup_date = max(EXAMDATE, na.rm = TRUE), .groups = "drop") %>% mutate(event = !is.na(first_mci_date) & (is.na(first_dem_date) | first_mci_date <= first_dem_date), analysis_end_date = case_when(event ~ first_mci_date, !is.na(first_dem_date) ~ first_dem_date, TRUE ~ last_followup_date), censoring_reason = case_when(event ~ "incident_MCI", !is.na(first_dem_date) ~ "censored_at_dementia", TRUE ~ "censored_no_MCI"))
  subject_df %>% left_join(ev, by = "RID") %>% filter(!is.na(last_followup_date)) %>% mutate(outcome_id = "DXSUM_2_MCI_only", time_to_event_years = pmax(as.numeric(analysis_end_date - anchor_date) / 365.25, 0.01))
}

build_cdr_single <- function(subject_df, cdr_df) {
  pre <- cdr_df %>% inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>% filter(!is.na(VISDATE), VISDATE <= anchor_date, !is.na(CDGLOBAL)) %>% arrange(RID, VISDATE)
  eligible <- pre %>% group_by(RID) %>% summarise(last_pre_cdr = dplyr::last(CDGLOBAL), .groups = "drop") %>% filter(last_pre_cdr == 0)
  post <- cdr_df %>% semi_join(eligible, by = "RID") %>% inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>% filter(!is.na(VISDATE), VISDATE > anchor_date, !is.na(CDGLOBAL)) %>% arrange(RID, VISDATE)
  ev <- post %>% group_by(RID) %>% summarise(first_event_date = if (any(CDGLOBAL >= 0.5, na.rm = TRUE)) VISDATE[which(CDGLOBAL >= 0.5)[1]] else as.Date(NA), last_followup_date = max(VISDATE, na.rm = TRUE), .groups = "drop") %>% mutate(event = !is.na(first_event_date), analysis_end_date = if_else(event, first_event_date, last_followup_date))
  subject_df %>% semi_join(eligible, by = "RID") %>% left_join(ev, by = "RID") %>% filter(!is.na(last_followup_date)) %>% mutate(outcome_id = "CDR_single", time_to_event_years = pmax(as.numeric(analysis_end_date - anchor_date) / 365.25, 0.01))
}

dxsum <- read_csv(file.path(data_dir, "DXSUM_10Apr2026.csv"), show_col_types = FALSE) %>% mutate(EXAMDATE = pdate(EXAMDATE))
cdr <- read_csv(file.path(data_dir, "CDR_10Apr2026.csv"), show_col_types = FALSE) %>% mutate(VISDATE = pdate(VISDATE))
ptdemog <- read_csv(file.path(data_dir, "PTDEMOG_10Apr2026.csv"), show_col_types = FALSE) %>% mutate(VISDATE = pdate(VISDATE))
apoe_raw <- read_csv(file.path(data_dir, "APOERES_10Apr2026.csv"), show_col_types = FALSE) %>% mutate(APTESTDT = pdate(APTESTDT))
phc_cogn <- read_csv(file.path(data_dir, "ADSP_PHC_COGN_10Apr2026.csv"), show_col_types = FALSE) %>% mutate(EXAMDATE = pdate(EXAMDATE))
uwpsych <- read_csv(file.path(data_dir, "UWNPSYCHSUM_10Apr2026.csv"), show_col_types = FALSE) %>% mutate(EXAMDATE = pdate(EXAMDATE))
embic_qcp <- read_csv(file.path(data_dir, "EMBICqCP_10Apr2026.csv"), show_col_types = FALSE) %>% mutate(EXAMDATE = pdate(EXAMDATE))
plasma_fuji <- read_csv(file.path(data_dir, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_10Apr2026.csv"), show_col_types = FALSE) %>% mutate(EXAMDATE = pdate(EXAMDATE))
amyloid_pet <- read_csv(file.path(data_dir, "ADSP_PHC_PET_Amyloid_Simple_10Apr2026.csv"), show_col_types = FALSE) %>% mutate(PHC_SCANDATE = pdate(PHC_SCANDATE))
csf_elecsys <- read_csv(file.path(data_dir, "UPENNBIOMK_ROCHE_ELECSYS_10Apr2026.csv"), show_col_types = FALSE) %>% mutate(EXAMDATE = pdate(EXAMDATE))

baseline_cu <- dxsum %>% filter(VISCODE2 == "bl", DIAGNOSIS == 1) %>% distinct(RID)

ptdemog_subject <- ptdemog %>% arrange(RID, VISDATE) %>% group_by(RID) %>% summarise(PTGENDER = first_nm(as.numeric(PTGENDER[PTGENDER >= 0])), education_years = first_nm(as.numeric(PTEDUCAT[PTEDUCAT >= 0])), .groups = "drop")
age_subject <- phc_cogn %>% filter(VISCODE2 == "bl") %>% group_by(RID) %>% summarise(baseline_age_years = first_nm(as.numeric(PHC_Age_Cognition)), .groups = "drop")
apoe_subject <- apoe_raw %>% arrange(RID, APTESTDT) %>% group_by(RID) %>% summarise(APOE4 = ifelse(str_detect(first_nm(GENOTYPE), "4"), 1, 0), .groups = "drop")

subject_cov <- baseline_cu %>% left_join(age_subject, by = "RID") %>% left_join(ptdemog_subject, by = "RID") %>% left_join(apoe_subject, by = "RID") %>% mutate(Age_z = sscale(baseline_age_years), Sex_F = factor(ifelse(PTGENDER == 2, "Female", "Male"), levels = c("Male", "Female")), Education_z = sscale(education_years))

phc_mem_long <- phc_cogn %>% semi_join(baseline_cu, by = "RID") %>% transmute(RID, VISCODE2, source_date = EXAMDATE, score = if_else(as.numeric(PHC_MEM_PreciseFilter) == 1, as.numeric(PHC_MEM), NA_real_))
uw_mem_long <- uwpsych %>% semi_join(baseline_cu, by = "RID") %>% transmute(RID, VISCODE2, source_date = EXAMDATE, score = as.numeric(ADNI_MEM))
embic_qcp_long <- embic_qcp %>% semi_join(baseline_cu, by = "RID") %>% transmute(RID, VISCODE2, source_date = EXAMDATE, score = rowMeans(cbind(as.numeric(N1M), as.numeric(N2M), as.numeric(N3M), as.numeric(N4M), as.numeric(R1M), as.numeric(R2M), as.numeric(R3M), as.numeric(MIFR1M), as.numeric(MIFR2M), as.numeric(MIFR3M)), na.rm = TRUE))

plasma_prepped <- plasma_fuji %>% transmute(RID, biomarker_date = EXAMDATE, pT217_F = as.numeric(pT217_F), pT217_F_z = sscale(as.numeric(pT217_F)))
amy_prepped <- amyloid_pet %>% transmute(RID, biomarker_date = PHC_SCANDATE, amyloid_centiloids = as.numeric(PHC_CENTILOIDS), amyloid_centiloids_z = sscale(as.numeric(PHC_CENTILOIDS)))
csf_prepped <- csf_elecsys %>% transmute(RID, biomarker_date = EXAMDATE, CSF_PTAU = as.numeric(PTAU), CSF_PTAU_z = sscale(as.numeric(PTAU)), CSF_ABETA42 = as.numeric(ABETA42), CSF_ABETA42_z = sscale(as.numeric(ABETA42)))

add_layer <- function(df, layer) {
  out <- df
  if (layer == "csf_elecsys") {
    csf_match <- nearest_bio(csf_prepped, "biomarker_date", df, c("CSF_PTAU_z", "CSF_ABETA42_z")) %>% transmute(RID, CSF_PTAU, CSF_PTAU_z, CSF_ABETA42, CSF_ABETA42_z)
    out <- out %>% left_join(csf_match, by = "RID")
  }
  if (layer == "main_biomarker") {
    p_match <- nearest_bio(plasma_prepped, "biomarker_date", df, "pT217_F_z") %>% transmute(RID, pT217_F, pT217_F_z)
    a_match <- nearest_bio(amy_prepped, "biomarker_date", df, "amyloid_centiloids_z") %>% transmute(RID, amyloid_centiloids, amyloid_centiloids_z)
    out <- out %>% left_join(p_match, by = "RID") %>% left_join(a_match, by = "RID")
  }
  out
}

add_all_biomarkers <- function(df) {
  df %>%
    left_join(nearest_bio(csf_prepped, "biomarker_date", df, c("CSF_PTAU_z", "CSF_ABETA42_z")) %>% transmute(RID, CSF_PTAU, CSF_ABETA42), by = "RID") %>%
    left_join(nearest_bio(plasma_prepped, "biomarker_date", df, "pT217_F_z") %>% transmute(RID, pT217_F), by = "RID") %>%
    left_join(nearest_bio(amy_prepped, "biomarker_date", df, "amyloid_centiloids_z") %>% transmute(RID, amyloid_centiloids), by = "RID")
}

base_rhs <- rhs("Age_z", "Sex_F", "Education_z", "APOE4", "baseline_score_z")
csf_rhs <- rhs(base_rhs, "CSF_PTAU_z", "CSF_ABETA42_z")
mainbio_rhs <- rhs(base_rhs, "pT217_F_z", "amyloid_centiloids_z")

fmt_mean_sd <- function(x, d = 1) if (sum(!is.na(x)) == 0) "NA" else sprintf(paste0("%.", d, "f (%.", d, "f)"), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
fmt_med_iqr <- function(x, d = 2, n_tag = FALSE) if (sum(!is.na(x)) == 0) "NA" else paste0(sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), median(x, na.rm = TRUE), quantile(x, 0.25, na.rm = TRUE), quantile(x, 0.75, na.rm = TRUE)), if (n_tag) paste0(" (n=", sum(!is.na(x)), ")") else "")
fmt_npct <- function(x) if (sum(!is.na(x)) == 0) "NA" else sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
fmt_p <- function(p) if (!is.finite(p)) "" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
cont_p <- function(df, var, kruskal = FALSE) { dat <- df %>% select(PE_group, value = all_of(var)) %>% filter(!is.na(PE_group), !is.na(value)); if (nrow(dat) == 0 || n_distinct(dat$PE_group) < 2) return(NA_real_); tryCatch(if (kruskal) kruskal.test(value ~ PE_group, dat)$p.value else summary(aov(value ~ PE_group, dat))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_) }
bin_p <- function(df, var) { dat <- df %>% select(PE_group, value = all_of(var)) %>% filter(!is.na(PE_group), !is.na(value)); if (nrow(dat) == 0 || n_distinct(dat$value) < 2) return(NA_real_); tab <- table(dat$PE_group, dat$value); tryCatch(if (any(tab < 5)) fisher.test(tab)$p.value else chisq.test(tab)$p.value, error = function(e) NA_real_) }
row_out <- function(df, label, fun, p) { grp <- c("Overall", "High_PE", "Moderate_PE", "Low_PE"); vals <- sapply(grp, function(g) fun(if (g == "Overall") df else dplyr::filter(df, PE_group == g))); tibble(characteristic = label, Overall = vals[1], High_PE = vals[2], Moderate_PE = vals[3], Low_PE = vals[4], p_value = fmt_p(p)) }
event_row <- function(df, outcome_id) tibble(outcome_id = outcome_id, n = nrow(df), n_events = sum(df$event), median_event_time_years = median(df$time_to_event_years[df$event], na.rm = TRUE), q25_event_time_years = quantile(df$time_to_event_years[df$event], 0.25, na.rm = TRUE), q75_event_time_years = quantile(df$time_to_event_years[df$event], 0.75, na.rm = TRUE), pct_events_within_1y = 100 * mean(df$time_to_event_years[df$event] <= 1, na.rm = TRUE), pct_events_within_2y = 100 * mean(df$time_to_event_years[df$event] <= 2, na.rm = TRUE), pct_events_within_5y = 100 * mean(df$time_to_event_years[df$event] <= 5, na.rm = TRUE))
follow_row <- function(df, outcome_id) { cens <- df %>% filter(!event); tibble(outcome_id = outcome_id, n = nrow(df), n_events = sum(df$event), n_censored = nrow(cens), median_followup_censored_years = median(cens$time_to_event_years, na.rm = TRUE), q25_followup_censored_years = quantile(cens$time_to_event_years, 0.25, na.rm = TRUE), q75_followup_censored_years = quantile(cens$time_to_event_years, 0.75, na.rm = TRUE), median_followup_all_years = median(df$time_to_event_years, na.rm = TRUE), max_followup_years = max(df$time_to_event_years, na.rm = TRUE), total_person_years = sum(df$time_to_event_years, na.rm = TRUE)) }

specs <- tribble(
  ~exposure_id, ~window_mode, ~long_df,
  "PHC_MEM", "flex", "phc_mem_long",
  "UW_MEM", "flex", "uw_mem_long",
  "EMBIC_QCP", "flex", "embic_qcp_long",
  "PHC_MEM", "exact", "phc_mem_long",
  "UW_MEM", "exact", "uw_mem_long",
  "EMBIC_QCP", "exact", "embic_qcp_long"
)
lookup <- list(phc_mem_long = phc_mem_long, uw_mem_long = uw_mem_long, embic_qcp_long = embic_qcp_long)

mci_rows <- list()
mci_n_rows <- list()
phc_mem_flex_mci_csf <- NULL

for (i in seq_len(nrow(specs))) {
  subj <- compute_pe(lookup[[specs$long_df[[i]]]], specs$exposure_id[[i]], specs$window_mode[[i]]) %>% left_join(subject_cov, by = "RID")
  mci <- build_mci_only(subj, dxsum)
  mci_n_rows[[length(mci_n_rows) + 1]] <- mci %>% summarise(exposure_id = first(exposure_id), window_mode = first(window_mode), n = n(), n_events = sum(event), n_censored_at_dementia = sum(censoring_reason == "censored_at_dementia"), n_censored_no_mci = sum(censoring_reason == "censored_no_MCI"))
  for (layer in c("base", "csf_elecsys", "main_biomarker")) {
    dat <- add_layer(mci, layer)
    vars <- c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z", "APOE4", "baseline_score_z", "pe_resid_z_decr", "PE_group")
    if (layer == "csf_elecsys") vars <- c(vars, "CSF_PTAU_z", "CSF_ABETA42_z")
    if (layer == "main_biomarker") vars <- c(vars, "pT217_F_z", "amyloid_centiloids_z")
    dat <- keep_complete(dat, vars)
    if (nrow(dat) >= 100 && sum(dat$event) >= 30) {
      rhs_here <- if (layer == "base") base_rhs else if (layer == "csf_elecsys") csf_rhs else mainbio_rhs
      mci_rows[[length(mci_rows) + 1]] <- extract_term(fit_cox(dat, rhs("pe_resid_z_decr", rhs_here)), "pe_resid_z_decr", list(exposure_id = specs$exposure_id[[i]], window_mode = specs$window_mode[[i]], outcome_id = "DXSUM_2_MCI_only", biomarker_layer = layer, pe_model = "continuous", n = nrow(dat), n_events = sum(dat$event)))
      mci_rows[[length(mci_rows) + 1]] <- extract_term(fit_cox(dat, rhs("PE_group", rhs_here)), "PE_groupLow_PE", list(exposure_id = specs$exposure_id[[i]], window_mode = specs$window_mode[[i]], outcome_id = "DXSUM_2_MCI_only", biomarker_layer = layer, pe_model = "grouped", n = nrow(dat), n_events = sum(dat$event)))
      if (specs$exposure_id[[i]] == "PHC_MEM" && specs$window_mode[[i]] == "flex" && layer == "csf_elecsys") phc_mem_flex_mci_csf <- dat
    }
  }
}

mci_results <- bind_rows(mci_rows)
write_csv(bind_rows(mci_n_rows), file.path(out_dir, "mci_censoring_summary.csv"))
write_csv(mci_results, file.path(out_dir, "mci_clinical_results.csv"))
if (!is.null(phc_mem_flex_mci_csf)) write_csv(phc_mem_flex_mci_csf, file.path(out_dir, "PHC_MEM_flex_MCI_only_csf_dataset.csv"))

phc_mem_flex <- compute_pe(phc_mem_long, "PHC_MEM", "flex") %>% left_join(subject_cov, by = "RID")
dx23 <- build_dx23(phc_mem_flex, dxsum)
cdr1 <- build_cdr_single(phc_mem_flex, cdr)
dx23_csf <- keep_complete(add_layer(dx23, "csf_elecsys"), c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z", "APOE4", "baseline_score_z", "pe_resid_z_decr", "PE_group", "CSF_PTAU_z", "CSF_ABETA42_z"))
cdr1_csf <- keep_complete(add_layer(cdr1, "csf_elecsys"), c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z", "APOE4", "baseline_score_z", "pe_resid_z_decr", "PE_group", "CSF_PTAU_z", "CSF_ABETA42_z"))
dx23_fit <- fit_cox(dx23_csf, rhs("pe_resid_z_decr", csf_rhs))
cdr1_fit <- fit_cox(cdr1_csf, rhs("pe_resid_z_decr", csf_rhs))

primary_results <- bind_rows(
  extract_term(dx23_fit, "pe_resid_z_decr", list(exposure_id = "PHC_MEM", window_mode = "flex", outcome_id = "DXSUM_23", biomarker_layer = "csf_elecsys", pe_model = "continuous", n = nrow(dx23_csf), n_events = sum(dx23_csf$event))),
  extract_term(cdr1_fit, "pe_resid_z_decr", list(exposure_id = "PHC_MEM", window_mode = "flex", outcome_id = "CDR_single", biomarker_layer = "csf_elecsys", pe_model = "continuous", n = nrow(cdr1_csf), n_events = sum(cdr1_csf$event)))
)
write_csv(primary_results, file.path(out_dir, "phase1_primary_results.csv"))

table1_df <- add_all_biomarkers(dx23) %>% keep_complete(c("PE_group", "baseline_age_years", "Sex_F", "education_years", "APOE4", "baseline_score", "time_to_event_years"))
table1 <- bind_rows(
  row_out(table1_df, "N", function(x) as.character(nrow(x)), NA_real_),
  row_out(table1_df, "Age, years", function(x) fmt_mean_sd(x$baseline_age_years, 1), cont_p(table1_df, "baseline_age_years")),
  row_out(table1_df, "Female sex", function(x) fmt_npct(x$Sex_F == "Female"), bin_p(table1_df %>% mutate(is_female = Sex_F == "Female"), "is_female")),
  row_out(table1_df, "Education, years", function(x) fmt_mean_sd(x$education_years, 1), cont_p(table1_df, "education_years")),
  row_out(table1_df, "APOE4 carrier", function(x) fmt_npct(x$APOE4 == 1), bin_p(table1_df %>% mutate(apoe4 = APOE4 == 1), "apoe4")),
  row_out(table1_df, "Baseline PHC_MEM score", function(x) fmt_mean_sd(x$baseline_score, 3), cont_p(table1_df, "baseline_score")),
  row_out(table1_df, "Follow-up time, years", function(x) fmt_med_iqr(x$time_to_event_years, 2), cont_p(table1_df, "time_to_event_years", TRUE)),
  row_out(table1_df, "CSF pTau", function(x) fmt_med_iqr(x$CSF_PTAU, 2, TRUE), cont_p(table1_df, "CSF_PTAU", TRUE)),
  row_out(table1_df, "CSF Abeta42", function(x) fmt_med_iqr(x$CSF_ABETA42, 2, TRUE), cont_p(table1_df, "CSF_ABETA42", TRUE)),
  row_out(table1_df, "Plasma pTau217", function(x) fmt_med_iqr(x$pT217_F, 2, TRUE), cont_p(table1_df, "pT217_F", TRUE)),
  row_out(table1_df, "Amyloid centiloids", function(x) fmt_med_iqr(x$amyloid_centiloids, 2, TRUE), cont_p(table1_df, "amyloid_centiloids", TRUE))
)
write_csv(table1_df, file.path(out_dir, "Table1_primary_cohort_dataset.csv"))
write_csv(table1, file.path(out_dir, "Table1_primary_cohort.csv"))

event_summary <- bind_rows(event_row(dx23_csf, "DXSUM_23"), event_row(cdr1_csf, "CDR_single"), if (!is.null(phc_mem_flex_mci_csf)) event_row(phc_mem_flex_mci_csf, "DXSUM_2_MCI_only") else NULL)
follow_summary <- bind_rows(follow_row(dx23_csf, "DXSUM_23"), follow_row(cdr1_csf, "CDR_single"), if (!is.null(phc_mem_flex_mci_csf)) follow_row(phc_mem_flex_mci_csf, "DXSUM_2_MCI_only") else NULL)
write_csv(event_summary, file.path(out_dir, "event_time_summary.csv"))
write_csv(follow_summary, file.path(out_dir, "followup_summary.csv"))

tidy_zph <- function(zph, outcome_id) as_tibble(as.data.frame(zph$table), rownames = "term") %>% mutate(outcome_id = outcome_id, .before = 1) %>% rename(chisq = `chisq`, df = `df`, p.value = `p`)
ph_details <- bind_rows(tidy_zph(cox.zph(dx23_fit), "DXSUM_23"), tidy_zph(cox.zph(cdr1_fit), "CDR_single"))
ph_summary <- ph_details %>% filter(term %in% c("pe_resid_z_decr", "GLOBAL"))
write_csv(ph_details, file.path(out_dir, "ph_assumption_details.csv"))
write_csv(ph_summary, file.path(out_dir, "ph_assumption_summary.csv"))

incremental <- bind_rows(lrt_row(dx23_csf, csf_rhs, "pe_resid_z_decr", "DXSUM_23", "continuous"), lrt_row(dx23_csf, csf_rhs, "PE_group", "DXSUM_23", "grouped"), lrt_row(cdr1_csf, csf_rhs, "pe_resid_z_decr", "CDR_single", "continuous"), lrt_row(cdr1_csf, csf_rhs, "PE_group", "CDR_single", "grouped"))
write_csv(incremental, file.path(out_dir, "incremental_discrimination_summary.csv"))

key_results <- bind_rows(primary_results, mci_results %>% filter(exposure_id == "PHC_MEM", window_mode == "flex", biomarker_layer == "csf_elecsys", pe_model == "continuous"))
write_csv(key_results, file.path(out_dir, "phase1_key_results.csv"))
capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))
