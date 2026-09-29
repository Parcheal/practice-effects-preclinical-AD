

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
# ADNI replication v3 for the PE manuscript
# - Re-selection pass designed for ADNI-native validation
# - Primary ideas:
#   * memory-centered exposure lines
#   * flexible early window: bl / m06 / (m12 or m24)
#   * DXSUM plus single-CDR progression as practical clinical outcomes
#   * CSF-forward biomarker strategy, with plasma+amyloid as A4-like sensitivity
# - Does not overwrite ADNI v1 or v2 outputs
# ============================================================================

suppressPackageStartupMessages({
  library(broom)
  library(dplyr)
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

out_dir <- file.path(pkg_root, "verify", "adni", "v3")
data_dir <- file.path(data_root, "adni")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260424)
bootstrap_reps <- 1000
cv_k <- 5
cv_repeats <- 10

clinical_n_min <- 100
clinical_events_min <- 30
lmm_subjects_min <- 100
lmm_rows_min <- 300
lmm_repeat_subjects_min <- 100

parse_date_safe <- function(x) {
  as.Date(x)
}

parse_viscode2_months <- function(x) {
  x <- tolower(as.character(x))
  out <- suppressWarnings(as.numeric(str_extract(x, "\\d+")))
  out[!is.na(x) & x %in% c("bl", "sc", "f")] <- 0
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

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }
  x[[1]]
}

pe_group_from_z <- function(z) {
  factor(
    case_when(
      is.na(z) ~ NA_character_,
      z > 0.2 ~ "High_PE",
      z < -0.2 ~ "Low_PE",
      TRUE ~ "Moderate_PE"
    ),
    levels = c("High_PE", "Moderate_PE", "Low_PE")
  )
}

combine_rhs <- function(...) {
  terms <- c(...)
  terms <- terms[!is.na(terms) & nzchar(terms)]
  paste(terms, collapse = " + ")
}

filter_complete_cases <- function(df, vars) {
  vars <- unique(vars)
  vars <- vars[vars %in% names(df)]
  if (length(vars) == 0) {
    return(df)
  }
  df %>% filter(if_all(all_of(vars), ~ !is.na(.)))
}

fit_cox_model <- function(data, rhs) {
  tryCatch(
    coxph(as.formula(paste0("Surv(time_to_event_years, event) ~ ", rhs)), data = data),
    error = function(e) NULL
  )
}

fit_lmm_model <- function(data, rhs) {
  tryCatch(
    lmerTest::lmer(as.formula(paste0("outcome_value ~ ", rhs, " + (1 + time_years | RID)")), data = data, REML = FALSE),
    error = function(e) NULL
  )
}

collect_cox_term <- function(model, exposure_id, window_mode, outcome_id, biomarker_layer, pe_model, term_name, n, n_events) {
  if (is.null(model)) {
    return(tibble(
      exposure_id = exposure_id,
      window_mode = window_mode,
      outcome_id = outcome_id,
      biomarker_layer = biomarker_layer,
      pe_model = pe_model,
      term = term_name,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      p.value = NA_real_,
      n = n,
      n_events = n_events
    ))
  }

  out <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == term_name) %>%
    select(term, estimate, conf.low, conf.high, p.value)

  if (nrow(out) == 0) {
    out <- tibble(term = term_name, estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_)
  }

  out %>%
    mutate(
      exposure_id = exposure_id,
      window_mode = window_mode,
      outcome_id = outcome_id,
      biomarker_layer = biomarker_layer,
      pe_model = pe_model,
      n = n,
      n_events = n_events,
      .before = 1
    )
}

collect_lmm_term <- function(model, exposure_id, window_mode, target_id, biomarker_layer, term_name, n_subjects, n_rows, n_repeat_subjects) {
  if (is.null(model)) {
    return(tibble(
      exposure_id = exposure_id,
      window_mode = window_mode,
      target_id = target_id,
      biomarker_layer = biomarker_layer,
      term = term_name,
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      p.value = NA_real_,
      n_subjects = n_subjects,
      n_rows = n_rows,
      n_repeat_subjects = n_repeat_subjects
    ))
  }

  out <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    filter(term == term_name) %>%
    select(term, estimate, conf.low, conf.high, p.value)

  if (nrow(out) == 0) {
    out <- tibble(term = term_name, estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_)
  }

  out %>%
    mutate(
      exposure_id = exposure_id,
      window_mode = window_mode,
      target_id = target_id,
      biomarker_layer = biomarker_layer,
      n_subjects = n_subjects,
      n_rows = n_rows,
      n_repeat_subjects = n_repeat_subjects,
      .before = 1
    )
}

harrell_c <- function(model) {
  out <- tryCatch(summary(model)$concordance[1], error = function(e) NA_real_)
  as.numeric(out)
}

harrell_c_from_lp <- function(time, status, lp) {
  cc <- tryCatch(survival::concordance(Surv(time, status) ~ lp), error = function(e) NULL)
  if (is.null(cc) || is.null(cc$concordance)) {
    return(NA_real_)
  }
  unname(1 - cc$concordance)
}

assign_stratified_folds <- function(status, k) {
  fold_id <- integer(length(status))
  for (level in unique(status)) {
    idx <- which(status == level)
    fold_id[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  fold_id
}

fit_blup_slope_model <- function(long_use) {
  lmer(
    score ~ time_years + (1 + time_years | RID),
    data = long_use,
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
}

predict_blup_slopes_from_fit <- function(fit, target_long) {
  beta <- lme4::fixef(fit)
  random_cov <- as.matrix(lme4::VarCorr(fit)$RID)
  sigma2 <- sigma(fit)^2

  bind_rows(lapply(split(target_long, target_long$RID), function(dat) {
    dat <- dat %>% arrange(time_years)
    x_mat <- cbind(`(Intercept)` = 1, time_years = dat$time_years)
    y_vec <- dat$score
    resid_vec <- as.numeric(y_vec - x_mat %*% beta[c("(Intercept)", "time_years")])
    v_mat <- x_mat %*% random_cov %*% t(x_mat) + diag(sigma2, nrow(x_mat))
    b_hat <- random_cov %*% t(x_mat) %*% solve(v_mat, resid_vec)
    tibble(RID = dat$RID[[1]], raw_slope_blup = unname(beta[["time_years"]] + b_hat[2, 1]))
  }))
}

estimate_blup_slopes <- function(long_use) {
  fit <- fit_blup_slope_model(long_use)
  predict_blup_slopes_from_fit(fit, long_use)
}

add_fold_specific_pe <- function(train, test, early_long = NULL, slope_col = "raw_slope_blup", baseline_col = "baseline_score") {
  if (slope_col == "raw_slope_blup") {
    if (is.null(early_long)) {
      stop("early_long is required for fold-specific BLUP PE estimation.")
    }

    train_long <- early_long %>% semi_join(train %>% select(RID), by = "RID")
    test_long <- early_long %>% semi_join(test %>% select(RID), by = "RID")
    blup_fit <- fit_blup_slope_model(train_long)

    train_slopes <- predict_blup_slopes_from_fit(blup_fit, train_long)
    test_slopes <- predict_blup_slopes_from_fit(blup_fit, test_long)

    train <- train %>%
      select(-any_of(c("pe_cv_slope", "pe_cv_z", "pe_cv_decr"))) %>%
      left_join(train_slopes %>% rename(pe_cv_slope = raw_slope_blup), by = "RID")
    test <- test %>%
      select(-any_of(c("pe_cv_slope", "pe_cv_z", "pe_cv_decr"))) %>%
      left_join(test_slopes %>% rename(pe_cv_slope = raw_slope_blup), by = "RID")
    slope_col <- "pe_cv_slope"
  }

  resid_fit <- lm(as.formula(paste0(slope_col, " ~ ", baseline_col)), data = train)
  train_resid <- resid(resid_fit)
  train_mean <- mean(train_resid, na.rm = TRUE)
  train_sd <- sd(train_resid, na.rm = TRUE)
  if (!is.finite(train_sd) || train_sd == 0) {
    train_sd <- 1
  }

  train <- train %>%
    mutate(
      pe_cv_z = (train_resid - train_mean) / train_sd,
      pe_cv_decr = -pe_cv_z
    )

  test_resid <- test[[slope_col]] - predict(resid_fit, newdata = test)
  test <- test %>%
    mutate(
      pe_cv_z = (test_resid - train_mean) / train_sd,
      pe_cv_decr = -pe_cv_z
    )

  list(train = train, test = test)
}

bootstrap_delta_c <- function(data, base_rhs, add_rhs, n_boot = 500) {
  out <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(nrow(data), size = nrow(data), replace = TRUE)
    boot_dat <- data[idx, , drop = FALSE]
    m_base <- fit_cox_model(boot_dat, base_rhs)
    m_add <- fit_cox_model(boot_dat, add_rhs)
    if (!is.null(m_base) && !is.null(m_add)) {
      out[i] <- harrell_c(m_add) - harrell_c(m_base)
    }
  }
  out
}

run_nested_cv_delta_c <- function(data, base_rhs, k = 5, repeats = 10, early_long = NULL) {
  rows <- list()
  idx_row <- 1

  for (rep_i in seq_len(repeats)) {
    fold_id <- assign_stratified_folds(data$event, k)
    for (fold_i in seq_len(k)) {
      train <- data[fold_id != fold_i, , drop = FALSE]
      test <- data[fold_id == fold_i, , drop = FALSE]

      if (nrow(test) == 0 || sum(test$event) == 0) {
        next
      }

      fold_dat <- add_fold_specific_pe(train, test, early_long = early_long)
      train <- fold_dat$train
      test <- fold_dat$test

      m_base <- fit_cox_model(train, base_rhs)
      m_add <- fit_cox_model(train, combine_rhs("pe_cv_decr", base_rhs))
      if (is.null(m_base) || is.null(m_add)) {
        next
      }

      lp_base <- tryCatch(predict(m_base, newdata = test, type = "lp"), error = function(e) NULL)
      lp_add <- tryCatch(predict(m_add, newdata = test, type = "lp"), error = function(e) NULL)
      if (is.null(lp_base) || is.null(lp_add)) {
        next
      }

      rows[[idx_row]] <- tibble(
        repeat_id = rep_i,
        fold = fold_i,
        n_train = nrow(train),
        n_test = nrow(test),
        events_test = sum(test$event),
        c_base = harrell_c_from_lp(test$time_to_event_years, test$event, lp_base),
        c_add = harrell_c_from_lp(test$time_to_event_years, test$event, lp_add)
      ) %>%
        mutate(delta_c = c_add - c_base)

      idx_row <- idx_row + 1
    }
  }

  if (length(rows) == 0) {
    return(tibble())
  }

  bind_rows(rows)
}

nearest_biomarker <- function(df, date_col, analysis_df, required_cols) {
  df %>%
    filter(if_all(all_of(required_cols), ~ !is.na(.))) %>%
    inner_join(analysis_df %>% select(RID, anchor_date, analysis_end_date), by = "RID") %>%
    filter(.data[[date_col]] <= analysis_end_date) %>%
    mutate(
      diff_days = as.numeric(.data[[date_col]] - anchor_date),
      abs_diff_days = abs(diff_days)
    ) %>%
    arrange(RID, abs_diff_days, diff_days) %>%
    group_by(RID) %>%
    slice(1) %>%
    ungroup()
}

select_pe_window <- function(long_df, window_mode) {
  long_df <- long_df %>%
    mutate(visit_months = parse_viscode2_months(VISCODE2))

  if (window_mode == "exact") {
    selected_ids <- long_df %>%
      filter(VISCODE2 %in% c("bl", "m06", "m12"), !is.na(score)) %>%
      distinct(RID, VISCODE2) %>%
      count(RID) %>%
      filter(n == 3) %>%
      select(RID)

    long_use <- long_df %>%
      semi_join(selected_ids, by = "RID") %>%
      filter(VISCODE2 %in% c("bl", "m06", "m12"))
  } else {
    elig <- long_df %>%
      filter(VISCODE2 %in% c("bl", "m06", "m12", "m24"), !is.na(score))

    selected_ids <- elig %>%
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

    long_use <- elig %>%
      inner_join(selected_ids, by = "RID") %>%
      filter(VISCODE2 %in% c("bl", "m06") | VISCODE2 == chosen_end)
  }

  long_use %>%
    filter(!is.na(score), !is.na(visit_months)) %>%
    mutate(time_years = visit_months / 12)
}

compute_pe_subjects <- function(long_df, exposure_id, exposure_label, source_role, window_mode) {
  long_use <- select_pe_window(long_df, window_mode)

  subject_df <- long_use %>%
    arrange(RID, visit_months) %>%
    group_by(RID) %>%
    summarise(
      PTID = first_non_missing(PTID),
      anchor_visit = dplyr::last(VISCODE2),
      anchor_date = dplyr::last(source_date),
      baseline_score = score[which.min(visit_months)],
      raw_slope_ols = coef(lm(score ~ I(visit_months / 12)))[["I(visit_months/12)"]],
      .groups = "drop"
    )

  blup_slopes <- estimate_blup_slopes(long_use)
  subject_df <- subject_df %>% left_join(blup_slopes, by = "RID")

  resid_fit <- lm(raw_slope_blup ~ baseline_score, data = subject_df)
  ols_resid_fit <- lm(raw_slope_ols ~ baseline_score, data = subject_df)

  out <- subject_df %>%
    mutate(
      exposure_id = exposure_id,
      exposure_label = exposure_label,
      source_role = source_role,
      window_mode = window_mode,
      pe_estimator = "BLUP",
      pe_resid_raw = resid(resid_fit),
      pe_resid_z = safe_scale(pe_resid_raw),
      pe_resid_z_decr = -pe_resid_z,
      PE_group = pe_group_from_z(pe_resid_z),
      ols_pe_resid_raw = resid(ols_resid_fit),
      ols_pe_resid_z = safe_scale(ols_pe_resid_raw),
      ols_pe_resid_z_decr = -ols_pe_resid_z,
      PE_group_ols = pe_group_from_z(ols_pe_resid_z),
      baseline_score_z = safe_scale(baseline_score)
    )

  attr(out, "early_long") <- long_use
  out
}

build_dxsum_outcome <- function(subject_df, dxsum_df, exclude_anchor_prevalent = TRUE) {
  if (isTRUE(exclude_anchor_prevalent)) {
    # Drop participants with a DXSUM MCI/dementia diagnosis on or before
    # the PE anchor date. Without this filter the post-anchor risk-set
    # silently includes anchor-prevalent cases, inflating events.
    prevalent_ids <- dxsum_df %>%
      inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>%
      filter(
        !is.na(EXAMDATE), !is.na(DIAGNOSIS),
        EXAMDATE <= anchor_date,
        DIAGNOSIS %in% c(2, 3)
      ) %>%
      distinct(RID)
    attr(subject_df, "anchor_prevalent_excluded") <- nrow(prevalent_ids)
    subject_df <- subject_df %>% anti_join(prevalent_ids, by = "RID")
  }

  dx_post <- dxsum_df %>%
    inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>%
    filter(!is.na(EXAMDATE), EXAMDATE > anchor_date) %>%
    arrange(RID, EXAMDATE)

  dx_post %>%
    group_by(RID) %>%
    summarise(
      n_post = n(),
      first_event_date = if (any(DIAGNOSIS %in% c(2, 3), na.rm = TRUE)) {
        EXAMDATE[which(DIAGNOSIS %in% c(2, 3))[1]]
      } else {
        as.Date(NA)
      },
      last_followup_date = max(EXAMDATE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      event = !is.na(first_event_date),
      analysis_end_date = if_else(event, first_event_date, last_followup_date),
      time_to_event_years = pmax(as.numeric(analysis_end_date - subject_df$anchor_date[match(RID, subject_df$RID)]) / 365.25, 0.01)
    ) %>%
    right_join(subject_df, by = "RID") %>%
    filter(!is.na(last_followup_date)) %>%
    mutate(outcome_id = "DXSUM_23")
}

build_cdr_single_outcome <- function(subject_df, cdr_df) {
  cdr_pre <- cdr_df %>%
    inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>%
    filter(!is.na(VISDATE), VISDATE <= anchor_date, !is.na(CDGLOBAL)) %>%
    arrange(RID, VISDATE)

  eligible_pre <- cdr_pre %>%
    group_by(RID) %>%
    summarise(last_pre_cdr = dplyr::last(CDGLOBAL), .groups = "drop") %>%
    filter(last_pre_cdr == 0)

  cdr_post <- cdr_df %>%
    semi_join(eligible_pre, by = "RID") %>%
    inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>%
    filter(!is.na(VISDATE), VISDATE > anchor_date, !is.na(CDGLOBAL)) %>%
    arrange(RID, VISDATE)

  cdr_post %>%
    group_by(RID) %>%
    summarise(
      n_post = n(),
      first_event_date = if (any(CDGLOBAL >= 0.5, na.rm = TRUE)) {
        VISDATE[which(CDGLOBAL >= 0.5)[1]]
      } else {
        as.Date(NA)
      },
      last_followup_date = max(VISDATE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      event = !is.na(first_event_date),
      analysis_end_date = if_else(event, first_event_date, last_followup_date),
      time_to_event_years = pmax(as.numeric(analysis_end_date - subject_df$anchor_date[match(RID, subject_df$RID)]) / 365.25, 0.01)
    ) %>%
    right_join(subject_df %>% semi_join(eligible_pre, by = "RID"), by = "RID") %>%
    filter(!is.na(last_followup_date)) %>%
    mutate(outcome_id = "CDR_single")
}

build_cdr_confirmed_outcome <- function(subject_df, cdr_df) {
  cdr_pre <- cdr_df %>%
    inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>%
    filter(!is.na(VISDATE), VISDATE <= anchor_date, !is.na(CDGLOBAL)) %>%
    arrange(RID, VISDATE)

  eligible_pre <- cdr_pre %>%
    group_by(RID) %>%
    summarise(last_pre_cdr = dplyr::last(CDGLOBAL), .groups = "drop") %>%
    filter(last_pre_cdr == 0)

  cdr_post <- cdr_df %>%
    semi_join(eligible_pre, by = "RID") %>%
    inner_join(subject_df %>% select(RID, anchor_date), by = "RID") %>%
    filter(!is.na(VISDATE), VISDATE > anchor_date, !is.na(CDGLOBAL)) %>%
    arrange(RID, VISDATE)

  cdr_post %>%
    group_by(RID) %>%
    summarise(
      n_post = n(),
      first_event_date = {
        idx <- which(CDGLOBAL[-length(CDGLOBAL)] >= 0.5 & CDGLOBAL[-1] >= 0.5)
        if (length(idx) == 0) as.Date(NA) else VISDATE[idx[1]]
      },
      last_followup_date = max(VISDATE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      event = !is.na(first_event_date),
      analysis_end_date = if_else(event, first_event_date, last_followup_date),
      time_to_event_years = pmax(as.numeric(analysis_end_date - subject_df$anchor_date[match(RID, subject_df$RID)]) / 365.25, 0.01)
    ) %>%
    right_join(subject_df %>% semi_join(eligible_pre, by = "RID"), by = "RID") %>%
    filter(!is.na(last_followup_date)) %>%
    mutate(outcome_id = "CDR_confirmed")
}

build_longitudinal_post_anchor <- function(subject_df, long_df, target_var, target_id) {
  long_df %>%
    inner_join(
      subject_df %>%
        select(
          RID, exposure_id, exposure_label, source_role, window_mode,
          anchor_date, raw_slope_ols, baseline_score, baseline_score_z,
          pe_resid_z_decr, PE_group
        ),
      by = "RID"
    ) %>%
    filter(!is.na(EXAMDATE), EXAMDATE > anchor_date) %>%
    mutate(
      target_id = target_id,
      outcome_value = .data[[target_var]],
      time_years = as.numeric(EXAMDATE - anchor_date) / 365.25
    ) %>%
    filter(!is.na(outcome_value), is.finite(time_years), time_years > 0)
}

# ----------------------------
# Load source data
# ----------------------------
dxsum <- read_csv(file.path(data_dir, "DXSUM_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(EXAMDATE = parse_date_safe(EXAMDATE))

cdr <- read_csv(file.path(data_dir, "CDR_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(VISDATE = parse_date_safe(VISDATE))

ptdemog <- read_csv(file.path(data_dir, "PTDEMOG_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(VISDATE = parse_date_safe(VISDATE))

apoe_raw <- read_csv(file.path(data_dir, "APOERES_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(APTESTDT = parse_date_safe(APTESTDT))

phc_cogn <- read_csv(file.path(data_dir, "ADSP_PHC_COGN_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(EXAMDATE = parse_date_safe(EXAMDATE))

uwpsych <- read_csv(file.path(data_dir, "UWNPSYCHSUM_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(EXAMDATE = parse_date_safe(EXAMDATE))

embic_qcp <- read_csv(file.path(data_dir, "EMBICqCP_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(EXAMDATE = parse_date_safe(EXAMDATE))

embic_dcb <- read_csv(file.path(data_dir, "ADNI_EMBICDCB_10Apr2026.csv"), show_col_types = FALSE)

plasma_fuji <- read_csv(file.path(data_dir, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(EXAMDATE = parse_date_safe(EXAMDATE))

amyloid_pet <- read_csv(file.path(data_dir, "ADSP_PHC_PET_Amyloid_Simple_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(PHC_SCANDATE = parse_date_safe(PHC_SCANDATE))

csf_elecsys <- read_csv(file.path(data_dir, "UPENNBIOMK_ROCHE_ELECSYS_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(EXAMDATE = parse_date_safe(EXAMDATE))

phc_biomarker <- read_csv(file.path(data_dir, "ADSP_PHC_BIOMARKER_10Apr2026.csv"), show_col_types = FALSE) %>%
  mutate(DRAWDATE = parse_date_safe(DRAWDATE))

baseline_cu <- dxsum %>%
  filter(VISCODE2 == "bl", DIAGNOSIS == 1) %>%
  distinct(RID)

visit_date_lookup <- dxsum %>%
  filter(!is.na(EXAMDATE)) %>%
  group_by(RID, VISCODE2) %>%
  summarise(source_date = min(EXAMDATE), .groups = "drop")

ptdemog_subject <- ptdemog %>%
  arrange(RID, VISDATE) %>%
  group_by(RID) %>%
  summarise(
    PTGENDER = first_non_missing(as.numeric(PTGENDER[PTGENDER >= 0])),
    PTEDUCAT = first_non_missing(as.numeric(PTEDUCAT[PTEDUCAT >= 0])),
    .groups = "drop"
  )

phc_baseline_demo <- phc_cogn %>%
  filter(VISCODE2 == "bl") %>%
  group_by(RID) %>%
  summarise(
    baseline_age_years = first_non_missing(as.numeric(PHC_Age_Cognition)),
    .groups = "drop"
  )

apoe_subject <- apoe_raw %>%
  arrange(RID, APTESTDT) %>%
  group_by(RID) %>%
  summarise(
    APOE_genotype = first_non_missing(GENOTYPE),
    APOE4 = case_when(
      is.na(APOE_genotype) ~ NA_real_,
      str_detect(APOE_genotype, "4") ~ 1,
      TRUE ~ 0
    ),
    .groups = "drop"
  )

subject_covariates <- baseline_cu %>%
  left_join(phc_baseline_demo, by = "RID") %>%
  left_join(ptdemog_subject, by = "RID") %>%
  left_join(apoe_subject %>% select(RID, APOE4), by = "RID") %>%
  mutate(
    Age_z = safe_scale(baseline_age_years),
    Sex_F = factor(ifelse(PTGENDER == 2, "Female", "Male"), levels = c("Male", "Female")),
    Education_z = safe_scale(PTEDUCAT)
  ) %>%
  select(RID, Age_z, Sex_F, Education_z, APOE4)

phc_global_long <- phc_cogn %>%
  semi_join(baseline_cu, by = "RID") %>%
  mutate(
    PHC_MEM_use = if_else(as.numeric(PHC_MEM_PreciseFilter) == 1, as.numeric(PHC_MEM), NA_real_),
    PHC_EXF_use = if_else(as.numeric(PHC_EXF_PreciseFilter) == 1, as.numeric(PHC_EXF), NA_real_),
    PHC_LAN_use = if_else(as.numeric(PHC_LAN_PreciseFilter) == 1, as.numeric(PHC_LAN), NA_real_),
    PHC_VSP_use = if_else(as.numeric(PHC_VSP_PreciseFilter) == 1, as.numeric(PHC_VSP), NA_real_),
    precise_n = rowSums(!is.na(cbind(PHC_MEM_use, PHC_EXF_use, PHC_LAN_use, PHC_VSP_use))),
    score = if_else(precise_n >= 3, rowMeans(cbind(PHC_MEM_use, PHC_EXF_use, PHC_LAN_use, PHC_VSP_use), na.rm = TRUE), NA_real_),
    source_date = EXAMDATE
  ) %>%
  select(RID, PTID, VISCODE2, source_date, score)

phc_mem_long <- phc_cogn %>%
  semi_join(baseline_cu, by = "RID") %>%
  transmute(
    RID,
    PTID,
    VISCODE2,
    source_date = EXAMDATE,
    score = if_else(as.numeric(PHC_MEM_PreciseFilter) == 1, as.numeric(PHC_MEM), NA_real_)
  )

uw_global_long <- uwpsych %>%
  semi_join(baseline_cu, by = "RID") %>%
  transmute(
    RID,
    PTID = NA_character_,
    VISCODE2,
    source_date = EXAMDATE,
    score = rowMeans(cbind(as.numeric(ADNI_MEM), as.numeric(ADNI_EF), as.numeric(ADNI_LAN), as.numeric(ADNI_VS)), na.rm = TRUE)
  )

uw_mem_long <- uwpsych %>%
  semi_join(baseline_cu, by = "RID") %>%
  transmute(
    RID,
    PTID = NA_character_,
    VISCODE2,
    source_date = EXAMDATE,
    score = as.numeric(ADNI_MEM)
  )

embic_qcp_long <- embic_qcp %>%
  semi_join(baseline_cu, by = "RID") %>%
  transmute(
    RID,
    PTID = NA_character_,
    VISCODE2,
    source_date = EXAMDATE,
    score = rowMeans(
      cbind(
        as.numeric(N1M), as.numeric(N2M), as.numeric(N3M), as.numeric(N4M),
        as.numeric(R1M), as.numeric(R2M), as.numeric(R3M),
        as.numeric(MIFR1M), as.numeric(MIFR2M), as.numeric(MIFR3M)
      ),
      na.rm = TRUE
    )
  )

embic_dcb_long <- embic_dcb %>%
  semi_join(baseline_cu, by = "RID") %>%
  left_join(visit_date_lookup, by = c("RID", "VISCODE2")) %>%
  transmute(
    RID,
    PTID = NA_character_,
    VISCODE2,
    source_date,
    score = rowMeans(
      cbind(
        as.numeric(N1), as.numeric(N2), as.numeric(N3), as.numeric(N4),
        as.numeric(R1), as.numeric(R2), as.numeric(R3),
        as.numeric(M1), as.numeric(M2), as.numeric(M3)
      ),
      na.rm = TRUE
    )
  )

plasma_prepped <- plasma_fuji %>%
  transmute(RID, biomarker_date = EXAMDATE, pT217_F = as.numeric(pT217_F)) %>%
  mutate(pT217_F_z = safe_scale(pT217_F))

amyloid_prepped <- amyloid_pet %>%
  transmute(RID, biomarker_date = PHC_SCANDATE, amyloid_centiloids = as.numeric(PHC_CENTILOIDS)) %>%
  mutate(amyloid_centiloids_z = safe_scale(amyloid_centiloids))

csf_prepped <- csf_elecsys %>%
  transmute(
    RID,
    biomarker_date = EXAMDATE,
    CSF_ABETA42 = as.numeric(ABETA42),
    CSF_PTAU = as.numeric(PTAU)
  ) %>%
  mutate(
    CSF_ABETA42_z = safe_scale(CSF_ABETA42),
    CSF_PTAU_z = safe_scale(CSF_PTAU)
  )

harmonized_csf_prepped <- phc_biomarker %>%
  transmute(
    RID,
    biomarker_date = DRAWDATE,
    PHC_AB42 = as.numeric(PHC_AB42),
    PHC_pTau = as.numeric(PHC_pTau)
  ) %>%
  mutate(
    PHC_AB42_z = safe_scale(PHC_AB42),
    PHC_pTau_z = safe_scale(PHC_pTau)
  )

uw_longitudinal <- uwpsych %>%
  semi_join(baseline_cu, by = "RID") %>%
  transmute(
    RID,
    EXAMDATE,
    UW_GLOBAL = rowMeans(cbind(as.numeric(ADNI_MEM), as.numeric(ADNI_EF), as.numeric(ADNI_LAN), as.numeric(ADNI_VS)), na.rm = TRUE),
    ADNI_MEM = as.numeric(ADNI_MEM)
  )

add_biomarker_layer <- function(analysis_df, biomarker_layer) {
  out <- analysis_df

  if (biomarker_layer == "main_biomarker") {
    plasma_match <- nearest_biomarker(plasma_prepped, "biomarker_date", analysis_df, "pT217_F_z") %>%
      transmute(RID, pT217_F, pT217_F_z)
    amyloid_match <- nearest_biomarker(amyloid_prepped, "biomarker_date", analysis_df, "amyloid_centiloids_z") %>%
      transmute(RID, amyloid_centiloids, amyloid_centiloids_z)
    out <- out %>% left_join(plasma_match, by = "RID") %>% left_join(amyloid_match, by = "RID")
  }

  if (biomarker_layer == "csf_elecsys") {
    csf_match <- nearest_biomarker(csf_prepped, "biomarker_date", analysis_df, c("CSF_PTAU_z", "CSF_ABETA42_z")) %>%
      transmute(RID, CSF_PTAU, CSF_PTAU_z, CSF_ABETA42, CSF_ABETA42_z)
    out <- out %>% left_join(csf_match, by = "RID")
  }

  if (biomarker_layer == "harmonized_csf") {
    harm_match <- nearest_biomarker(harmonized_csf_prepped, "biomarker_date", analysis_df, c("PHC_pTau_z", "PHC_AB42_z")) %>%
      transmute(RID, PHC_pTau, PHC_pTau_z, PHC_AB42, PHC_AB42_z)
    out <- out %>% left_join(harm_match, by = "RID")
  }

  out %>% mutate(biomarker_layer = biomarker_layer)
}

exposure_specs <- tribble(
  ~exposure_id, ~exposure_label, ~source_role, ~window_mode, ~long_obj,
  "PHC_STRICT", "ADSP_PHC strict global", "comparator", "flex", "phc_global_long",
  "PHC_MEM", "ADSP_PHC memory", "primary", "flex", "phc_mem_long",
  "UW_GLOBAL", "UW global mean", "comparator", "flex", "uw_global_long",
  "UW_MEM", "UW memory", "primary", "flex", "uw_mem_long",
  "EMBIC_DCB", "EMBIC DCB 10-item mean", "primary_backup", "flex", "embic_dcb_long",
  "EMBIC_QCP", "EMBIC qCP 10-item mean", "primary", "flex", "embic_qcp_long",
  "PHC_STRICT", "ADSP_PHC strict global", "sensitivity", "exact", "phc_global_long",
  "PHC_MEM", "ADSP_PHC memory", "sensitivity", "exact", "phc_mem_long",
  "UW_MEM", "UW memory", "sensitivity", "exact", "uw_mem_long",
  "EMBIC_QCP", "EMBIC qCP 10-item mean", "sensitivity", "exact", "embic_qcp_long"
)

long_lookup <- list(
  phc_global_long = phc_global_long,
  phc_mem_long = phc_mem_long,
  uw_global_long = uw_global_long,
  uw_mem_long = uw_mem_long,
  embic_dcb_long = embic_dcb_long,
  embic_qcp_long = embic_qcp_long
)

exposure_objects <- lapply(seq_len(nrow(exposure_specs)), function(i) {
  spec <- exposure_specs[i, ]
  compute_pe_subjects(
    long_lookup[[spec$long_obj[[1]]]],
    exposure_id = spec$exposure_id[[1]],
    exposure_label = spec$exposure_label[[1]],
    source_role = spec$source_role[[1]],
    window_mode = spec$window_mode[[1]]
  )
})

exposure_subject_level <- bind_rows(exposure_objects)

window_gain_summary <- exposure_subject_level %>%
  count(exposure_id, exposure_label, window_mode, name = "n_subjects") %>%
  tidyr::pivot_wider(names_from = window_mode, values_from = n_subjects) %>%
  mutate(gain_flex_vs_exact = flex - exact)

cohort_rows <- list()
model_rows <- list()
clinical_rows <- list()
longitudinal_rows <- list()
validation_draw_rows <- list()
validation_fold_rows <- list()
validation_summary_rows <- list()
validation_input_list <- list()
ols_sensitivity_rows <- list()

base_rhs <- combine_rhs("Age_z", "Sex_F", "Education_z", "APOE4", "baseline_score_z")
layer_rhs <- list(
  base = base_rhs,
  csf_elecsys = combine_rhs(base_rhs, "CSF_PTAU_z", "CSF_ABETA42_z"),
  harmonized_csf = combine_rhs(base_rhs, "PHC_pTau_z", "PHC_AB42_z"),
  main_biomarker = combine_rhs(base_rhs, "pT217_F_z", "amyloid_centiloids_z")
)
layer_vars <- list(
  base = character(0),
  csf_elecsys = c("CSF_PTAU_z", "CSF_ABETA42_z"),
  harmonized_csf = c("PHC_pTau_z", "PHC_AB42_z"),
  main_biomarker = c("pT217_F_z", "amyloid_centiloids_z")
)

for (subject_df in exposure_objects) {
  early_long_here <- attr(subject_df, "early_long")

  subject_df <- subject_df %>%
    left_join(subject_covariates, by = "RID") %>%
    mutate(APOE4 = as.numeric(APOE4))

  cohort_rows[[length(cohort_rows) + 1]] <- subject_df %>%
    summarise(
      exposure_id = first(exposure_id),
      exposure_label = first(exposure_label),
      source_role = first(source_role),
      window_mode = first(window_mode),
      outcome_id = "PE_subjects",
      biomarker_layer = "none",
      n_subjects = n(),
      n_events = NA_real_
    )

  outcomes <- list(
    build_dxsum_outcome(subject_df, dxsum),
    build_cdr_single_outcome(subject_df, cdr),
    build_cdr_confirmed_outcome(subject_df, cdr)
  )

  for (outcome_df in outcomes) {
    outcome_id <- first(outcome_df$outcome_id)
    if (nrow(outcome_df) == 0) {
      next
    }

    cohort_rows[[length(cohort_rows) + 1]] <- outcome_df %>%
      summarise(
        exposure_id = first(exposure_id),
        exposure_label = first(exposure_label),
        source_role = first(source_role),
        window_mode = first(window_mode),
        outcome_id = first(outcome_id),
        biomarker_layer = "screening",
        n_subjects = n(),
        n_events = sum(event, na.rm = TRUE)
      )

    for (layer in c("base", "csf_elecsys", "harmonized_csf", "main_biomarker")) {
      model_df <- add_biomarker_layer(outcome_df, layer)
      model_df <- filter_complete_cases(
        model_df,
        c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z", "APOE4", "baseline_score_z", "pe_resid_z_decr", layer_vars[[layer]])
      )

      status_here <- if (nrow(model_df) < clinical_n_min || sum(model_df$event, na.rm = TRUE) < clinical_events_min) {
        "below_threshold"
      } else {
        "ready"
      }

      fit_cont <- NULL
      fit_group <- NULL
      if (status_here == "ready") {
        fit_cont <- fit_cox_model(model_df, combine_rhs("pe_resid_z_decr", layer_rhs[[layer]]))
        fit_group <- fit_cox_model(model_df, combine_rhs("PE_group", layer_rhs[[layer]]))
        status_here <- if (is.null(fit_cont)) "fit_failed" else "ran"
      }

      model_rows[[length(model_rows) + 1]] <- tibble(
        exposure_id = first(subject_df$exposure_id),
        exposure_label = first(subject_df$exposure_label),
        source_role = first(subject_df$source_role),
        window_mode = first(subject_df$window_mode),
        outcome_id = outcome_id,
        biomarker_layer = layer,
        status = status_here,
        n_subjects = nrow(model_df),
        n_events = sum(model_df$event, na.rm = TRUE)
      )

      clinical_rows[[length(clinical_rows) + 1]] <- collect_cox_term(
        fit_cont,
        first(subject_df$exposure_id),
        first(subject_df$window_mode),
        outcome_id,
        layer,
        "continuous",
        "pe_resid_z_decr",
        nrow(model_df),
        sum(model_df$event, na.rm = TRUE)
      ) %>%
        mutate(exposure_label = first(subject_df$exposure_label), source_role = first(subject_df$source_role))

      clinical_rows[[length(clinical_rows) + 1]] <- collect_cox_term(
        fit_group,
        first(subject_df$exposure_id),
        first(subject_df$window_mode),
        outcome_id,
        layer,
        "grouped",
        "PE_groupLow_PE",
        nrow(model_df),
        sum(model_df$event, na.rm = TRUE)
      ) %>%
        mutate(exposure_label = first(subject_df$exposure_label), source_role = first(subject_df$source_role))

      if (
        first(subject_df$exposure_id) == "PHC_MEM" &&
          first(subject_df$window_mode) == "flex" &&
          outcome_id == "DXSUM_23" &&
          layer %in% c("base", "csf_elecsys", "main_biomarker")
      ) {
        ols_model_df <- filter_complete_cases(
          model_df,
          c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z", "APOE4", "baseline_score_z", "ols_pe_resid_z_decr", "PE_group_ols", layer_vars[[layer]])
        )
        ols_fit_cont <- NULL
        ols_fit_group <- NULL
        if (nrow(ols_model_df) >= clinical_n_min && sum(ols_model_df$event, na.rm = TRUE) >= clinical_events_min) {
          ols_fit_cont <- fit_cox_model(ols_model_df, combine_rhs("ols_pe_resid_z_decr", layer_rhs[[layer]]))
          ols_fit_group <- fit_cox_model(ols_model_df, combine_rhs("PE_group_ols", layer_rhs[[layer]]))
        }
        ols_sensitivity_rows[[length(ols_sensitivity_rows) + 1]] <- collect_cox_term(
          ols_fit_cont,
          first(subject_df$exposure_id),
          first(subject_df$window_mode),
          outcome_id,
          layer,
          "continuous_ols",
          "ols_pe_resid_z_decr",
          nrow(ols_model_df),
          sum(ols_model_df$event, na.rm = TRUE)
        ) %>%
          mutate(exposure_label = first(subject_df$exposure_label), source_role = first(subject_df$source_role))
        ols_sensitivity_rows[[length(ols_sensitivity_rows) + 1]] <- collect_cox_term(
          ols_fit_group,
          first(subject_df$exposure_id),
          first(subject_df$window_mode),
          outcome_id,
          layer,
          "grouped_ols",
          "PE_group_olsLow_PE",
          nrow(ols_model_df),
          sum(ols_model_df$event, na.rm = TRUE)
        ) %>%
          mutate(exposure_label = first(subject_df$exposure_label), source_role = first(subject_df$source_role))
      }

      if (
        first(subject_df$window_mode) == "flex" &&
          first(subject_df$exposure_id) %in% c("PHC_MEM", "UW_MEM", "EMBIC_DCB", "EMBIC_QCP") &&
          outcome_id == "DXSUM_23" &&
          layer %in% c("csf_elecsys", "main_biomarker") &&
          status_here == "ran"
      ) {
        validation_input_list[[paste(first(subject_df$exposure_id), outcome_id, layer, sep = "_")]] <- list(
          exposure_id = first(subject_df$exposure_id),
          exposure_label = first(subject_df$exposure_label),
          window_mode = first(subject_df$window_mode),
          biomarker_layer = layer,
          outcome_id = outcome_id,
          data = model_df,
          base_rhs = layer_rhs[[layer]],
          early_long = early_long_here
        )
      }
    }
  }

  if (first(subject_df$window_mode) == "flex") {
    for (target_id in c("ADNI_MEM", "UW_GLOBAL")) {
      long_df <- build_longitudinal_post_anchor(subject_df, uw_longitudinal, target_id, target_id) %>%
        left_join(subject_covariates, by = "RID")

      for (layer in c("base", "csf_elecsys")) {
        long_model_df <- if (layer == "base") {
          long_df
        } else {
          subject_level <- long_df %>%
            group_by(RID) %>%
            summarise(anchor_date = first(anchor_date), analysis_end_date = max(EXAMDATE), .groups = "drop") %>%
            add_biomarker_layer("csf_elecsys")

          long_df %>%
            left_join(subject_level %>% select(RID, CSF_PTAU_z, CSF_ABETA42_z), by = "RID")
        }

        long_model_df <- filter_complete_cases(
          long_model_df,
          c("RID", "outcome_value", "time_years", "Age_z", "Sex_F", "Education_z", "APOE4", "baseline_score_z", "pe_resid_z_decr", if (layer == "csf_elecsys") c("CSF_PTAU_z", "CSF_ABETA42_z") else character(0))
        )

        n_subjects_here <- n_distinct(long_model_df$RID)
        n_rows_here <- nrow(long_model_df)
        n_repeat_subjects_here <- long_model_df %>%
          count(RID) %>%
          summarise(n_repeat = sum(n >= 2), .groups = "drop") %>%
          pull(n_repeat)
        if (length(n_repeat_subjects_here) == 0) {
          n_repeat_subjects_here <- 0
        }

        status_here <- if (
          n_subjects_here < lmm_subjects_min ||
            n_rows_here < lmm_rows_min ||
            n_repeat_subjects_here < lmm_repeat_subjects_min
        ) {
          "below_threshold"
        } else {
          "ready"
        }

        rhs_here <- if (layer == "base") {
          combine_rhs("time_years * pe_resid_z_decr", base_rhs)
        } else {
          combine_rhs("time_years * pe_resid_z_decr", layer_rhs[[layer]])
        }

        fit_lmm <- NULL
        if (status_here == "ready") {
          fit_lmm <- fit_lmm_model(long_model_df, rhs_here)
          status_here <- if (is.null(fit_lmm)) "fit_failed" else "ran"
        }

        model_rows[[length(model_rows) + 1]] <- tibble(
          exposure_id = first(subject_df$exposure_id),
          exposure_label = first(subject_df$exposure_label),
          source_role = first(subject_df$source_role),
          window_mode = first(subject_df$window_mode),
          outcome_id = target_id,
          biomarker_layer = paste0("LMM_", layer),
          status = status_here,
          n_subjects = n_subjects_here,
          n_events = NA_real_
        )

        longitudinal_rows[[length(longitudinal_rows) + 1]] <- collect_lmm_term(
          fit_lmm,
          first(subject_df$exposure_id),
          first(subject_df$window_mode),
          target_id,
          layer,
          "time_years:pe_resid_z_decr",
          n_subjects_here,
          n_rows_here,
          n_repeat_subjects_here
        ) %>%
          mutate(exposure_label = first(subject_df$exposure_label), source_role = first(subject_df$source_role))
      }
    }
  }
}

for (nm in names(validation_input_list)) {
  item <- validation_input_list[[nm]]
  boot <- bootstrap_delta_c(item$data, item$base_rhs, combine_rhs("pe_resid_z_decr", item$base_rhs), bootstrap_reps)
  cv <- run_nested_cv_delta_c(item$data, item$base_rhs, cv_k, cv_repeats, early_long = item$early_long)

  validation_draw_rows[[length(validation_draw_rows) + 1]] <- tibble(
    exposure_id = item$exposure_id,
    outcome_id = item$outcome_id,
    biomarker_layer = item$biomarker_layer,
    validation_type = "bootstrap",
    draw_id = seq_along(boot),
    delta_c = boot
  )

  validation_summary_rows[[length(validation_summary_rows) + 1]] <- tibble(
    exposure_id = item$exposure_id,
    outcome_id = item$outcome_id,
    biomarker_layer = item$biomarker_layer,
    validation_type = "bootstrap",
    n = sum(is.finite(boot)),
    mean_delta_c = mean(boot, na.rm = TRUE),
    ci_lower_2.5 = quantile(boot, 0.025, na.rm = TRUE),
    ci_upper_97.5 = quantile(boot, 0.975, na.rm = TRUE),
    pct_positive = mean(boot > 0, na.rm = TRUE) * 100
  )

  if (nrow(cv) > 0) {
    validation_fold_rows[[length(validation_fold_rows) + 1]] <- cv %>%
      mutate(exposure_id = item$exposure_id, outcome_id = item$outcome_id, biomarker_layer = item$biomarker_layer, validation_type = "nested_cv", .before = 1)

    validation_summary_rows[[length(validation_summary_rows) + 1]] <- tibble(
      exposure_id = item$exposure_id,
      outcome_id = item$outcome_id,
      biomarker_layer = item$biomarker_layer,
      validation_type = "nested_cv",
      n = nrow(cv),
      mean_delta_c = mean(cv$delta_c, na.rm = TRUE),
      ci_lower_2.5 = quantile(cv$delta_c, 0.025, na.rm = TRUE),
      ci_upper_97.5 = quantile(cv$delta_c, 0.975, na.rm = TRUE),
      pct_positive = mean(cv$delta_c > 0, na.rm = TRUE) * 100
    )
  }
}

cohort_matrix <- bind_rows(cohort_rows)
model_registry <- bind_rows(model_rows)
clinical_results <- bind_rows(clinical_rows)
longitudinal_results <- bind_rows(longitudinal_rows)
validation_bootstrap_draws <- bind_rows(validation_draw_rows)
validation_cv_fold_results <- bind_rows(validation_fold_rows)
validation_summary <- bind_rows(validation_summary_rows)
clinical_results_ols_sensitivity <- bind_rows(ols_sensitivity_rows)

ranking_flex <- exposure_specs %>%
  filter(window_mode == "flex") %>%
  distinct(exposure_id, exposure_label, source_role) %>%
  left_join(
    clinical_results %>%
      filter(window_mode == "flex", outcome_id == "DXSUM_23", biomarker_layer == "csf_elecsys", pe_model == "continuous") %>%
      transmute(exposure_id, dx_csf_hr = estimate, dx_csf_p = p.value, dx_csf_n = n),
    by = "exposure_id"
  ) %>%
  left_join(
    clinical_results %>%
      filter(window_mode == "flex", outcome_id == "CDR_single", biomarker_layer == "csf_elecsys", pe_model == "continuous") %>%
      transmute(exposure_id, cdr_single_csf_hr = estimate, cdr_single_csf_p = p.value, cdr_single_csf_n = n),
    by = "exposure_id"
  ) %>%
  left_join(
    clinical_results %>%
      filter(window_mode == "flex", outcome_id == "DXSUM_23", biomarker_layer == "main_biomarker", pe_model == "continuous") %>%
      transmute(exposure_id, dx_mainbio_hr = estimate, dx_mainbio_p = p.value),
    by = "exposure_id"
  ) %>%
  left_join(
    clinical_results %>%
      filter(window_mode == "exact", outcome_id == "DXSUM_23", biomarker_layer == "csf_elecsys", pe_model == "continuous") %>%
      transmute(exposure_id, exact_dx_csf_hr = estimate, exact_dx_csf_p = p.value),
    by = "exposure_id"
  ) %>%
  left_join(
    longitudinal_results %>%
      filter(window_mode == "flex", target_id == "ADNI_MEM", biomarker_layer == "base") %>%
      transmute(exposure_id, mem_lmm_beta = estimate, mem_lmm_p = p.value),
    by = "exposure_id"
  ) %>%
  left_join(
    validation_summary %>%
      filter(validation_type == "bootstrap", outcome_id == "DXSUM_23", biomarker_layer == "csf_elecsys") %>%
      transmute(exposure_id, bootstrap_mean_delta_c = mean_delta_c, bootstrap_pct_positive = pct_positive),
    by = "exposure_id"
  ) %>%
  left_join(
    validation_summary %>%
      filter(validation_type == "nested_cv", outcome_id == "DXSUM_23", biomarker_layer == "csf_elecsys") %>%
      transmute(exposure_id, cv_mean_delta_c = mean_delta_c, cv_pct_positive = pct_positive),
    by = "exposure_id"
  ) %>%
  mutate(
    dx_csf_supported = !is.na(dx_csf_hr) & dx_csf_hr > 1 & dx_csf_p < 0.05,
    cdr_single_csf_supported = !is.na(cdr_single_csf_hr) & cdr_single_csf_hr > 1 & cdr_single_csf_p < 0.05,
    mainbio_supported = !is.na(dx_mainbio_hr) & dx_mainbio_hr > 1 & dx_mainbio_p < 0.05,
    exact_supported = !is.na(exact_dx_csf_hr) & exact_dx_csf_hr > 1 & exact_dx_csf_p < 0.05,
    downstream_supported = !is.na(mem_lmm_beta) & mem_lmm_beta < 0 & mem_lmm_p < 0.10,
    rank_score = 100 * as.numeric(dx_csf_supported) +
      80 * as.numeric(cdr_single_csf_supported) +
      30 * as.numeric(mainbio_supported) +
      20 * as.numeric(exact_supported) +
      15 * as.numeric(downstream_supported) +
      coalesce(bootstrap_mean_delta_c, 0) * 1000 +
      coalesce(cv_mean_delta_c, 0) * 1000,
    recommended_status = case_when(
      dx_csf_supported & cdr_single_csf_supported & (mainbio_supported | exact_supported) ~ "advance",
      dx_csf_supported & (cdr_single_csf_supported | mainbio_supported | exact_supported) ~ "keep_supplement",
      TRUE ~ "stop"
    )
  ) %>%
  arrange(desc(rank_score), exposure_id) %>%
  mutate(rank = row_number()) %>%
  select(
    rank, exposure_id, exposure_label, source_role,
    dx_csf_hr, dx_csf_p, cdr_single_csf_hr, cdr_single_csf_p,
    dx_mainbio_hr, dx_mainbio_p, exact_dx_csf_hr, exact_dx_csf_p,
    mem_lmm_beta, mem_lmm_p, bootstrap_mean_delta_c, cv_mean_delta_c,
    recommended_status
  )

write_csv(exposure_subject_level, file.path(out_dir, "exposure_subject_level.csv"))
write_csv(window_gain_summary, file.path(out_dir, "window_gain_summary.csv"))
write_csv(cohort_matrix, file.path(out_dir, "cohort_matrix.csv"))
write_csv(model_registry, file.path(out_dir, "model_registry.csv"))
write_csv(clinical_results, file.path(out_dir, "clinical_results.csv"))
write_csv(clinical_results_ols_sensitivity, file.path(out_dir, "clinical_results_ols_sensitivity.csv"))
write_csv(longitudinal_results, file.path(out_dir, "longitudinal_results.csv"))
write_csv(validation_bootstrap_draws, file.path(out_dir, "validation_bootstrap_draws.csv"))
write_csv(validation_cv_fold_results, file.path(out_dir, "validation_cv_fold_results.csv"))
write_csv(validation_summary, file.path(out_dir, "validation_summary.csv"))
write_csv(ranking_flex, file.path(out_dir, "scenario_ranking.csv"))

# ----------------------------------------------------------------------------
# Sensitivity: include anchor-prevalent DXSUM cases for the PHC_MEM / flex
# DXSUM_23 cell only. The strict default (excluding anchor-prevalent) is
# the inferential primary; this block reproduces the older risk-set so the
# directional comparison can be reported as a supplement note. All other
# cells are not re-run.
# ----------------------------------------------------------------------------
sens_subject_df <- exposure_objects[[
  which(exposure_specs$exposure_id == "PHC_MEM" & exposure_specs$window_mode == "flex")[1]
]] %>%
  left_join(subject_covariates, by = "RID") %>%
  mutate(APOE4 = as.numeric(APOE4))

sens_outcome <- build_dxsum_outcome(sens_subject_df, dxsum, exclude_anchor_prevalent = FALSE)

sens_clinical_rows <- list()
for (layer in c("base", "csf_elecsys", "main_biomarker")) {
  m <- add_biomarker_layer(sens_outcome, layer)
  m <- filter_complete_cases(
    m,
    c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z", "APOE4",
      "baseline_score_z", "pe_resid_z_decr", layer_vars[[layer]])
  )
  if (nrow(m) < clinical_n_min || sum(m$event, na.rm = TRUE) < clinical_events_min) next
  fit_cont <- fit_cox_model(m, combine_rhs("pe_resid_z_decr", layer_rhs[[layer]]))
  fit_group <- fit_cox_model(m, combine_rhs("PE_group", layer_rhs[[layer]]))
  sens_clinical_rows[[length(sens_clinical_rows) + 1]] <- collect_cox_term(
    fit_cont, "PHC_MEM", "flex", "DXSUM_23", layer, "continuous",
    "pe_resid_z_decr", nrow(m), sum(m$event, na.rm = TRUE)
  )
  sens_clinical_rows[[length(sens_clinical_rows) + 1]] <- collect_cox_term(
    fit_group, "PHC_MEM", "flex", "DXSUM_23", layer, "grouped",
    "PE_groupLow_PE", nrow(m), sum(m$event, na.rm = TRUE)
  )
}
sens_clinical <- bind_rows(sens_clinical_rows) %>%
  mutate(risk_set = "include_anchor_prevalent")
write_csv(sens_clinical, file.path(out_dir, "clinical_results_inclprevalent_sensitivity.csv"))

# Also record how many participants were excluded by the strict primary rule.
phc_mem_flex_idx <- which(
  exposure_specs$exposure_id == "PHC_MEM" & exposure_specs$window_mode == "flex"
)[1]
sens_excluded <- exposure_objects[[phc_mem_flex_idx]] %>%
  inner_join(dxsum %>% select(RID, EXAMDATE, DIAGNOSIS), by = "RID") %>%
  filter(
    !is.na(EXAMDATE), !is.na(DIAGNOSIS),
    EXAMDATE <= anchor_date,
    DIAGNOSIS %in% c(2, 3)
  ) %>%
  distinct(RID)
write_csv(
  tibble(
    rule = "anchor_prevalent_DXSUM_exclusion_PHC_MEM_flex",
    n_excluded = nrow(sens_excluded)
  ),
  file.path(out_dir, "anchor_prevalent_exclusion_count.csv")
)

capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))
