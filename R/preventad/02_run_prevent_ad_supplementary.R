

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
# PREVENT-AD supplementary validation analyses
# - Event-time distribution for MCI conversion
# - Bootstrap delta-C
# - Nested CV delta-C with fold-specific PE re-estimation
# - Non-nested CV comparison
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(survival)
  library(tibble)
})

out_dir <- file.path(pkg_root, "verify", "preventad")
data_dir <- file.path(data_root, "preventad", "Clinical_data")

parse_visit_months <- function(x) {
  out <- suppressWarnings(as.numeric(str_extract(x, "\\d+")))
  out[!is.na(x) & x == "BL00"] <- 0
  out
}

harrell_c_from_lp <- function(time, status, lp) {
  cc <- survival::concordance(Surv(time, status) ~ lp)
  unname(1 - cc$concordance)
}

summarise_delta_c <- function(x, approach, note = NA_character_) {
  x <- x[is.finite(x)]
  tibble(
    approach = approach,
    note = note,
    n_folds = length(x),
    mean = round(mean(x), 4),
    median = round(median(x), 4),
    sd = round(sd(x), 4),
    ci_lower_2.5 = round(unname(quantile(x, 0.025)), 4),
    ci_upper_97.5 = round(unname(quantile(x, 0.975)), 4),
    pct_positive = round(100 * mean(x > 0), 1)
  )
}

assign_stratified_folds <- function(status, k) {
  fold_id <- integer(length(status))
  for (level in unique(status)) {
    idx <- which(status == level)
    fold_id[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  fold_id
}

add_fold_specific_pe <- function(train, test, slope_col = "raw_slope_ols", baseline_col = "baseline_RBANS_total") {
  resid_fit <- lm(
    as.formula(paste0(slope_col, " ~ ", baseline_col)),
    data = train
  )

  train_resid <- resid(resid_fit)
  train_mean <- mean(train_resid, na.rm = TRUE)
  train_sd <- sd(train_resid, na.rm = TRUE)
  if (!is.finite(train_sd) || train_sd == 0) {
    train_sd <- 1
  }

  train <- train %>%
    mutate(
      pe_cv_z = (train_resid - train_mean) / train_sd,
      PE_group_cv = case_when(
        pe_cv_z > 0.2 ~ "High_PE",
        pe_cv_z < -0.2 ~ "Low_PE",
        TRUE ~ "Moderate_PE"
      )
    )

  test_resid <- test[[slope_col]] - predict(resid_fit, newdata = test)
  test <- test %>%
    mutate(
      pe_cv_z = (test_resid - train_mean) / train_sd,
      PE_group_cv = case_when(
        pe_cv_z > 0.2 ~ "High_PE",
        pe_cv_z < -0.2 ~ "Low_PE",
        TRUE ~ "Moderate_PE"
      )
    )

  train$PE_group_cv <- factor(train$PE_group_cv, levels = c("High_PE", "Moderate_PE", "Low_PE"))
  test$PE_group_cv <- factor(test$PE_group_cv, levels = c("High_PE", "Moderate_PE", "Low_PE"))

  list(train = train, test = test)
}

# ----------------------------
# Load data
# ----------------------------
surv_primary <- read_csv(file.path(out_dir, "survival_dataset_primary.csv"), show_col_types = FALSE)

surv_complete <- surv_primary %>%
  mutate(
    PE_group_lmer = factor(PE_group_lmer, levels = c("High_PE", "Moderate_PE", "Low_PE"))
  ) %>%
  filter(
    !is.na(Age_z),
    !is.na(Education_z),
    !is.na(Sex_F),
    !is.na(APOE4_carrier),
    !is.na(baseline_RBANS_z),
    !is.na(PE_group_lmer)
  )

rbans_early <- read_csv(file.path(data_dir, "RBANS_Registered_PREVENTAD_2024-08-01.csv"), show_col_types = FALSE) %>%
  mutate(
    visit_months = parse_visit_months(Visit_label),
    time_years = visit_months / 12,
    total_scale_index_score = as.numeric(total_scale_index_score)
  ) %>%
  filter(visit_months %in% c(0, 12, 24), !is.na(total_scale_index_score))

raw_slopes <- rbans_early %>%
  group_by(CONP_ID) %>%
  filter(any(visit_months == 0), n_distinct(visit_months) >= 2) %>%
  summarise(
    raw_slope_ols = coef(lm(total_scale_index_score ~ time_years))[["time_years"]],
    .groups = "drop"
  )

cv_dat <- surv_complete %>%
  left_join(raw_slopes, by = "CONP_ID") %>%
  filter(!is.na(raw_slope_ols), !is.na(baseline_RBANS_total))

# ----------------------------
# 1. Event-time distribution
# ----------------------------
events <- surv_primary %>% filter(converted)

event_summary <- tibble(
  statistic = c(
    "N_events", "mean", "sd", "min", "Q25", "median", "Q75", "max",
    "pct_within_1yr", "pct_within_2yr", "pct_after_3yr"
  ),
  value = c(
    nrow(events),
    round(mean(events$time_to_event_years), 3),
    round(sd(events$time_to_event_years), 3),
    round(min(events$time_to_event_years), 3),
    round(quantile(events$time_to_event_years, 0.25), 3),
    round(median(events$time_to_event_years), 3),
    round(quantile(events$time_to_event_years, 0.75), 3),
    round(max(events$time_to_event_years), 3),
    round(100 * mean(events$time_to_event_years <= 1), 1),
    round(100 * mean(events$time_to_event_years <= 2), 1),
    round(100 * mean(events$time_to_event_years > 3), 1)
  )
)
write_csv(event_summary, file.path(out_dir, "MCI_event_time_summary.csv"))

p_event <- ggplot(events, aes(x = time_to_event_years)) +
  geom_histogram(bins = 15, fill = "#B2182B", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(events$time_to_event_years), linetype = "dashed") +
  labs(
    title = "PREVENT-AD MCI event-time distribution",
    subtitle = "Time from end of PE window (FU24)",
    x = "Years to incident MCI",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "Supp_MCI_event_time_distribution.pdf"), p_event, width = 7, height = 5)
ggsave(file.path(out_dir, "Supp_MCI_event_time_distribution.png"), p_event, width = 7, height = 5, dpi = 300)

# ----------------------------
# 2. Bootstrap delta-C
# ----------------------------
bootstrap_delta_c <- function(data, n_boot = 1000) {
  out <- rep(NA_real_, n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(nrow(data), size = nrow(data), replace = TRUE)
    boot_dat <- data[idx, , drop = FALSE]

    m_base <- tryCatch(
      coxph(Surv(time_to_event_years, converted) ~ Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z, data = boot_dat),
      error = function(e) NULL
    )
    m_pe <- tryCatch(
      coxph(Surv(time_to_event_years, converted) ~ PE_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z, data = boot_dat),
      error = function(e) NULL
    )

    if (!is.null(m_base) && !is.null(m_pe)) {
      out[i] <- summary(m_pe)$concordance[1] - summary(m_base)$concordance[1]
    }
  }
  out
}

set.seed(42)
boot_dc <- bootstrap_delta_c(surv_complete, 1000)
boot_summary <- tibble(
  statistic = c("mean", "median", "ci_lower_2.5", "ci_upper_97.5", "pct_positive", "n_boot"),
  value = c(
    round(mean(boot_dc, na.rm = TRUE), 4),
    round(median(boot_dc, na.rm = TRUE), 4),
    round(quantile(boot_dc, 0.025, na.rm = TRUE), 4),
    round(quantile(boot_dc, 0.975, na.rm = TRUE), 4),
    round(100 * mean(boot_dc > 0, na.rm = TRUE), 1),
    sum(is.finite(boot_dc))
  )
)
write_csv(boot_summary, file.path(out_dir, "bootstrap_delta_C_summary.csv"))

p_boot <- ggplot(tibble(delta_c = boot_dc), aes(x = delta_c)) +
  geom_histogram(bins = 40, fill = "#4393C3", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Bootstrap distribution of delta-C",
    x = "Delta C-index (PE model - base model)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "Supp_bootstrap_deltaC.pdf"), p_boot, width = 7, height = 5)
ggsave(file.path(out_dir, "Supp_bootstrap_deltaC.png"), p_boot, width = 7, height = 5, dpi = 300)

# ----------------------------
# 3. Nested and non-nested CV delta-C
# ----------------------------
run_cv_delta_c <- function(data, nested = TRUE, k = 5, repeats = 10) {
  rows <- list()
  idx_row <- 1

  for (rep_i in seq_len(repeats)) {
    fold_id <- assign_stratified_folds(data$converted, k)
    for (fold_i in seq_len(k)) {
      train <- data[fold_id != fold_i, , drop = FALSE]
      test <- data[fold_id == fold_i, , drop = FALSE]

      if (nested) {
        fold_dat <- add_fold_specific_pe(train, test)
        train <- fold_dat$train
        test <- fold_dat$test
        pe_term <- "PE_group_cv"
      } else {
        train <- train %>% mutate(PE_group_cv = PE_group_lmer)
        test <- test %>% mutate(PE_group_cv = PE_group_lmer)
        pe_term <- "PE_group_cv"
      }

      if (sum(test$converted) == 0) {
        next
      }

      m_base <- tryCatch(
        coxph(Surv(time_to_event_years, converted) ~ Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z, data = train),
        error = function(e) NULL
      )
      m_pe <- tryCatch(
        coxph(
          as.formula(paste0("Surv(time_to_event_years, converted) ~ ", pe_term, " + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_RBANS_z")),
          data = train
        ),
        error = function(e) NULL
      )

      if (is.null(m_base) || is.null(m_pe)) {
        next
      }

      lp_base <- predict(m_base, newdata = test, type = "lp")
      lp_pe <- predict(m_pe, newdata = test, type = "lp")

      c_base <- harrell_c_from_lp(test$time_to_event_years, test$converted, lp_base)
      c_pe <- harrell_c_from_lp(test$time_to_event_years, test$converted, lp_pe)

      rows[[idx_row]] <- tibble(
        replicate = rep_i,
        fold = fold_i,
        approach = ifelse(nested, "nested", "non_nested"),
        c_base = c_base,
        c_pe = c_pe,
        delta_c = c_pe - c_base,
        n_train = nrow(train),
        n_test = nrow(test),
        n_test_events = sum(test$converted)
      )
      idx_row <- idx_row + 1
    }
  }

  bind_rows(rows)
}

set.seed(123)
cv_nested <- run_cv_delta_c(cv_dat, nested = TRUE, k = 5, repeats = 10)
cv_non_nested <- run_cv_delta_c(cv_dat, nested = FALSE, k = 5, repeats = 10)

write_csv(cv_nested, file.path(out_dir, "cv_delta_C_fold_results.csv"))
write_csv(cv_non_nested, file.path(out_dir, "cv_delta_C_non_nested_fold_results.csv"))

cv_nested_summary <- summarise_delta_c(
  cv_nested$delta_c,
  approach = "nested",
  note = "Fold-specific PE re-estimation used subject-level early RBANS slopes residualized within the training fold."
)
cv_non_nested_summary <- summarise_delta_c(
  cv_non_nested$delta_c,
  approach = "non_nested",
  note = "Comparison analysis using the saved full-sample PE grouping."
)

write_csv(cv_nested_summary, file.path(out_dir, "cv_delta_C_summary.csv"))
write_csv(cv_non_nested_summary, file.path(out_dir, "cv_delta_C_non_nested_summary.csv"))
write_csv(bind_rows(cv_nested_summary, cv_non_nested_summary), file.path(out_dir, "cv_delta_C_comparison.csv"))

# ----------------------------
# 4. Adjustment-set comparison
# ----------------------------
adjustment_comparison <- tibble(
  Covariate = c(
    "Age", "Sex", "Education", "APOE4", "Baseline cognition",
    "Amyloid PET SUVR", "Plasma p-tau217", "Baseline MoCA", "Censoring definition", "Internal validation"
  ),
  A4 = c(
    "Yes (z)", "Yes", "Yes (z)", "Yes",
    "Baseline Cogstate", "Yes", "Yes", "No",
    "Confirmed CDR follow-up", "Bootstrap + CV"
  ),
  PREVENT_AD = c(
    "Yes (z)", "Yes", "Yes (z)", "Yes",
    "Baseline RBANS", "PET subset only", "Not available", "Sensitivity analysis only",
    "RBANS/contact proxy after FU24", "Bootstrap + nested CV"
  )
)
write_csv(adjustment_comparison, file.path(out_dir, "adjustment_set_comparison.csv"))

cat("PREVENT-AD supplementary analyses complete. Outputs written to:", out_dir, "\n")
