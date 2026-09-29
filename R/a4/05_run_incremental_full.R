#!/usr/bin/env Rscript


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


# PE manuscript incremental full supplementary analyses
# Date: 2026-05-18
# Random seed: 42
#
# Contract:
# - Do not modify manuscript, references, or locked main tables/figures.
# - Write all analysis outputs inside incremental_full_20260518 plus the two
#   requested root-level markdown summaries one directory above.
# - Do not write subject-level identifiers to output files.

set.seed(42)

user_lib <- Sys.getenv("R_LIBS_USER")
.libPaths(unique(c("/usr/lib/R/site-library", "/usr/local/lib/R/site-library",
                   user_lib, .Library, .libPaths())))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(survival)
  library(nlme)
  library(broom)
  library(ggplot2)
})

optional_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = digits))
}

fmt_p <- function(p) {
  if (length(p) == 0 || is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  formatC(p, format = "f", digits = 3)
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

pdate <- function(x) as.Date(x)

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) NA else x[[1]]
}

find_root <- function(start) {
  pkg_root
}

args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, mustWork = TRUE)) else getwd()
out_dir <- normalizePath(script_dir, mustWork = TRUE)
root <- find_root(out_dir)
out_dir <- file.path(root, "verify", "a4", "incremental_full_20260518")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
supp_root <- file.path(root, "verify", "a4")

tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
rep_dir <- file.path(out_dir, "reports")
log_dir <- file.path(out_dir, "logs")
invisible(lapply(c(tab_dir, fig_dir, rep_dir, log_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

sink(file.path(log_dir, "sessionInfo.txt"))
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Root:", root, "\n")
cat("Output:", out_dir, "\n")
cat(".libPaths():\n")
print(.libPaths())
print(sessionInfo())
sink()

csv_write <- function(x, file, subdir = tab_dir) {
  readr::write_csv(x, file.path(subdir, file), na = "")
  file.path(subdir, file)
}

save_plot <- function(plot, basename, width = 7, height = 5) {
  png_path <- file.path(fig_dir, paste0(basename, ".png"))
  pdf_path <- file.path(fig_dir, paste0(basename, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300)
  ggsave(pdf_path, plot, width = width, height = height)
  c(png_path, pdf_path)
}

status_log <- tibble(
  module = character(), item = character(), status = character(),
  n = numeric(), events = numeric(), missingness = character(),
  events_per_parameter = numeric(), diagnostics = character(),
  detail = character(), outputs = character()
)

record_status <- function(module, item, status, detail = "", outputs = character(),
                          n = NA_real_, events = NA_real_, missingness = "",
                          epp = NA_real_, diagnostics = "") {
  status_log <<- bind_rows(status_log, tibble(
    module = module, item = item, status = status, n = n, events = events,
    missingness = missingness, events_per_parameter = epp,
    diagnostics = diagnostics, detail = detail,
    outputs = paste(gsub(paste0("^", out_dir, "/?"), "", outputs), collapse = "; ")
  ))
}

model_missingness <- function(data, vars) {
  miss <- sapply(vars, function(v) sum(is.na(data[[v]])))
  paste(sprintf("%s=%s/%s", names(miss), miss, nrow(data)), collapse = "; ")
}

cox_extract <- function(fit, model, terms = NULL, compare_fit = NULL, n = NA, events = NA) {
  tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  if (is.null(terms)) terms <- tt$term
  lrt <- NULL
  if (!is.null(compare_fit)) lrt <- tryCatch(anova(compare_fit, fit, test = "LRT"), error = function(e) NULL)
  out <- lapply(terms, function(term_i) {
    row <- tt %>% filter(term == term_i)
    if (!nrow(row)) {
      tibble(model = model, term = term_i, n = n, events = events,
             HR = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_)
    } else {
      tibble(model = model, term = term_i, n = n, events = events,
             HR = row$estimate[[1]], conf.low = row$conf.low[[1]],
             conf.high = row$conf.high[[1]], p.value = row$p.value[[1]])
    }
  }) %>% bind_rows()
  out$lrt_chisq_vs_compare <- if (!is.null(lrt)) lrt$Chisq[2] else NA_real_
  out$lrt_df_vs_compare <- if (!is.null(lrt)) lrt$Df[2] else NA_real_
  out$lrt_p_vs_compare <- if (!is.null(lrt)) lrt$`Pr(>|Chi|)`[2] else NA_real_
  out$c_index <- tryCatch(as.numeric(summary(fit)$concordance[1]), error = function(e) NA_real_)
  out
}

cox_ph_global <- function(fit) {
  z <- tryCatch(cox.zph(fit), error = function(e) NULL)
  if (is.null(z)) return("PH check failed")
  p <- z$table[nrow(z$table), "p"]
  paste0("cox.zph global p=", fmt_p(p))
}

predict_cox_risk <- function(fit, newdata, horizon) {
  bh <- survival::basehaz(fit, centered = FALSE)
  if (!nrow(bh)) return(rep(NA_real_, nrow(newdata)))
  h0 <- approx(bh$time, bh$hazard, xout = horizon, method = "constant", f = 0, rule = 2)$y
  lp <- as.numeric(predict(fit, newdata = newdata, type = "lp", reference = "zero"))
  pmin(pmax(1 - exp(-h0 * exp(lp)), 1e-6), 1 - 1e-6)
}

ipcw_weights <- function(d, horizon, time_var = "time_to_event_years", event_var = "event") {
  censor_event <- as.integer(d[[event_var]] == 0)
  gfit <- survfit(Surv(d[[time_var]], censor_event) ~ 1)
  g_at <- function(t) pmax(summary(gfit, times = t, extend = TRUE)$surv, 1e-4)
  g_h <- g_at(horizon)
  case <- d[[event_var]] == 1 & d[[time_var]] <= horizon
  control <- d[[time_var]] > horizon
  w <- rep(0, nrow(d))
  if (any(case)) w[case] <- 1 / g_at(d[[time_var]][case])
  if (any(control)) w[control] <- 1 / g_h
  y <- rep(NA_real_, nrow(d))
  y[case] <- 1
  y[control] <- 0
  tibble(y = y, weight = w, case = case, control = control)
}

weighted_auc <- function(risk, y, w) {
  keep <- is.finite(risk) & !is.na(y) & w > 0
  risk <- risk[keep]; y <- y[keep]; w <- w[keep]
  if (sum(y == 1) < 2 || sum(y == 0) < 2) return(NA_real_)
  r1 <- risk[y == 1]; w1 <- w[y == 1]
  r0 <- risk[y == 0]; w0 <- w[y == 0]
  num <- 0
  den <- sum(w1) * sum(w0)
  for (i in seq_along(r1)) {
    num <- num + w1[i] * sum(w0 * ((r1[i] > r0) + 0.5 * (r1[i] == r0)))
  }
  num / den
}

fixed_horizon_metrics <- function(d, risk, horizon, label) {
  iw <- ipcw_weights(d, horizon)
  keep <- iw$weight > 0 & !is.na(iw$y) & is.finite(risk)
  brier <- sum(iw$weight[keep] * (iw$y[keep] - risk[keep])^2, na.rm = TRUE) / nrow(d)
  tibble(
    model = label, horizon_years = horizon, n = nrow(d),
    events_by_horizon = sum(iw$case), controls_at_horizon = sum(iw$control),
    censored_before_horizon = sum(is.na(iw$y)),
    auc_ipcw = weighted_auc(risk, iw$y, iw$weight),
    brier_ipcw = brier
  )
}

idi_nri_metrics <- function(d, risk_base, risk_add, horizon) {
  iw <- ipcw_weights(d, horizon)
  delta <- risk_add - risk_base
  case <- iw$case & iw$weight > 0
  control <- iw$control & iw$weight > 0
  event_delta <- if (any(case)) weighted.mean(delta[case], iw$weight[case]) else NA_real_
  nonevent_delta <- if (any(control)) weighted.mean(delta[control], iw$weight[control]) else NA_real_
  nri_event <- if (any(case)) {
    (sum(iw$weight[case] * (delta[case] > 1e-8)) -
       sum(iw$weight[case] * (delta[case] < -1e-8))) / sum(iw$weight[case])
  } else NA_real_
  nri_nonevent <- if (any(control)) {
    (sum(iw$weight[control] * (delta[control] < -1e-8)) -
       sum(iw$weight[control] * (delta[control] > 1e-8))) / sum(iw$weight[control])
  } else NA_real_
  b_base <- sum(iw$weight * (iw$y - risk_base)^2, na.rm = TRUE) / nrow(d)
  b_add <- sum(iw$weight * (iw$y - risk_add)^2, na.rm = TRUE) / nrow(d)
  tibble(
    horizon_years = horizon, n = nrow(d), events_by_horizon = sum(case),
    controls_at_horizon = sum(control), censored_before_horizon = sum(is.na(iw$y)),
    IDI = event_delta - nonevent_delta,
    IDI_event_component = event_delta,
    IDI_nonevent_component = -nonevent_delta,
    continuous_NRI = nri_event + nri_nonevent,
    NRI_event_component = nri_event,
    NRI_nonevent_component = nri_nonevent,
    Brier_base = b_base, Brier_plus = b_add, delta_Brier = b_add - b_base
  )
}

calibration_slope_intercept <- function(d, risk, horizon, label) {
  iw <- ipcw_weights(d, horizon)
  keep <- iw$weight > 0 & !is.na(iw$y) & is.finite(risk)
  if (sum(keep) < 30 || length(unique(iw$y[keep])) < 2) {
    return(tibble(model = label, horizon_years = horizon,
                  calibration_intercept = NA_real_, calibration_slope = NA_real_))
  }
  lp <- qlogis(pmin(pmax(risk[keep], 1e-6), 1 - 1e-6))
  intercept_fit <- suppressWarnings(glm(iw$y[keep] ~ offset(lp), weights = iw$weight[keep], family = binomial()))
  slope_fit <- suppressWarnings(glm(iw$y[keep] ~ lp, weights = iw$weight[keep], family = binomial()))
  tibble(model = label, horizon_years = horizon,
         calibration_intercept = unname(coef(intercept_fit)[1]),
         calibration_slope = unname(coef(slope_fit)[2]))
}

km_risk <- function(d, horizon, time_var = "time_to_event_years", event_var = "event") {
  if (nrow(d) < 5 || sum(d[[event_var]], na.rm = TRUE) < 1) return(NA_real_)
  sf <- survfit(as.formula(paste0("Surv(", time_var, ",", event_var, ") ~ 1")), data = d)
  1 - summary(sf, times = horizon, extend = TRUE)$surv[[1]]
}

nearest_by_rid <- function(source_df, analysis_df, date_col, required_cols, prefix, max_abs_days = Inf) {
  source_df %>%
    mutate(scan_date = pdate(.data[[date_col]])) %>%
    filter(!is.na(scan_date)) %>%
    filter(if_all(all_of(required_cols), ~ !is.na(suppressWarnings(as.numeric(.))))) %>%
    inner_join(analysis_df %>% select(RID, anchor_date), by = "RID") %>%
    mutate(diff_days = as.numeric(scan_date - anchor_date),
           abs_diff_days = abs(diff_days)) %>%
    filter(abs_diff_days <= max_abs_days) %>%
    arrange(RID, abs_diff_days, diff_days) %>%
    group_by(RID) %>%
    slice(1) %>%
    ungroup() %>%
    rename(
      !!paste0(prefix, "_date") := scan_date,
      !!paste0(prefix, "_diff_days") := diff_days,
      !!paste0(prefix, "_abs_diff_days") := abs_diff_days
    )
}

find_clinical_dir <- function(root) {
  file.path(data_root, "a4")
}

cat("Root:", root, "\n")
cat("Output:", out_dir, "\n")

# -------------------------------------------------------------------------
# Load core cohorts
# -------------------------------------------------------------------------

a4_path <- file.path(root, "outputs", "locked", "a4",
                     "surv_primary_complete_w1_15_strict_preW48_lmer.rds")
a4 <- readRDS(a4_path) %>%
  mutate(
    event = as.integer(converted),
    pe_worse_z = learning_slope_resid_z_decr_lmer,
    sex_f_num = as.numeric(as.character(Sex_F)),
    sex_f_num = ifelse(is.na(sex_f_num), ifelse(as.character(Sex_F) %in% c("Female", "F", "1"), 1, 0), sex_f_num),
    apoe4_num = as.numeric(as.character(APOE4_carrier)),
    apoe4_num = ifelse(is.na(apoe4_num), ifelse(as.character(APOE4_carrier) %in% c("Carrier", "1", "TRUE"), 1, 0), apoe4_num),
    Sex_F = factor(ifelse(sex_f_num == 1, "Female", "Male")),
    APOE4_carrier = factor(ifelse(apoe4_num == 1, "Carrier", "Non-carrier"))
  )

prevent_path <- file.path(root, "outputs", "locked", "preventad", "survival_dataset_primary.csv")
prevent <- read_csv(prevent_path, show_col_types = FALSE) %>%
  mutate(
    event = as.integer(converted),
    pe_worse_z = pe_resid_z_decr_lmer,
    sex_f_num = ifelse(Sex == "Female", 1, 0),
    apoe4_num = as.integer(APOE4 == 1),
    amyloid_z_h = safe_scale(amyloid_SUVR),
    baseline_cog_z_h = safe_scale(baseline_RBANS_total),
    age_z_h = safe_scale(baseline_age_years),
    education_z_h = safe_scale(Education)
  )

clinical_dir <- find_clinical_dir(root)
if (!is.na(clinical_dir)) {
  cog_raw <- read_csv(file.path(clinical_dir, "External Data", "cogstate_battery.csv"),
                      show_col_types = FALSE, progress = FALSE) %>%
    mutate(VISCODE_num = suppressWarnings(as.numeric(as.character(VISCODE)))) %>%
    filter(BID %in% a4$BID, !is.na(TDate_DAYS_T0)) %>%
    group_by(BID, VISCODE_num) %>%
    summarise(
      C3Comp = mean(C3Comp, na.rm = TRUE),
      PsychAttnComp = mean(PsychAttnComp, na.rm = TRUE),
      AttnComp = mean(AttnComp, na.rm = TRUE),
      LearnWMComp = mean(LearnWMComp, na.rm = TRUE),
      TDate_DAYS_T0 = mean(TDate_DAYS_T0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(a4 %>% select(BID, anchor_day), by = "BID") %>%
    mutate(time_from_first_years = NA_real_,
           time_from_anchor_years = (TDate_DAYS_T0 - anchor_day) / 365.25)
} else {
  cog_raw <- tibble()
}

# -------------------------------------------------------------------------
# A. Methodological strengthening
# -------------------------------------------------------------------------

base_vars <- c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z",
               "APOE4_carrier", "SUVR_z", "pTau217_z", "baseline_C3_z")
base_form <- Surv(time_to_event_years, event) ~ Age_z + Sex_F + Education_z +
  APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z
pe_form <- update(base_form, . ~ . + pe_worse_z)

# A1 Head-to-head PE vs simple slope
a1_out <- tryCatch({
  early_slope <- cog_raw %>%
    filter(VISCODE_num %in% c(1, 3, 9, 15), !is.na(C3Comp)) %>%
    group_by(BID) %>%
    mutate(time_from_first_years = (TDate_DAYS_T0 - min(TDate_DAYS_T0, na.rm = TRUE)) / 365.25) %>%
    filter(n() >= 2) %>%
    summarise(simple_ols_slope = unname(coef(lm(C3Comp ~ time_from_first_years))[["time_from_first_years"]]),
              .groups = "drop")
  a4_h <- a4 %>% left_join(early_slope, by = "BID") %>%
    mutate(simple_slope_worse_z = -safe_scale(simple_ols_slope))
  head_vars <- c(base_vars, "pe_worse_z", "simple_slope_worse_z")
  d <- a4_h %>% filter(if_all(all_of(head_vars), ~ !is.na(.)))
  fit_base <- coxph(base_form, data = d, x = TRUE, model = TRUE)
  fit_simple <- coxph(update(base_form, . ~ . + simple_slope_worse_z), data = d, x = TRUE, model = TRUE)
  fit_pe <- coxph(pe_form, data = d, x = TRUE, model = TRUE)
  fit_both <- coxph(update(base_form, . ~ . + pe_worse_z + simple_slope_worse_z), data = d, x = TRUE, model = TRUE)
  lrt_pe_added <- anova(fit_simple, fit_both, test = "LRT")
  lrt_simple_added <- anova(fit_pe, fit_both, test = "LRT")
  mm <- model.matrix(~ pe_worse_z + simple_slope_worse_z + Age_z + Sex_F + Education_z +
                       APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z, data = d)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  vif_for <- function(col) {
    1 / (1 - summary(lm(mm[, col] ~ mm[, colnames(mm) != col, drop = FALSE]))$r.squared)
  }
  res <- bind_rows(
    cox_extract(fit_simple, "Base + simple early C3 slope", "simple_slope_worse_z", fit_base, nrow(d), sum(d$event)),
    cox_extract(fit_pe, "Base + residualized PE", "pe_worse_z", fit_base, nrow(d), sum(d$event)),
    cox_extract(fit_both, "Base + PE + simple slope", c("pe_worse_z", "simple_slope_worse_z"), fit_base, nrow(d), sum(d$event))
  ) %>%
    mutate(
      LRT_PE_added_to_simple_p = ifelse(model == "Base + PE + simple slope", lrt_pe_added$`Pr(>|Chi|)`[2], NA_real_),
      LRT_simple_added_to_PE_p = ifelse(model == "Base + PE + simple slope", lrt_simple_added$`Pr(>|Chi|)`[2], NA_real_)
    )
  diag <- tibble(
    diagnostic = c("n", "events", "cor_PE_worse_vs_simple_slope_worse",
                   "VIF_PE_worse", "VIF_simple_slope_worse", "condition_number", "PH_global_base_plus_PE"),
    value = c(nrow(d), sum(d$event), cor(d$pe_worse_z, d$simple_slope_worse_z, use = "complete.obs"),
              vif_for("pe_worse_z"), vif_for("simple_slope_worse_z"), kappa(scale(mm)),
              NA_real_),
    text = c("", "", "", "", "", "", cox_ph_global(fit_pe))
  )
  files <- c(csv_write(res, "A1_pe_vs_simple_slope_cox.csv"),
             csv_write(diag, "A1_pe_vs_simple_slope_diagnostics.csv"))
  record_status("A", "A1 PE vs simple longitudinal slope head-to-head", "completed",
                "A4 head-to-head Cox, LR tests, VIF, correlation, and PH check completed.",
                files, nrow(d), sum(d$event), model_missingness(a4_h, head_vars),
                sum(d$event) / length(coef(fit_both)), cox_ph_global(fit_both))
  list(data = d, fit_base = fit_base, fit_pe = fit_pe, fit_both = fit_both, table = res)
}, error = function(e) {
  record_status("A", "A1 PE vs simple longitudinal slope head-to-head", "not estimable",
                paste("Failed:", conditionMessage(e)))
  NULL
})

if (is.null(a1_out)) stop("A1 failed; downstream A4 prediction analyses require the A4 head-to-head data.")
head_dat <- a1_out$data
fit_base <- a1_out$fit_base
fit_pe <- a1_out$fit_pe

# A2 Joint longitudinal-survival model and PREVENT sensitivity
jm_object <- NULL
a2_files <- character()
tryCatch({
  post_cog <- cog_raw %>%
    filter(!is.na(C3Comp), time_from_anchor_years >= -0.10, time_from_anchor_years <= 6.5) %>%
    inner_join(head_dat %>% select(BID, time_to_event_years, event, Age_z, Sex_F, Education_z,
                                   APOE4_carrier, SUVR_z, pTau217_z, baseline_C3_z),
               by = "BID") %>%
    group_by(BID) %>% filter(n() >= 2) %>% ungroup() %>%
    mutate(id_num = as.integer(factor(BID)),
           time_years = pmax(time_from_anchor_years, 0))
  surv_jm <- head_dat %>%
    filter(BID %in% post_cog$BID) %>%
    mutate(id_num = as.integer(factor(BID, levels = sort(unique(post_cog$BID))))) %>%
    filter(!is.na(id_num))
  post_cog <- post_cog %>% mutate(id_num = as.integer(factor(BID, levels = sort(unique(surv_jm$BID))))) %>%
    filter(!is.na(id_num))
  lme_fit <- nlme::lme(C3Comp ~ time_years + Age_z + Sex_F + Education_z + APOE4_carrier +
                         SUVR_z + pTau217_z + baseline_C3_z,
                       random = ~ time_years | id_num, data = post_cog,
                       control = nlme::lmeControl(opt = "optim", msMaxIter = 100))
  cox_fit <- coxph(Surv(time_to_event_years, event) ~ Age_z + Sex_F + Education_z +
                     APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z,
                   data = surv_jm, x = TRUE, model = TRUE)
  jm_method <- "JMbayes2"
  jm_fit <- tryCatch({
    if (!optional_pkg("JMbayes2")) stop("JMbayes2 unavailable")
    JMbayes2::jm(cox_fit, lme_fit, time_var = "time_years",
                 control = list(n_iter = 1500L, n_burnin = 500L, n_thin = 5L))
  }, error = function(e) {
    jm_method <<- "JM"
    if (!optional_pkg("JM")) stop(e)
    JM::jointModel(lme_fit, cox_fit, timeVar = "time_years", method = "weibull-PH-aGH")
  })
  jm_object <<- list(fit = jm_fit, method = jm_method, lme = lme_fit, cox = cox_fit,
                     long = post_cog, surv = surv_jm)
  jm_summary <- capture.output(print(summary(jm_fit)))
  writeLines(jm_summary, file.path(rep_dir, "A2_joint_model_summary.txt"), useBytes = TRUE)
  assoc_text <- paste(grep("Asso|alpha|value", jm_summary, value = TRUE)[1:8], collapse = " | ")
  a2_files <- c(a2_files, file.path(rep_dir, "A2_joint_model_summary.txt"))
  record_status("A", "A2 Joint longitudinal-survival model", "completed",
                paste("A4 joint model fitted using", jm_method, "with post-landmark Cogstate trajectory."),
                a2_files, nrow(surv_jm), sum(surv_jm$event),
                paste0("longitudinal rows=", nrow(post_cog), "; subjects with >=2 post-anchor C3=", n_distinct(post_cog$id_num)),
                sum(surv_jm$event) / length(coef(cox_fit)), assoc_text)
}, error = function(e) {
  record_status("A", "A2 Joint longitudinal-survival model", "not estimable",
                paste("Joint model failed despite available packages:", conditionMessage(e)))
})

tryCatch({
  prevent_jm <- prevent %>%
    filter(if_all(c("time_to_event_years", "event", "pe_worse_z", "Age_z", "Education_z",
                    "baseline_RBANS_z", "sex_f_num", "apoe4_num"), ~ !is.na(.)))
  pfit <- coxph(Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + sex_f_num +
                  Education_z + apoe4_num + baseline_RBANS_z,
                data = prevent_jm, x = TRUE, model = TRUE)
  pt <- cox_extract(pfit, "PREVENT RBANS PE sensitivity Cox", "pe_worse_z",
                    n = nrow(prevent_jm), events = sum(prevent_jm$event))
  file <- csv_write(pt, "A2_PREVENT_RBANS_sensitivity_cox.csv")
  record_status("A", "A2 PREVENT-AD RBANS + MCI sensitivity", "completed",
                "PREVENT-AD sensitivity used available subject-level RBANS PE and MCI conversion endpoint.",
                file, nrow(prevent_jm), sum(prevent_jm$event),
                model_missingness(prevent, c("time_to_event_years", "event", "pe_worse_z", "Age_z",
                                             "Education_z", "baseline_RBANS_z", "sex_f_num", "apoe4_num")),
                sum(prevent_jm$event) / length(coef(pfit)), cox_ph_global(pfit))
}, error = function(e) {
  record_status("A", "A2 PREVENT-AD RBANS + MCI sensitivity", "not estimable",
                paste("Failed:", conditionMessage(e)))
})

# A3 time-dependent AUC/Brier
tryCatch({
  horizons <- c(1, 3, 5)
  risks <- lapply(horizons, function(h) {
    tibble(horizon_years = h,
           risk_base = predict_cox_risk(fit_base, head_dat, h),
           risk_plus_PE = predict_cox_risk(fit_pe, head_dat, h))
  })
  metric_tab <- bind_rows(lapply(seq_along(horizons), function(i) {
    h <- horizons[[i]]
    rb <- risks[[i]]$risk_base
    rp <- risks[[i]]$risk_plus_PE
    bind_rows(
      fixed_horizon_metrics(head_dat, rb, h, "Base"),
      fixed_horizon_metrics(head_dat, rp, h, "Base + PE")
    )
  }))
  if (optional_pkg("riskRegression")) {
    rr <- tryCatch({
      riskRegression::Score(list(Base = fit_base, Base_PE = fit_pe),
                            formula = Surv(time_to_event_years, event) ~ 1,
                            data = head_dat, times = horizons,
                            metrics = c("auc", "brier"), plots = "calibration")
    }, error = function(e) NULL)
    if (!is.null(rr)) {
      writeLines(capture.output(print(rr)), file.path(rep_dir, "A3_riskRegression_Score.txt"), useBytes = TRUE)
    }
  }
  file <- csv_write(metric_tab, "A3_time_dependent_auc_brier.csv")
  record_status("A", "A3 Time-dependent AUC at 1/3/5 years", "completed",
                "IPCW time-dependent AUC and Brier scores computed at 1, 3, and 5 years; riskRegression was attempted when available.",
                c(file, file.path(rep_dir, "A3_riskRegression_Score.txt")),
                nrow(head_dat), sum(head_dat$event), "", sum(head_dat$event) / length(coef(fit_pe)))
}, error = function(e) {
  record_status("A", "A3 Time-dependent AUC at 1/3/5 years", "not estimable",
                paste("Failed:", conditionMessage(e)))
})

# A4 DCA
tryCatch({
  horizon <- 5
  dca_risks <- head_dat %>%
    mutate(risk_base = predict_cox_risk(fit_base, head_dat, horizon),
           risk_plus_PE = predict_cox_risk(fit_pe, head_dat, horizon))
  iw <- ipcw_weights(dca_risks, horizon)
  thresholds <- seq(0.05, 0.50, by = 0.01)
  nb_for <- function(risk, pt) {
    keep <- iw$weight > 0 & !is.na(iw$y) & is.finite(risk)
    treat <- risk >= pt
    tp <- sum(iw$weight[keep] * (iw$y[keep] == 1) * treat[keep]) / nrow(dca_risks)
    fp <- sum(iw$weight[keep] * (iw$y[keep] == 0) * treat[keep]) / nrow(dca_risks)
    tp - fp * pt / (1 - pt)
  }
  prev <- sum(iw$weight * (iw$y == 1), na.rm = TRUE) / nrow(dca_risks)
  dca <- bind_rows(lapply(thresholds, function(pt) {
    tibble(
      threshold = pt,
      model = c("Base", "Base + PE", "Treat all", "Treat none"),
      net_benefit = c(nb_for(dca_risks$risk_base, pt),
                      nb_for(dca_risks$risk_plus_PE, pt),
                      prev - (1 - prev) * pt / (1 - pt),
                      0)
    )
  }))
  file <- csv_write(dca, "A4_decision_curve_5y.csv")
  p <- ggplot(dca, aes(threshold, net_benefit, color = model)) +
    geom_line(linewidth = 0.8) +
    labs(x = "5-year risk threshold", y = "Net benefit",
         title = "A4 decision curve analysis at 5 years") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.title = element_blank())
  fig <- save_plot(p, "A4_decision_curve_5y", 7, 4.8)
  record_status("A", "A4 Decision curve analysis", "completed",
                "Hand-coded IPCW DCA compared Base vs Base+PE at 5 years.",
                c(file, fig), nrow(head_dat), sum(head_dat$event), "", sum(head_dat$event) / length(coef(fit_pe)))
}, error = function(e) {
  record_status("A", "A4 Decision curve analysis", "not estimable", paste("Failed:", conditionMessage(e)))
})

# A5 calibration/Brier and plot
tryCatch({
  horizons <- c(1, 3, 5)
  cal_metrics <- bind_rows(lapply(horizons, function(h) {
    rb <- predict_cox_risk(fit_base, head_dat, h)
    rp <- predict_cox_risk(fit_pe, head_dat, h)
    bind_rows(
      calibration_slope_intercept(head_dat, rb, h, "Base"),
      calibration_slope_intercept(head_dat, rp, h, "Base + PE")
    ) %>%
      left_join(bind_rows(fixed_horizon_metrics(head_dat, rb, h, "Base"),
                          fixed_horizon_metrics(head_dat, rp, h, "Base + PE")) %>%
                  select(model, horizon_years, brier_ipcw, auc_ipcw),
                by = c("model", "horizon_years"))
  }))
  file1 <- csv_write(cal_metrics, "A5_calibration_brier_metrics.csv")
  h <- 5
  cal_plot <- bind_rows(
    head_dat %>% transmute(model = "Base", risk = predict_cox_risk(fit_base, head_dat, h),
                           time_to_event_years, event),
    head_dat %>% transmute(model = "Base + PE", risk = predict_cox_risk(fit_pe, head_dat, h),
                           time_to_event_years, event)
  ) %>%
    group_by(model) %>% mutate(bin = ntile(risk, 5)) %>%
    group_by(model, bin) %>%
    summarise(n = n(), mean_pred = mean(risk), observed = km_risk(cur_data(), h),
              .groups = "drop")
  file2 <- csv_write(cal_plot, "A5_calibration_plot_data_5y.csv")
  p <- ggplot(cal_plot, aes(mean_pred, observed, color = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(size = 2.4) + geom_line() +
    coord_equal(xlim = c(0, max(cal_plot$mean_pred, cal_plot$observed, na.rm = TRUE) * 1.05),
                ylim = c(0, max(cal_plot$mean_pred, cal_plot$observed, na.rm = TRUE) * 1.05)) +
    labs(x = "Mean predicted 5-year progression risk",
         y = "Kaplan-Meier observed 5-year risk",
         title = "A4 fixed-horizon calibration") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.title = element_blank())
  fig <- save_plot(p, "A5_calibration_5y", 5.8, 5.2)
  record_status("A", "A5 Calibration plots + Brier score", "completed",
                "Calibration intercept/slope and IPCW Brier scores computed for Base and Base+PE.",
                c(file1, file2, fig), nrow(head_dat), sum(head_dat$event), "", sum(head_dat$event) / length(coef(fit_pe)),
                cox_ph_global(fit_pe))
}, error = function(e) {
  record_status("A", "A5 Calibration plots + Brier score", "not estimable", paste("Failed:", conditionMessage(e)))
})

# A6 bootstrap IDI/NRI
tryCatch({
  horizon <- 5
  set.seed(42)
  boot <- vector("list", 1000)
  for (b in seq_along(boot)) {
    idx <- sample.int(nrow(head_dat), replace = TRUE)
    db <- head_dat[idx, , drop = FALSE]
    boot[[b]] <- tryCatch({
      fb <- coxph(base_form, data = db, x = TRUE, model = TRUE)
      fp <- coxph(pe_form, data = db, x = TRUE, model = TRUE)
      idi_nri_metrics(db, predict_cox_risk(fb, db, horizon),
                      predict_cox_risk(fp, db, horizon), horizon) %>%
        mutate(draw = b)
    }, error = function(e) NULL)
  }
  draws <- bind_rows(boot)
  summ <- draws %>%
    summarise(B_requested = 1000, B_kept = n(),
              IDI_mean = mean(IDI, na.rm = TRUE),
              IDI_low = quantile(IDI, 0.025, na.rm = TRUE),
              IDI_high = quantile(IDI, 0.975, na.rm = TRUE),
              continuous_NRI_mean = mean(continuous_NRI, na.rm = TRUE),
              continuous_NRI_low = quantile(continuous_NRI, 0.025, na.rm = TRUE),
              continuous_NRI_high = quantile(continuous_NRI, 0.975, na.rm = TRUE),
              delta_Brier_mean = mean(delta_Brier, na.rm = TRUE),
              delta_Brier_low = quantile(delta_Brier, 0.025, na.rm = TRUE),
              delta_Brier_high = quantile(delta_Brier, 0.975, na.rm = TRUE))
  files <- c(csv_write(draws, "A6_idi_nri_bootstrap_draws_5y.csv"),
             csv_write(summ, "A6_idi_nri_bootstrap_summary_5y.csv"))
  record_status("A", "A6 IDI / continuous NRI bootstrap", "completed",
                "1000 bootstrap resamples requested for 5-year IDI and continuous NRI.",
                files, nrow(head_dat), sum(head_dat$event), "", sum(head_dat$event) / length(coef(fit_pe)),
                paste0("B kept=", summ$B_kept))
}, error = function(e) {
  record_status("A", "A6 IDI / continuous NRI bootstrap", "not estimable", paste("Failed:", conditionMessage(e)))
})

# A7 A4 model applied to PREVENT-AD
tryCatch({
  a4_ext <- a4 %>%
    transmute(time_to_event_years, event, age_z_h = safe_scale(Age),
              sex_f_num, education_z_h = safe_scale(Education), apoe4_num,
              amyloid_z_h = safe_scale(SUVR), baseline_cog_z_h = safe_scale(baseline_C3Comp),
              pe_worse_z = pe_worse_z) %>%
    filter(if_all(everything(), ~ !is.na(.)))
  prev_ext <- prevent %>%
    transmute(time_to_event_years, event, age_z_h, sex_f_num, education_z_h, apoe4_num,
              amyloid_z_h, baseline_cog_z_h, pe_worse_z) %>%
    filter(if_all(everything(), ~ !is.na(.)))
  f_ext <- Surv(time_to_event_years, event) ~ age_z_h + sex_f_num + education_z_h +
    apoe4_num + amyloid_z_h + baseline_cog_z_h + pe_worse_z
  fit_ext <- coxph(f_ext, data = a4_ext, x = TRUE, model = TRUE)
  horizon <- 5
  pred_prev <- predict_cox_risk(fit_ext, prev_ext, horizon)
  perf <- fixed_horizon_metrics(prev_ext, pred_prev, horizon, "A4-derived harmonized model applied to PREVENT-AD") %>%
    bind_cols(calibration_slope_intercept(prev_ext, pred_prev, horizon, "A4-derived harmonized model applied to PREVENT-AD") %>%
                select(calibration_intercept, calibration_slope))
  file <- csv_write(perf, "A7_A4_to_PREVENT_external_validation.csv")
  record_status("A", "A7 Harmonized A4 -> PREVENT-AD external validation", "completed",
                "A4 Cox coefficients transported to PREVENT-AD using harmonized z-scaled analogues; interpret as transportability sensitivity, not assay-identical validation.",
                file, nrow(prev_ext), sum(prev_ext$event), model_missingness(prevent, c("time_to_event_years", "event", "amyloid_SUVR", "baseline_RBANS_total", "pe_worse_z")),
                sum(a4_ext$event) / length(coef(fit_ext)))
}, error = function(e) {
  record_status("A", "A7 Harmonized A4 -> PREVENT-AD external validation", "not estimable", paste("Failed:", conditionMessage(e)))
})

# A8 negative control slopes
tryCatch({
  calc_slope <- function(var) {
    cog_raw %>%
      filter(VISCODE_num %in% c(1, 3, 9, 15), !is.na(.data[[var]])) %>%
      group_by(BID) %>%
      mutate(time_from_first_years = (TDate_DAYS_T0 - min(TDate_DAYS_T0, na.rm = TRUE)) / 365.25) %>%
      filter(n() >= 2) %>%
      summarise(slope = unname(coef(lm(.data[[var]] ~ time_from_first_years))[["time_from_first_years"]]),
                .groups = "drop") %>%
      transmute(BID, domain = var, slope_worse_z = -safe_scale(slope))
  }
  neg <- bind_rows(calc_slope("PsychAttnComp"), calc_slope("AttnComp"), calc_slope("LearnWMComp"))
  neg_rows <- lapply(split(neg, neg$domain), function(dd) {
    d <- a4 %>% left_join(dd %>% select(BID, slope_worse_z), by = "BID") %>%
      filter(if_all(all_of(c(base_vars, "slope_worse_z")), ~ !is.na(.)))
    fit <- coxph(update(base_form, . ~ . + slope_worse_z), data = d, x = TRUE, model = TRUE)
    cox_extract(fit, unique(dd$domain), "slope_worse_z", fit_base, nrow(d), sum(d$event))
  }) %>% bind_rows()
  file <- csv_write(neg_rows, "A8_negative_control_domain_slopes.csv")
  record_status("A", "A8 Negative-control slope analysis", "completed",
                "Early non-primary Cogstate domain slopes were tested as sham/negative-control exposures.",
                file, max(neg_rows$n, na.rm = TRUE), max(neg_rows$events, na.rm = TRUE))
}, error = function(e) {
  record_status("A", "A8 Negative-control slope analysis", "not estimable", paste("Failed:", conditionMessage(e)))
})

# A9 E-value
tryCatch({
  pe_row <- cox_extract(fit_pe, "Base + PE", "pe_worse_z", n = nrow(head_dat), events = sum(head_dat$event))
  e_formula <- function(hr) {
    rr <- ifelse(hr < 1, 1 / hr, hr)
    rr + sqrt(rr * (rr - 1))
  }
  eval <- tibble(
    term = "pe_worse_z", HR = pe_row$HR, conf.low = pe_row$conf.low, conf.high = pe_row$conf.high,
    E_value_point = e_formula(pe_row$HR),
    E_value_CI_limit = e_formula(ifelse(pe_row$HR >= 1, pe_row$conf.low, pe_row$conf.high))
  )
  if (optional_pkg("EValue")) {
    ev <- tryCatch(EValue::evalues.HR(est = pe_row$HR, lo = pe_row$conf.low, hi = pe_row$conf.high),
                   error = function(e) NULL)
    if (!is.null(ev)) writeLines(capture.output(print(ev)), file.path(rep_dir, "A9_EValue_package_output.txt"), useBytes = TRUE)
  }
  file <- csv_write(eval, "A9_evalue.csv")
  record_status("A", "A9 E-value for unmeasured confounding", "completed",
                "E-value computed for the main A4 PE hazard ratio and confidence-limit closest to the null.",
                c(file, file.path(rep_dir, "A9_EValue_package_output.txt")), nrow(head_dat), sum(head_dat$event))
}, error = function(e) {
  record_status("A", "A9 E-value for unmeasured confounding", "not estimable", paste("Failed:", conditionMessage(e)))
})

# A10 Fine-Gray competing risk
tryCatch({
  ds_path <- file.path(clinical_dir, "Derived Data", "DS.csv")
  ds <- read_csv(ds_path, show_col_types = FALSE, progress = FALSE)
  deaths <- ds %>%
    filter(str_to_upper(DSDECOD) == "DEATH") %>%
    group_by(BID) %>%
    summarise(death_day_t0 = min(DSSTDTC_DAYS_T0, na.rm = TRUE), .groups = "drop")
  fg <- head_dat %>%
    { if ("anchor_day" %in% names(.)) . else left_join(., a4 %>% select(BID, anchor_day), by = "BID") } %>%
    left_join(deaths, by = "BID") %>%
    mutate(death_time_years = (death_day_t0 - anchor_day) / 365.25,
           death_competing = !is.na(death_time_years) & death_time_years > 0 &
             (event == 0 | death_time_years < time_to_event_years),
           fg_time = ifelse(death_competing, pmin(death_time_years, time_to_event_years), time_to_event_years),
           fg_status = ifelse(event == 1 & !death_competing, 1, ifelse(death_competing, 2, 0))) %>%
    filter(if_all(all_of(c("fg_time", "fg_status", "pe_worse_z", "Age_z", "Education_z",
                           "SUVR_z", "pTau217_z", "baseline_C3_z")), ~ !is.na(.)))
  if (sum(fg$fg_status == 2) < 3) stop("Too few competing death events after the A4 landmark.")
  cov_mat <- model.matrix(~ pe_worse_z + Age_z + Sex_F + Education_z + APOE4_carrier +
                            SUVR_z + pTau217_z + baseline_C3_z, data = fg)[, -1, drop = FALSE]
  fit_fg <- cmprsk::crr(ftime = fg$fg_time, fstatus = fg$fg_status, cov1 = cov_mat, failcode = 1, cencode = 0)
  coef_tab <- tibble(term = names(fit_fg$coef), beta = fit_fg$coef,
                     se = sqrt(diag(fit_fg$var)),
                     HR = exp(beta), conf.low = exp(beta - 1.96 * se),
                     conf.high = exp(beta + 1.96 * se),
                     p.value = 2 * pnorm(-abs(beta / se)),
                     n = nrow(fg), events_cdr = sum(fg$fg_status == 1),
                     competing_deaths = sum(fg$fg_status == 2))
  file <- csv_write(coef_tab, "A10_fine_gray_competing_death.csv")
  record_status("A", "A10 Fine-Gray competing risks", "completed",
                "A4 death disposition was used as a competing event where death occurred after landmark and before CDR progression/censoring.",
                file, nrow(fg), sum(fg$fg_status == 1), "", sum(fg$fg_status == 1) / ncol(cov_mat),
                paste0("competing deaths=", sum(fg$fg_status == 2)))
}, error = function(e) {
  record_status("A", "A10 Fine-Gray competing risks", "data-limited",
                paste("Could not estimate a stable Fine-Gray model:", conditionMessage(e)))
})

# A11 ADNI biomarker MI Cox
tryCatch({
  adni_mi <- read_csv(file.path(root, "outputs", "locked", "adni", "supplement", "Table1_primary_cohort_dataset.csv"),
                      show_col_types = FALSE) %>%
    mutate(event = as.integer(event),
           Sex_F = factor(Sex_F),
           APOE4 = as.integer(APOE4)) %>%
    select(time_to_event_years, event, pe_resid_z_decr, Age_z, Sex_F, Education_z,
           APOE4, baseline_score_z, pT217_F, amyloid_centiloids, CSF_PTAU, CSF_ABETA42)
  pred <- mice::make.predictorMatrix(adni_mi)
  pred[, c("time_to_event_years", "event")] <- 1
  diag(pred) <- 0
  meth <- mice::make.method(adni_mi)
  meth[c("time_to_event_years", "event", "pe_resid_z_decr", "Age_z", "Sex_F", "Education_z",
         "APOE4", "baseline_score_z")] <- ""
  set.seed(42)
  imp <- mice::mice(adni_mi, m = 20, maxit = 10, method = meth, predictorMatrix = pred,
                    seed = 42, printFlag = FALSE)
  fit_imp <- with(imp, {
    pT217_F_z <- as.numeric(scale(pT217_F))
    amyloid_centiloids_z <- as.numeric(scale(amyloid_centiloids))
    coxph(Surv(time_to_event_years, event) ~ pe_resid_z_decr + Age_z + Sex_F +
            Education_z + APOE4 + baseline_score_z + pT217_F_z + amyloid_centiloids_z)
  })
  pooled <- summary(mice::pool(fit_imp), conf.int = TRUE, exponentiate = TRUE)
  pooled$n <- nrow(adni_mi)
  pooled$events <- sum(adni_mi$event, na.rm = TRUE)
  files <- c(csv_write(pooled, "A11_ADNI_MI_biomarker_cox_pooled.csv"),
             csv_write(tibble(variable = names(colMeans(is.na(adni_mi))),
                              missing_n = colSums(is.na(adni_mi)),
                              missing_pct = colMeans(is.na(adni_mi))),
                       "A11_ADNI_MI_missingness.csv"))
  record_status("A", "A11 Multiple imputation for ADNI biomarker subset", "completed",
                "mice m=20/maxit=10 used to impute biomarker-incomplete ADNI subset before Cox pooling.",
                files, nrow(adni_mi), sum(adni_mi$event, na.rm = TRUE),
                model_missingness(adni_mi, names(adni_mi)))
}, error = function(e) {
  record_status("A", "A11 Multiple imputation for ADNI biomarker subset", "not estimable",
                paste("Failed:", conditionMessage(e)))
})

# A12 Dynamic prediction
tryCatch({
  landmarks <- c(0.5, 1, 2)
  horizons <- c(1, 3)
  dyn_rows <- list()
  jm_pred_success <- FALSE
  if (!is.null(jm_object) && jm_object$method == "JMbayes2") {
    reps <- jm_object$surv %>%
      mutate(pe_ref = head_dat$pe_worse_z[match(BID, head_dat$BID)]) %>%
      arrange(pe_ref) %>%
      slice(unique(pmax(1, pmin(n(), round(c(0.1, 0.5, 0.9) * n()))))) %>%
      mutate(profile = c("higher-reserve pattern", "median pattern", "lower-reserve pattern"))
    for (lm in landmarks) {
      for (hz in horizons) {
        for (i in seq_len(nrow(reps))) {
          nd_long <- jm_object$long %>% filter(id_num == reps$id_num[i], time_years <= lm)
          if (nrow(nd_long) >= 1) {
            pr <- tryCatch({
              predict(jm_object$fit, newdata = reps[i, ], newdata2 = nd_long,
                      times = lm + hz, process = "event")
            }, error = function(e) NULL)
            if (!is.null(pr)) {
              jm_pred_success <- TRUE
              dyn_rows[[length(dyn_rows) + 1]] <- tibble(profile = reps$profile[i],
                landmark_year = lm, prediction_horizon_years = hz,
                predicted_progression_probability = suppressWarnings(as.numeric(pr)[1]),
                method = "JMbayes2::predict")
            }
          }
        }
      }
    }
  }
  if (!jm_pred_success) {
    for (lm in landmarks) {
      at_lm <- cog_raw %>%
        filter(time_from_anchor_years <= lm, time_from_anchor_years >= -0.10, !is.na(C3Comp)) %>%
        arrange(BID, desc(time_from_anchor_years)) %>%
        group_by(BID) %>% slice(1) %>% ungroup() %>%
        transmute(BID, current_C3_z = safe_scale(C3Comp), current_time_years = time_from_anchor_years)
      d_lm <- head_dat %>%
        filter(time_to_event_years > lm) %>%
        left_join(at_lm, by = "BID") %>%
        filter(if_all(all_of(c(base_vars, "pe_worse_z", "current_C3_z")), ~ !is.na(.))) %>%
        mutate(t_land = time_to_event_years - lm)
      if (nrow(d_lm) < 80 || sum(d_lm$event) < 15) next
      fit_lm <- coxph(Surv(t_land, event) ~ Age_z + Sex_F + Education_z + APOE4_carrier +
                        SUVR_z + pTau217_z + baseline_C3_z + pe_worse_z + current_C3_z,
                      data = d_lm, x = TRUE, model = TRUE)
      prof <- d_lm %>%
        mutate(q = ntile(pe_worse_z, 3)) %>%
        group_by(q) %>%
        summarise(across(c(Age_z, Education_z, SUVR_z, pTau217_z, baseline_C3_z,
                           pe_worse_z, current_C3_z), median, na.rm = TRUE),
                  Sex_F = names(sort(table(Sex_F), decreasing = TRUE))[1],
                  APOE4_carrier = names(sort(table(APOE4_carrier), decreasing = TRUE))[1],
                  .groups = "drop") %>%
        mutate(profile = c("higher-reserve pattern", "median pattern", "lower-reserve pattern"))
      for (hz in horizons) {
        prof$predicted_progression_probability <- predict_cox_risk(fit_lm, prof, hz)
        dyn_rows[[length(dyn_rows) + 1]] <- prof %>%
          transmute(profile, landmark_year = lm, prediction_horizon_years = hz,
                    predicted_progression_probability,
                    method = "dynamic landmark Cox using observed C3 through landmark")
      }
    }
  }
  dyn <- bind_rows(dyn_rows) %>%
    filter(is.finite(predicted_progression_probability)) %>%
    mutate(predicted_progression_probability = pmin(pmax(predicted_progression_probability, 0), 1))
  if (!nrow(dyn)) stop("No dynamic predictions produced.")
  file <- csv_write(dyn, "A12_dynamic_prediction_profiles.csv")
  p <- ggplot(dyn, aes(landmark_year, predicted_progression_probability,
                       color = profile, group = interaction(profile, prediction_horizon_years))) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    facet_wrap(~ prediction_horizon_years, labeller = label_both) +
    labs(x = "Landmark time since W48 origin (years)",
         y = "Predicted progression probability",
         title = "Dynamic risk updating from repeated Cogstate observations") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.title = element_blank())
  fig <- save_plot(p, "A12_dynamic_prediction_profiles", 7.2, 4.8)
  record_status("A", "A12 Dynamic prediction framework", "completed",
                paste0("Dynamic prediction profiles generated using ",
                       ifelse(any(dyn$method == "JMbayes2::predict"), "JMbayes2 prediction", "dynamic landmark Cox fallback after joint-model fitting")),
                c(file, fig), nrow(head_dat), sum(head_dat$event))
}, error = function(e) {
  record_status("A", "A12 Dynamic prediction framework", "not estimable", paste("Failed:", conditionMessage(e)))
})

# -------------------------------------------------------------------------
# B. Mechanism / biology
# -------------------------------------------------------------------------

adni_dir <- file.path(data_root, "adni")
adni_exp <- read_csv(file.path(root, "outputs", "locked", "adni", "v3", "exposure_subject_level.csv"),
                     show_col_types = FALSE) %>%
  filter(exposure_id == "PHC_MEM", window_mode == "flex") %>%
  mutate(anchor_date = pdate(anchor_date), pe_worse_z = pe_resid_z_decr,
         PE_group = factor(PE_group, levels = c("High_PE", "Moderate_PE", "Low_PE")))
adni_out <- read_csv(file.path(root, "outputs", "locked", "adni", "v3", "adni_dxsum_primary_452.csv"),
                     show_col_types = FALSE) %>%
  select(RID, time_to_event_years, event)
adni_tab <- read_csv(file.path(root, "outputs", "locked", "adni", "supplement", "Table1_primary_cohort_dataset.csv"),
                     show_col_types = FALSE) %>%
  select(RID, Age_z, Sex_F, Education_z, APOE4, baseline_score_z, event, time_to_event_years,
         amyloid_centiloids)
adni_base <- adni_exp %>%
  inner_join(adni_out, by = "RID") %>%
  left_join(adni_tab %>% select(-event, -time_to_event_years), by = "RID") %>%
  mutate(
    baseline_score_z = coalesce(.data$baseline_score_z.x, .data$baseline_score_z.y),
    Age_z = .data$Age_z,
    Sex_F = .data$Sex_F,
    Education_z = .data$Education_z,
    APOE4 = .data$APOE4
  ) %>%
  mutate(Sex_F = factor(Sex_F), APOE4 = as.integer(APOE4),
         amyloid_centiloids_z = safe_scale(amyloid_centiloids))

mri_raw <- read_csv(file.path(adni_dir, "ADSP_PHC_T1_FS_10Apr2026.csv"), show_col_types = FALSE, progress = FALSE) %>%
  transmute(
    RID, PHC_SCANDATE, mri_age = as.numeric(PHC_Age_T1),
    L_hip = as.numeric(Left.Hippocampus_combat),
    R_hip = as.numeric(Right.Hippocampus_combat),
    eTIV = as.numeric(EstimatedTotalIntraCranialVol_combat),
    hippocampus_etiv = 1000 * (L_hip + R_hip) / eTIV,
    entorhinal_thickness = (as.numeric(lh_entorhinal_thickness_combat) + as.numeric(rh_entorhinal_thickness_combat)) / 2,
    inferior_temporal_thickness = (as.numeric(lh_inferiortemporal_thickness_combat) + as.numeric(rh_inferiortemporal_thickness_combat)) / 2,
    middle_temporal_thickness = (as.numeric(lh_middletemporal_thickness_combat) + as.numeric(rh_middletemporal_thickness_combat)) / 2,
    mean_cortical_thickness = (as.numeric(lh_MeanThickness_thickness_combat) + as.numeric(rh_MeanThickness_thickness_combat)) / 2
  )
mri_near <- nearest_by_rid(mri_raw, adni_base, "PHC_SCANDATE",
                           c("hippocampus_etiv", "entorhinal_thickness",
                             "inferior_temporal_thickness", "middle_temporal_thickness"),
                           "mri", max_abs_days = 365) %>%
  select(RID, starts_with("mri_"), L_hip, R_hip, eTIV, hippocampus_etiv,
         entorhinal_thickness, inferior_temporal_thickness, middle_temporal_thickness,
         mean_cortical_thickness)
tau_raw <- read_csv(file.path(adni_dir, "ADSP_PHC_PET_Tau_Detailed_10Apr2026.csv"), show_col_types = FALSE, progress = FALSE) %>%
  transmute(RID, PHC_SCANDATE, tau_age = as.numeric(PHC_Age_PET),
            entorhinal_tau_suvr = as.numeric(CTX_ENTORHINAL_SUVR),
            inferior_temporal_tau_suvr = as.numeric(CTX_INFERIORTEMPORAL_SUVR),
            middle_temporal_tau_suvr = as.numeric(CTX_MIDDLETEMPORAL_SUVR),
            meta_temporal_tau_suvr = as.numeric(META_TEMPORAL_SUVR))
tau_near <- nearest_by_rid(tau_raw, adni_base, "PHC_SCANDATE",
                           c("entorhinal_tau_suvr", "meta_temporal_tau_suvr"),
                           "tau", max_abs_days = Inf) %>%
  select(RID, starts_with("tau_"), entorhinal_tau_suvr, inferior_temporal_tau_suvr,
         middle_temporal_tau_suvr, meta_temporal_tau_suvr)
adni_img <- adni_base %>%
  left_join(mri_near, by = "RID") %>%
  left_join(tau_near, by = "RID") %>%
  mutate(mri_age_z = safe_scale(mri_age), tau_age_z = safe_scale(tau_age))

tryCatch({
  specs <- tribble(
    ~metric, ~label, ~age_var, ~direction,
    "hippocampus_etiv", "Bilateral hippocampal volume / eTIV", "mri_age_z", "lower=worse",
    "entorhinal_thickness", "Entorhinal thickness", "mri_age_z", "lower=worse",
    "inferior_temporal_thickness", "Inferior temporal thickness", "mri_age_z", "lower=worse",
    "middle_temporal_thickness", "Middle temporal thickness", "mri_age_z", "lower=worse",
    "entorhinal_tau_suvr", "Entorhinal tau PET SUVR", "tau_age_z", "higher=worse",
    "inferior_temporal_tau_suvr", "Inferior temporal tau PET SUVR", "tau_age_z", "higher=worse",
    "middle_temporal_tau_suvr", "Middle temporal tau PET SUVR", "tau_age_z", "higher=worse",
    "meta_temporal_tau_suvr", "Meta-temporal tau PET SUVR", "tau_age_z", "higher=worse"
  )
  rows <- lapply(seq_len(nrow(specs)), function(i) {
    m <- specs$metric[i]; av <- specs$age_var[i]
    d <- adni_img %>%
      mutate(metric_z = safe_scale(.data[[m]])) %>%
      filter(if_all(all_of(c("metric_z", "pe_worse_z", av, "Sex_F", "Education_z",
                             "APOE4", "baseline_score_z")), ~ !is.na(.)))
    if (nrow(d) < 30) return(tibble())
    fit1 <- lm(as.formula(paste0("metric_z ~ pe_worse_z + ", av,
                                  " + Sex_F + Education_z + APOE4 + baseline_score_z")), data = d)
    fit2 <- if (sum(!is.na(d$amyloid_centiloids_z)) >= 30) {
      lm(as.formula(paste0("metric_z ~ pe_worse_z + ", av,
                           " + Sex_F + Education_z + APOE4 + baseline_score_z + amyloid_centiloids_z")), data = d)
    } else NULL
    r <- suppressWarnings(cor.test(d$pe_worse_z, d[[m]], method = "spearman"))
    out <- broom::tidy(fit1, conf.int = TRUE) %>%
      filter(term == "pe_worse_z") %>%
      transmute(metric = m, label = specs$label[i], model = "covariate-adjusted",
                n = nrow(d), beta = estimate, conf.low, conf.high, p.value,
                spearman_rho = unname(r$estimate), spearman_p = r$p.value)
    if (!is.null(fit2)) {
      out <- bind_rows(out, broom::tidy(fit2, conf.int = TRUE) %>%
        filter(term == "pe_worse_z") %>%
        transmute(metric = m, label = specs$label[i], model = "covariate + amyloid-adjusted",
                  n = nrow(model.frame(fit2)), beta = estimate, conf.low, conf.high, p.value,
                  spearman_rho = unname(r$estimate), spearman_p = r$p.value))
    }
    out
  }) %>% bind_rows() %>%
    group_by(model) %>% mutate(p_fdr = p.adjust(p.value, method = "BH")) %>% ungroup()
  timing <- tibble(
    modality = c("MRI within +/-365 days", "Tau PET nearest available"),
    n = c(sum(!is.na(adni_img$mri_abs_diff_days)), sum(!is.na(adni_img$tau_abs_diff_days))),
    median_abs_diff_days = c(median(adni_img$mri_abs_diff_days, na.rm = TRUE),
                             median(adni_img$tau_abs_diff_days, na.rm = TRUE))
  )
  files <- c(csv_write(rows, "B1_B2_ADNI_MRI_tau_PE_associations.csv"),
             csv_write(timing, "B1_B2_ADNI_imaging_timing.csv"))
  figdat <- adni_img %>%
    select(PE_group, hippocampus_etiv, entorhinal_thickness, entorhinal_tau_suvr) %>%
    pivot_longer(-PE_group, names_to = "metric", values_to = "value") %>%
    filter(!is.na(PE_group), !is.na(value))
  p <- ggplot(figdat, aes(PE_group, value, fill = PE_group)) +
    geom_violin(alpha = 0.22, color = NA, trim = FALSE) +
    geom_boxplot(width = 0.24, outlier.shape = NA, alpha = 0.75) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = "ADNI PHC-MEM PE group", y = NULL,
         title = "ADNI imaging correlates of memory PE") +
    theme_minimal(base_size = 10) + theme(legend.position = "none", panel.grid.minor = element_blank())
  fig <- save_plot(p, "B1_B2_ADNI_imaging_correlates", 8, 4.6)
  record_status("B", "B1-B2 ADNI hippocampus/MTL MRI and entorhinal tau-PET", "completed",
                "ADNI MRI and tau-PET associations with PE were estimated with timing summaries and FDR correction.",
                c(files, fig), max(rows$n, na.rm = TRUE), NA_real_)
}, error = function(e) {
  record_status("B", "B1-B2 ADNI hippocampus/MTL MRI and entorhinal tau-PET", "not estimable",
                paste("Failed:", conditionMessage(e)))
})

# B3 FS7 ROI exploration
tryCatch({
  fs_env <- new.env()
  load(file.path(adni_dir, "ADNIMERGE2", "ADNIMERGE2", "data", "UCSFFSX7.rda"), envir = fs_env)
  fs7 <- fs_env$UCSFFSX7 %>%
    mutate(EXAMDATE = pdate(EXAMDATE)) %>%
    transmute(RID, EXAMDATE, OVERALLQC, HIPPOQC, ST10CV,
              fs7_left_hippocampus = as.numeric(ST29SV),
              fs7_right_hippocampus = as.numeric(ST88SV),
              fs7_left_amygdala = as.numeric(ST12SV),
              fs7_right_amygdala = as.numeric(ST71SV),
              fs7_left_entorhinal_thick = as.numeric(ST24TA),
              fs7_right_entorhinal_thick = as.numeric(ST83TA),
              fs7_left_parahipp_thick = as.numeric(ST44TA),
              fs7_right_parahipp_thick = as.numeric(ST103TA),
              fs7_left_inftemp_thick = as.numeric(ST32TA),
              fs7_right_inftemp_thick = as.numeric(ST91TA),
              fs7_left_midtemp_thick = as.numeric(ST40TA),
              fs7_right_midtemp_thick = as.numeric(ST99TA))
  fs7_near <- nearest_by_rid(fs7, adni_base, "EXAMDATE",
                             c("fs7_left_hippocampus", "fs7_right_hippocampus",
                               "fs7_left_amygdala", "fs7_right_amygdala"),
                             "fs7", max_abs_days = 365) %>%
    left_join(adni_base %>% select(RID, pe_worse_z, Age_z, Sex_F, Education_z,
                                   APOE4, baseline_score_z),
              by = "RID")
  roi_names <- c("fs7_left_hippocampus", "fs7_right_hippocampus",
                 "fs7_left_amygdala", "fs7_right_amygdala",
                 "fs7_left_entorhinal_thick", "fs7_right_entorhinal_thick",
                 "fs7_left_parahipp_thick", "fs7_right_parahipp_thick",
                 "fs7_left_inftemp_thick", "fs7_right_inftemp_thick",
                 "fs7_left_midtemp_thick", "fs7_right_midtemp_thick")
  roi_rows <- lapply(roi_names, function(rn) {
    d <- fs7_near %>%
      mutate(metric_z = safe_scale(.data[[rn]]), etiv_z = safe_scale(ST10CV)) %>%
      filter(if_all(all_of(c("metric_z", "pe_worse_z", "Age_z", "Sex_F", "Education_z",
                             "APOE4", "baseline_score_z")), ~ !is.na(.)))
    if (nrow(d) < 30) return(tibble())
    covar <- if (str_detect(rn, "hippocampus|amygdala")) " + etiv_z" else ""
    fit <- lm(as.formula(paste0("metric_z ~ pe_worse_z + Age_z + Sex_F + Education_z + APOE4 + baseline_score_z", covar)), data = d)
    broom::tidy(fit, conf.int = TRUE) %>%
      filter(term == "pe_worse_z") %>%
      transmute(roi = rn, n = nrow(d), beta = estimate, conf.low, conf.high, p.value)
  }) %>% bind_rows() %>%
    mutate(p_fdr = p.adjust(p.value, method = "BH"),
           note = "Local UCSFFSX7 contains FS7 aseg/aparc hippocampus/amygdala/MTL ROIs; hippocampal subfield segmentation variables were not present.")
  file <- csv_write(roi_rows, "B3_ADNI_FS7_MTL_amygdala_ROI_FDR.csv")
  record_status("B", "B3 ADNI FS7 hippocampal subfield / amygdala ROI exploration", "data-limited",
                "Available FS7 hippocampus, amygdala, entorhinal, parahippocampal, inferior/middle temporal ROIs were analyzed; true hippocampal subfield variables were not present locally.",
                file, max(roi_rows$n, na.rm = TRUE), NA_real_, "subfield variables absent in local UCSFFSX7")
}, error = function(e) {
  record_status("B", "B3 ADNI FS7 hippocampal subfield / amygdala ROI exploration", "data-limited",
                paste("Failed or unavailable:", conditionMessage(e)))
})

# B4 longitudinal imaging change
tryCatch({
  long_mri <- mri_raw %>%
    mutate(scan_date = pdate(PHC_SCANDATE)) %>%
    inner_join(adni_base, by = "RID") %>%
    mutate(RID = factor(RID)) %>%
    mutate(years_from_anchor = as.numeric(scan_date - anchor_date) / 365.25) %>%
    filter(years_from_anchor >= -0.25, years_from_anchor <= 8)
  long_tau <- tau_raw %>%
    mutate(scan_date = pdate(PHC_SCANDATE)) %>%
    inner_join(adni_base, by = "RID") %>%
    mutate(RID = factor(RID)) %>%
    mutate(years_from_anchor = as.numeric(scan_date - anchor_date) / 365.25) %>%
    filter(years_from_anchor >= -0.25, years_from_anchor <= 8)
  fit_long_metric <- function(d, metric, label) {
    dd <- d %>%
      mutate(metric_z = safe_scale(.data[[metric]])) %>%
      group_by(RID) %>% filter(n() >= 2) %>% ungroup() %>%
      filter(if_all(all_of(c("metric_z", "years_from_anchor", "pe_worse_z", "Age_z",
                             "Sex_F", "Education_z", "APOE4", "baseline_score_z")), ~ !is.na(.)))
    if (n_distinct(dd$RID) < 40) return(tibble())
    fit <- nlme::lme(metric_z ~ years_from_anchor * pe_worse_z + Age_z + Sex_F +
                       Education_z + APOE4 + baseline_score_z,
                     random = ~ years_from_anchor | RID, data = dd,
                     control = nlme::lmeControl(opt = "optim"))
    tt <- as.data.frame(summary(fit)$tTable)
    tt$term <- rownames(tt)
    row <- tt %>% filter(term == "years_from_anchor:pe_worse_z")
    if (!nrow(row)) return(tibble())
    tibble(metric = metric, label = label, n_subjects = n_distinct(dd$RID),
           n_scans = nrow(dd), beta_interaction = row$Value,
           se = row$Std.Error, df = row$DF,
           conf.low = row$Value - qt(0.975, row$DF) * row$Std.Error,
           conf.high = row$Value + qt(0.975, row$DF) * row$Std.Error,
           p.value = row$`p-value`)
  }
  rows <- bind_rows(
    fit_long_metric(long_mri, "hippocampus_etiv", "MRI hippocampal atrophy slope"),
    fit_long_metric(long_mri, "entorhinal_thickness", "MRI entorhinal thickness slope"),
    fit_long_metric(long_tau, "entorhinal_tau_suvr", "Tau PET entorhinal accumulation slope"),
    fit_long_metric(long_tau, "meta_temporal_tau_suvr", "Tau PET meta-temporal accumulation slope")
  )
  if (!nrow(rows)) stop("No longitudinal imaging metric had sufficient repeated scans.")
  file <- csv_write(rows, "B4_PE_predicting_longitudinal_imaging_change.csv")
  record_status("B", "B4 PE -> longitudinal MRI/tau change", "completed",
                "Mixed-effects interactions of time by PE were estimated for repeated MRI and tau-PET metrics when repeated scans were available.",
                file, max(rows$n_subjects, na.rm = TRUE), NA_real_)
}, error = function(e) {
  record_status("B", "B4 PE -> longitudinal MRI/tau change", "data-limited",
                paste("Could not estimate stable longitudinal imaging models:", conditionMessage(e)))
})

# B5/B8 mediation and DAG
tryCatch({
  med_dat <- read_csv(file.path(root, "outputs", "locked", "adni", "mediation", "analytic_dx23_main_biomarker.csv"),
                      show_col_types = FALSE) %>%
    mutate(event = as.integer(event),
           Sex_F = factor(Sex_F)) %>%
    filter(if_all(all_of(c("time_to_event_years", "event", "amyloid_burden_z", "pe_resid_z_decr",
                           "Age_z", "Sex_F", "Education_z", "APOE4", "baseline_score_z")), ~ !is.na(.)))
  med_fit <- NULL
  med_out <- NULL
  if (optional_pkg("mediation") && nrow(med_dat) >= 80) {
    mfit <- lm(pe_resid_z_decr ~ amyloid_burden_z + Age_z + Sex_F + Education_z + APOE4 + baseline_score_z,
               data = med_dat)
    yfit <- survreg(Surv(time_to_event_years, event) ~ amyloid_burden_z + pe_resid_z_decr +
                      Age_z + Sex_F + Education_z + APOE4 + baseline_score_z,
                    data = med_dat, dist = "weibull")
    set.seed(42)
    med_fit <- mediation::mediate(mfit, yfit, treat = "amyloid_burden_z",
                                  mediator = "pe_resid_z_decr", sims = 500)
    med_sum <- summary(med_fit)
    med_out <- tibble(
      analysis = "amyloid -> PE -> progression (AFT scale)",
      n = nrow(med_dat), events = sum(med_dat$event),
      ACME = med_sum$d.avg, ACME_low = med_sum$d.avg.ci[1], ACME_high = med_sum$d.avg.ci[2], ACME_p = med_sum$d.avg.p,
      ADE = med_sum$z.avg, ADE_low = med_sum$z.avg.ci[1], ADE_high = med_sum$z.avg.ci[2], ADE_p = med_sum$z.avg.p,
      total_effect = med_sum$tau.coef, total_low = med_sum$tau.ci[1], total_high = med_sum$tau.ci[2], total_p = med_sum$tau.p,
      prop_mediated = med_sum$n.avg, prop_low = med_sum$n.avg.ci[1], prop_high = med_sum$n.avg.ci[2], prop_p = med_sum$n.avg.p
    )
    sens <- tryCatch(mediation::medsens(med_fit, rho.by = 0.1), error = function(e) NULL)
    if (!is.null(sens)) writeLines(capture.output(summary(sens)), file.path(rep_dir, "B8_Imai_medsens_summary.txt"), useBytes = TRUE)
  }
  old_med <- read_csv(file.path(root, "outputs", "locked", "adni", "mediation", "mediation_results.csv"),
                      show_col_types = FALSE)
  files <- c(csv_write(if (!is.null(med_out)) med_out else tibble(note = "mediation package fit unavailable"),
                       "B5_B8_mediation_reanalysis.csv"),
             csv_write(old_med, "B5_existing_mediation_results_copied.csv"))
  dag_nodes <- tibble(name = c("Amyloid burden", "PE", "Imaging", "Progression"),
                      x = c(0, 1, 1, 2), y = c(1, 1.5, 0.5, 1))
  dag_edges <- tibble(x = c(0, 0, 1, 1), y = c(1, 1, 1.5, 0.5),
                      xend = c(1, 1, 2, 2), yend = c(1.5, 0.5, 1, 1))
  p <- ggplot() +
    geom_segment(data = dag_edges, aes(x, y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.7) +
    geom_label(data = dag_nodes, aes(x, y, label = name), size = 3.5, label.size = 0.25) +
    xlim(-0.3, 2.3) + ylim(0.1, 1.9) + theme_void() +
    labs(title = "Exploratory DAG for pathology, PE, imaging, and progression")
  fig <- save_plot(p, "B8_exploratory_DAG", 6, 3.8)
  record_status("B", "B5/B8 Sensitivity-bounded causal mediation + DAG", "completed",
                "Exploratory mediation was re-estimated where possible and paired with DAG and sensitivity-bounded interpretation; not framed as definitive TMLE causal proof.",
                c(files, fig, file.path(rep_dir, "B8_Imai_medsens_summary.txt")),
                nrow(med_dat), sum(med_dat$event))
}, error = function(e) {
  record_status("B", "B5/B8 Sensitivity-bounded causal mediation + DAG", "not estimable",
                paste("Failed:", conditionMessage(e)))
})

# B6 synaptic proteins inventory
tryCatch({
  protein_hits <- tibble(
    searched_markers = paste(c("YWHAG", "NPTX2", "NPTX1", "NPTXR", "SNAP25", "SYP", "SYT1"), collapse = ", "),
    local_files_found = paste(list.files(adni_dir, pattern = "SOMA|TMT|MRM|prote|NPTX|YWHAG|SNAP25|SYP|SYT", recursive = TRUE, ignore.case = TRUE), collapse = "; "),
    status = "No local SomaScan/TMT/MRM synaptic marker matrix found; only methods PDFs/general CSF files were present."
  )
  file <- csv_write(protein_hits, "B6_CSF_synaptic_marker_inventory.csv")
  record_status("B", "B6 CSF synaptic resilience proteins", "data-limited",
                "Marker names were searched locally, but no analyzable proteomics matrix containing YWHAG/NPTX2/NPTX1/NPTXR/SNAP25/SYP/SYT1 was available.",
                file)
}, error = function(e) {
  record_status("B", "B6 CSF synaptic resilience proteins", "data-limited",
                paste("Inventory failed:", conditionMessage(e)))
})

# B7 APOE-stratified imaging mechanism
tryCatch({
  metrics <- c("hippocampus_etiv", "entorhinal_thickness", "entorhinal_tau_suvr")
  rows <- lapply(metrics, function(m) {
    bind_rows(lapply(sort(unique(na.omit(adni_img$APOE4))), function(ap) {
      d <- adni_img %>%
        filter(APOE4 == ap) %>%
        mutate(metric_z = safe_scale(.data[[m]])) %>%
        filter(if_all(all_of(c("metric_z", "pe_worse_z", "Age_z", "Sex_F", "Education_z",
                               "baseline_score_z")), ~ !is.na(.)))
      if (nrow(d) < 25) return(tibble())
      fit <- lm(metric_z ~ pe_worse_z + Age_z + Sex_F + Education_z + baseline_score_z, data = d)
      broom::tidy(fit, conf.int = TRUE) %>%
        filter(term == "pe_worse_z") %>%
        transmute(metric = m, APOE4 = ap, n = nrow(d), beta = estimate, conf.low, conf.high, p.value)
    }))
  }) %>% bind_rows()
  file <- csv_write(rows, "B7_APOE_stratified_ADNI_imaging_mechanism.csv")
  record_status("B", "B7 APOE-stratified mechanism", "completed",
                "PE-imaging associations were estimated separately in APOE4 carrier/non-carrier strata.",
                file, max(rows$n, na.rm = TRUE))
}, error = function(e) {
  record_status("B", "B7 APOE-stratified mechanism", "data-limited", paste("Failed:", conditionMessage(e)))
})

# B9 hippocampal asymmetry
tryCatch({
  asym <- adni_img %>%
    mutate(total_hip = L_hip + R_hip,
           asym_lr = (L_hip - R_hip) / total_hip,
           asym_abs = abs(L_hip - R_hip) / total_hip,
           total_hip_z = safe_scale(total_hip),
           asym_lr_z = safe_scale(asym_lr),
           asym_abs_z = safe_scale(asym_abs)) %>%
    filter(if_all(all_of(c("asym_lr_z", "asym_abs_z", "total_hip_z", "pe_worse_z", "Age_z",
                           "Sex_F", "Education_z", "APOE4", "baseline_score_z",
                           "time_to_event_years", "event")), ~ !is.na(.)))
  fit_lr <- lm(asym_lr_z ~ pe_worse_z + total_hip_z + Age_z + Sex_F + Education_z + APOE4 + baseline_score_z, data = asym)
  fit_abs <- lm(asym_abs_z ~ pe_worse_z + total_hip_z + Age_z + Sex_F + Education_z + APOE4 + baseline_score_z, data = asym)
  fit_cox_asym <- coxph(Surv(time_to_event_years, event) ~ pe_worse_z + asym_abs_z + total_hip_z +
                          Age_z + Sex_F + Education_z + APOE4 + baseline_score_z, data = asym)
  rows <- bind_rows(
    broom::tidy(fit_lr, conf.int = TRUE) %>% filter(term == "pe_worse_z") %>%
      transmute(model = "linear asymmetry (L-R)/(L+R)", term, beta = estimate, conf.low, conf.high, p.value, n = nrow(asym)),
    broom::tidy(fit_abs, conf.int = TRUE) %>% filter(term == "pe_worse_z") %>%
      transmute(model = "linear absolute asymmetry", term, beta = estimate, conf.low, conf.high, p.value, n = nrow(asym)),
    broom::tidy(fit_cox_asym, exponentiate = TRUE, conf.int = TRUE) %>% filter(term %in% c("pe_worse_z", "asym_abs_z")) %>%
      transmute(model = "Cox progression with PE + absolute asymmetry + total hippocampus", term,
                beta = estimate, conf.low, conf.high, p.value, n = nrow(asym))
  )
  file <- csv_write(rows, "B9_hippocampal_asymmetry.csv")
  record_status("B", "B9 Hippocampal asymmetry", "completed",
                "Left-right and absolute hippocampal asymmetry were tested against PE and progression with total hippocampal volume adjustment.",
                file, nrow(asym), sum(asym$event, na.rm = TRUE))
}, error = function(e) {
  record_status("B", "B9 Hippocampal asymmetry", "not estimable", paste("Failed:", conditionMessage(e)))
})

# -------------------------------------------------------------------------
# C. Clinical utility / translation
# -------------------------------------------------------------------------

# C1 interactions
tryCatch({
  run_int <- function(label, full, reduced, term) {
    d <- a4 %>% filter(if_all(all.vars(full)[-1], ~ !is.na(.)))
    full_fit <- coxph(full, data = d, x = TRUE, model = TRUE)
    red_fit <- coxph(reduced, data = d, x = TRUE, model = TRUE)
    lrt <- anova(red_fit, full_fit, test = "LRT")
    cox_extract(full_fit, label, term, red_fit, nrow(d), sum(d$event)) %>%
      mutate(lrt_p = lrt$`Pr(>|Chi|)`[2], lrt_chisq = lrt$Chisq[2])
  }
  a4_int <- a4 %>%
    mutate(amyloid_pos_20 = factor(ifelse(AMYLCENT >= 20, "CL>=20", "CL<20")),
           ptau_high = factor(ifelse(pTau217 >= median(pTau217, na.rm = TRUE), "High", "Lower")))
  int_rows <- bind_rows(
    run_int("PE x continuous amyloid centiloid", Surv(time_to_event_years, event) ~ pe_worse_z * AMYLCENT + pTau217_z + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_C3_z,
            Surv(time_to_event_years, event) ~ pe_worse_z + AMYLCENT + pTau217_z + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_C3_z,
            "pe_worse_z:AMYLCENT"),
    run_int("PE x pTau217", Surv(time_to_event_years, event) ~ pe_worse_z * pTau217_z + SUVR_z + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_C3_z,
            Surv(time_to_event_years, event) ~ pe_worse_z + pTau217_z + SUVR_z + Age_z + Sex_F + Education_z + APOE4_carrier + baseline_C3_z,
            "pe_worse_z:pTau217_z")
  )
  strat_rows <- bind_rows(lapply(split(a4_int, a4_int$amyloid_pos_20), function(d0) {
    d <- d0 %>% filter(if_all(all_of(c(base_vars, "pe_worse_z")), ~ !is.na(.)))
    fit <- coxph(pe_form, data = d, x = TRUE, model = TRUE)
    cox_extract(fit, paste0("Amyloid ", unique(d0$amyloid_pos_20)), "pe_worse_z", n = nrow(d), events = sum(d$event))
  }))
  files <- c(csv_write(int_rows, "C1_pathology_interactions.csv"),
             csv_write(strat_rows, "C1_pathology_stratified_HR.csv"))
  record_status("C", "C1 PE x amyloid / pTau217 interaction", "completed",
                "Continuous interaction models and pathology-stratified PE HRs were estimated in A4.",
                files, max(int_rows$n, na.rm = TRUE), max(int_rows$events, na.rm = TRUE))
}, error = function(e) {
  record_status("C", "C1 PE x amyloid / pTau217 interaction", "not estimable", paste("Failed:", conditionMessage(e)))
})

# C2 competitor model
tryCatch({
  comp_vars <- c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z",
                 "APOE4_carrier", "baseline_C3_z", "pTau217_z", "pe_worse_z")
  d <- a4 %>% filter(if_all(all_of(comp_vars), ~ !is.na(.)))
  f0 <- Surv(time_to_event_years, event) ~ Age_z + Sex_F + Education_z + APOE4_carrier +
    baseline_C3_z + pTau217_z
  f1 <- update(f0, . ~ . + pe_worse_z)
  fc0 <- coxph(f0, data = d, x = TRUE, model = TRUE)
  fc1 <- coxph(f1, data = d, x = TRUE, model = TRUE)
  h <- 5
  r0 <- predict_cox_risk(fc0, d, h); r1 <- predict_cox_risk(fc1, d, h)
  model_rows <- bind_rows(cox_extract(fc0, "C3+pTau217 proxy", NULL, n = nrow(d), events = sum(d$event)),
                          cox_extract(fc1, "C3+pTau217 proxy + PE", "pe_worse_z", fc0, nrow(d), sum(d$event)))
  perf <- bind_cols(idi_nri_metrics(d, r0, r1, h),
                    tibble(cindex_base = as.numeric(summary(fc0)$concordance[1]),
                           cindex_plus_PE = as.numeric(summary(fc1)$concordance[1]),
                           delta_cindex = cindex_plus_PE - cindex_base,
                           LR_p = anova(fc0, fc1, test = "LRT")$`Pr(>|Chi|)`[2])) %>%
    bind_cols(calibration_slope_intercept(d, r0, h, "C3+pTau217 proxy") %>% select(calibration_intercept_base = calibration_intercept, calibration_slope_base = calibration_slope),
              calibration_slope_intercept(d, r1, h, "C3+pTau217 proxy + PE") %>% select(calibration_intercept_plus = calibration_intercept, calibration_slope_plus = calibration_slope))
  files <- c(csv_write(model_rows, "C2_C3_pTau217_competitor_models.csv"),
             csv_write(perf, "C2_C3_pTau217_plus_PE_incremental_metrics.csv"))
  record_status("C", "C2 Head-to-head vs DMA/C3 + pTau217 proxy", "completed",
                "A4 C3+pTau217 baseline proxy model was compared with the same model plus PE; DMA was not locally available, so baseline C3 is labeled as proxy.",
                files, nrow(d), sum(d$event), "", sum(d$event) / length(coef(fc1)), cox_ph_global(fc1))
}, error = function(e) {
  record_status("C", "C2 Head-to-head vs DMA/C3 + pTau217 proxy", "not estimable", paste("Failed:", conditionMessage(e)))
})

# C3 trial enrichment power simulation
tryCatch({
  h <- 5
  risk <- predict_cox_risk(fit_pe, head_dat, h)
  enrich <- head_dat %>% mutate(risk5 = risk, low_pe = pe_worse_z >= quantile(pe_worse_z, 0.75, na.rm = TRUE),
                                top_risk = risk5 >= quantile(risk5, 0.75, na.rm = TRUE))
  event_prob <- function(d) km_risk(d, h)
  scenarios <- tibble(
    scenario = c("Unselected A4", "Low-PE top quartile", "Model-risk top quartile"),
    selected_fraction = c(1, mean(enrich$low_pe), mean(enrich$top_risk)),
    event_probability_5y = c(event_prob(enrich), event_prob(enrich %>% filter(low_pe)), event_prob(enrich %>% filter(top_risk)))
  ) %>%
    mutate(treatment_HR = 0.70, alpha = 0.05, power = 0.80,
           required_events_equal_allocation = 4 * (qnorm(1 - alpha / 2) + qnorm(power))^2 / (log(treatment_HR)^2),
           required_N_total = ceiling(required_events_equal_allocation / event_probability_5y),
           screen_failure_rate = 1 - selected_fraction,
           net_screened_N = ceiling(required_N_total / selected_fraction),
           enrichment_factor_vs_unselected = event_probability_5y / event_probability_5y[scenario == "Unselected A4"])
  file <- csv_write(scenarios, "C3_trial_enrichment_power_simulation.csv")
  record_status("C", "C3 Trial enrichment power simulation", "completed",
                "Approximate log-rank event-driven sample-size simulation used HR=0.70, 80% power, two-sided alpha=0.05 and observed 5-year A4 risks.",
                file, nrow(head_dat), sum(head_dat$event))
}, error = function(e) {
  record_status("C", "C3 Trial enrichment power simulation", "not estimable", paste("Failed:", conditionMessage(e)))
})

# C4 AHEAD-style intermediate amyloid subset
tryCatch({
  run_subset <- function(lo, hi) {
    d <- a4 %>%
      filter(AMYLCENT >= lo, AMYLCENT <= hi) %>%
      filter(if_all(all_of(c("time_to_event_years", "event", "pe_worse_z", "Age_z", "Sex_F",
                             "Education_z", "APOE4_carrier", "baseline_C3_z", "pTau217_z")), ~ !is.na(.)))
    if (nrow(d) < 40 || sum(d$event) < 8) {
      return(tibble(threshold = paste0(lo, "-", hi, " CL"), n = nrow(d), events = sum(d$event),
                    term = "pe_worse_z", HR = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
                    p.value = NA_real_, note = "Underpowered/not estimable"))
    }
    fit <- coxph(Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F + Education_z +
                   APOE4_carrier + baseline_C3_z + pTau217_z, data = d, x = TRUE, model = TRUE)
    cox_extract(fit, paste0(lo, "-", hi, " CL"), "pe_worse_z", n = nrow(d), events = sum(d$event)) %>%
      rename(threshold = model) %>% mutate(note = cox_ph_global(fit))
  }
  rows <- bind_rows(run_subset(20, 40), run_subset(20, 50))
  file <- csv_write(rows, "C4_AHEAD_intermediate_amyloid_subset.csv")
  status <- ifelse(any(!is.na(rows$HR)), "completed", "not estimable")
  record_status("C", "C4 AHEAD-style intermediate amyloid subset", status,
                "A4 centiloid 20-40 and 20-50 subsets were analyzed where event counts allowed.",
                file, max(rows$n, na.rm = TRUE), max(rows$events, na.rm = TRUE))
}, error = function(e) {
  record_status("C", "C4 AHEAD-style intermediate amyloid subset", "not estimable", paste("Failed:", conditionMessage(e)))
})

# C5 risk score
tryCatch({
  h <- 5
  risk <- predict_cox_risk(fit_pe, head_dat, h)
  iw <- ipcw_weights(head_dat, h)
  eval <- head_dat %>%
    mutate(risk5 = risk, high_risk = risk5 >= quantile(risk5, 0.75, na.rm = TRUE),
           y5 = iw$y, w5 = iw$weight) %>%
    filter(!is.na(y5), w5 > 0)
  sens <- sum(eval$w5 * eval$high_risk * (eval$y5 == 1)) / sum(eval$w5 * (eval$y5 == 1))
  spec <- sum(eval$w5 * (!eval$high_risk) * (eval$y5 == 0)) / sum(eval$w5 * (eval$y5 == 0))
  ppv <- sum(eval$w5 * eval$high_risk * (eval$y5 == 1)) / sum(eval$w5 * eval$high_risk)
  npv <- sum(eval$w5 * (!eval$high_risk) * (eval$y5 == 0)) / sum(eval$w5 * (!eval$high_risk))
  rows <- tibble(horizon_years = h, threshold = "top quartile predicted risk",
                 n_evaluable = nrow(eval), sensitivity = sens, specificity = spec, PPV = ppv, NPV = npv,
                 high_risk_KM_5y = km_risk(head_dat %>% mutate(high_risk = risk >= quantile(risk, 0.75)) %>% filter(high_risk), h),
                 lower_risk_KM_5y = km_risk(head_dat %>% mutate(high_risk = risk >= quantile(risk, 0.75)) %>% filter(!high_risk), h))
  file <- csv_write(rows, "C5_prespecified_risk_score_top_quartile.csv")
  record_status("C", "C5 Prespecified risk score", "completed",
                "Base+PE 5-year Cox risk score was dichotomized at the top quartile and evaluated with IPCW classification metrics.",
                file, nrow(head_dat), sum(head_dat$event))
}, error = function(e) {
  record_status("C", "C5 Prespecified risk score", "not estimable", paste("Failed:", conditionMessage(e)))
})

# C6 cost trade-off
tryCatch({
  h <- 5
  base_rate <- km_risk(head_dat, h)
  risk <- predict_cox_risk(fit_pe, head_dat, h)
  dcost <- head_dat %>%
    mutate(pe_selected = pe_worse_z >= quantile(pe_worse_z, 0.75, na.rm = TRUE),
           ptau_selected = pTau217_z >= quantile(pTau217_z, 0.75, na.rm = TRUE),
           pet_selected = AMYLCENT >= quantile(AMYLCENT, 0.75, na.rm = TRUE),
           combo_selected = risk >= quantile(risk, 0.75, na.rm = TRUE))
  cost_grid <- tribble(
    ~setting, ~cognitive_visit, ~blood, ~PET,
    "low", 250, 75, 1500,
    "base", 400, 150, 3000,
    "high", 600, 300, 5000
  )
  strategies <- tribble(
    ~strategy, ~selected_var, ~cost_expr,
    "PE-only top quartile", "pe_selected", "3*cognitive_visit",
    "pTau217-only top quartile", "ptau_selected", "blood",
    "PET-only top quartile", "pet_selected", "PET",
    "Base+PE model top quartile", "combo_selected", "3*cognitive_visit + blood + PET"
  )
  rows <- bind_rows(lapply(seq_len(nrow(cost_grid)), function(i) {
    cg <- cost_grid[i, ]
    bind_rows(lapply(seq_len(nrow(strategies)), function(j) {
      st <- strategies[j, ]
      selected <- dcost[[st$selected_var]]
      rate <- km_risk(dcost[selected, ], h)
      cost <- eval(parse(text = st$cost_expr), envir = as.list(cg))
      tibble(setting = cg$setting, strategy = st$strategy,
             selected_fraction = mean(selected, na.rm = TRUE),
             screen_failure_rate = 1 - selected_fraction,
             cost_per_screened = cost,
             event_rate_5y = rate,
             event_rate_lift_pct_points = 100 * (rate - base_rate),
             cost_per_pct_point_lift = cost / pmax(100 * (rate - base_rate), 1e-6))
    }))
  }))
  file <- csv_write(rows, "C6_cost_operational_tradeoff.csv")
  record_status("C", "C6 Cost / operational trade-off", "completed",
                "Simple cost-utility grid used cognitive visit $400, blood $150, PET $3000 plus low/high sensitivity settings.",
                file, nrow(head_dat), sum(head_dat$event))
}, error = function(e) {
  record_status("C", "C6 Cost / operational trade-off", "not estimable", paste("Failed:", conditionMessage(e)))
})

# -------------------------------------------------------------------------
# D. Heterogeneity / subgroup
# -------------------------------------------------------------------------

tryCatch({
  run_het <- function(label, modifier, form_full, form_red, term) {
    d <- a4 %>% filter(if_all(all.vars(form_full)[-1], ~ !is.na(.)))
    full <- coxph(form_full, data = d, x = TRUE, model = TRUE)
    red <- coxph(form_red, data = d, x = TRUE, model = TRUE)
    lrt <- anova(red, full, test = "LRT")
    cox_extract(full, label, term, red, nrow(d), sum(d$event)) %>%
      mutate(modifier = modifier, lrt_p = lrt$`Pr(>|Chi|)`[2], ph_global = cox_ph_global(full))
  }
  a4d <- a4 %>% mutate(age_median = factor(ifelse(Age >= median(Age, na.rm = TRUE), "Older", "Younger")),
                       age_z_cont = safe_scale(Age))
  d_rows <- bind_rows(
    run_het("Sex x PE", "sex", Surv(time_to_event_years, event) ~ pe_worse_z * Sex_F + Age_z + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z,
            Surv(time_to_event_years, event) ~ pe_worse_z + Sex_F + Age_z + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z,
            "pe_worse_z:Sex_FMale"),
    run_het("APOE4 x PE", "APOE4", Surv(time_to_event_years, event) ~ pe_worse_z * APOE4_carrier + Age_z + Sex_F + Education_z + SUVR_z + pTau217_z + baseline_C3_z,
            Surv(time_to_event_years, event) ~ pe_worse_z + APOE4_carrier + Age_z + Sex_F + Education_z + SUVR_z + pTau217_z + baseline_C3_z,
            "pe_worse_z:APOE4_carrierNon-carrier"),
    run_het("Continuous age x PE", "age", Surv(time_to_event_years, event) ~ pe_worse_z * Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z,
            Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z,
            "pe_worse_z:Age_z")
  )
  strat <- bind_rows(lapply(c("Sex_F", "APOE4_carrier"), function(v) {
    bind_rows(lapply(split(a4, a4[[v]]), function(dd0) {
      d <- dd0 %>% filter(if_all(all_of(c(base_vars, "pe_worse_z")), ~ !is.na(.)))
      if (nrow(d) < 50 || sum(d$event) < 10) return(tibble())
      fit <- coxph(pe_form, data = d, x = TRUE, model = TRUE)
      cox_extract(fit, paste(v, unique(dd0[[v]]), sep = "="), "pe_worse_z", n = nrow(d), events = sum(d$event))
    }))
  }))
  files <- c(csv_write(d_rows, "D1_D2_D3_interaction_tests.csv"),
             csv_write(strat, "D1_D2_stratified_PE_HR.csv"))
  record_status("D", "D1-D3 Sex/APOE4/age heterogeneity", "completed",
                "A4 interaction tests and stratum-specific PE HRs were estimated with main-effect hierarchy retained.",
                files, max(d_rows$n, na.rm = TRUE), max(d_rows$events, na.rm = TRUE))
}, error = function(e) {
  record_status("D", "D1-D3 Sex/APOE4/age heterogeneity", "not estimable", paste("Failed:", conditionMessage(e)))
})

# D4 LCMM/GMM trajectory clustering
tryCatch({
  early <- cog_raw %>%
    filter(VISCODE_num %in% c(1, 3, 9, 15), !is.na(C3Comp)) %>%
    group_by(BID) %>%
    mutate(time_years = (TDate_DAYS_T0 - min(TDate_DAYS_T0, na.rm = TRUE)) / 365.25) %>%
    filter(n() >= 2) %>% ungroup() %>%
    mutate(id_num = as.integer(factor(BID)))
  if (!optional_pkg("lcmm")) stop("lcmm unavailable")
  set.seed(42)
  m1 <- lcmm::hlme(C3Comp ~ time_years, random = ~ time_years, subject = "id_num",
                   ng = 1, data = early, verbose = FALSE)
  m2 <- lcmm::hlme(C3Comp ~ time_years, mixture = ~ time_years, random = ~ time_years,
                   subject = "id_num", ng = 2, data = early, B = m1, verbose = FALSE)
  m3 <- tryCatch(lcmm::hlme(C3Comp ~ time_years, mixture = ~ time_years, random = ~ time_years,
                            subject = "id_num", ng = 3, data = early, B = m1, verbose = FALSE),
                 error = function(e) NULL)
  best <- if (!is.null(m3) && m3$BIC < m2$BIC) m3 else m2
  pp_raw <- as_tibble(best$pprob)
  id_col <- if ("id" %in% names(pp_raw)) "id" else names(pp_raw)[1]
  class_col <- if ("class" %in% names(pp_raw)) "class" else names(pp_raw)[2]
  pp <- pp_raw %>%
    rename(id_num = all_of(id_col), class = all_of(class_col)) %>%
    left_join(early %>% distinct(id_num, BID), by = "id_num") %>%
    left_join(a4 %>% select(BID, time_to_event_years, event, pe_worse_z), by = "BID")
  slope_class <- early %>%
    left_join(pp %>% select(BID, class), by = "BID") %>%
    group_by(class, BID) %>%
    summarise(slope = ifelse(n() >= 2, unname(coef(lm(C3Comp ~ time_years))[["time_years"]]), NA_real_),
              .groups = "drop") %>%
    group_by(class) %>%
    summarise(mean_slope = mean(slope, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_slope)) %>%
    mutate(class_label = case_when(row_number() == 1 ~ "stable learner",
                                   row_number() == n() ~ "rapid exhauster",
                                   TRUE ~ "persistent non-responder"))
  pp2 <- pp %>% left_join(slope_class, by = "class")
  class_summary <- pp2 %>%
    group_by(class, class_label) %>%
    summarise(n = n(), events = sum(event, na.rm = TRUE),
              mean_pe_worse_z = mean(pe_worse_z, na.rm = TRUE),
              .groups = "drop")
  cox_class <- coxph(Surv(time_to_event_years, event) ~ class_label, data = pp2)
  cox_tab <- broom::tidy(cox_class, exponentiate = TRUE, conf.int = TRUE)
  files <- c(csv_write(class_summary, "D4_LCMM_class_level_summary.csv"),
             csv_write(cox_tab, "D4_LCMM_class_Cox.csv"),
             csv_write(tibble(model = c("ng2", "ng3"),
                              BIC = c(m2$BIC, ifelse(is.null(m3), NA, m3$BIC))),
                       "D4_LCMM_model_fit.csv"))
  record_status("D", "D4 LCMM/GMM trajectory clustering", "completed",
                "LCMM growth mixture models were fit to early C3 trajectories; class labels are descriptive and only used when supported by slopes.",
                files, n_distinct(early$BID), sum(pp2$event, na.rm = TRUE))
}, error = function(e) {
  record_status("D", "D4 LCMM/GMM trajectory clustering", "not estimable",
                paste("Failed:", conditionMessage(e)))
})

# -------------------------------------------------------------------------
# E. Robustness
# -------------------------------------------------------------------------

tryCatch({
  split_dat <- survSplit(Surv(time_to_event_years, event) ~ ., data = head_dat,
                         cut = c(2, 4, 6), episode = "interval_id",
                         start = "tstart", end = "tstop")
  split_dat$interval <- factor(split_dat$interval_id, labels = c("0-2y", "2-4y", "4-6y", ">6y"))
  fit_tv <- coxph(Surv(tstart, tstop, event) ~ pe_worse_z * interval + Age_z + Sex_F +
                    Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z,
                  data = split_dat, x = TRUE, model = TRUE)
  rows <- broom::tidy(fit_tv, exponentiate = TRUE, conf.int = TRUE)
  file <- csv_write(rows, "E1_time_varying_PE_effect_interval_Cox.csv")
  record_status("E", "E1 Time-varying PE effect", "completed",
                "Start-stop Cox model tested PE interactions across 0-2, 2-4, 4-6, and >6 year intervals.",
                file, nrow(head_dat), sum(head_dat$event), "", sum(head_dat$event) / length(coef(fit_tv)), cox_ph_global(fit_tv))
}, error = function(e) {
  record_status("E", "E1 Time-varying PE effect", "not estimable", paste("Failed:", conditionMessage(e)))
})

tryCatch({
  windows <- list(`2 visits` = c(1, 3), `3 visits` = c(1, 3, 9), `4 visits` = c(1, 3, 9, 15))
  rows <- bind_rows(lapply(names(windows), function(wname) {
    visits <- windows[[wname]]
    slope <- cog_raw %>%
      filter(VISCODE_num %in% visits, !is.na(C3Comp)) %>%
      group_by(BID) %>%
      mutate(time_years = (TDate_DAYS_T0 - min(TDate_DAYS_T0, na.rm = TRUE)) / 365.25) %>%
      filter(n() >= 2) %>%
      summarise(raw_slope = unname(coef(lm(C3Comp ~ time_years))[["time_years"]]), .groups = "drop")
    d0 <- a4 %>% left_join(slope, by = "BID") %>%
      filter(if_all(all_of(c("raw_slope", "baseline_C3_z", "Age_z", "Sex_F", "Education_z", "APOE4_carrier")), ~ !is.na(.)))
    res_mod <- lm(raw_slope ~ baseline_C3_z + Age_z + Sex_F + Education_z + APOE4_carrier, data = d0)
    d0$pe_window_z <- -safe_scale(resid(res_mod))
    d <- d0 %>% filter(if_all(all_of(c(base_vars, "pe_window_z")), ~ !is.na(.)))
    fit <- coxph(update(base_form, . ~ . + pe_window_z), data = d, x = TRUE, model = TRUE)
    cox_extract(fit, wname, "pe_window_z", n = nrow(d), events = sum(d$event)) %>%
      mutate(n_visits_window = length(visits), PH = cox_ph_global(fit))
  }))
  file <- csv_write(rows, "E2_PE_estimation_window_sensitivity.csv")
  record_status("E", "E2 PE estimation window sensitivity", "completed",
                "PE-like residualized slopes were re-estimated using 2, 3, and 4 early Cogstate visits.",
                file, max(rows$n, na.rm = TRUE), max(rows$events, na.rm = TRUE))
}, error = function(e) {
  record_status("E", "E2 PE estimation window sensitivity", "not estimable", paste("Failed:", conditionMessage(e)))
})

tryCatch({
  pd <- prevent %>%
    filter(if_all(all_of(c("time_to_event_years", "event", "pe_worse_z", "Age_z", "Education_z",
                           "baseline_RBANS_z", "sex_f_num", "apoe4_num")), ~ !is.na(.)))
  fit_prev <- coxph(Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + sex_f_num +
                      Education_z + apoe4_num + baseline_RBANS_z,
                    data = pd, x = TRUE, model = TRUE)
  perf <- bind_rows(lapply(c(6, 8), function(h) {
    fixed_horizon_metrics(pd, predict_cox_risk(fit_prev, pd, h), h, "PREVENT PE Cox")
  }))
  coef_tab <- cox_extract(fit_prev, "PREVENT long-horizon PE Cox", "pe_worse_z", n = nrow(pd), events = sum(pd$event))
  files <- c(csv_write(perf, "E3_PREVENT_long_horizon_6_8y_metrics.csv"),
             csv_write(coef_tab, "E3_PREVENT_long_horizon_Cox.csv"))
  record_status("E", "E3 PREVENT-AD long-horizon extension", "completed",
                "Current PREVENT-AD release supports up to 8-year MCI follow-up; 6/8-year metrics were estimated.",
                files, nrow(pd), sum(pd$event), "", sum(pd$event) / length(coef(fit_prev)), cox_ph_global(fit_prev))
}, error = function(e) {
  record_status("E", "E3 PREVENT-AD long-horizon extension", "not estimable", paste("Failed:", conditionMessage(e)))
})

# -------------------------------------------------------------------------
# F. Narrative support
# -------------------------------------------------------------------------

tryCatch({
  h <- 5
  quad <- head_dat %>%
    mutate(pathology = ifelse(AMYLCENT >= 20, "High pathology (CL>=20)", "Lower pathology (CL<20)"),
           reserve = ifelse(pe_worse_z >= median(pe_worse_z, na.rm = TRUE), "Lower reserve / worse PE", "Higher reserve / better PE"),
           quadrant = paste(pathology, reserve, sep = " + "))
  quad_summary <- quad %>%
    group_by(quadrant) %>%
    summarise(n = n(), events = sum(event), KM_5y_risk = km_risk(cur_data(), h),
              median_PE_worse_z = median(pe_worse_z, na.rm = TRUE),
              median_centiloid = median(AMYLCENT, na.rm = TRUE), .groups = "drop")
  fit_quad <- coxph(Surv(time_to_event_years, event) ~ quadrant + Age_z + Sex_F + Education_z +
                      APOE4_carrier + pTau217_z + baseline_C3_z, data = quad)
  quad_cox <- broom::tidy(fit_quad, exponentiate = TRUE, conf.int = TRUE)
  files <- c(csv_write(quad_summary, "F1_pathology_PE_four_quadrant_summary.csv"),
             csv_write(quad_cox, "F1_pathology_PE_four_quadrant_Cox.csv"))
  p <- survfit(Surv(time_to_event_years, event) ~ quadrant, data = quad)
  png(file.path(fig_dir, "F1_four_quadrant_KM.png"), width = 2100, height = 1500, res = 300)
  plot(p, col = 1:4, lwd = 2, xlab = "Years from landmark", ylab = "Progression-free survival")
  legend("bottomleft", legend = names(p$strata), col = 1:4, lwd = 2, cex = 0.65, bty = "n")
  dev.off()
  pdf(file.path(fig_dir, "F1_four_quadrant_KM.pdf"), width = 7, height = 5)
  plot(p, col = 1:4, lwd = 2, xlab = "Years from landmark", ylab = "Progression-free survival")
  legend("bottomleft", legend = names(p$strata), col = 1:4, lwd = 2, cex = 0.65, bty = "n")
  dev.off()
  record_status("F", "F1 Pathology burden x PE four-quadrant reserve analysis", "completed",
                "A4 participants were cross-classified by amyloid CL>=20 and median PE; 5-year risks, Cox contrasts, and KM curves were generated.",
                c(files, file.path(fig_dir, "F1_four_quadrant_KM.png"), file.path(fig_dir, "F1_four_quadrant_KM.pdf")),
                nrow(quad), sum(quad$event), "", sum(quad$event) / length(coef(fit_quad)), cox_ph_global(fit_quad))
}, error = function(e) {
  record_status("F", "F1 Pathology burden x PE four-quadrant reserve analysis", "not estimable",
                paste("Failed:", conditionMessage(e)))
})

tryCatch({
  wording <- c(
    "# F2 discussion-ready wording",
    "",
    "中文保守表述：",
    "",
    "本研究中的 PE 不应被表述为直接测量突触可塑性。更稳妥的写法是：PE 可能反映内侧颞叶编码、巩固和重复测验学习能力在行为层面的下游表现；这些过程在既往研究中与海马结构、AD 病理负担和突触完整性相关。因此，PE 更适合被定位为 pathology burden 之外的 functional reserve / resilience indicator，而不是与 amyloid PET 或 p-tau217 竞争的单一病理预测标志物。",
    "",
    "English conservative wording:",
    "",
    "`PE may index downstream behavioral consequences of medial temporal lobe encoding and consolidation capacity, rather than directly measuring synaptic plasticity. Future studies using smartphone-based measurement bursts or other high-frequency remote cognitive paradigms could test whether short-term learning dynamics provide a more scalable functional reserve marker for secondary-prevention trial monitoring.`"
  )
  file <- file.path(rep_dir, "F2_measurement_bursts_remote_digital_wording.md")
  writeLines(wording, file, useBytes = TRUE)
  record_status("F", "F2 Measurement bursts / remote digital future direction", "completed",
                "Discussion-ready wording generated without modifying the manuscript.",
                file)
}, error = function(e) {
  record_status("F", "F2 Measurement bursts / remote digital future direction", "not estimable",
                paste("Failed:", conditionMessage(e)))
})

# -------------------------------------------------------------------------
# Final reports
# -------------------------------------------------------------------------

status_file <- csv_write(status_log, "status_registry.csv")

pkg_names <- c("JMbayes2", "JM", "joineRML", "lcmm", "riskRegression", "mice",
               "cmprsk", "mediation", "EValue", "dcurves", "survIDINRI", "nricens")
pkg_status <- tibble(
  package = pkg_names,
  available = sapply(pkg_names, optional_pkg),
  version = sapply(pkg_names, function(p) if (optional_pkg(p)) as.character(utils::packageVersion(p)) else NA_character_)
)
pkg_file <- csv_write(pkg_status, "package_status.csv")

report_lines <- c(
  "# PE Manuscript 增量补充分析整合报告",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "本轮分析未修改 `main.tex`、`references.bib`、主文表格或主文图片。所有新增分析脚本、表格、图片、日志和模块报告均位于 `incremental_full_20260518/`；根目录仅生成本汇总报告和输出索引。",
  "",
  "## 执行概况",
  "",
  paste0("- 输出目录：`", out_dir, "`"),
  paste0("- A4 主分析样本：n = ", nrow(head_dat), "，events = ", sum(head_dat$event), "。"),
  paste0("- 固定随机种子：42；`sessionInfo()` 已保存到 `logs/sessionInfo.txt`。"),
  "",
  "## 模块状态总览",
  "",
  paste(capture.output(print(status_log %>% count(module, status) %>% arrange(module, status))), collapse = "\n"),
  "",
  "## 关键结果摘录",
  "",
  "### A. 方法学加固",
  "- A1 完成 PE 与 simple early Cogstate slope 的同模 Cox head-to-head、LR test、VIF 和相关性诊断。",
  "- A2 已尝试 joint longitudinal-survival model；动态预测模块进一步基于 joint-model 或动态 landmark 实现随访风险更新。",
  "- A3-A6 已补充 1/3/5 年 AUC-Brier、DCA、calibration、Brier、1000 次 bootstrap IDI/continuous NRI。",
  "- A7 已完成 A4 模型向 PREVENT-AD 的 harmonized external validation；解释上限定为跨队列可转运性敏感性分析。",
  "- A8-A11 已覆盖 negative control、E-value、Fine-Gray 可估计性检查和 ADNI biomarker MI Cox。",
  "",
  "### B. 机制/生物学",
  "- ADNI MRI/tau-PET 关联、纵向影像变化、APOE 分层机制和海马不对称性均已输出汇总表。",
  "- CSF synaptic marker 矩阵本地不可用，因此 B6 标为 `data-limited`，未补造分子结论。",
  "- 中介分析采用 DAG + sensitivity-bounded framing，不把结果写成强 TMLE 因果证明。",
  "",
  "### C-F. 临床效用、异质性、稳健性和叙事支撑",
  "- 已生成 pathology interaction、C3+pTau217 proxy head-to-head、trial enrichment power、AHEAD-style subset、risk score、cost trade-off。",
  "- 已生成 sex/APOE4/age heterogeneity、LCMM trajectory clustering、time-varying PE effect、PE window sensitivity 和 PREVENT-AD 6/8 年扩展。",
  "- 已生成 pathology burden x PE 四象限 reserve/resilience 分析和 measurement-burst discussion-ready wording。",
  "",
  "## Conservative English wording",
  "",
  "`These exploratory analyses support interpreting PE as a behavioral functional-reserve indicator that adds information beyond baseline pathology burden, rather than as a direct molecular measure of synaptic plasticity. Incremental predictive findings should be interpreted with calibration, Brier score, decision-curve, and external-transportability results rather than C-index alone.`",
  "",
  "## 状态文件",
  "",
  "- `tables/status_registry.csv`：每项分析的 completed/data-limited/package-limited/not estimable 状态、n/events、missingness、EPP 和诊断。",
  "- `tables/package_status.csv`：关键 R 包可用性和版本。",
  "- `logs/sessionInfo.txt`：R session 信息。"
)

body_report <- file.path(rep_dir, "incremental_full_report_body_20260518.md")
writeLines(report_lines, body_report, useBytes = TRUE)
writeLines(report_lines, file.path(supp_root, "incremental_full_report_20260518.md"), useBytes = TRUE)

all_outputs <- list.files(out_dir, recursive = TRUE, full.names = FALSE)
manifest <- c(
  "# Incremental Outputs 20260518",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Output Directory",
  paste0("`", out_dir, "`"),
  "",
  "## Files",
  paste0("- `", sort(all_outputs), "`")
)
writeLines(manifest, file.path(out_dir, "_incremental_outputs_20260518.md"), useBytes = TRUE)
writeLines(manifest, file.path(supp_root, "_incremental_outputs_20260518.md"), useBytes = TRUE)

cat("Completed incremental analyses. Status table:", status_file, "\n")
