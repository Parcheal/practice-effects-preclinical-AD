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


# Analysis: PE manuscript upgrade analyses D1-D7
# Date: 2026-05-20
# Random seed: 42
#
# Contract:
# - Do not modify manuscript .tex files.
# - Do not export subject-level identifiers.
# - Write all outputs inside this analysis bundle.

set.seed(42)

user_lib <- Sys.getenv("R_LIBS_USER")
.libPaths(unique(c("/usr/lib/R/site-library", "/usr/local/lib/R/site-library",
                   user_lib, .Library, .libPaths())))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(survival)
  library(broom)
  library(ggplot2)
  library(patchwork)
  library(survminer)
  library(cmprsk)
  library(metafor)
})

optional_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

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

find_root <- function(start) {
  pkg_root
}

find_clinical_dir <- function(root) {
  file.path(data_root, "a4")
}

args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, mustWork = TRUE)) else getwd()
out_dir <- normalizePath(script_dir, mustWork = TRUE)
root <- find_root(out_dir)
out_dir <- file.path(root, "verify", "common", "upgrade_D1_D7_20260520")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
clinical_dir <- find_clinical_dir(root)

tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
rep_dir <- file.path(out_dir, "reports")
log_dir <- file.path(out_dir, "logs")
src_dir <- file.path(out_dir, "source_data")
invisible(lapply(c(tab_dir, fig_dir, rep_dir, log_dir, src_dir),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

sink(file.path(log_dir, "sessionInfo.txt"))
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Root:", root, "\n")
cat("Output:", out_dir, "\n")
cat("Clinical dir:", clinical_dir, "\n")
cat(".libPaths():\n")
print(.libPaths())
print(sessionInfo())
sink()

csv_write <- function(x, file, subdir = tab_dir) {
  readr::write_csv(x, file.path(subdir, file), na = "")
  file.path(subdir, file)
}

source_write <- function(x, file) {
  csv_write(x, file, src_dir)
}

save_plot <- function(plot, basename, width = 7, height = 5, dpi = 600) {
  paths <- c(
    png = file.path(fig_dir, paste0(basename, ".png")),
    pdf = file.path(fig_dir, paste0(basename, ".pdf")),
    svg = file.path(fig_dir, paste0(basename, ".svg")),
    tiff = file.path(fig_dir, paste0(basename, ".tiff"))
  )
  ggsave(paths[["png"]], plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(paths[["pdf"]], plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  if (optional_pkg("svglite")) {
    ggsave(paths[["svg"]], plot, width = width, height = height, device = svglite::svglite, bg = "white")
  }
  if (optional_pkg("ragg")) {
    ggsave(paths[["tiff"]], plot, width = width, height = height, device = ragg::agg_tiff,
           dpi = dpi, compression = "lzw", bg = "white")
  }
  unname(paths[file.exists(paths)])
}

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.title = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = base_size + 1),
      strip.text = element_text(face = "bold", size = base_size),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

status_log <- tibble(
  task = character(), item = character(), status = character(), n = numeric(),
  events = numeric(), detail = character(), outputs = character()
)

record_status <- function(task, item, status, detail = "", outputs = character(),
                          n = NA_real_, events = NA_real_) {
  status_log <<- bind_rows(status_log, tibble(
    task = task, item = item, status = status, n = n, events = events,
    detail = detail,
    outputs = paste(gsub(paste0("^", out_dir, "/?"), "", outputs), collapse = "; ")
  ))
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
  den <- sum(w1) * sum(w0)
  num <- 0
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

km_risk <- function(d, horizon, time_var = "time_to_event_years", event_var = "event") {
  if (nrow(d) < 5 || sum(d[[event_var]], na.rm = TRUE) < 1) return(NA_real_)
  sf <- survfit(as.formula(paste0("Surv(", time_var, ",", event_var, ") ~ 1")), data = d)
  1 - summary(sf, times = horizon, extend = TRUE)$surv[[1]]
}

extract_loghr <- function(fit, term, cohort, n, events) {
  tt <- broom::tidy(fit, exponentiate = FALSE, conf.int = TRUE)
  row <- tt %>% filter(term == !!term)
  tibble(
    cohort = cohort, term = term, n = n, events = events,
    logHR = row$estimate[[1]], SE = row$std.error[[1]],
    HR = exp(row$estimate[[1]]),
    conf.low = exp(row$estimate[[1]] - 1.96 * row$std.error[[1]]),
    conf.high = exp(row$estimate[[1]] + 1.96 * row$std.error[[1]]),
    p.value = row$p.value[[1]]
  )
}

km_source <- function(fit, stratum_label) {
  ss <- summary(fit)
  tibble(
    stratum = stratum_label, time = ss$time, n.risk = ss$n.risk,
    n.event = ss$n.event, n.censor = ss$n.censor, surv = ss$surv,
    event_probability = 1 - ss$surv,
    curve = if (is.null(ss$strata)) "All" else as.character(ss$strata)
  )
}

cox_for_strata <- function(data, split_var, split_label, cox_formula, pe_term = "pe_worse_z") {
  bind_rows(lapply(split(data, data[[split_var]]), function(d0) {
    d <- d0 %>% filter(if_all(all.vars(cox_formula)[-1], ~ !is.na(.)))
    fit <- coxph(cox_formula, data = d, x = TRUE, model = TRUE)
    cox_extract(fit, unique(d[[split_label]]), pe_term, n = nrow(d), events = sum(d$event))
  }))
}

stratum_model_sensitivity <- function(data, split_var, split_label, formulas,
                                      pe_term = "pe_worse_z") {
  bind_rows(lapply(names(formulas), function(model_name) {
    bind_rows(lapply(split(data, data[[split_var]]), function(d0) {
      formula_i <- formulas[[model_name]]
      d <- d0 %>% filter(if_all(all.vars(formula_i)[-1], ~ !is.na(.)))
      fit <- coxph(formula_i, data = d, x = TRUE, model = TRUE)
      cox_extract(fit, unique(d[[split_label]]), pe_term,
                  n = nrow(d), events = sum(d$event)) %>%
        mutate(adjustment_model = model_name, .before = model)
    }))
  }))
}

make_km_panel <- function(d, title, label_text) {
  sf <- survfit(Surv(time_to_event_years, event) ~ PE_lowest_quartile, data = d)
  g <- survminer::ggsurvplot(
    sf, data = d, fun = "event", conf.int = FALSE, risk.table = FALSE,
    palette = c("Other PE" = "#6F7C85", "Lowest PE quartile" = "#B13F3B"),
    censor = FALSE, linewidth = 0.55
  )$plot
  g +
    labs(title = title, x = "Years after landmark", y = "Cumulative progression") +
    annotate("text", x = 0.15, y = 0.82, label = label_text, hjust = 0, vjust = 1, size = 2.2) +
    coord_cartesian(xlim = c(0, 8), ylim = c(0, 0.9), clip = "off") +
    theme_pub(7)
}

make_stratified_figure <- function(data, split_var, split_label, cox_rows, basename, forest_title) {
  forest_data <- cox_rows %>%
    mutate(
      model = factor(model, levels = rev(unique(model))),
      label = paste0("HR ", fmt(HR, 2), " (", fmt(conf.low, 2), "-", fmt(conf.high, 2),
                     "), P=", vapply(p.value, fmt_p, character(1)))
    )
  p_forest <- ggplot(forest_data, aes(HR, model)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, colour = "grey45") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.12, linewidth = 0.45) +
    geom_point(size = 1.8, colour = "#2D536A") +
    geom_text(aes(label = label), x = max(forest_data$conf.high, na.rm = TRUE) * 1.03,
              hjust = 0, size = 2.2) +
    scale_x_log10() +
    coord_cartesian(clip = "off") +
    labs(title = forest_title, x = "HR per SD lower PE", y = NULL) +
    theme_pub(7) +
    theme(legend.position = "none", plot.margin = margin(4, 42, 4, 4))

  km_plots <- lapply(levels(data[[split_var]]), function(level_i) {
    d0 <- data %>% filter(.data[[split_var]] == level_i)
    row <- cox_rows %>% filter(model == unique(d0[[split_label]]))
    label <- paste0("Continuous PE HR ", fmt(row$HR[1], 2), "\n95% CI ",
                    fmt(row$conf.low[1], 2), "-", fmt(row$conf.high[1], 2),
                    "\nP=", fmt_p(row$p.value[1]))
    make_km_panel(d0, unique(d0[[split_label]]), label)
  })

  plot <- p_forest + km_plots[[1]] + km_plots[[2]] +
    plot_layout(widths = c(1.05, 1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 8))
  save_plot(plot, basename, width = 10.5, height = 3.6)
}

parse_date_safe <- function(x) {
  suppressWarnings(as.Date(x, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%m/%d/%y", "%Y/%m/%d")))
}

finegray_table <- function(data, formula_covariates, cohort, pe_term = "pe_worse_z") {
  cov_mat <- model.matrix(formula_covariates, data = data)[, -1, drop = FALSE]
  fit <- cmprsk::crr(ftime = data$fg_time, fstatus = data$fg_status,
                     cov1 = cov_mat, failcode = 1, cencode = 0)
  beta <- fit$coef
  se <- sqrt(diag(fit$var))
  tibble(
    cohort = cohort, term = names(beta), beta = beta, se = se,
    HR = exp(beta), conf.low = exp(beta - 1.96 * se),
    conf.high = exp(beta + 1.96 * se),
    p.value = 2 * pnorm(-abs(beta / se)),
    n = nrow(data), clinical_events = sum(data$fg_status == 1),
    competing_deaths = sum(data$fg_status == 2)
  ) %>%
    mutate(primary_pe_term = term == pe_term)
}

cat("Root:", root, "\n")
cat("Output:", out_dir, "\n")

# -------------------------------------------------------------------------
# Load analysis cohorts
# -------------------------------------------------------------------------

a4 <- readRDS(file.path(root, "outputs", "locked", "a4",
                        "surv_primary_complete_w1_15_strict_preW48_lmer.rds")) %>%
  mutate(
    event = as.integer(converted),
    pe_worse_z = learning_slope_resid_z_decr_lmer,
    AMYLCENT_z = safe_scale(AMYLCENT),
    Sex_F = factor(ifelse(as.character(Sex_F) %in% c("1", "Female", "F", "TRUE"),
                          "Female", "Male")),
    APOE4_carrier = factor(ifelse(as.character(APOE4_carrier) %in% c("1", "Carrier", "TRUE"),
                                  "Carrier", "Non-carrier"))
  )

prevent <- read_csv(file.path(root, "outputs", "locked", "preventad", "survival_dataset_primary.csv"),
                    show_col_types = FALSE) %>%
  mutate(
    event = as.integer(converted),
    pe_worse_z = pe_resid_z_decr_lmer,
    Sex_F = factor(ifelse(as.character(Sex_F) %in% c("1", "Female", "F", "TRUE"),
                          "Female", "Male")),
    APOE4_carrier = factor(ifelse(as.character(APOE4_carrier) %in% c("1", "Carrier", "TRUE"),
                                  "Carrier", "Non-carrier"))
  )

adni_tab <- read_csv(file.path(root, "outputs", "locked", "adni", "supplement",
                               "Table1_primary_cohort_dataset.csv"),
                     show_col_types = FALSE) %>%
  mutate(
    event = as.integer(event),
    pe_worse_z = pe_resid_z_decr,
    CSF_PTAU_z = safe_scale(CSF_PTAU),
    CSF_ABETA42_z = safe_scale(CSF_ABETA42),
    amyloid_centiloids_z = safe_scale(amyloid_centiloids),
    Sex_F = factor(ifelse(as.character(Sex_F) %in% c("1", "Female", "F", "TRUE"),
                          "Female", "Male")),
    APOE4 = as.numeric(APOE4)
  )

adni_exp <- read_csv(file.path(root, "outputs", "locked", "adni", "v3", "exposure_subject_level.csv"),
                     show_col_types = FALSE) %>%
  filter(exposure_id == "PHC_MEM", source_role == "primary", window_mode == "flex") %>%
  select(RID, anchor_visit)

adni <- adni_tab %>%
  left_join(adni_exp, by = "RID")

# -------------------------------------------------------------------------
# D1. A4 amyloid-centiloid median stratification
# -------------------------------------------------------------------------

d1_outputs <- tryCatch({
  amy_median <- median(a4$AMYLCENT, na.rm = TRUE)
  d1 <- a4 %>%
    mutate(
      amyloid_stratum = factor(ifelse(AMYLCENT <= amy_median, "Lower amyloid burden", "Higher amyloid burden"),
                               levels = c("Lower amyloid burden", "Higher amyloid burden")),
      amyloid_stratum_label = paste0(as.character(amyloid_stratum), "\nmedian split at CL ",
                                     fmt(amy_median, 1))
    ) %>%
    group_by(amyloid_stratum) %>%
    mutate(
      PE_lowest_quartile = factor(ifelse(pe_worse_z >= quantile(pe_worse_z, 0.75, na.rm = TRUE),
                                         "Lowest PE quartile", "Other PE"),
                                  levels = c("Other PE", "Lowest PE quartile"))
    ) %>%
    ungroup()

  d1_form <- Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F + Education_z +
    APOE4_carrier + AMYLCENT_z + pTau217_z + baseline_C3_z
  d1_cox <- cox_for_strata(d1, "amyloid_stratum", "amyloid_stratum_label", d1_form)
  d1_sensitivity <- stratum_model_sensitivity(
    d1, "amyloid_stratum", "amyloid_stratum_label",
    list(
      unadjusted = Surv(time_to_event_years, event) ~ pe_worse_z,
      demographics = Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F +
        Education_z + APOE4_carrier,
      demographics_baseline_cognition = Surv(time_to_event_years, event) ~ pe_worse_z +
        Age_z + Sex_F + Education_z + APOE4_carrier + baseline_C3_z,
      plus_ptau217 = Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F +
        Education_z + APOE4_carrier + baseline_C3_z + pTau217_z,
      plus_ptau217_centiloid = d1_form
    )
  )
  d1_rates <- d1 %>%
    group_by(amyloid_stratum, PE_lowest_quartile) %>%
    summarise(n = n(), events = sum(event),
              event_probability_5y = km_risk(pick(time_to_event_years, event), 5),
              .groups = "drop") %>%
    mutate(amyloid_median_centiloid = amy_median)
  d1_km <- bind_rows(lapply(levels(d1$amyloid_stratum), function(level_i) {
    d0 <- d1 %>% filter(amyloid_stratum == level_i)
    sf <- survfit(Surv(time_to_event_years, event) ~ PE_lowest_quartile, data = d0)
    km_source(sf, unique(d0$amyloid_stratum_label))
  }))

  files <- c(
    csv_write(d1_cox, "D1_A4_amyloid_median_stratified_PE_Cox.csv"),
    csv_write(d1_sensitivity, "D1_A4_amyloid_median_stratified_model_sensitivity.csv"),
    csv_write(d1_rates, "D1_A4_amyloid_median_stratified_event_rates.csv"),
    source_write(d1_km, "D1_A4_amyloid_median_stratified_KM_source.csv"),
    make_stratified_figure(d1, "amyloid_stratum", "amyloid_stratum_label", d1_cox,
                           "D1_A4_amyloid_median_stratified_PE",
                           "A4 median amyloid strata")
  )
  record_status("D1", "A4 pathology-stratified PE Cox and KM", "completed",
                paste0("A4 risk set split at amyloid centiloid median ", fmt(amy_median, 2),
                       "; main continuous PE Cox adjusted for demographics, APOE4, baseline cognition, centiloid and p-tau217; model-sensitivity table also exported."),
                files, nrow(d1), sum(d1$event))
  list(data = d1, cox = d1_cox, sensitivity = d1_sensitivity, rates = d1_rates, files = files)
}, error = function(e) {
  record_status("D1", "A4 pathology-stratified PE Cox and KM", "failed", conditionMessage(e))
  NULL
})

# -------------------------------------------------------------------------
# D2. Two-step trial-enrichment workflow
# -------------------------------------------------------------------------

d2_outputs <- tryCatch({
  horizon <- 5
  ptau_cut <- quantile(a4$pTau217, 0.75, na.rm = TRUE)
  pe_cut <- quantile(a4$pe_worse_z, 0.75, na.rm = TRUE)
  required_events <- 4 * (qnorm(1 - 0.05 / 2) + qnorm(0.80))^2 / (log(0.70)^2)

  d2 <- a4 %>%
    mutate(
      ptau217_high = pTau217 >= ptau_cut,
      pe_lowest_quartile = pe_worse_z >= pe_cut,
      stage_unselected = TRUE,
      stage_ptau = ptau217_high,
      stage_ptau_pe = ptau217_high & pe_lowest_quartile
    )

  summarize_stage <- function(var, label, previous_n = nrow(d2)) {
    dd <- d2 %>% filter(.data[[var]])
    risk <- km_risk(dd, horizon)
    selected_fraction <- nrow(dd) / nrow(d2)
    tibble(
      stage = label, selected_rule = var, n_selected = nrow(dd), events = sum(dd$event),
      selected_fraction_of_A4 = selected_fraction,
      selected_fraction_of_previous = nrow(dd) / previous_n,
      event_probability_5y = risk,
      required_events_equal_allocation = required_events,
      required_randomized_N = ceiling(required_events / risk),
      net_screened_N = ceiling(ceiling(required_events / risk) / selected_fraction),
      NNS_screened_per_selected = 1 / selected_fraction,
      NNS_screened_per_5y_event = 1 / (selected_fraction * risk)
    )
  }

  s0 <- summarize_stage("stage_unselected", "Unselected A4", previous_n = nrow(d2))
  s1 <- summarize_stage("stage_ptau", "p-tau217 top quartile", previous_n = s0$n_selected)
  s2 <- summarize_stage("stage_ptau_pe", "p-tau217 top quartile + PE lowest quartile",
                        previous_n = s1$n_selected)
  d2_summary <- bind_rows(s0, s1, s2) %>%
    mutate(
      pTau217_top_quartile_cutoff = ptau_cut,
      PE_lowest_quartile_cutoff = pe_cut,
      randomized_N_reduction_vs_unselected_pct =
        100 * (1 - required_randomized_N / required_randomized_N[stage == "Unselected A4"]),
      net_screened_N_change_vs_unselected_pct =
        100 * (net_screened_N / net_screened_N[stage == "Unselected A4"] - 1)
    )

  p_n <- ggplot(d2_summary, aes(factor(stage, levels = stage), n_selected)) +
    geom_col(width = 0.62, fill = "#6F7C85") +
    geom_text(aes(label = paste0("n=", n_selected, "\n", fmt(100 * selected_fraction_of_A4, 1), "%")),
              vjust = -0.2, size = 2.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Screen-retained sample", x = NULL, y = "Participants retained") +
    theme_pub(7) +
    theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "none")

  p_event <- ggplot(d2_summary, aes(factor(stage, levels = stage), 100 * event_probability_5y)) +
    geom_col(width = 0.62, fill = "#B13F3B") +
    geom_text(aes(label = paste0(fmt(100 * event_probability_5y, 1), "%\n",
                                 "NNS ", fmt(NNS_screened_per_5y_event, 1))),
              vjust = -0.2, size = 2.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Observed 5-year event probability", x = NULL, y = "Progression probability (%)") +
    theme_pub(7) +
    theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "none")

  p_rand <- ggplot(d2_summary, aes(factor(stage, levels = stage), required_randomized_N)) +
    geom_col(width = 0.62, fill = "#2D536A") +
    geom_text(aes(label = paste0("N=", required_randomized_N, "\n",
                                 fmt(randomized_N_reduction_vs_unselected_pct, 1), "% red.")),
              vjust = -0.2, size = 2.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "Randomized N needed for fixed event target", x = NULL, y = "Required randomized N") +
    theme_pub(7) +
    theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "none")

  d2_fig <- (p_n | p_event | p_rand) +
    plot_annotation(title = "Two-step p-tau217 then PE trial-enrichment scenario",
                    subtitle = "Internal A4 cutoffs: p-tau217 top quartile followed by PE lowest quartile") &
    theme(plot.title = element_text(face = "bold", size = 9),
          plot.subtitle = element_text(size = 7))

  files <- c(
    csv_write(d2_summary, "D2_A4_two_step_ptau217_PE_enrichment.csv"),
    save_plot(d2_fig, "D2_A4_two_step_ptau217_PE_enrichment", width = 10.5, height = 3.7)
  )
  record_status("D2", "Two-step p-tau217 then PE enrichment", "completed",
                "p-tau217 high was operationalized as the A4 top quartile because no external assay cutoff is encoded in the locked dataset; PE lowest quartile used the global A4 PE threshold.",
                files, nrow(d2), sum(d2$event))
  list(summary = d2_summary, files = files)
}, error = function(e) {
  record_status("D2", "Two-step p-tau217 then PE enrichment", "failed", conditionMessage(e))
  NULL
})

# -------------------------------------------------------------------------
# D3. Decision curve analysis at 5 years
# -------------------------------------------------------------------------

d3_outputs <- tryCatch({
  d3 <- a4 %>%
    filter(if_all(all_of(c("time_to_event_years", "event", "Age_z", "Sex_F", "Education_z",
                           "APOE4_carrier", "baseline_C3_z", "AMYLCENT_z",
                           "pTau217_z", "pe_worse_z")), ~ !is.na(.)))
  horizon <- 5
  f_base <- Surv(time_to_event_years, event) ~ Age_z + Sex_F + Education_z +
    APOE4_carrier + baseline_C3_z + AMYLCENT_z
  f_pe <- update(f_base, . ~ . + pe_worse_z)
  f_ptau <- update(f_base, . ~ . + pTau217_z)
  f_ptau_pe <- update(f_base, . ~ . + pTau217_z + pe_worse_z)
  fits <- list(
    "Base" = coxph(f_base, data = d3, x = TRUE, model = TRUE),
    "Base + PE" = coxph(f_pe, data = d3, x = TRUE, model = TRUE),
    "Base + p-tau217" = coxph(f_ptau, data = d3, x = TRUE, model = TRUE),
    "Base + p-tau217 + PE" = coxph(f_ptau_pe, data = d3, x = TRUE, model = TRUE)
  )
  risks <- as_tibble(lapply(fits, predict_cox_risk, newdata = d3, horizon = horizon))
  names(risks) <- names(fits)
  iw <- ipcw_weights(d3, horizon)
  thresholds <- seq(0.05, 0.60, by = 0.01)
  nb_for <- function(risk, pt) {
    keep <- iw$weight > 0 & !is.na(iw$y) & is.finite(risk)
    treat <- risk >= pt
    tp <- sum(iw$weight[keep] * (iw$y[keep] == 1) * treat[keep]) / nrow(d3)
    fp <- sum(iw$weight[keep] * (iw$y[keep] == 0) * treat[keep]) / nrow(d3)
    tp - fp * pt / (1 - pt)
  }
  prev <- sum(iw$weight * (iw$y == 1), na.rm = TRUE) / nrow(d3)
  dca <- bind_rows(lapply(thresholds, function(pt) {
    bind_rows(lapply(names(risks), function(model_name) {
      tibble(threshold = pt, model = model_name,
             net_benefit = nb_for(risks[[model_name]], pt))
    })) %>%
      bind_rows(tibble(
        threshold = pt, model = c("Treat all", "Treat none"),
        net_benefit = c(prev - (1 - prev) * pt / (1 - pt), 0)
      ))
  }))
  perf <- bind_rows(lapply(names(fits), function(nm) fixed_horizon_metrics(d3, risks[[nm]], horizon, nm))) %>%
    mutate(c_index = vapply(fits, function(f) as.numeric(summary(f)$concordance[1]), numeric(1))[model])
  coef_tab <- bind_rows(lapply(names(fits), function(nm) {
    broom::tidy(fits[[nm]], exponentiate = TRUE, conf.int = TRUE) %>%
      transmute(model = nm, term, HR = estimate, conf.low, conf.high, p.value,
                n = nrow(d3), events = sum(d3$event))
  }))

  dca_plot <- dca %>%
    mutate(model = factor(model, levels = c("Base", "Base + PE", "Base + p-tau217",
                                            "Base + p-tau217 + PE", "Treat all", "Treat none")))
  p <- ggplot(dca_plot, aes(threshold, net_benefit, color = model, linetype = model)) +
    geom_line(linewidth = 0.65) +
    scale_color_manual(values = c(
      "Base" = "#6F7C85", "Base + PE" = "#B13F3B",
      "Base + p-tau217" = "#2D536A", "Base + p-tau217 + PE" = "#1E7A62",
      "Treat all" = "grey55", "Treat none" = "black"
    )) +
    scale_linetype_manual(values = c(
      "Base" = "solid", "Base + PE" = "solid", "Base + p-tau217" = "solid",
      "Base + p-tau217 + PE" = "solid", "Treat all" = "dashed", "Treat none" = "dotted"
    )) +
    coord_cartesian(ylim = c(-0.05, max(dca_plot$net_benefit, na.rm = TRUE) * 1.05)) +
    labs(title = "A4 decision curve analysis at 5 years",
         x = "5-year progression-risk threshold", y = "Net benefit") +
    theme_pub(7)

  files <- c(
    csv_write(dca, "D3_A4_decision_curve_5y.csv"),
    csv_write(perf, "D3_A4_decision_curve_model_performance.csv"),
    csv_write(coef_tab, "D3_A4_DCA_model_coefficients.csv"),
    save_plot(p, "D3_A4_decision_curve_5y", width = 6.9, height = 4.2)
  )
  record_status("D3", "A4 5-year DCA with base, PE, p-tau217 and combined models",
                "completed",
                "Base model included demographics, APOE4, baseline cognition and continuous amyloid centiloid; PE and p-tau217 were then added separately and jointly.",
                files, nrow(d3), sum(d3$event))
  list(dca = dca, perf = perf, files = files)
}, error = function(e) {
  record_status("D3", "A4 5-year DCA with base, PE, p-tau217 and combined models",
                "failed", conditionMessage(e))
  NULL
})

# -------------------------------------------------------------------------
# D4. Random-effects meta-analysis across cohorts
# -------------------------------------------------------------------------

d4_outputs <- tryCatch({
  a4_base <- a4 %>%
    filter(if_all(all_of(c("time_to_event_years", "event", "pe_worse_z", "Age_z",
                           "Sex_F", "Education_z", "APOE4_carrier", "baseline_C3_z")), ~ !is.na(.)))
  fit_a4 <- coxph(Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F +
                    Education_z + APOE4_carrier + baseline_C3_z,
                  data = a4_base, x = TRUE, model = TRUE)

  prev_base <- prevent %>%
    filter(if_all(all_of(c("time_to_event_years", "event", "pe_worse_z", "Age_z",
                           "Sex_F", "Education_z", "APOE4_carrier", "baseline_RBANS_z")), ~ !is.na(.)))
  fit_prev <- coxph(Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F +
                      Education_z + APOE4_carrier + baseline_RBANS_z,
                    data = prev_base, x = TRUE, model = TRUE)

  adni_base <- adni %>%
    filter(if_all(all_of(c("time_to_event_years", "event", "pe_worse_z", "Age_z",
                           "Sex_F", "Education_z", "APOE4", "baseline_score_z")), ~ !is.na(.)))
  fit_adni <- coxph(Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F +
                      Education_z + APOE4 + baseline_score_z,
                    data = adni_base, x = TRUE, model = TRUE)

  meta_input <- bind_rows(
    extract_loghr(fit_a4, "pe_worse_z", "A4", nrow(a4_base), sum(a4_base$event)),
    extract_loghr(fit_prev, "pe_worse_z", "PREVENT-AD", nrow(prev_base), sum(prev_base$event)),
    extract_loghr(fit_adni, "pe_worse_z", "ADNI", nrow(adni_base), sum(adni_base$event))
  )
  meta_fit <- metafor::rma(yi = logHR, sei = SE, data = meta_input, method = "REML")
  meta_summary <- tibble(
    k = meta_fit$k, pooled_logHR = as.numeric(meta_fit$b), pooled_SE = meta_fit$se,
    pooled_HR = exp(as.numeric(meta_fit$b)),
    conf.low = exp(meta_fit$ci.lb), conf.high = exp(meta_fit$ci.ub),
    p.value = meta_fit$pval, tau2 = meta_fit$tau2, I2 = meta_fit$I2,
    H2 = meta_fit$H2, Q = meta_fit$QE, Q_p.value = meta_fit$QEp,
    method = "REML random-effects meta-analysis of log HRs"
  )

  forest_df <- meta_input %>%
    mutate(cohort = factor(cohort, levels = rev(cohort)),
           label = paste0(fmt(HR, 2), " (", fmt(conf.low, 2), "-", fmt(conf.high, 2), ")"))
  pooled_row <- tibble(
    cohort = factor("Pooled", levels = c(levels(forest_df$cohort), "Pooled")),
    HR = meta_summary$pooled_HR, conf.low = meta_summary$conf.low,
    conf.high = meta_summary$conf.high,
    label = paste0(fmt(meta_summary$pooled_HR, 2), " (",
                   fmt(meta_summary$conf.low, 2), "-", fmt(meta_summary$conf.high, 2), ")")
  )
  plot_df <- bind_rows(
    forest_df %>% transmute(cohort = as.character(cohort), HR, conf.low, conf.high, label),
    pooled_row %>% transmute(cohort = as.character(cohort), HR, conf.low, conf.high, label)
  ) %>%
    mutate(cohort = factor(cohort, levels = c("Pooled", "ADNI", "PREVENT-AD", "A4")))
  p <- ggplot(plot_df, aes(HR, cohort)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, colour = "grey45") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.13, linewidth = 0.45) +
    geom_point(aes(shape = cohort == "Pooled"), size = 2.0, colour = "#2D536A") +
    geom_text(aes(label = label), x = max(plot_df$conf.high) * 1.04, hjust = 0, size = 2.4) +
    scale_x_log10() +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
    coord_cartesian(clip = "off") +
    labs(title = "Cross-cohort random-effects meta-analysis",
         subtitle = paste0("Pooled HR ", fmt(meta_summary$pooled_HR, 2),
                           "; I2=", fmt(meta_summary$I2, 1), "%"),
         x = "HR per SD lower PE", y = NULL) +
    theme_pub(7) +
    theme(plot.margin = margin(4, 54, 4, 4), legend.position = "none")

  files <- c(
    csv_write(meta_input, "D4_cross_cohort_continuous_PE_meta_inputs.csv"),
    csv_write(meta_summary, "D4_cross_cohort_random_effects_meta_summary.csv"),
    save_plot(p, "D4_cross_cohort_random_effects_meta_forest", width = 5.8, height = 3.4)
  )
  record_status("D4", "Cross-cohort random-effects meta-analysis", "completed",
                "Continuous PE HRs were refit from the three local risk sets using the base covariate layer described in the manuscript.",
                files, sum(meta_input$n), sum(meta_input$events))
  list(input = meta_input, summary = meta_summary, files = files)
}, error = function(e) {
  record_status("D4", "Cross-cohort random-effects meta-analysis", "failed", conditionMessage(e))
  NULL
})

# -------------------------------------------------------------------------
# D5. Death counts and Fine-Gray sensitivity where records permit
# -------------------------------------------------------------------------

d5_outputs <- tryCatch({
  death_summary <- list()
  finegray_rows <- list()

  # A4 death records
  if (!is.na(clinical_dir)) {
    ds_path <- file.path(clinical_dir, "Derived Data", "DS.csv")
    if (file.exists(ds_path)) {
      ds <- read_csv(ds_path, show_col_types = FALSE, progress = FALSE)
      deaths <- ds %>%
        filter(str_to_upper(DSDECOD) == "DEATH") %>%
        group_by(BID) %>%
        summarise(death_day_t0 = min(DSSTDTC_DAYS_T0, na.rm = TRUE), .groups = "drop")
      a4_fg <- a4 %>%
        left_join(deaths, by = "BID") %>%
        mutate(
          death_time_years = (death_day_t0 - anchor_day) / 365.25,
          death_after_landmark = !is.na(death_time_years) & death_time_years > 0,
          death_competing = death_after_landmark & (event == 0 | death_time_years < time_to_event_years),
          fg_time = ifelse(death_competing, death_time_years, time_to_event_years),
          fg_status = ifelse(event == 1 & !death_competing, 1, ifelse(death_competing, 2, 0))
        ) %>%
        filter(if_all(all_of(c("fg_time", "fg_status", "pe_worse_z", "Age_z", "Sex_F",
                               "Education_z", "APOE4_carrier", "AMYLCENT_z",
                               "pTau217_z", "baseline_C3_z")), ~ !is.na(.)))
      death_summary[["A4"]] <- tibble(
        cohort = "A4", n = nrow(a4_fg), clinical_events = sum(a4_fg$fg_status == 1),
        death_records_available = TRUE,
        deaths_after_landmark_any_timing = sum(a4_fg$death_after_landmark, na.rm = TRUE),
        competing_deaths_before_event_or_censor = sum(a4_fg$fg_status == 2),
        competing_death_pct = 100 * competing_deaths_before_event_or_censor / n,
        fine_gray_triggered_by_gt5pct_rule = competing_death_pct > 5,
        note = "A4 death disposition from Derived Data/DS.csv."
      )
      if (sum(a4_fg$fg_status == 2) >= 3) {
        finegray_rows[["A4"]] <- finegray_table(
          a4_fg,
          ~ pe_worse_z + Age_z + Sex_F + Education_z + APOE4_carrier +
            AMYLCENT_z + pTau217_z + baseline_C3_z,
          "A4", "pe_worse_z"
        )
      }
    }
  }

  # PREVENT-AD local files do not encode death/disposition in the survival dataset.
  death_summary[["PREVENT-AD"]] <- tibble(
    cohort = "PREVENT-AD", n = nrow(prevent), clinical_events = sum(prevent$event, na.rm = TRUE),
    death_records_available = FALSE,
    deaths_after_landmark_any_timing = NA_real_,
    competing_deaths_before_event_or_censor = NA_real_,
    competing_death_pct = NA_real_,
    fine_gray_triggered_by_gt5pct_rule = FALSE,
    note = "No death/deceased/mortality field was found in the local PREVENT-AD replication survival or clinical CSV files used by this project."
  )

  # ADNI death records from adverse-event files.
  adv_path <- file.path(data_root, "adni", "ADVERSE_10Apr2026.csv")
  recadv_path <- file.path(data_root, "adni", "RECADV_10Apr2026.csv")
  if (file.exists(adv_path) && file.exists(recadv_path)) {
    adv <- suppressWarnings(read_csv(adv_path, col_types = cols(.default = col_character()), progress = FALSE))
    recadv <- suppressWarnings(read_csv(recadv_path, col_types = cols(.default = col_character()), progress = FALSE))
    deaths_adv <- adv %>%
      filter(SAEDEATH %in% c("Yes", "1", "TRUE") | !is.na(AEHDTHDT)) %>%
      transmute(RID = as.numeric(RID), death_date = parse_date_safe(AEHDTHDT))
    deaths_rec <- recadv %>%
      filter(AEHDEATH %in% c("Yes", "1", "TRUE") | !is.na(AEHDTHDT)) %>%
      transmute(RID = as.numeric(RID), death_date = parse_date_safe(AEHDTHDT))
    adni_deaths <- bind_rows(deaths_adv, deaths_rec) %>%
      filter(!is.na(RID), !is.na(death_date)) %>%
      group_by(RID) %>%
      summarise(death_date = min(death_date), .groups = "drop")
    adni_fg <- adni %>%
      left_join(adni_deaths, by = "RID") %>%
      mutate(
        anchor_date = as.Date(anchor_date),
        first_event_date = as.Date(first_event_date),
        analysis_end_date = as.Date(analysis_end_date),
        death_after_landmark = !is.na(death_date) & death_date > anchor_date,
        death_competing = death_after_landmark & death_date <= analysis_end_date &
          (event == 0 | is.na(first_event_date) | death_date < first_event_date),
        death_time_years = as.numeric(death_date - anchor_date) / 365.25,
        fg_time = ifelse(death_competing, death_time_years, time_to_event_years),
        fg_status = ifelse(event == 1 & !death_competing, 1, ifelse(death_competing, 2, 0))
      ) %>%
      filter(if_all(all_of(c("fg_time", "fg_status", "pe_worse_z", "Age_z", "Sex_F",
                             "Education_z", "APOE4", "baseline_score_z")), ~ !is.na(.)))
    death_summary[["ADNI"]] <- tibble(
      cohort = "ADNI", n = nrow(adni_fg), clinical_events = sum(adni_fg$fg_status == 1),
      death_records_available = TRUE,
      deaths_after_landmark_any_timing = sum(adni_fg$death_after_landmark, na.rm = TRUE),
      competing_deaths_before_event_or_censor = sum(adni_fg$fg_status == 2),
      competing_death_pct = 100 * competing_deaths_before_event_or_censor / n,
      fine_gray_triggered_by_gt5pct_rule = competing_death_pct > 5,
      note = "ADNI deaths extracted from ADVERSE and RECADV adverse-event death fields; only deaths before clinical event/censoring are competing events."
    )
    if (sum(adni_fg$fg_status == 2) >= 3) {
      finegray_rows[["ADNI"]] <- finegray_table(
        adni_fg,
        ~ pe_worse_z + Age_z + Sex_F + Education_z + APOE4 + baseline_score_z,
        "ADNI", "pe_worse_z"
      )
    }
  }

  death_tab <- bind_rows(death_summary)
  fg_tab <- bind_rows(finegray_rows)
  text <- c(
    "# D5 competing-risk sensitivity",
    "",
    paste0("A4 competing deaths before event/censoring were ",
           death_tab$competing_deaths_before_event_or_censor[death_tab$cohort == "A4"],
           "/", death_tab$n[death_tab$cohort == "A4"], " (",
           fmt(death_tab$competing_death_pct[death_tab$cohort == "A4"], 1), "%)."),
    paste0("ADNI competing deaths before event/censoring were ",
           death_tab$competing_deaths_before_event_or_censor[death_tab$cohort == "ADNI"],
           "/", death_tab$n[death_tab$cohort == "ADNI"], " (",
           fmt(death_tab$competing_death_pct[death_tab$cohort == "ADNI"], 1), "%)."),
    "PREVENT-AD could not be evaluated for death competing risk because the local replication files used here do not contain a death/disposition field.",
    "Because the analyzable competing-death proportions were below 5% in cohorts with records, standard Cox censoring is unlikely to be driven by competing mortality. Exploratory Fine-Gray estimates are provided where death records were available and at least three competing deaths were observed."
  )
  writeLines(text, file.path(rep_dir, "D5_competing_risk_text.md"), useBytes = TRUE)

  files <- c(
    csv_write(death_tab, "D5_competing_death_counts_by_cohort.csv"),
    if (nrow(fg_tab)) csv_write(fg_tab, "D5_fine_gray_sensitivity.csv") else character(),
    file.path(rep_dir, "D5_competing_risk_text.md")
  )
  record_status("D5", "Death counts and Fine-Gray competing risk sensitivity",
                "completed",
                "Death proportions were summarized for all cohorts; Fine-Gray rows are provided where local death records permitted estimation.",
                files, sum(death_tab$n, na.rm = TRUE), sum(death_tab$clinical_events, na.rm = TRUE))
  list(death = death_tab, finegray = fg_tab, files = files)
}, error = function(e) {
  record_status("D5", "Death counts and Fine-Gray competing risk sensitivity",
                "failed", conditionMessage(e))
  NULL
})

# -------------------------------------------------------------------------
# D6. ADNI terminal-window sensitivity
# -------------------------------------------------------------------------

d6_outputs <- tryCatch({
  d6 <- adni %>%
    filter(!is.na(anchor_visit)) %>%
    mutate(anchor_visit = factor(anchor_visit, levels = c("m12", "m24")))
  summ <- d6 %>%
    group_by(anchor_visit) %>%
    summarise(
      n = n(), pct = 100 * n / nrow(d6), events = sum(event),
      PE_mean = mean(pe_worse_z, na.rm = TRUE), PE_sd = sd(pe_worse_z, na.rm = TRUE),
      PE_median = median(pe_worse_z, na.rm = TRUE),
      PE_q1 = quantile(pe_worse_z, 0.25, na.rm = TRUE),
      PE_q3 = quantile(pe_worse_z, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      wilcox_p_m12_vs_m24 = wilcox.test(pe_worse_z ~ anchor_visit, data = d6)$p.value,
      ttest_p_m12_vs_m24 = t.test(pe_worse_z ~ anchor_visit, data = d6)$p.value
    )

  fit_all_data <- d6 %>%
    filter(if_all(all_of(c("time_to_event_years", "event", "pe_worse_z", "Age_z",
                           "Sex_F", "Education_z", "APOE4", "baseline_score_z")), ~ !is.na(.)))
  fit_m12_data <- fit_all_data %>% filter(anchor_visit == "m12")
  fit_all <- coxph(Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F +
                     Education_z + APOE4 + baseline_score_z,
                   data = fit_all_data, x = TRUE, model = TRUE)
  fit_m12 <- coxph(Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F +
                     Education_z + APOE4 + baseline_score_z,
                   data = fit_m12_data, x = TRUE, model = TRUE)
  cox_tab <- bind_rows(
    cox_extract(fit_all, "All ADNI primary risk set", "pe_worse_z",
                n = nrow(fit_all_data), events = sum(fit_all_data$event)),
    cox_extract(fit_m12, "Exclude M24 terminal-window participants", "pe_worse_z",
                n = nrow(fit_m12_data), events = sum(fit_m12_data$event))
  )

  p <- ggplot(d6, aes(anchor_visit, pe_worse_z, fill = anchor_visit)) +
    geom_boxplot(width = 0.52, outlier.shape = NA, alpha = 0.8, linewidth = 0.35) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.28, size = 0.6, colour = "grey25") +
    scale_fill_manual(values = c("m12" = "#6F7C85", "m24" = "#B13F3B")) +
    labs(title = "ADNI PE distribution by terminal early-window visit",
         subtitle = paste0("M12 n=", summ$n[summ$anchor_visit == "m12"],
                           "; M24 n=", summ$n[summ$anchor_visit == "m24"],
                           "; Wilcoxon P=", fmt_p(summ$wilcox_p_m12_vs_m24[1])),
         x = "Terminal visit used for PE estimation", y = "PE worse z-score") +
    theme_pub(7) +
    theme(legend.position = "none")

  files <- c(
    csv_write(summ, "D6_ADNI_terminal_window_PE_distribution.csv"),
    csv_write(cox_tab, "D6_ADNI_exclude_M24_Cox_sensitivity.csv"),
    save_plot(p, "D6_ADNI_terminal_window_PE_distribution", width = 4.3, height = 3.5)
  )
  record_status("D6", "ADNI M12 versus M24 terminal-window sensitivity",
                "completed",
                "PHC_MEM primary flexible-window exposure file supplied anchor_visit; Cox sensitivity excluded M24 terminal-window participants.",
                files, nrow(d6), sum(d6$event))
  list(summary = summ, cox = cox_tab, files = files)
}, error = function(e) {
  record_status("D6", "ADNI M12 versus M24 terminal-window sensitivity",
                "failed", conditionMessage(e))
  NULL
})

# -------------------------------------------------------------------------
# D7. ADNI CSF p-tau181 median stratification
# -------------------------------------------------------------------------

d7_outputs <- tryCatch({
  ptau_median <- median(adni$CSF_PTAU, na.rm = TRUE)
  d7 <- adni %>%
    filter(!is.na(CSF_PTAU), !is.na(CSF_ABETA42)) %>%
    mutate(
      ptau181_stratum = factor(ifelse(CSF_PTAU <= ptau_median, "Lower CSF p-tau181", "Higher CSF p-tau181"),
                               levels = c("Lower CSF p-tau181", "Higher CSF p-tau181")),
      ptau181_stratum_label = paste0(as.character(ptau181_stratum), "\nmedian split at ",
                                     fmt(ptau_median, 2), " pg/mL")
    ) %>%
    group_by(ptau181_stratum) %>%
    mutate(
      PE_lowest_quartile = factor(ifelse(pe_worse_z >= quantile(pe_worse_z, 0.75, na.rm = TRUE),
                                         "Lowest PE quartile", "Other PE"),
                                  levels = c("Other PE", "Lowest PE quartile"))
    ) %>%
    ungroup()
  d7_form <- Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F +
    Education_z + APOE4 + baseline_score_z + CSF_PTAU_z + CSF_ABETA42_z
  d7_cox <- cox_for_strata(d7, "ptau181_stratum", "ptau181_stratum_label", d7_form)
  d7_sensitivity <- stratum_model_sensitivity(
    d7, "ptau181_stratum", "ptau181_stratum_label",
    list(
      unadjusted = Surv(time_to_event_years, event) ~ pe_worse_z,
      demographics = Surv(time_to_event_years, event) ~ pe_worse_z + Age_z + Sex_F +
        Education_z + APOE4,
      demographics_baseline_cognition = Surv(time_to_event_years, event) ~ pe_worse_z +
        Age_z + Sex_F + Education_z + APOE4 + baseline_score_z,
      plus_CSF_ptau181_abeta42 = d7_form
    )
  )
  d7_rates <- d7 %>%
    group_by(ptau181_stratum, PE_lowest_quartile) %>%
    summarise(n = n(), events = sum(event),
              event_probability_5y = km_risk(pick(time_to_event_years, event), 5),
              .groups = "drop") %>%
    mutate(CSF_PTAU_median = ptau_median)
  d7_km <- bind_rows(lapply(levels(d7$ptau181_stratum), function(level_i) {
    d0 <- d7 %>% filter(ptau181_stratum == level_i)
    sf <- survfit(Surv(time_to_event_years, event) ~ PE_lowest_quartile, data = d0)
    km_source(sf, unique(d0$ptau181_stratum_label))
  }))

  files <- c(
    csv_write(d7_cox, "D7_ADNI_CSF_ptau181_median_stratified_PE_Cox.csv"),
    csv_write(d7_sensitivity, "D7_ADNI_CSF_ptau181_median_stratified_model_sensitivity.csv"),
    csv_write(d7_rates, "D7_ADNI_CSF_ptau181_median_stratified_event_rates.csv"),
    source_write(d7_km, "D7_ADNI_CSF_ptau181_median_stratified_KM_source.csv"),
    make_stratified_figure(d7, "ptau181_stratum", "ptau181_stratum_label", d7_cox,
                           "D7_ADNI_CSF_ptau181_median_stratified_PE",
                           "ADNI CSF p-tau181 strata")
  )
  record_status("D7", "ADNI CSF p-tau181 median-stratified PE Cox and KM",
                "completed",
                paste0("ADNI CSF subset split at p-tau181 median ", fmt(ptau_median, 2),
                       "; small-stratum estimates are treated as directionality checks."),
                files, nrow(d7), sum(d7$event))
  list(data = d7, cox = d7_cox, sensitivity = d7_sensitivity, rates = d7_rates, files = files)
}, error = function(e) {
  record_status("D7", "ADNI CSF p-tau181 median-stratified PE Cox and KM",
                "failed", conditionMessage(e))
  NULL
})

# -------------------------------------------------------------------------
# Write status, manifest, and concise report
# -------------------------------------------------------------------------

status_file <- csv_write(status_log, "status_registry.csv")

pkg_status <- tibble(
  package = c("survival", "broom", "ggplot2", "patchwork", "survminer",
              "cmprsk", "metafor", "svglite", "ragg"),
  available = vapply(package, optional_pkg, logical(1)),
  version = vapply(package, function(p) if (optional_pkg(p)) as.character(packageVersion(p)) else NA_character_, character(1))
)
pkg_file <- csv_write(pkg_status, "package_status.csv")

all_files <- list.files(out_dir, recursive = TRUE, full.names = FALSE)
manifest <- c(
  "# Upgrade D1-D7 Analysis Outputs",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Output folder: `", out_dir, "`"),
  "",
  "## Contract",
  "- Manuscript files were not modified.",
  "- Outputs are aggregate tables, figure source curves, figures, logs, and reports.",
  "- Subject-level identifiers were not exported.",
  "",
  "## Main Entry Points",
  "- `run_upgrade_D1_D7.R` -- reproducible R script.",
  "- `tables/status_registry.csv` -- task status and output routing.",
  "- `reports/upgrade_D1_D7_report_20260520.md` -- concise numeric report.",
  "- `logs/sessionInfo.txt` -- R session information.",
  "",
  "## Files",
  paste0("- `", all_files, "`")
)
writeLines(manifest, file.path(out_dir, "_analysis_outputs.md"), useBytes = TRUE)

get_row <- function(obj, idx = 1) if (is.null(obj) || !nrow(obj)) NULL else obj[idx, , drop = FALSE]

report <- c(
  "# D1-D7 upgrade analysis report",
  "",
  "This bundle is analysis-only. It does not modify the manuscript source.",
  "",
  "## D1 A4 amyloid-centiloid median stratification",
  if (!is.null(d1_outputs)) {
    paste0("- Median amyloid centiloid cutoff: ", fmt(median(a4$AMYLCENT, na.rm = TRUE), 2),
           ". Cox rows are in `tables/D1_A4_amyloid_median_stratified_PE_Cox.csv`; figure in `figures/D1_A4_amyloid_median_stratified_PE.*`.")
  } else "- D1 failed; see status registry.",
  if (!is.null(d1_outputs)) paste0(capture.output(print(d1_outputs$cox)), collapse = "\n") else "",
  if (!is.null(d1_outputs)) "Model-sensitivity rows:" else "",
  if (!is.null(d1_outputs)) paste0(capture.output(print(d1_outputs$sensitivity)), collapse = "\n") else "",
  "",
  "## D2 two-step p-tau217 then PE enrichment",
  if (!is.null(d2_outputs)) {
    paste0("- p-tau217 top-quartile cutoff: ", fmt(d2_outputs$summary$pTau217_top_quartile_cutoff[1], 3),
           "; PE lowest-quartile cutoff: ", fmt(d2_outputs$summary$PE_lowest_quartile_cutoff[1], 3),
           ". Summary is in `tables/D2_A4_two_step_ptau217_PE_enrichment.csv`; figure in `figures/D2_A4_two_step_ptau217_PE_enrichment.*`.")
  } else "- D2 failed; see status registry.",
  if (!is.null(d2_outputs)) paste0(capture.output(print(d2_outputs$summary)), collapse = "\n") else "",
  "",
  "## D3 five-year decision curve analysis",
  if (!is.null(d3_outputs)) "- DCA curves and model performance were generated for Base, Base+PE, Base+p-tau217, and Base+p-tau217+PE. See `tables/D3_A4_decision_curve_5y.csv` and `figures/D3_A4_decision_curve_5y.*`." else "- D3 failed; see status registry.",
  if (!is.null(d3_outputs)) paste0(capture.output(print(d3_outputs$perf)), collapse = "\n") else "",
  "",
  "## D4 cross-cohort random-effects meta-analysis",
  if (!is.null(d4_outputs)) {
    paste0("- Pooled HR=", fmt(d4_outputs$summary$pooled_HR, 2),
           " (95% CI ", fmt(d4_outputs$summary$conf.low, 2), "-",
           fmt(d4_outputs$summary$conf.high, 2), "), I2=",
           fmt(d4_outputs$summary$I2, 1), "%. See `tables/D4_cross_cohort_random_effects_meta_summary.csv` and `figures/D4_cross_cohort_random_effects_meta_forest.*`.")
  } else "- D4 failed; see status registry.",
  if (!is.null(d4_outputs)) paste0(capture.output(print(d4_outputs$input)), collapse = "\n") else "",
  "",
  "## D5 competing-risk sensitivity",
  if (!is.null(d5_outputs)) "- Death-count table and available Fine-Gray sensitivity rows are in `tables/D5_competing_death_counts_by_cohort.csv` and `tables/D5_fine_gray_sensitivity.csv`; prose is in `reports/D5_competing_risk_text.md`." else "- D5 failed; see status registry.",
  if (!is.null(d5_outputs)) paste0(capture.output(print(d5_outputs$death)), collapse = "\n") else "",
  "",
  "## D6 ADNI terminal-window sensitivity",
  if (!is.null(d6_outputs)) "- M12/M24 proportions, PE distribution comparison, and M24-excluded Cox sensitivity were generated. See `tables/D6_ADNI_terminal_window_PE_distribution.csv` and `tables/D6_ADNI_exclude_M24_Cox_sensitivity.csv`." else "- D6 failed; see status registry.",
  if (!is.null(d6_outputs)) paste0(capture.output(print(d6_outputs$summary)), collapse = "\n") else "",
  if (!is.null(d6_outputs)) paste0(capture.output(print(d6_outputs$cox)), collapse = "\n") else "",
  "",
  "## D7 ADNI CSF p-tau181 median stratification",
  if (!is.null(d7_outputs)) {
    paste0("- CSF p-tau181 median cutoff: ", fmt(median(adni$CSF_PTAU, na.rm = TRUE), 2),
           ". Cox rows are in `tables/D7_ADNI_CSF_ptau181_median_stratified_PE_Cox.csv`; figure in `figures/D7_ADNI_CSF_ptau181_median_stratified_PE.*`.")
  } else "- D7 failed; see status registry.",
  if (!is.null(d7_outputs)) paste0(capture.output(print(d7_outputs$cox)), collapse = "\n") else "",
  if (!is.null(d7_outputs)) "Model-sensitivity rows:" else "",
  if (!is.null(d7_outputs)) paste0(capture.output(print(d7_outputs$sensitivity)), collapse = "\n") else "",
  "",
  "## Status",
  paste0(capture.output(print(status_log)), collapse = "\n")
)
writeLines(report, file.path(rep_dir, "upgrade_D1_D7_report_20260520.md"), useBytes = TRUE)

cat("Completed D1-D7 bundle.\n")
cat("Status:", status_file, "\n")
cat("Manifest:", file.path(out_dir, "_analysis_outputs.md"), "\n")
