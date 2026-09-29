

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


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(broom)
})

pkg_dir <- pkg_root
clinical_dir <- file.path(data_root, "a4")
locked_dir <- file.path(pkg_dir, "outputs", "locked", "a4")
out_dir <- file.path(pkg_dir, "verify", "a4", "tables_and_sensitivity")
tex_dir <- file.path(pkg_dir, "verify", "a4", "tex")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tex_dir, recursive = TRUE, showWarnings = FALSE)

pe_levels <- c("High_PE", "Moderate_PE", "Low_PE")
pe_labels <- c(High_PE = "High PE", Moderate_PE = "Moderate PE", Low_PE = "Low PE")
tex_row <- "\\\\"

analysis <- read_csv(file.path(locked_dir, "analysis_baseline_w1_15_strict_preW48_lmer.csv"), show_col_types = FALSE) %>%
  mutate(LC_group_lmer = factor(LC_group_lmer, levels = pe_levels))
learning <- read_csv(file.path(locked_dir, "learningcurve_PE_lmer_w1_15_strict_preW48.csv"), show_col_types = FALSE) %>%
  select(BID, learning_slope_resid_lmer)
analysis <- analysis %>% left_join(learning, by = "BID")
surv_primary <- readRDS(file.path(locked_dir, "surv_primary_complete_w1_15_strict_preW48_lmer.rds")) %>%
  mutate(LC_group_lmer = factor(LC_group_lmer, levels = pe_levels))
surv_group <- readRDS(file.path(locked_dir, "surv_group_complete_w1_15_strict_preW48_lmer.rds")) %>%
  mutate(LC_group_lmer = factor(LC_group_lmer, levels = pe_levels))

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}

fmt_p <- function(p) {
  if (is.na(p)) return("--")
  if (p < 0.001) return("$<0.001$")
  formatC(p, format = "f", digits = 3)
}

fmt_mean_sd <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  if (!length(x)) return("--")
  paste0(fmt_num(mean(x), digits), " (", fmt_num(sd(x), digits), ")")
}

fmt_pct <- function(x, denom = length(x), digits = 1) {
  x <- x[!is.na(x)]
  if (!denom || !length(x)) return("--")
  fmt_num(mean(x) * 100, digits)
}

anova_p <- function(dat, value, group = "LC_group_lmer") {
  dat <- dat %>% filter(!is.na(.data[[value]]), !is.na(.data[[group]]))
  if (n_distinct(dat[[group]]) < 2) return(NA_real_)
  summary(aov(dat[[value]] ~ dat[[group]]))[[1]][["Pr(>F)"]][1]
}

chisq_p <- function(dat, value, group = "LC_group_lmer") {
  dat <- dat %>% filter(!is.na(.data[[value]]), !is.na(.data[[group]]))
  if (n_distinct(dat[[group]]) < 2 || n_distinct(dat[[value]]) < 2) return(NA_real_)
  suppressWarnings(chisq.test(table(dat[[value]], dat[[group]]))$p.value)
}

group_counts <- function(dat) {
  table(factor(dat$LC_group_lmer, levels = pe_levels))
}

summ_by_group <- function(dat, value, type, digits = 2) {
  sapply(pe_levels, function(g) {
    x <- dat[[value]][dat$LC_group_lmer == g]
    if (type == "mean_sd") fmt_mean_sd(x, digits)
    else if (type == "female") fmt_pct(x == "Female", length(x), 1)
    else if (type == "apoe") fmt_pct(x == 1, length(x), 1)
    else stop("unknown type")
  })
}

make_group_table <- function(dat) {
  dat <- dat %>% filter(!is.na(LC_group_lmer))
  rows <- tibble(
    characteristic = c("Age, years", "Female, \\%", "Education, years", "APOE $\\varepsilon4$+, \\%",
                       "Amyloid PET SUVR", "Plasma \\textit{p}-tau217$^\\dag$, pg/mL",
                       "Baseline Cogstate ($z$-score)", "Baseline PACC"),
    high = c(summ_by_group(dat, "Age", "mean_sd", 2)[1], summ_by_group(dat, "Sex", "female", 1)[1],
             summ_by_group(dat, "Education", "mean_sd", 2)[1], summ_by_group(dat, "APOE4", "apoe", 1)[1],
             summ_by_group(dat, "SUVR", "mean_sd", 3)[1], summ_by_group(dat, "pTau217", "mean_sd", 3)[1],
             summ_by_group(dat, "baseline_C3_z", "mean_sd", 2)[1], summ_by_group(dat, "baseline_pacc", "mean_sd", 2)[1]),
    moderate = c(summ_by_group(dat, "Age", "mean_sd", 2)[2], summ_by_group(dat, "Sex", "female", 1)[2],
                 summ_by_group(dat, "Education", "mean_sd", 2)[2], summ_by_group(dat, "APOE4", "apoe", 1)[2],
                 summ_by_group(dat, "SUVR", "mean_sd", 3)[2], summ_by_group(dat, "pTau217", "mean_sd", 3)[2],
                 summ_by_group(dat, "baseline_C3_z", "mean_sd", 2)[2], summ_by_group(dat, "baseline_pacc", "mean_sd", 2)[2]),
    low = c(summ_by_group(dat, "Age", "mean_sd", 2)[3], summ_by_group(dat, "Sex", "female", 1)[3],
            summ_by_group(dat, "Education", "mean_sd", 2)[3], summ_by_group(dat, "APOE4", "apoe", 1)[3],
            summ_by_group(dat, "SUVR", "mean_sd", 3)[3], summ_by_group(dat, "pTau217", "mean_sd", 3)[3],
            summ_by_group(dat, "baseline_C3_z", "mean_sd", 2)[3], summ_by_group(dat, "baseline_pacc", "mean_sd", 2)[3]),
    p = c(anova_p(dat, "Age"), chisq_p(dat, "Sex"), anova_p(dat, "Education"), chisq_p(dat, "APOE4"),
          anova_p(dat, "SUVR"), anova_p(dat, "pTau217"), anova_p(dat, "baseline_C3_z"), anova_p(dat, "baseline_pacc"))
  ) %>% mutate(p_display = vapply(p, fmt_p, character(1)))
  rows
}

write_table1_tex <- function(rows, dat) {
  counts <- group_counts(dat %>% filter(!is.na(LC_group_lmer)))
  ptau_counts <- sapply(pe_levels, function(g) sum(!is.na(dat$pTau217) & dat$LC_group_lmer == g, na.rm = TRUE))
  cut_pct <- quantile(dat$learning_slope_resid_z_lmer[!is.na(dat$learning_slope_resid_z_lmer)], probs = c(0.2, 0.8))
  con <- file(file.path(tex_dir, "a4_table1_corrected.tex"), "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "\\begin{tabularx}{\\textwidth}{>{\\raggedright\\arraybackslash}p{0.28\\textwidth} >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}p{0.10\\textwidth}}",
    "\\toprule",
    paste0(sprintf("Characteristic & High PE ($n = %d$) & Moderate PE ($n = %d$) & Low PE ($n = %d$) & $p$-value", counts["High_PE"], counts["Moderate_PE"], counts["Low_PE"]), " ", tex_row),
    "\\midrule"
  ), con)
  for (i in seq_len(nrow(rows))) {
    writeLines(paste0(sprintf("%s & %s & %s & %s & %s", rows$characteristic[i], rows$high[i], rows$moderate[i], rows$low[i], rows$p_display[i]), " ", tex_row), con)
  }
  writeLines(c(
    "\\bottomrule",
    "\\end{tabularx}",
    "\\captionsetup{font=small}",
    sprintf("\\caption*{\\textit{Note:} Values are mean (SD) or \\%%. PE groups defined by residualized $z$-score: High ($z > 0.2$), Moderate ($-0.2 \\leq z \\leq 0.2$), Low ($z < -0.2$). In the corrected V1/V3/V9/V15 PE distribution, these cutpoints correspond to the %.1fth and %.1fth empirical percentiles, respectively. $p$-values are from one-way ANOVA for continuous variables and Pearson $\\chi^2$ tests for categorical variables. $^\\dag$Plasma \\textit{p}-tau217 was available in $n = %d$, %d, and %d for High, Moderate, and Low PE groups, respectively.}", mean(dat$learning_slope_resid_z_lmer <= -0.2, na.rm = TRUE) * 100, mean(dat$learning_slope_resid_z_lmer <= 0.2, na.rm = TRUE) * 100, ptau_counts[1], ptau_counts[2], ptau_counts[3])
  ), con)
}

make_two_col_table <- function(dat, included_ids, caption, label, file, note) {
  dat <- dat %>% filter(!is.na(LC_group_lmer)) %>% mutate(included = BID %in% included_ids)
  groups <- list(Included = dat %>% filter(included), Excluded = dat %>% filter(!included))
  row_val <- function(d, value, type, digits = 2) {
    if (type == "mean_sd") fmt_mean_sd(d[[value]], digits)
    else if (type == "female") fmt_pct(d[[value]] == "Female", nrow(d), 1)
    else if (type == "apoe") fmt_pct(d[[value]] == 1, nrow(d), 1)
  }
  p_for <- function(value, type) if (type == "mean_sd") {
    summary(aov(dat[[value]] ~ dat$included))[[1]][["Pr(>F)"]][1]
  } else {
    suppressWarnings(chisq.test(table(dat[[value]], dat$included))$p.value)
  }
  specs <- tibble(
    characteristic = c("Age, years", "Female, \\%", "Education, years", "APOE $\\varepsilon4$+, \\%",
                       "Amyloid PET SUVR", "Plasma \\textit{p}-tau217$^\\dag$, pg/mL",
                       "Baseline Cogstate ($z$-score)", "Baseline PACC"),
    value = c("Age", "Sex", "Education", "APOE4", "SUVR", "pTau217", "baseline_C3_z", "baseline_pacc"),
    type = c("mean_sd", "female", "mean_sd", "apoe", "mean_sd", "mean_sd", "mean_sd", "mean_sd"),
    digits = c(2, 1, 2, 1, 3, 3, 2, 2)
  )
  rows <- specs %>% rowwise() %>% mutate(
    included_val = row_val(groups$Included, value, type, digits),
    excluded_val = row_val(groups$Excluded, value, type, digits),
    p_display = fmt_p(p_for(value, type))
  ) %>% ungroup()
  write_csv(rows, file.path(out_dir, sub("\\.tex$", ".csv", file)))
  con <- file(file.path(tex_dir, file), "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    sprintf("\\caption*{\\textbf{%s}}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{tabularx}{\\textwidth}{>{\\raggedright\\arraybackslash}p{0.34\\textwidth} >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}p{0.12\\textwidth}}",
    "\\toprule",
    paste0(sprintf("Characteristic & Included ($n = %d$) & Excluded ($n = %d$) & $p$-value", nrow(groups$Included), nrow(groups$Excluded)), " ", tex_row),
    "\\midrule"
  ), con)
  for (i in seq_len(nrow(rows))) writeLines(paste0(sprintf("%s & %s & %s & %s", rows$characteristic[i], rows$included_val[i], rows$excluded_val[i], rows$p_display[i]), " ", tex_row), con)
  writeLines(c("\\bottomrule", "\\end{tabularx}", "\\captionsetup{font=small}", sprintf("\\caption*{\\textit{Note:} %s}", note), "\\end{table}"), con)
}

model_group <- coxph(Surv(time_to_event_years, converted) ~ LC_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z, data = surv_primary)
model_group_pacc <- coxph(Surv(time_to_event_years, converted) ~ LC_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z + scale(baseline_pacc), data = surv_primary)
model_cont <- coxph(Surv(time_to_event_years, converted) ~ learning_slope_resid_z_decr_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z, data = surv_primary)
model_cont_pacc <- coxph(Surv(time_to_event_years, converted) ~ learning_slope_resid_z_decr_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z + scale(baseline_pacc), data = surv_primary)
base_model <- coxph(Surv(time_to_event_years, converted) ~ Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z, data = surv_primary)

extract_terms <- function(fit, terms, label) {
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(term %in% terms) %>%
    transmute(analysis = label, term, HR = estimate, conf.low, conf.high, p.value)
}

s3 <- bind_rows(
  extract_terms(model_group_pacc, c("LC_group_lmerModerate_PE", "LC_group_lmerLow_PE"), "Grouped PE + baseline PACC"),
  extract_terms(model_cont_pacc, "learning_slope_resid_z_decr_lmer", "Continuous PE + baseline PACC")
)
write_csv(s3, file.path(out_dir, "a4_s3_pacc_adjusted_cox.csv"))

tertile_dat <- surv_primary %>%
  mutate(LC_tertile = ntile(learning_slope_resid_z_lmer, 3),
         LC_tertile = factor(case_when(LC_tertile == 1 ~ "Low tertile", LC_tertile == 2 ~ "Middle tertile", TRUE ~ "High tertile"),
                             levels = c("High tertile", "Middle tertile", "Low tertile")))
model_tertile <- coxph(Surv(time_to_event_years, converted) ~ LC_tertile + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z, data = tertile_dat)
s4 <- bind_rows(
  extract_terms(model_group, c("LC_group_lmerModerate_PE", "LC_group_lmerLow_PE"), "Primary threshold groups"),
  tidy(model_tertile, exponentiate = TRUE, conf.int = TRUE) %>% filter(grepl("LC_tertile", term)) %>% transmute(analysis = "Data-driven tertiles", term, HR = estimate, conf.low, conf.high, p.value)
)
write_csv(s4, file.path(out_dir, "a4_s4_threshold_vs_tertile.csv"))

cog <- read_csv(file.path(clinical_dir, "External Data", "cogstate_battery.csv"), show_col_types = FALSE, progress = FALSE) %>%
  mutate(VISCODE_num = suppressWarnings(as.numeric(as.character(VISCODE)))) %>%
  filter(VISCODE_num %in% c(1, 3, 9, 15), !is.na(C3Comp), !is.na(TDate_DAYS_T0)) %>%
  semi_join(analysis %>% filter(!is.na(LC_group_lmer)), by = "BID") %>%
  group_by(BID, VISCODE_num) %>%
  summarise(C3Comp = mean(C3Comp, na.rm = TRUE), TDate_DAYS_T0 = mean(TDate_DAYS_T0, na.rm = TRUE), .groups = "drop") %>%
  group_by(BID) %>%
  mutate(time_years = (TDate_DAYS_T0 - min(TDate_DAYS_T0, na.rm = TRUE)) / 365.25) %>%
  filter(n() >= 2) %>%
  summarise(ols_slope = coef(lm(C3Comp ~ time_years))[2], .groups = "drop")
ols_pe <- analysis %>%
  filter(!is.na(LC_group_lmer)) %>%
  inner_join(cog, by = "BID") %>%
  { mutate(.,
           ols_resid = resid(lm(ols_slope ~ baseline_C3Comp, data = .)),
           ols_resid_z = as.numeric(scale(ols_resid)),
           ols_resid_z_decr = -ols_resid_z,
           ols_group = factor(case_when(
             ols_resid_z > 0.2 ~ "High_PE",
             ols_resid_z < -0.2 ~ "Low_PE",
             TRUE ~ "Moderate_PE"
           ), levels = pe_levels)) }
ols_surv <- surv_primary %>% inner_join(ols_pe %>% select(BID, ols_resid_z, ols_resid_z_decr, ols_group), by = "BID")
ols_group_fit <- coxph(Surv(time_to_event_years, converted) ~ ols_group + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z, data = ols_surv)
ols_cont_fit <- coxph(Surv(time_to_event_years, converted) ~ ols_resid_z_decr + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z, data = ols_surv)
s5 <- tibble(
  metric = c("PE-estimable subjects with OLS and BLUP slopes", "SD(BLUP residualized slope) / SD(OLS residualized slope)", "Correlation between OLS and BLUP residualized slopes"),
  value = c(nrow(ols_pe), sd(ols_pe$learning_slope_resid_lmer, na.rm = TRUE) / sd(ols_pe$ols_resid, na.rm = TRUE), cor(ols_pe$ols_resid_z, ols_pe$learning_slope_resid_z_lmer, use = "complete.obs"))
)
s5_cox <- bind_rows(
  tidy(ols_group_fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(term %in% c("ols_groupModerate_PE", "ols_groupLow_PE")) %>% transmute(analysis = "OLS grouped PE", term, HR = estimate, conf.low, conf.high, p.value),
  tidy(ols_cont_fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "ols_resid_z_decr") %>% transmute(analysis = "OLS continuous PE", term, HR = estimate, conf.low, conf.high, p.value)
)
write_csv(s5, file.path(out_dir, "a4_s5_ols_blup_metrics.csv"))
write_csv(s5_cox, file.path(out_dir, "a4_s5_ols_blup_cox.csv"))

ph_group <- cox.zph(model_group)
ph_cont <- cox.zph(model_cont)
s_ph <- tibble(
  model = c("Grouped PE", "Continuous PE"),
  pe_p = c(ph_group$table["LC_group_lmer", "p"], ph_cont$table["learning_slope_resid_z_decr_lmer", "p"]),
  global_p = c(ph_group$table["GLOBAL", "p"], ph_cont$table["GLOBAL", "p"])
)
write_csv(s_ph, file.path(out_dir, "a4_global_schoenfeld.csv"))

bootstrap <- read_csv(file.path(locked_dir, "survival_bootstrap_delta_c_w1_15_strict_preW48_lmer.csv"), show_col_types = FALSE)
cv <- read_csv(file.path(locked_dir, "survival_cv_delta_c_w1_15_strict_preW48_lmer.csv"), show_col_types = FALSE)
lrt_p <- anova(base_model, model_cont, test = "LRT")$`Pr(>|Chi|)`[2]
s6 <- tibble(
  validation = c("Apparent", "Bootstrap paired resampling", "Nested 5-fold CV, 10 repeats"),
  delta_c = c(concordance(model_cont)$concordance - concordance(base_model)$concordance, bootstrap$mean_dC[1], cv$mean_dC[1]),
  lower = c(NA_real_, bootstrap$ci_lo[1], cv$min_dC[1]),
  upper = c(NA_real_, bootstrap$ci_hi[1], cv$max_dC[1]),
  p.value = c(lrt_p, NA_real_, NA_real_)
)
write_csv(s6, file.path(out_dir, "a4_s6_incremental_discrimination.csv"))

sens_exclude_0_5 <- surv_primary %>% filter(!(converted & time_to_event_years <= 0.5))
sens_exclude_1_0 <- surv_primary %>% filter(!(converted & time_to_event_years <= 1.0))
fit_sens <- function(dat, label) {
  fit <- coxph(Surv(time_to_event_years, converted) ~ LC_group_lmer + Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z, data = dat)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "LC_group_lmerLow_PE") %>%
    transmute(analysis = label, n = nrow(dat), events = sum(dat$converted), HR = estimate, conf.low, conf.high, p.value)
}
s9 <- bind_rows(
  fit_sens(surv_primary, "Primary post-W48 risk set"),
  fit_sens(sens_exclude_0_5, "Exclude events within 0.5 y"),
  fit_sens(sens_exclude_1_0, "Exclude events within 1.0 y")
)
write_csv(s9, file.path(out_dir, "a4_s9_event_timing_sensitivity.csv"))

write_model_tex <- function(df, file, caption, label, note) {
  con <- file(file.path(tex_dir, file), "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(c("\\begin{table}[H]", "\\centering", "\\small", sprintf("\\caption*{\\textbf{%s}}", caption), sprintf("\\label{%s}", label),
               "\\begin{tabularx}{\\textwidth}{>{\\raggedright\\arraybackslash}p{0.32\\textwidth} >{\\raggedright\\arraybackslash}p{0.26\\textwidth} >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}p{0.12\\textwidth}}",
               "\\toprule", paste0("Analysis & Term & HR & 95\\% CI & $p$-value ", tex_row), "\\midrule"), con)
  for (i in seq_len(nrow(df))) {
    term <- gsub("LC_group_lmer", "", df$term[i]); term <- gsub("LC_tertile", "", term); term <- gsub("learning_slope_resid_z_decr_lmer", "Per-SD decrease", term)
    term <- gsub("_", " ", term)
    writeLines(paste0(sprintf("%s & %s & %.2f & %.2f to %.2f & %s", df$analysis[i], term, df$HR[i], df$conf.low[i], df$conf.high[i], fmt_p(df$p.value[i])), " ", tex_row), con)
  }
  writeLines(c("\\bottomrule", "\\end{tabularx}", "\\captionsetup{font=small}", sprintf("\\caption*{\\textit{Note:} %s}", note), "\\end{table}"), con)
}

write_model_tex(s3, "a4_s3_corrected.tex", "Supplementary Table S3. Cox model results for confirmed CDR progression with additional adjustment for baseline PACC.", "tab:pacc_adjusted", "Models use the corrected V1/V3/V9/V15 PE derivation and real-time post-Week-48 CDR endpoint. The primary covariate set is additionally adjusted for baseline PACC.")
write_model_tex(s4, "a4_s4_corrected.tex", "Supplementary Table S4. Cox proportional hazards results comparing the primary threshold-based PE groups with data-driven tertile groups.", "tab:threshold_tertile", "All models use the corrected biomarker-complete primary Cox sample. High PE is the reference group.")

con <- file(file.path(tex_dir, "a4_s5_corrected.tex"), "w", encoding = "UTF-8")
writeLines(c("\\begin{table}[H]", "\\centering", "\\small", "\\caption*{\\textbf{Supplementary Table S5.} Comparison of per-subject OLS and BLUP (empirical Bayes) PE slope estimates after correction of the PE window.}", "\\label{tab:ols_blup}",
             "\\begin{tabularx}{\\textwidth}{>{\\raggedright\\arraybackslash}p{0.58\\textwidth} >{\\centering\\arraybackslash}X}", "\\toprule", paste0("Metric & Value ", tex_row), "\\midrule"), con)
for (i in seq_len(nrow(s5))) {
  value_display <- if (grepl("subjects", s5$metric[i], ignore.case = TRUE)) {
    sprintf("%d", as.integer(round(s5$value[i])))
  } else {
    sprintf("%.3f", s5$value[i])
  }
  writeLines(paste0(sprintf("%s & %s", s5$metric[i], value_display), " ", tex_row), con)
}
writeLines(c("\\midrule", paste0("Analysis & HR (95\\% CI), $p$-value ", tex_row)), con)
for (i in seq_len(nrow(s5_cox))) {
  term <- gsub("ols_group", "", s5_cox$term[i]); term <- gsub("ols_resid_z_decr", "Per-SD decrease", term)
  term <- gsub("_", " ", term)
  writeLines(paste0(sprintf("%s: %s & %.2f (%.2f to %.2f), %s", s5_cox$analysis[i], term, s5_cox$HR[i], s5_cox$conf.low[i], s5_cox$conf.high[i], fmt_p(s5_cox$p.value[i])), " ", tex_row), con)
}
writeLines(c("\\bottomrule", "\\end{tabularx}", "\\captionsetup{font=small}", "\\caption*{\\textit{Note:} OLS PE slopes were estimated within subject across V1/V3/V9/V15, residualized against baseline Cogstate, and standardized using the same direction as the BLUP metric. Cox models use the corrected biomarker-complete primary risk set and the same covariate set as the primary A4 model.}", "\\end{table}"), con)
close(con)

con <- file(file.path(tex_dir, "a4_s6_corrected.tex"), "w", encoding = "UTF-8")
writeLines(c("\\begin{table}[H]", "\\centering", "\\small", "\\caption*{\\textbf{Supplementary Table S6.} Incremental discrimination for the A4 clinical endpoint after correction of the time axis.}", "\\label{tab:nested_cv}",
             "\\begin{tabularx}{\\textwidth}{>{\\raggedright\\arraybackslash}p{0.42\\textwidth} >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}p{0.12\\textwidth}}", "\\toprule", paste0("Validation approach & $\\Delta$C-index & Interval & $p$-value ", tex_row), "\\midrule"), con)
fmt_s6_interval <- function(validation, lower, upper) {
  if (is.na(lower)) return("--")
  if (grepl("Bootstrap", validation)) return(sprintf("%.5f to %.3f", lower, upper))
  if (grepl("Nested", validation)) return(sprintf("%.3f to %.3f", lower, ceiling(upper * 1000) / 1000))
  sprintf("%.3f to %.3f", lower, upper)
}
for (i in seq_len(nrow(s6))) writeLines(paste0(sprintf("%s & %.3f & %s & %s", s6$validation[i], s6$delta_c[i], fmt_s6_interval(s6$validation[i], s6$lower[i], s6$upper[i]), fmt_p(s6$p.value[i])), " ", tex_row), con)
writeLines(c("\\bottomrule", "\\end{tabularx}", "\\captionsetup{font=small}", "\\caption*{\\textit{Note:} The apparent row compares the base model with and without continuous PE using the corrected primary Cox sample. Bootstrap values are from 1,000 paired resamples. The nested-CV row re-standardizes PE within each training fold.}", "\\end{table}"), con)
close(con)

con <- file(file.path(tex_dir, "a4_s9_corrected.tex"), "w", encoding = "UTF-8")
writeLines(c("\\begin{table}[H]", "\\centering", "\\small", "\\caption*{\\textbf{Supplementary Table S9.} Event-timing sensitivity analyses for the A4 clinical endpoint.}", "\\label{tab:a4_event_timing}",
             "\\begin{tabularx}{\\textwidth}{>{\\raggedright\\arraybackslash}p{0.36\\textwidth} >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}X >{\\centering\\arraybackslash}p{0.12\\textwidth}}", "\\toprule", paste0("Analysis & $n$ & Events & HR & 95\\% CI & $p$-value ", tex_row), "\\midrule"), con)
for (i in seq_len(nrow(s9))) writeLines(paste0(sprintf("%s & %d & %d & %.2f & %.2f to %.2f & %s", s9$analysis[i], s9$n[i], s9$events[i], s9$HR[i], s9$conf.low[i], s9$conf.high[i], fmt_p(s9$p.value[i])), " ", tex_row), con)
writeLines(c("\\bottomrule", "\\end{tabularx}", "\\captionsetup{font=small}", "\\caption*{\\textit{Note:} All models report the Low vs High PE contrast from the corrected biomarker-complete Cox model. Early-event exclusions are based on real time from the per-subject Week-48 landmark. Results are identical across rows because no confirmed events in the primary risk set occurred within 0.5 or 1.0 years of the Week-48 landmark.}", "\\end{table}"), con)
close(con)

table1 <- make_group_table(analysis)
write_csv(table1, file.path(out_dir, "a4_table1_corrected.csv"))
write_table1_tex(table1, analysis)

long_ids <- analysis %>% filter(!is.na(LC_group_lmer), !is.na(SUVR_z), !is.na(pTau217_z)) %>% pull(BID)
make_two_col_table(analysis, long_ids,
                   "Supplementary Table S1. Baseline characteristics: biomarker-complete longitudinal sample ($n = 871$) versus excluded PE-estimable participants ($n = 127$).",
                   "tab:s1_longitudinal_baseline", "a4_s1_corrected.tex",
                   "Values are mean (SD) or \\%. The included group is the corrected biomarker-complete longitudinal PACC analysis sample; the excluded group consists of PE-estimable participants lacking complete biomarker covariates. $^\\dag$Plasma \\textit{p}-tau217 is unavailable by definition in the excluded group because complete biomarker covariates were required for inclusion.")

surv_ids <- surv_primary$BID
make_two_col_table(analysis, surv_ids,
                   "Supplementary Table S2. Baseline characteristics: biomarker-complete primary Cox sample ($n = 719$) versus excluded PE-estimable participants ($n = 279$).",
                   "tab:s2_survival_baseline", "a4_s2_corrected.tex",
                   "Values are mean (SD) or \\%. The included group is the corrected primary Cox risk set after excluding V18-prevalent participants and requiring complete biomarker covariates; V18-prevalent participants are described separately in Supplementary Table S12. $^\\dag$The excluded-group plasma \\textit{p}-tau217 value is calculated only among excluded participants with available plasma data, primarily the V18-prevalent biomarker-complete subset, and should not be interpreted as representing all 279 excluded participants.")

write_csv(s_ph, file.path(out_dir, "a4_ph_summary_for_text.csv"))
cat("Generated corrected tables and sensitivity outputs in:\n", out_dir, "\n", tex_dir, "\n")

