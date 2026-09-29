

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


# C2: A4 PE x covariate interaction tests (CDR Cox, primary biomarker-adjusted sample)
# Uses the saved primary Cox risk set (n = 719, 175 events).
# Tests PE x {APOE4 carrier, Sex, Age, baseline Cogstate} via LRT.

suppressPackageStartupMessages({
  library(dplyr); library(survival); library(broom); library(readr); library(tibble)
})

root <- pkg_root
rds  <- file.path(root, "outputs", "locked", "a4", "surv_primary_complete_w1_15_strict_preW48_lmer.rds")
out  <- file.path(root, "verify", "a4", "c2_pe_covariate_interactions.csv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

d <- readRDS(rds)
cat("Loaded A4 primary Cox sample: n =", nrow(d), ", events =", sum(d$converted), "\n")

base_form <- Surv(time_to_event_years, converted) ~ learning_slope_resid_z_decr_lmer +
  Age_z + Sex_F + Education_z + APOE4_carrier + SUVR_z + pTau217_z + baseline_C3_z

fit_base <- coxph(base_form, data = d)

run_lrt <- function(label, int_var) {
  # Build a formula that adds PE:int_var as an interaction term
  fml_int <- update(base_form, as.formula(paste0(". ~ . + learning_slope_resid_z_decr_lmer:", int_var)))
  fit_int <- coxph(fml_int, data = d)
  lrt <- anova(fit_base, fit_int)
  term_row <- tidy(fit_int, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl(paste0("learning_slope_resid_z_decr_lmer:", int_var), term))
  tibble(
    interaction = label,
    n = nrow(d),
    events = sum(d$converted),
    lrt_chisq = lrt[["Chisq"]][2],
    lrt_df    = lrt[["Df"]][2],
    lrt_p     = lrt[["Pr(>|Chi|)"]][2],
    hr_interaction       = if (nrow(term_row)) term_row$estimate[1] else NA,
    hr_interaction_lo    = if (nrow(term_row)) term_row$conf.low[1] else NA,
    hr_interaction_hi    = if (nrow(term_row)) term_row$conf.high[1] else NA,
    wald_p_interaction   = if (nrow(term_row)) term_row$p.value[1] else NA
  )
}

res <- bind_rows(
  run_lrt("PE x APOE4_carrier", "APOE4_carrier"),
  run_lrt("PE x Sex_F",         "Sex_F"),
  run_lrt("PE x Age_z",         "Age_z"),
  run_lrt("PE x baseline_C3_z", "baseline_C3_z")
)

write_csv(res, out)
cat("\n=== PE x covariate interaction results (A4 primary Cox, n = 719, 175 events) ===\n")
print(res, width = 200)
cat("\nWritten to:", out, "\n")
