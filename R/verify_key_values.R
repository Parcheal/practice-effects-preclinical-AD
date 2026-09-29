#!/usr/bin/env Rscript
# Lightweight key-value checklist for the locked package outputs.

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg)) script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
root <- if (!is.na(script_path)) normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)

expected <- data.frame(
  metric = c(
    "combined_n",
    "combined_age",
    "female_pct",
    "pooled_base_hr",
    "base_hr_by_cohort",
    "a4_biomarker_adjusted_hr",
    "adni_csf_adjusted_hr",
    "sequential_enrichment"
  ),
  expected = c(
    "1791 = 998(A4)+322(PREVENT-AD)+471(ADNI)",
    "71.0 [6.3]",
    "48.7%",
    "1.40 (95% CI 1.25-1.57), I2=0%",
    "A4 1.32; PREVENT-AD 1.61; ADNI 1.46",
    "1.24 (1.07-1.43)",
    "1.57 (1.22-2.02)",
    "p-tau217 top quartile 48.2% -> +PE 63.5%"
  ),
  stringsAsFactors = FALSE
)

dir.create(file.path(root, "verify"), recursive = TRUE, showWarnings = FALSE)
write.csv(expected, file.path(root, "verify", "expected_key_values.csv"), row.names = FALSE, fileEncoding = "UTF-8")
print(expected)
cat("\nKey-value checklist written to verify/expected_key_values.csv.\n")
