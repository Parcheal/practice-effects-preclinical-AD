#!/usr/bin/env Rscript
# Package orchestration script.
#
# This script documents the dependency order and can source the retained
# scripts when PE_RUN_ANALYSES=true. Default mode is a dry run because several
# stages are long-running. Rerun outputs are configured to land under verify/.

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)
args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg)) script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
if (!is.na(script_path)) {
  setwd(normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE))
}
source(file.path("R", "config.R"))
setwd(pkg_root)

message("PE analysis code package")
message("Package root: ", pkg_root)
message("Locked outputs are preserved under outputs/locked/")
message("Use verify/ for rerun outputs and lightweight checks. See README.md and PACKAGING_NOTES.md.")

steps <- data.frame(
  label = c(
    "A4 primary corrected analysis",
    "A4 table and sensitivity outputs",
    "A4 covariate interactions",
    "PREVENT-AD primary replication",
    "PREVENT-AD bootstrap/CV supplement",
    "ADNI v3 primary replication",
    "ADNI supplement phase 1",
    "ADNI supplement phase 2",
    "A4 incremental analyses",
    "Upgrade/meta/DCA/enrichment analyses"
  ),
  path = c(
    "R/a4/01_run_corrected_analysis.R",
    "R/a4/02_generate_corrected_tables_and_sensitivity.R",
    "R/a4/03_run_c2_pe_covariate_interactions.R",
    "R/preventad/01_run_prevent_ad_replication.R",
    "R/preventad/02_run_prevent_ad_supplementary.R",
    "R/adni/01_run_adni_replication_v3.R",
    "R/adni/02_run_adni_replication_v3_supplement.R",
    "R/adni/03_run_adni_replication_v3_phase2.R",
    "R/a4/05_run_incremental_full.R",
    "R/common/06_run_upgrade_D1_D7.R"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(steps))) {
  cat(sprintf("%d. %s: %s\n", i, steps$label[i], steps$path[i]))
}

run_analyses <- identical(tolower(Sys.getenv("PE_RUN_ANALYSES", "false")), "true")
if (!run_analyses) {
  cat("\nDry run only. Set PE_RUN_ANALYSES=true to source these R scripts in order; keep outputs under verify/.\n")
} else {
  for (i in seq_len(nrow(steps))) {
    cat("\n--- Running ", steps$label[i], " ---\n", sep = "")
    source(steps$path[i], local = new.env(parent = globalenv()))
  }
}
