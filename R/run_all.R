# Run from the repository root:  Rscript R/run_all.R
for (s in c("R/01_prepare_results.R", "R/02_qc_plots.R", "R/03_geneset_analysis.R")) {
  message("== ", s); source(s, echo = FALSE)
}
if (file.exists("data/raw/counts.csv")) source("R/edgeR_from_counts.R") else
  message("No data/raw/counts.csv: skipping the edgeR count-level pipeline.")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message("Done. See results/tables and results/figures.")
