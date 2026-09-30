# Install everything the pipeline needs (run once).
cran <- c("readxl")
bioc <- c("edgeR", "limma")

for (p in cran) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
for (p in bioc) if (!requireNamespace(p, quietly = TRUE))
  BiocManager::install(p, update = FALSE, ask = FALSE)

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message("All packages available.")
