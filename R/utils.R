# Shared configuration and helper functions.
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(readxl)
})

cfg <- list(
  deseq_csv = "data/raw/DESeq2_results_4control_4mutant.csv",
  wnt_xlsx  = "data/raw/wnt_targets_EMBOJ2018.xlsx",
  genesets  = "data/genesets/curated_skeletal_craniofacial_genesets.tsv",
  counts    = "data/raw/counts.csv",        # only needed for edgeR_from_counts.R
  samples   = "data/raw/samples.csv",       # only needed for edgeR_from_counts.R
  tab       = "results/tables",
  fig       = "results/figures",
  fdr       = 0.05,
  inter_gene_cor = 0.01                     # cameraPR default
)
for (d in c(cfg$tab, cfg$fig)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

save_png <- function(name, expr, width = 900, height = 750, res = 120) {
  png(file.path(cfg$fig, name), width = width, height = height, res = res)
  on.exit(dev.off())
  force(expr)
}

# Human (HGNC) -> mouse (MGI) symbol matching by case-insensitive identity.
# This recovers most 1:1 orthologs (AXIN2 -> Axin2) but misses lncRNAs,
# primate-specific genes and renamed orthologs. For publication, replace with
# a proper ortholog table (e.g. MGI HOM_MouseHumanSequence.rpt or biomaRt).
human_to_mouse <- function(human_symbols, mouse_universe) {
  key <- setNames(mouse_universe, toupper(mouse_universe))
  out <- unname(key[toupper(human_symbols)])
  out[!is.na(out)]
}
