# Craniofacial / skeletal RNA-seq: Wnt/β-catenin target analysis

Analysis of a mouse RNA-seq experiment (4 control vs 4 mutant) for craniofacial,
skeletal and osteogenic programmes, and for Wnt/β-catenin targets defined in
Doumpas et al., *EMBO J* 2019 (TCF/LEF-dependent and -independent targets).

The full write-up is in [`report/report.md`](report/report.md).

## Repository layout

```
R/
  00_install_packages.R   install edgeR, limma, readxl
  utils.R                 paths, settings, helpers (edit here)
  01_prepare_results.R    clean DESeq2 table, BH FDR, signed z-statistic
  02_qc_plots.R           p-value histogram, MA plot, volcano
  03_geneset_analysis.R   cameraPR tests: paper Wnt sets + skeletal sets, concordance plots
  edgeR_from_counts.R     full edgeR QL pipeline (needs raw counts)
  run_all.R               runs everything in order
data/
  raw/                    input files (git-ignored; see data/raw/README.md)
  genesets/               curated craniofacial / cartilage / bone / control gene sets
results/                  generated tables and figures (git-ignored, regenerate)
report/                   report.md and its figures
setup_github.fish         git init + push to a private GitHub repo
```

## Running

From the repository root:

```sh
Rscript R/00_install_packages.R   # once
Rscript R/run_all.R
```

`run_all.R` runs the results-level analysis on the DESeq2 table. If
`data/raw/counts.csv` and `data/raw/samples.csv` exist, it also runs the edgeR
quasi-likelihood pipeline in `R/edgeR_from_counts.R`.

## Important limitation

The supplied file contains DESeq2 **results** (log2FC and p-values), not raw
counts. edgeR models counts directly, so the edgeR differential expression step
cannot be run until the count matrix is added. The current analysis therefore
uses limma's `cameraPR`, which tests gene sets from a per-gene statistic.

## Reference

Doumpas N, Lampart F, Robinson MD, et al. TCF/LEF dependent and independent
transcriptional regulation of Wnt/β-catenin target genes. *EMBO J.*
2019;38:e98873. doi:10.15252/embj.201798873
