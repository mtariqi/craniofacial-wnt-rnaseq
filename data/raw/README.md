# Raw inputs (not tracked by git)

| File | Source | Used by |
|---|---|---|
| `DESeq2_results_4control_4mutant.csv` | DESeq2 output, 4 control vs 4 mutant (mouse). Columns: gene_id, symbol, baseMean, log2FoldChange, pvalue | scripts 01–03 |
| `wnt_targets_EMBOJ2018.xlsx` | Table EV1, Doumpas et al., EMBO J 2019;38:e98873 | script 03 |
| `counts.csv` *(to add)* | Raw integer counts, genes × samples, first column gene_id | `edgeR_from_counts.R` |
| `samples.csv` *(to add)* | Columns `sample,group` (+ optional `sex,litter,batch`), same order as counts columns | `edgeR_from_counts.R` |

These files are excluded by `.gitignore` so unpublished data is not pushed to GitHub.
