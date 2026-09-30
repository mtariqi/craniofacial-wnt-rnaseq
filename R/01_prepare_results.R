# Load the DESeq2 results table, clean it, and derive a signed statistic.
source("R/utils.R")

res <- read.csv(cfg$deseq_csv, check.names = FALSE)
stopifnot(all(c("gene_id", "symbol", "baseMean", "log2FoldChange", "pvalue") %in% names(res)))
cat("Rows in file:", nrow(res), "\n")

res <- res[!is.na(res$symbol) & res$symbol != "" & !is.na(res$pvalue), ]
res <- res[order(res$pvalue), ]
res <- res[!duplicated(res$symbol), ]          # keep most significant row per symbol
cat("Genes with symbol and p-value (deduplicated):", nrow(res), "\n")

# The file has no padj column, so BH is recomputed over all tested genes.
# This will differ slightly from DESeq2's padj, which uses independent filtering.
res$padj_BH <- p.adjust(res$pvalue, method = "BH")

# Signed z-statistic: direction from log2FC, magnitude from the two-sided p-value.
p <- pmax(res$pvalue, .Machine$double.xmin)
res$z <- sign(res$log2FoldChange) * qnorm(p / 2, lower.tail = FALSE)

res$status <- ifelse(res$padj_BH < cfg$fdr,
                     ifelse(res$log2FoldChange > 0, "Up", "Down"), "NotSig")
print(table(res$status))

rownames(res) <- res$symbol
saveRDS(res, file.path(cfg$tab, "deseq2_clean.rds"))
write.csv(res[res$status != "NotSig", ], file.path(cfg$tab, "DE_genes_BH05.csv"), row.names = FALSE)
