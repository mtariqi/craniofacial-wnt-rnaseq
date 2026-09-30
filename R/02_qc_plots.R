# Distribution checks and overview plots for the mutant vs control comparison.
source("R/utils.R")
res <- readRDS(file.path(cfg$tab, "deseq2_clean.rds"))
col <- c(Up = "firebrick", Down = "steelblue", NotSig = "grey70")[res$status]

save_png("01_pvalue_histogram.png",
  hist(res$pvalue, breaks = 50, col = "grey80", border = "white",
       main = "p-value distribution (mutant vs control)", xlab = "p-value"))

save_png("02_MA_plot.png", {
  plot(log10(res$baseMean + 1), res$log2FoldChange, pch = 16, cex = 0.4, col = col,
       xlab = "log10(baseMean + 1)", ylab = "log2 fold change (mutant / control)",
       main = "MA plot")
  abline(h = 0, col = "black")
})

save_png("03_volcano.png", {
  plot(res$log2FoldChange, -log10(res$pvalue), pch = 16, cex = 0.4, col = col,
       xlab = "log2 fold change (mutant / control)", ylab = "-log10 p-value",
       main = "Volcano plot (BH FDR < 0.05 coloured)")
  top <- head(res[res$status != "NotSig", ], 20)
  text(top$log2FoldChange, -log10(top$pvalue), top$symbol, cex = 0.6, pos = 3)
})
