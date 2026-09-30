# Competitive gene set tests (limma::cameraPR) for:
#   (a) Wnt/b-catenin target sets from Doumpas et al., EMBO J 2019 (Table EV1)
#   (b) curated craniofacial, cartilage, bone and control gene sets
source("R/utils.R")
res <- readRDS(file.path(cfg$tab, "deseq2_clean.rds"))
universe <- res$symbol

## (a) Paper-derived sets --------------------------------------------------
sheets <- c(bCat_indep = "bCat indep (65 out of 231)",
            bCat_dep   = "bCat dep (166 out of 231)",
            TCF_indep  = "TCF indep (27 out of 166)",
            TCF_dep    = "TCF dep (139 out of 166)")
paper <- lapply(sheets, function(s) {
  x <- as.data.frame(read_excel(cfg$wnt_xlsx, sheet = s, .name_repair = "minimal"))[, 1:5]
  names(x) <- c("EntrezID", "symbol", "gene", "logFC_WT", "pValue_WT")
  x <- x[!is.na(x$symbol), ]
  x$logFC_WT <- as.numeric(x$logFC_WT)
  x
})

sets <- list(); mapping <- data.frame()
for (k in names(paper)) {
  x <- paper[[k]]
  for (dirn in c("up", "down")) {
    h <- x$symbol[if (dirn == "up") x$logFC_WT > 0 else x$logFC_WT < 0]
    m <- human_to_mouse(h, universe)
    nm <- paste0("EMBOJ_", k, "_CHIR", dirn)
    sets[[nm]] <- m
    mapping <- rbind(mapping, data.frame(set = nm, n_human = length(h), n_mouse_tested = length(m)))
  }
}

## (b) Curated sets --------------------------------------------------------
cur <- read.delim(cfg$genesets)
for (s in unique(cur$set)) {
  g <- cur$gene[cur$set == s]
  sets[[s]] <- intersect(g, universe)
  mapping <- rbind(mapping, data.frame(set = s, n_human = length(g), n_mouse_tested = length(sets[[s]])))
}
write.csv(mapping, file.path(cfg$tab, "geneset_mapping.csv"), row.names = FALSE)

## Tests -------------------------------------------------------------------
idx <- lapply(sets, function(g) which(universe %in% g))
idx <- idx[lengths(idx) >= 3]
cam <- cameraPR(res$z, idx, inter.gene.cor = cfg$inter_gene_cor, sort = FALSE)
cam$set <- rownames(cam)
cam$median_log2FC <- sapply(idx, function(i) median(res$log2FoldChange[i]))
cam$frac_up       <- sapply(idx, function(i) mean(res$log2FoldChange[i] > 0))
cam$n_nominal_p05 <- sapply(idx, function(i) sum(res$pvalue[i] < 0.05))
cam$genes_nominal_p05 <- sapply(idx, function(i) {
  s <- res[i, ]; s <- s[s$pvalue < 0.05, ]
  paste(sprintf("%s(%.2f)", s$symbol, s$log2FoldChange), collapse = "; ")
})
cam <- cam[order(cam$PValue), c("set", "NGenes", "Direction", "PValue", "FDR",
                                "median_log2FC", "frac_up", "n_nominal_p05", "genes_nominal_p05")]
write.csv(cam, file.path(cfg$tab, "geneset_cameraPR.csv"), row.names = FALSE)
print(cam[, 1:7], digits = 3)

## Per-gene tables for each set -------------------------------------------
for (s in names(idx)) {
  t <- res[idx[[s]], c("symbol", "baseMean", "log2FoldChange", "pvalue", "padj_BH")]
  write.csv(t[order(t$pvalue), ], file.path(cfg$tab, paste0("set_", s, ".csv")), row.names = FALSE)
}

## Figures -----------------------------------------------------------------
save_png("04_geneset_summary.png", width = 1100, height = 800, {
  sc <- -log10(cam$PValue) * ifelse(cam$Direction == "Up", 1, -1)
  o <- order(sc)
  par(mar = c(5, 17, 3, 2))
  barplot(sc[o], names.arg = cam$set[o], horiz = TRUE, las = 1, cex.names = 0.7,
          col = ifelse(cam$FDR[o] < 0.05, ifelse(sc[o] > 0, "firebrick", "steelblue"), "grey75"),
          xlab = "signed -log10 p (cameraPR); right = up in mutant",
          main = "Gene set tests (coloured: FDR < 0.05)")
  abline(v = 0)
})

save_png("05_geneset_log2FC_boxplot.png", width = 1100, height = 800, {
  par(mar = c(5, 17, 3, 2))
  L <- lapply(idx[rev(cam$set)], function(i) res$log2FoldChange[i])
  boxplot(L, horizontal = TRUE, las = 1, cex.axis = 0.7, outline = FALSE,
          ylim = c(-3, 3), col = "grey90", xlab = "log2 fold change (mutant / control)")
  stripchart(L, add = TRUE, method = "jitter", pch = 16, cex = 0.5, col = "grey30")
  abline(v = 0, lty = 2)
})

key <- intersect(c("Osteoblast_osteogenesis", "Canonical_Wnt_readouts",
                   "EMBOJ_TCF_dep_CHIRup", "EMBOJ_TCF_dep_CHIRdown"), names(idx))
save_png("06_barcodeplots.png", width = 1100, height = 900, {
  par(mfrow = c(2, 2))
  for (s in key) barcodeplot(res$z, index = idx[[s]], main = s, xlab = "signed z (mutant vs control)")
})

## Concordance with the paper's CHIR response ------------------------------
save_png("07_concordance_with_EMBOJ.png", width = 1100, height = 900, {
  par(mfrow = c(2, 2))
  for (k in names(paper)) {
    x <- paper[[k]]; x$mouse <- NA
    key2 <- setNames(universe, toupper(universe)); x$mouse <- key2[toupper(x$symbol)]
    x <- x[!is.na(x$mouse), ]
    y <- res[x$mouse, "log2FoldChange"]
    r <- cor.test(x$logFC_WT, y, method = "spearman", exact = FALSE)
    plot(x$logFC_WT, y, pch = 16, col = "grey30",
         xlab = "log2FC, CHIR in WT HEK293T (paper)", ylab = "log2FC, mutant / control",
         main = sprintf("%s: rho = %.2f, p = %.2g (n = %d)", k, r$estimate, r$p.value, nrow(x)))
    abline(h = 0, v = 0, lty = 2)
  }
})
