# Full edgeR quasi-likelihood pipeline. REQUIRES RAW COUNTS, which are not in
# the DESeq2 results file. Provide:
#   data/raw/counts.csv   genes x samples integer matrix, first column = gene_id
#   data/raw/samples.csv  columns: sample, group  (group = control / mutant)
#                         optional: sex, litter, batch (added to the design)
source("R/utils.R")
if (!file.exists(cfg$counts)) stop("Place raw counts at ", cfg$counts, " to run this script.")

counts  <- as.matrix(read.csv(cfg$counts, row.names = 1, check.names = FALSE))
samples <- read.csv(cfg$samples)
stopifnot(identical(colnames(counts), samples$sample))
samples$group <- relevel(factor(samples$group), ref = "control")

ann <- read.csv(cfg$deseq_csv)[, c("gene_id", "symbol")]
y <- DGEList(counts, samples = samples, group = samples$group,
             genes = ann[match(rownames(counts), ann$gene_id), ])

covars <- intersect(c("sex", "litter", "batch"), names(samples))
design <- model.matrix(as.formula(paste("~", paste(c(covars, "group"), collapse = " + "))), data = samples)

keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y)
cat("Genes kept:", nrow(y), "\n"); print(y$samples)

save_png("edgeR_01_MDS.png", {
  plotMDS(y, col = as.numeric(y$samples$group), labels = colnames(y))
})

y   <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
save_png("edgeR_02_BCV.png", plotBCV(y))
save_png("edgeR_03_QLdisp.png", plotQLDisp(fit))

coef <- "groupmutant"
qlf <- glmQLFTest(fit, coef = coef)
print(summary(decideTests(qlf)))
tt <- topTags(qlf, n = Inf)$table
write.csv(tt, file.path(cfg$tab, "edgeR_QLF_mutant_vs_control.csv"), row.names = FALSE)

# Fold-change threshold test (|log2FC| > log2(1.5)), as used in the paper
tr <- glmTreat(fit, coef = coef, lfc = log2(1.5))
write.csv(topTags(tr, n = Inf)$table, file.path(cfg$tab, "edgeR_TREAT_lfc1.5.csv"), row.names = FALSE)

save_png("edgeR_04_MD.png", { plotMD(qlf); abline(h = c(-1, 1), col = "blue") })

# Gene set tests on sample-level data: camera (competitive) and fry (self-contained)
cur <- read.delim(cfg$genesets)
gs  <- split(cur$gene, cur$set)
ids <- lapply(gs, function(g) which(y$genes$symbol %in% g))
ids <- ids[lengths(ids) >= 3]
write.csv(camera(y, ids, design, contrast = coef), file.path(cfg$tab, "edgeR_camera.csv"))
write.csv(fry(y, ids, design, contrast = coef),    file.path(cfg$tab, "edgeR_fry.csv"))

logcpm <- cpm(y, log = TRUE)
top <- head(order(tt$PValue), 40)
save_png("edgeR_05_heatmap_top40.png", width = 900, height = 1100, {
  m <- logcpm[rownames(tt)[top], ]; rownames(m) <- tt$symbol[top]
  heatmap(m, scale = "row", col = hcl.colors(50, "Blue-Red 3"))
})
