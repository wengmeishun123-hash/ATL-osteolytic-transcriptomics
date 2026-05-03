# TASK G v2: MuSiC cross-validation vs BayesPrism
# Fix: Create pseudo-samples for MuSiC (requires multiple samples for variance estimation)
suppressPackageStartupMessages({
  library(MuSiC); library(Seurat); library(SingleCellExperiment)
  library(GEOquery); library(limma); library(ggplot2); library(dplyr)
})

OUT <- "D:/ATL research/analysis/results/deconv"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

cat(">>> Loading Seurat reference...\n")
seu <- readRDS("D:/ATL research/analysis/results/scRNA/seurat_GSE195674_refined.rds")
cat(sprintf("    Seurat object: %d cells, %d genes\n", ncol(seu), nrow(seu)))

# Subset to high-quality cell types
ct <- as.character(seu$cell_type_refined)
ct[is.na(ct)] <- "Other"
keep <- !ct %in% c("Mitochondrial_lowQ","Unknown","Other")
seu2 <- subset(seu, cells = Cells(seu)[keep])
ct2 <- ct[keep]
cat(sprintf("    After filtering: %d cells, cell types: %s\n", ncol(seu2),
            paste(unique(ct2), collapse=", ")))

# Create pseudo-samples (5 random splits) for MuSiC variance estimation
n_pseudo <- 5
cells_use <- colnames(seu2)
pseudo_id <- sample(paste0("pseudo_", 1:n_pseudo), length(cells_use), replace = TRUE)
names(pseudo_id) <- cells_use
cat(sprintf("    Created %d pseudo-samples for MuSiC\n", n_pseudo))
print(table(pseudo_id))

# Build SingleCellExperiment for MuSiC
cat(">>> Building SingleCellExperiment...\n")
counts_mat <- as.matrix(GetAssayData(seu2, assay = "RNA", layer = "counts"))
sce <- SingleCellExperiment(
  assays = list(counts = counts_mat),
  colData = data.frame(cellType = ct2,
                       sampleID = pseudo_id[colnames(counts_mat)],
                       row.names = colnames(seu2))
)
cat(sprintf("    SCE: %d genes x %d cells\n", nrow(sce), ncol(sce)))

# Load bulk GSE33615
cat(">>> Loading bulk GSE33615...\n")
gse <- readRDS("D:/ATL research/analysis/data/geo/GSE33615.rds")
if (is.list(gse) && !inherits(gse, "ExpressionSet")) gse <- gse[[1]]
cat(sprintf("    ExpressionSet: %d genes x %d samples\n", nrow(exprs(gse)), ncol(exprs(gse))))

expr <- exprs(gse)
fdat <- fData(gse)
sym_col <- intersect(c("Gene symbol", "Gene Symbol", "GENE_SYMBOL", "Symbol"), colnames(fdat))[1]
cat(sprintf("    Using symbol column: %s\n", sym_col))

pid2sym <- fdat[[sym_col]]
names(pid2sym) <- rownames(fdat)
expr2 <- expr
rownames(expr2) <- pid2sym[rownames(expr2)]
expr2 <- expr2[!is.na(rownames(expr2)) & rownames(expr2) != "" & !grepl("///", rownames(expr2)), ]
expr2 <- limma::avereps(expr2, ID = rownames(expr2))
bulk_linear <- 2^expr2
# Clean: remove genes with NA, Inf, or all-zero rows
bulk_linear <- bulk_linear[complete.cases(bulk_linear), ]
bulk_linear <- bulk_linear[is.finite(rowSums(bulk_linear)), ]
gene_sums <- rowSums(bulk_linear, na.rm = TRUE)
bulk_linear <- bulk_linear[gene_sums > 0, ]
cat(sprintf("    Bulk linear (cleaned): %d genes x %d samples\n", nrow(bulk_linear), ncol(bulk_linear)))

# Run MuSiC
common <- intersect(rownames(bulk_linear), rownames(sce))
cat(sprintf(">>> MuSiC common genes: %d\n", length(common)))
sce_use <- sce[common, ]
bulk_use <- bulk_linear[common, ]

# Final check: ensure no NA in bulk_use
if (any(is.na(bulk_use))) {
  cat("    Warning: NA values found in bulk_use, removing affected genes\n")
  good_genes <- rownames(bulk_use)[complete.cases(bulk_use)]
  bulk_use <- bulk_use[good_genes, ]
  sce_use <- sce_use[good_genes, ]
  cat(sprintf("    After NA removal: %d genes\n", nrow(bulk_use)))
}

cat(">>> Running MuSiC deconvolution (this may take a few minutes)...\n")
result <- music_prop(bulk.mtx = bulk_use, sc.sce = sce_use,
                     clusters = "cellType", samples = "sampleID",
                     verbose = TRUE)
fractions <- result$Est.prop.weighted
cat(sprintf("    MuSiC fractions: %d samples x %d cell types\n", nrow(fractions), ncol(fractions)))
cat("    MuSiC mean fractions:\n")
print(round(colMeans(fractions), 4))

write.csv(fractions, file.path(OUT, "MuSiC_GSE33615_fractions.csv"))
cat(sprintf(">>> Saved %s\n", file.path(OUT, "MuSiC_GSE33615_fractions.csv")))

# Compare to BayesPrism
cat(">>> Comparing MuSiC vs BayesPrism...\n")
bp_frac <- read.csv(file.path(OUT, "GSE33615_celltype_fractions_refined.csv"), row.names = 1)
common_samples <- intersect(rownames(fractions), rownames(bp_frac))
common_types <- intersect(colnames(fractions), colnames(bp_frac))
cat(sprintf("    Common samples: %d, Common cell types: %d\n", length(common_samples), length(common_types)))
cat(sprintf("    Common types: %s\n", paste(common_types, collapse = ", ")))

concord <- sapply(common_types, function(ct) {
  bp_v <- as.numeric(bp_frac[common_samples, ct])
  m_v  <- as.numeric(fractions[common_samples, ct])
  if (sum(m_v) == 0 || sum(bp_v) == 0) return(c(r = NA, rho = NA, p = NA))
  r <- cor(bp_v, m_v, method = "pearson")
  rho <- cor(bp_v, m_v, method = "spearman")
  p <- tryCatch(cor.test(bp_v, m_v, method = "pearson")$p.value, error = function(e) NA)
  c(r = r, rho = rho, p = p)
})
df <- data.frame(CellType = colnames(concord),
                 Pearson_r = concord["r", ],
                 Spearman_rho = concord["rho", ],
                 P_value = concord["p", ])
df <- df[order(-df$Pearson_r, na.last = TRUE), ]
write.csv(df, file.path(OUT, "cross_method_concordance_MuSiC_vs_BayesPrism.csv"), row.names = FALSE)
cat(">>> MuSiC vs BayesPrism concordance:\n")
print(df)

# Plot key cell types
key_ct <- intersect(c("Tumor_T_osteolytic", "Tumor_T_proliferating",
                       "T_naive_or_activated", "Myeloid"), common_types)
if (length(key_ct) == 0) {
  key_ct <- common_types[1:min(4, length(common_types))]
}
cat(sprintf("    Plotting cell types: %s\n", paste(key_ct, collapse = ", ")))

plots <- lapply(key_ct, function(ct) {
  d <- data.frame(BayesPrism = as.numeric(bp_frac[common_samples, ct]),
                  MuSiC = as.numeric(fractions[common_samples, ct]))
  r_val <- cor(d$BayesPrism, d$MuSiC)
  ggplot(d, aes(BayesPrism, MuSiC)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    theme_classic() +
    labs(title = sprintf("%s\nPearson r = %.3f", ct, r_val))
})
library(patchwork)
combined <- wrap_plots(plots, ncol = 2)
ggsave(file.path(OUT, "Fig3e_v2_MuSiC_vs_BayesPrism.pdf"), combined, width = 10, height = 8)
ggsave(file.path(OUT, "Fig3e_v2_MuSiC_vs_BayesPrism.png"), combined, width = 10, height = 8, dpi = 150)
cat(">>> Saved Fig3e_v2_MuSiC_vs_BayesPrism.pdf/png\n")

cat("\n=== TASK_G v2 完成 ===\n")
