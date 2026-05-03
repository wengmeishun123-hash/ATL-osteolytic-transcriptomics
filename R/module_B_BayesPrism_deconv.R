# =============================================================
# FIX 06 v2: 精细 cluster 注释 + 重跑 BayesPrism
# 修复: cell name 匹配问题
# =============================================================
suppressPackageStartupMessages({
  library(Seurat); library(BayesPrism); library(GEOquery)
  library(limma); library(ggplot2); library(reshape2); library(pheatmap); library(dplyr)
})
set.seed(42)

OUT_DIR   <- "D:/ATL research/analysis/results/deconv"
SCRNA_RDS <- "D:/ATL research/analysis/results/scRNA/seurat_GSE195674_annotated.rds"
GEO_DIR   <- "D:/ATL research/analysis/data/geo"
PANEL_FILE <- "D:/ATL research/analysis/gene_panels/osteolytic_core_panel.txt"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1: Load Seurat and check cluster IDs ----
cat(">>> Loading Seurat object\n")
seu <- readRDS(SCRNA_RDS)
res_use <- "SCT_snn_res.0.5"

# Check what cluster IDs exist
cat(sprintf(">>> Available identity classes: %s\n", paste(levels(Idents(seu)), collapse = ", ")))
cat(sprintf(">>> Available meta columns: %s\n", paste(colnames(seu@meta.data), collapse = ", ")))

# Set identity to res 0.5
Idents(seu) <- res_use
cluster_ids <- as.character(Idents(seu))
cat(sprintf(">>> Unique clusters: %s\n", paste(sort(unique(cluster_ids)), collapse = ", ")))

# ---- Step 2: Manual cluster -> cell_type mapping ----
# Based on cluster_markers_all.csv analysis from Claude's Day 4 report
manual_celltype <- c(
  "0"  = "T_naive_or_activated",
  "1"  = "B_or_GerminalCenter",
  "2"  = "Proliferating",
  "3"  = "Endothelial",
  "4"  = "T_resting",
  "5"  = "Myeloid",
  "6"  = "Fibroblast_osteolytic",
  "7"  = "Tumor_T_osteolytic",
  "8"  = "SmoothMuscle_or_pericyte",
  "9"  = "Mitochondrial_lowQ",
  "10" = "Other_epithelial",
  "11" = "Keratinocyte",
  "12" = "Tumor_T_proliferating",
  "13" = "NK_CTL",
  "14" = "Apocrine_glandular",
  "15" = "Lymphatic_endothelial"
)

# Apply mapping using the cluster IDs directly
cell_type_vec <- unname(manual_celltype[cluster_ids])
cell_type_vec[is.na(cell_type_vec)] <- "Unknown"
seu$cell_type_refined <- cell_type_vec

cat(">>> Refined cell_type distribution:\n")
print(table(seu$cell_type_refined))

# Filter out low-quality clusters
exclude_types <- c("Mitochondrial_lowQ", "Unknown")
keep <- !(seu$cell_type_refined %in% exclude_types)
seu_ref <- seu[, keep]
cat(sprintf(">>> Reference after filtering: %d cells\n", ncol(seu_ref)))

# ---- Step 3: Extract counts + cell type ----
sc_counts <- as.matrix(GetAssayData(seu_ref, assay = "RNA", layer = "counts"))
sc_celltype <- as.character(seu_ref$cell_type_refined)
sc_state <- sc_celltype

# Save refined seurat
saveRDS(seu, file.path("D:/ATL research/analysis/results/scRNA", "seurat_GSE195674_refined.rds"))

# ---- Step 4: Load + map GSE33615 bulk ----
cat("\n>>> Loading + mapping GSE33615\n")
gse <- readRDS(file.path(GEO_DIR, "GSE33615.rds"))
if (is.list(gse) && !inherits(gse, "ExpressionSet")) gse <- gse[[1]]
expr_matrix <- exprs(gse)

fdat <- fData(gse)
sym_col <- intersect(c("Gene symbol", "Gene Symbol"), colnames(fdat))[1]
if (is.na(sym_col)) sym_col <- colnames(fdat)[grepl("symbol", colnames(fdat), ignore.case = TRUE)][1]
cat(sprintf("    Using symbol column: %s\n", sym_col))

pid_to_sym <- fdat[[sym_col]]; names(pid_to_sym) <- rownames(fdat)
expr_t <- expr_matrix
rownames(expr_t) <- pid_to_sym[rownames(expr_t)]
expr_t <- expr_t[!is.na(rownames(expr_t)) & rownames(expr_t) != "" & !grepl("///", rownames(expr_t)), ]
expr_t <- limma::avereps(expr_t, ID = rownames(expr_t))

# Keep in log2 scale for BayesPrism
bk_dat <- t(expr_t)
bk_dat[!is.finite(bk_dat)] <- 0
cat(sprintf("    Bulk after mapping: %d samples x %d genes\n", nrow(bk_dat), ncol(bk_dat)))

# Intersect with sc reference
common <- intersect(colnames(bk_dat), rownames(sc_counts))
cat(sprintf("    Common genes: %d\n", length(common)))
sc_counts_use <- sc_counts[common, ]
bk_dat_use <- bk_dat[, common]

# ---- Step 5: Run BayesPrism with refined cell types ----
cat("\n>>> Running BayesPrism (refined types, ~10-30 min)\n")
prism <- new.prism(
  reference = t(sc_counts_use),
  mixture = bk_dat_use,
  input.type = "log2",
  cell.type.labels = sc_celltype,
  cell.state.labels = sc_state,
  outlier.cut = 0.01,
  outlier.fraction = 0.1,
  key = "Tumor_T_osteolytic"
)
bp_res <- run.prism(prism = prism, n.cores = max(1, parallel::detectCores() - 1))
saveRDS(bp_res, file.path(OUT_DIR, "GSE33615_bayesprism_refined.rds"))

# ---- Step 6: Extract refined fractions ----
theta <- get.fraction(bp = bp_res, which.theta = "final", state.or.type = "type")
write.csv(theta, file.path(OUT_DIR, "GSE33615_celltype_fractions_refined.csv"), row.names = TRUE)
cat(sprintf(">>> Refined fractions: %d x %d\n", nrow(theta), ncol(theta)))
cat("\n>>> Mean fraction by cell type:\n")
print(round(colMeans(theta), 3))

# ---- Step 7: Z matrix ----
Z <- get.exp(bp = bp_res, state.or.type = "type", cell.name = NULL)
saveRDS(Z, file.path(OUT_DIR, "GSE33615_celltype_Z_refined.rds"))

# ---- Step 8: Compare ATL vs Control ----
pheno <- pData(gse)
status_col <- grep("source|disease|tissue|condition|title", colnames(pheno), value = TRUE, ignore.case = TRUE)[1]
cat(sprintf("    Using phenotype column: %s\n", status_col))
group <- ifelse(grepl("ATL|patient|tumor", pheno[[status_col]], ignore.case = TRUE), "ATL", "Control")
names(group) <- rownames(pheno)
group <- group[rownames(theta)]

df_long <- as.data.frame(theta)
df_long$Sample <- rownames(df_long)
df_long$Group <- group[df_long$Sample]
df_long <- reshape2::melt(df_long, id.vars = c("Sample","Group"), variable.name = "CellType", value.name = "Fraction")

p_compare <- ggplot(df_long, aes(Group, Fraction, fill = Group)) +
  geom_boxplot() + facet_wrap(~CellType, scales = "free_y", ncol = 4) +
  theme_classic() +
  scale_fill_manual(values = c("ATL"="#d9534f","Control"="#5bc0de")) +
  labs(title = "GSE33615 cell type fractions: ATL vs Control (refined)")
ggsave(file.path(OUT_DIR, "Fig3c_ATL_vs_Control_fractions.pdf"), p_compare, width = 14, height = 10)

# Wilcoxon per cell type
wilcox_rows <- list()
for (ct in colnames(theta)) {
  a <- theta[group=="ATL", ct]; b <- theta[group=="Control", ct]
  if (length(a) >= 3 && length(b) >= 3) {
    w <- wilcox.test(a, b)
    wilcox_rows[[ct]] <- data.frame(CellType = ct,
                                    ATL_median = median(a), Control_median = median(b),
                                    P = w$p.value)
  }
}
wilcox_df <- do.call(rbind, wilcox_rows)
wilcox_df$padj <- p.adjust(wilcox_df$P, method = "BH")
write.csv(wilcox_df, file.path(OUT_DIR, "GSE33615_fractions_ATL_vs_Control_wilcox.csv"), row.names = FALSE)
cat("\n>>> Wilcoxon (ATL vs Control fractions):\n")
print(wilcox_df)

# ---- Step 9: Bone genes by refined cell type ----
GENE_PANEL <- readLines(PANEL_FILE)
GENE_PANEL <- GENE_PANEL[!grepl("^#", GENE_PANEL) & nchar(GENE_PANEL) > 0]
ct_means <- sapply(names(Z), function(ct) {
  z_ct <- Z[[ct]]
  bg <- intersect(GENE_PANEL, colnames(z_ct))
  if (length(bg) == 0) return(rep(NA, length(GENE_PANEL)))
  result <- rep(NA, length(GENE_PANEL)); names(result) <- GENE_PANEL
  result[bg] <- colMeans(z_ct[, bg, drop = FALSE])
  result
})
ct_means_log <- log2(ct_means + 1)
ct_means_log <- ct_means_log[rowSums(!is.na(ct_means_log)) > 0, ]
write.csv(ct_means_log, file.path(OUT_DIR, "GSE33615_bone_genes_by_celltype_refined.csv"), row.names = TRUE)
pheatmap(ct_means_log, scale = "row", cluster_cols = FALSE,
         filename = file.path(OUT_DIR, "Fig3b_bone_gene_celltype_heatmap_refined.pdf"),
         width = 8, height = 7,
         main = "Bone-axis cell-of-origin (refined cell types)")

cat("\n=========================================================\n")
cat("Module B refined COMPLETE\n")
cat(sprintf("   - Tumor_T_osteolytic in result: %s\n",
            "Tumor_T_osteolytic" %in% colnames(theta)))
cat(sprintf("   - ATL vs Control significant: %s\n",
            ifelse(any(wilcox_df$padj < 0.05, na.rm = TRUE), "YES", "NO")))
cat("=========================================================\n")
