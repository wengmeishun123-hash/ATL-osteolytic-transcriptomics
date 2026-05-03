# =============================================================
# FIX 04 v4: Module C - CCLE DE + CMap input (all fixes)
# =============================================================
suppressPackageStartupMessages({
  library(DESeq2); library(tibble); library(dplyr); library(tidyr)
  library(EnhancedVolcano); library(ggplot2)
})

OUT_DIR <- "D:/ATL research/analysis/results/celline"
CCLE_DIR <- "D:/ATL research/analysis/data/ccle"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1: Load ATL mapping ----
cat(">>> Loading ATL cell line mapping\n")
atl_map <- read.csv(file.path(CCLE_DIR, "atl_cell_line_mapping_final.csv"), stringsAsFactors = FALSE)

# ---- Step 2: Load counts file ----
cat("\n>>> Loading raw counts file\n")
expr <- read.csv(file.path(CCLE_DIR, "OmicsExpressionGenesExpectedCountProfile.csv"),
                 row.names = 1, check.names = FALSE)
cat(sprintf("    Expression: %d samples x %d genes\n", nrow(expr), ncol(expr)))

# ---- Step 3: Subset to ATL cell lines ----
atl_pids <- atl_map$ProfileID
atl_pids <- atl_pids[atl_pids %in% rownames(expr)]
cat(sprintf("    Found %d / %d ATL cell lines\n", length(atl_pids), nrow(atl_map)))

expr_subset <- expr[atl_pids, ]
# Replace ProfileID rownames with CellLine names
cellline_names <- atl_map$CellLine[match(rownames(expr_subset), atl_map$ProfileID)]
rownames(expr_subset) <- cellline_names

# Clean gene names
gene_ids <- colnames(expr_subset)
gene_symbols <- sub(" \\(ENSG[0-9]+\\)$", "", gene_ids)
keep_cols <- !duplicated(gene_symbols)
expr_subset <- expr_subset[, keep_cols]
colnames(expr_subset) <- gene_symbols[keep_cols]
cat(sprintf("    Subset: %d samples x %d unique genes\n", nrow(expr_subset), ncol(expr_subset)))

# ---- Step 4: Build coldata properly ----
# Samples are in ROWS, genes are in COLUMNS
# DESeq2 needs: counts matrix with genes as rows, samples as columns
cm <- t(round(as.matrix(expr_subset)))  # transpose: genes x samples
mode(cm) <- "integer"

coldata <- data.frame(
  group = factor(ifelse(colnames(cm) %in% atl_map$CellLine[atl_map$Group == "ATLpos"], "ATLpos", "ATLneg"),
                 levels = c("ATLneg", "ATLpos")),
  row.names = colnames(cm)
)
cat(sprintf("    coldata:\n"))
print(coldata)

# Filter low-count genes (at least 10 reads in at least 1 sample)
keep <- rowSums(cm >= 10) >= 1
cm <- cm[keep, ]
cat(sprintf("    After filtering: %d genes\n", nrow(cm)))

# ---- Step 5: DESeq2 ----
cat("\n>>> Running DESeq2\n")
dds <- DESeqDataSetFromMatrix(cm, colData = coldata, design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "ATLpos", "ATLneg"))

res_df <- as.data.frame(res) %>% rownames_to_column("gene")
write.csv(res_df, file.path(OUT_DIR, "DE_ATLpos_vs_ATLneg_full.csv"), row.names = FALSE)

n_sig <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1)
cat(sprintf("    DE results: %d significant (padj<0.05, |LFC|>1)\n", n_sig))

# ---- Step 6: CMap input ----
sig <- res_df %>% filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
if (nrow(sig) < 10) {
  # Relax criteria if too few significant genes
  cat("    Relaxing significance criteria...\n")
  sig <- res_df %>% filter(!is.na(pvalue) & pvalue < 0.05 & abs(log2FoldChange) > 0.5)
}
up150 <- sig %>% arrange(desc(log2FoldChange)) %>% head(150)
dn150 <- sig %>% arrange(log2FoldChange) %>% head(150)
writeLines(up150$gene, file.path(OUT_DIR, "cmap_input_up.txt"))
writeLines(dn150$gene, file.path(OUT_DIR, "cmap_input_down.txt"))
cat(sprintf("    CMap input: up=%d, down=%d\n", nrow(up150), nrow(dn150)))

# ---- Step 7: Volcano plot ----
p_volcano <- EnhancedVolcano(res_df,
  lab = res_df$gene, x = "log2FoldChange", y = "padj",
  pCutoff = 0.05, FCcutoff = 1,
  title = "ATL+ vs ATL- T cell lines (CCLE)",
  subtitle = sprintf("ATL+: %s | ATL-: %s",
    paste(atl_map$CellLine[atl_map$Group == "ATLpos"], collapse = ","),
    paste(atl_map$CellLine[atl_map$Group == "ATLneg"], collapse = ",")))
ggsave(file.path(OUT_DIR, "Fig4_volcano.pdf"), p_volcano, width = 9, height = 9)

cat("\n=========================================================\n")
cat("Module C COMPLETE\n")
cat(sprintf("  DE table: %s\n", file.path(OUT_DIR, "DE_ATLpos_vs_ATLneg_full.csv")))
cat(sprintf("  CMap up: %s (%d genes)\n", file.path(OUT_DIR, "cmap_input_up.txt"), nrow(up150)))
cat(sprintf("  CMap down: %s (%d genes)\n", file.path(OUT_DIR, "cmap_input_down.txt"), nrow(dn150)))
cat("=========================================================\n")
