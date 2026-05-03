# =============================================================
# FIX 02 v2: Module A - scRNA-seq (using verified 10X data)
# N=1 sample mode
# =============================================================
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(scDblFinder)
  library(SingleCellExperiment); library(UCell)
  library(dplyr); library(ggplot2); library(patchwork)
})
set.seed(42)

# ---- I/O ----
# Use the VERIFIED 10X data from the original RAW.tar extraction
DATA_DIR <- "D:/ATL research/analysis/data/geo/GSE195674/GSM5847946"
OUT_DIR  <- "D:/ATL research/analysis/results/scRNA"
PANEL_FILE <- "D:/ATL research/analysis/gene_panels/osteolytic_core_panel.txt"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(dir.exists(DATA_DIR), file.exists(PANEL_FILE))

GENE_PANEL <- readLines(PANEL_FILE)
GENE_PANEL <- GENE_PANEL[!grepl("^#", GENE_PANEL) & nchar(GENE_PANEL) > 0]
cat(sprintf(">>> Osteolytic panel: %d genes\n", length(GENE_PANEL)))

# ---- Step 1: Load data ----
cat(">>> Loading GSM5847946 (single GEX sample)\n")
m <- Read10X(DATA_DIR)
obj <- CreateSeuratObject(counts = m, project = "GSM5847946", min.cells = 3, min.features = 200)
obj$sample <- "GSM5847946"
N_SAMPLES <- 1
cat(sprintf(">>> Initial: %d cells x %d genes\n", ncol(obj), nrow(obj)))

# ---- Step 2: QC ----
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

# Doublet removal
sce <- as.SingleCellExperiment(obj)
sce <- scDblFinder(sce)
obj$doublet <- sce$scDblFinder.class
obj <- subset(obj, subset = doublet == "singlet")
cat(sprintf(">>> After QC + doublet removal: %d cells\n", ncol(obj)))

# ---- Step 3: Normalize + Reduction ----
obj <- SCTransform(obj, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = c(0.3, 0.5, 0.8))

# ---- Step 4: Lineage scores ----
DefaultAssay(obj) <- "SCT"
obj <- AddModuleScore(obj, features = list(c("CD3D","CD3E","CD4","CD8A","TRAC")), name = "T_score")
obj <- AddModuleScore(obj, features = list(c("CD68","CD14","C1QA","C1QB","CD163")), name = "Mac_score")
obj <- AddModuleScore(obj, features = list(c("KRT5","KRT14","KRT1","KRT10")), name = "Keratinocyte_score")
obj <- AddModuleScore(obj, features = list(c("FOXP3","IL2RA","IKZF2","CTLA4")), name = "Treg_score")

# ---- Step 5: Osteolytic module score ----
genes_in_data <- intersect(GENE_PANEL, rownames(obj))
cat(sprintf(">>> Osteolytic genes in data: %d / %d\n", length(genes_in_data), length(GENE_PANEL)))
obj <- AddModuleScore(obj, features = list(Osteo = genes_in_data), name = "OsteoScore")
obj <- AddModuleScore_UCell(obj, features = list(Osteo_UCell = genes_in_data))

# ---- Step 6: Plots ----
res_use <- "SCT_snn_res.0.5"
Idents(obj) <- res_use
p1 <- DimPlot(obj, group.by = res_use, label = TRUE, repel = TRUE) +
  ggtitle(sprintf("Clusters (%d cells, 1 sample, res 0.5)", ncol(obj)))
p2 <- FeaturePlot(obj, features = "OsteoScore1", min.cutoff = "q5", max.cutoff = "q95") +
  ggtitle("Osteolytic module score")
p3 <- VlnPlot(obj, features = "OsteoScore1", group.by = res_use, pt.size = 0) +
  NoLegend() + ggtitle("Osteolytic score per cluster")
p4 <- DotPlot(obj, features = genes_in_data, group.by = res_use) +
  RotatedAxis() + ggtitle("Bone genes by cluster")
p5 <- FeaturePlot(obj, features = c("T_score1","Mac_score1","Keratinocyte_score1","Treg_score1"),
                   min.cutoff = "q5", max.cutoff = "q95") & ggtitle("Lineage scores")

ggsave(file.path(OUT_DIR, "Fig1_UMAP_clusters.pdf"), p1, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "Fig2_osteolytic_score.pdf"), p2, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "Fig2c_violin_osteolytic.pdf"), p3, width = 8, height = 5)
ggsave(file.path(OUT_DIR, "Fig2b_dotplot_bone_genes.pdf"), p4, width = 14, height = 6)
ggsave(file.path(OUT_DIR, "Fig1b_lineage_scores.pdf"), p5, width = 12, height = 9)

# ---- Step 7: Cluster markers ----
markers_top <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_top, file.path(OUT_DIR, "cluster_markers_all.csv"), row.names = FALSE)

# ---- Step 8: Osteolytic score stats ----
osteo_stats <- data.frame()
for (cl in unique(Idents(obj))) {
  in_cl <- WhichCells(obj, idents = cl)
  out_cl <- setdiff(Cells(obj), in_cl)
  if (length(in_cl) >= 10 && length(out_cl) >= 10) {
    w <- wilcox.test(obj$OsteoScore1[in_cl], obj$OsteoScore1[out_cl], alternative = "greater")
    osteo_stats <- rbind(osteo_stats, data.frame(
      cluster = cl, n_cells = length(in_cl),
      median_score = median(obj$OsteoScore1[in_cl]),
      p_wilcox = w$p.value
    ))
  }
}
if (nrow(osteo_stats) > 0) {
  osteo_stats$padj <- p.adjust(osteo_stats$p_wilcox, method = "BH")
  osteo_stats <- osteo_stats[order(-osteo_stats$median_score), ]
  write.csv(osteo_stats, file.path(OUT_DIR, "OsteoScore_by_cluster.csv"), row.names = FALSE)
  cat("\n>>> Osteolytic score by cluster:\n")
  print(head(osteo_stats, 10))
}

# ---- Step 9: Auto cell type annotation ----
obj$cell_type <- "Other"
for (cl in unique(Idents(obj))) {
  m <- markers_top[markers_top$cluster == cl, ]
  m <- m[order(m$avg_log2FC, decreasing = TRUE), ]
  top_genes <- head(m$gene, 20)
  if (any(top_genes %in% c("CD3D","CD3E","CD4","TRAC","CD8A"))) {
    obj$cell_type[Idents(obj) == cl] <- if (any(top_genes %in% c("CD8A","CD8B","GZMK","GZMB"))) "CD8_T" else "CD4_T"
  } else if (any(top_genes %in% c("CD68","CD14","C1QA","CD163","LYZ"))) {
    obj$cell_type[Idents(obj) == cl] <- "Myeloid"
  } else if (any(top_genes %in% c("KRT5","KRT14","KRT1","KRT10"))) {
    obj$cell_type[Idents(obj) == cl] <- "Keratinocyte"
  } else if (any(top_genes %in% c("MS4A1","CD79A","CD19"))) {
    obj$cell_type[Idents(obj) == cl] <- "B"
  } else if (any(top_genes %in% c("FOXP3","IL2RA","IKZF2"))) {
    obj$cell_type[Idents(obj) == cl] <- "Treg_like"
  }
}

# Mark high-osteolytic CD4_T cluster as Tumor_T_osteolytic
if (nrow(osteo_stats) > 0) {
  top_osteo_cl <- osteo_stats$cluster[osteo_stats$padj < 0.01 & osteo_stats$median_score > median(osteo_stats$median_score)][1]
  if (!is.na(top_osteo_cl)) {
    cells_osteo <- WhichCells(obj, idents = top_osteo_cl)
    if (any(obj$cell_type[cells_osteo] %in% c("CD4_T","Treg_like"))) {
      obj$cell_type[cells_osteo] <- "Tumor_T_osteolytic"
    }
  }
}

cat("\n>>> Cell type distribution:\n")
print(table(obj$cell_type))

# ---- Step 10: Save ----
saveRDS(obj, file.path(OUT_DIR, "seurat_GSE195674_annotated.rds"))

cat("\n=========================================================\n")
cat("Module A COMPLETE\n")
cat(sprintf("  Seurat object: %s\n", file.path(OUT_DIR, "seurat_GSE195674_annotated.rds")))
cat(sprintf("  Cell types: %s\n", paste(names(table(obj$cell_type)), collapse = ", ")))
cat("\nWARNING: Only 1 GEX sample - cannot claim 'cell-of-origin'\n")
cat("  Module A must be presented as hypothesis-generating / corroboration\n")
cat("  Main narrative should use Bulk deconvolution (Module B)\n")
cat("=========================================================\n")
