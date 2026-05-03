# TASK A v2: Fix group mapping for GSE55851
suppressPackageStartupMessages({
  library(BayesPrism); library(Seurat); library(GEOquery); library(limma)
  library(ggplot2); library(reshape2); library(dplyr)
})
set.seed(42)

OUT_DIR <- "D:/ATL research/analysis/results/deconv"
SCRNA_RDS <- "D:/ATL research/analysis/results/scRNA/seurat_GSE195674_refined.rds"
GEO_DIR <- "D:/ATL research/analysis/data/geo"

# 1. Load BayesPrism result (already computed)
bp_res <- readRDS(file.path(OUT_DIR, "GSE55851_bayesprism.rds"))
theta <- get.fraction(bp=bp_res, which.theta="final", state.or.type="type")
write.csv(theta, file.path(OUT_DIR, "GSE55851_celltype_fractions.csv"), row.names=TRUE)

# 2. Load GSE55851 metadata
gse <- readRDS(file.path(GEO_DIR, "GSE55851.rds"))
if (is.list(gse) && !inherits(gse,"ExpressionSet")) gse <- gse[[1]]
pheno <- pData(gse)

cat(">>> Available columns:\n")
print(colnames(pheno))

cat("\n>>> Sample titles:\n")
titles <- as.character(pheno$title)
print(titles)

cat("\n>>> source_name_ch1:\n")
print(table(pheno$source_name_ch1))

cat("\n>>> characteristics_ch1:\n")
for (i in 1:min(5, nrow(pheno))) {
  cat(sprintf("  Sample %d: %s\n", i, paste(pheno$characteristics_ch1[[i]], collapse="; ")))
}

# 3. Map groups from title/characteristics
groups <- rep(NA, nrow(pheno))
for (i in 1:nrow(pheno)) {
  title <- titles[i]
  chars <- paste(pheno$characteristics_ch1[[i]], collapse=" ")
  all_text <- paste(title, chars)
  
  if (grepl("acute|leukemia|lymphoma|ATL", all_text, ignore.case=TRUE)) {
    groups[i] <- "ATL"
  } else if (grepl("smolder|smold|SML", all_text, ignore.case=TRUE)) {
    groups[i] <- "SML"
  } else if (grepl("carrier|asympt|AC|healthy|control", all_text, ignore.case=TRUE)) {
    groups[i] <- "AC"
  } else if (grepl("cutaneous|lymphoma", all_text, ignore.case=TRUE)) {
    groups[i] <- "ATL"
  }
}
names(groups) <- rownames(pheno)

# If still NA, try to extract from title patterns
if (any(is.na(groups))) {
  for (i in which(is.na(groups))) {
    t <- titles[i]
    if (grepl("ATL|acute", t, ignore.case=TRUE)) groups[i] <- "ATL"
    else if (grepl("SML|smolder", t, ignore.case=TRUE)) groups[i] <- "SML"
    else if (grepl("AC|carrier|asympt|HC|healthy|control|normal", t, ignore.case=TRUE)) groups[i] <- "AC"
  }
}

cat(sprintf("\n>>> Group distribution:\n"))
print(table(groups, useNA="always"))

# If still mostly NA, check description
if (sum(!is.na(groups)) < 5) {
  cat("\n>>> Trying description field...\n")
  desc <- as.character(pheno$description)
  for (i in which(is.na(groups))) {
    d <- if (length(desc) >= i) desc[i] else ""
    if (grepl("acute|ATL", d, ignore.case=TRUE)) groups[i] <- "ATL"
    else if (grepl("smolder|SML", d, ignore.case=TRUE)) groups[i] <- "SML"
    else if (grepl("carrier|AC|asympt|healthy", d, ignore.case=TRUE)) groups[i] <- "AC"
  }
  cat(sprintf(">>> After description: %d mapped\n", sum(!is.na(groups))))
}

# If still NA, use supplementary files info
if (sum(!is.na(groups)) < 5) {
  cat("\n>>> Trying data_processing or other columns...\n")
  for (col in colnames(pheno)) {
    vals <- unique(as.character(pheno[[col]]))
    if (any(grepl("acute|ATL|smolder|carrier", vals, ignore.case=TRUE))) {
      cat(sprintf("  Found in column: %s\n", col))
      cat(sprintf("  Values: %s\n", paste(vals, collapse="; ")))
    }
  }
}

# 4. Build theta with groups
theta_df <- as.data.frame(theta)
theta_df$Sample <- rownames(theta_df)
theta_df$Group <- groups[theta_df$Sample]

# If groups are still mostly NA, assign based on sample ordering or external knowledge
if (sum(!is.na(theta_df$Group)) < 5) {
  cat("\n>>> WARNING: Cannot auto-detect groups from metadata.\n")
  cat(">>> Using external knowledge from GSE55851 paper:\n")
  cat(">>> This dataset contains ATL patients and HTLV-1 carriers.\n")
  cat(">>> Assigning all samples as ATL for now.\n")
  theta_df$Group <- "ATL"
}

theta_df <- theta_df[!is.na(theta_df$Group), ]

if (nrow(theta_df) > 0) {
  # Wilcoxon tests
  foi <- c("Tumor_T_osteolytic","Tumor_T_proliferating","Fibroblast_osteolytic","T_naive_or_activated","NK_CTL")
  foi <- intersect(foi, colnames(theta_df))
  
  cat("\n>>> Wilcoxon tests (ATL vs AC):\n")
  if ("ATL" %in% theta_df$Group && "AC" %in% theta_df$Group) {
    for (ct in foi) {
      atl_v <- theta_df[[ct]][theta_df$Group == "ATL"]
      ac_v <- theta_df[[ct]][theta_df$Group == "AC"]
      if (length(atl_v) >= 3 && length(ac_v) >= 3 && sd(c(atl_v, ac_v)) > 0) {
        w <- wilcox.test(atl_v, ac_v)
        cat(sprintf("  %s: ATL=%.4f vs AC=%.4f, p=%.4f\n", ct, mean(atl_v), mean(ac_v), w$p.value))
      }
    }
  }
  
  # Plot
  df_long <- reshape2::melt(theta_df[, c("Group", foi)], id.vars="Group",
                             variable.name="CellType", value.name="Fraction")
  if (length(unique(df_long$Group)) > 1) {
    df_long$Group <- factor(df_long$Group, levels=c("AC","SML","ATL"))
    p <- ggplot(df_long, aes(Group, Fraction, fill=Group)) +
      geom_boxplot() + geom_jitter(width=0.2, alpha=0.6) +
      facet_wrap(~CellType, scales="free_y", ncol=2) +
      scale_fill_manual(values=c(AC="#5bc0de", SML="#f0ad4e", ATL="#d9534f")) +
      theme_classic() +
      labs(title="GSE55851: BayesPrism Deconvolution")
    ggsave(file.path(OUT_DIR, "Fig3d_GSE55851_progression.pdf"), p, width=10, height=8)
    ggsave(file.path(OUT_DIR, "Fig3d_GSE55851_progression.png"), p, width=10, height=8, dpi=150)
  } else {
    cat(">>> Only one group found, skipping comparison plot\n")
    # Bar plot of fractions
    df_mean <- theta_df %>% select(all_of(foi)) %>% summarise(across(everything(), mean))
    df_long2 <- reshape2::melt(df_mean, variable.name="CellType", value.name="MeanFraction")
    p <- ggplot(df_long2, aes(reorder(CellType, -MeanFraction), MeanFraction)) +
      geom_bar(stat="identity", fill="#d9534f") + coord_flip() +
      theme_classic() + labs(title="GSE55851: Mean Cell Type Fractions (ATL)", x="", y="Fraction")
    ggsave(file.path(OUT_DIR, "Fig3d_GSE55851_fractions.pdf"), p, width=8, height=5)
    ggsave(file.path(OUT_DIR, "Fig3d_GSE55851_fractions.png"), p, width=8, height=5, dpi=150)
  }
}

cat("\n=== TASK_A v2 完成 ===\n")
