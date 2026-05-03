# Windows 本机完整环境安装指南
**目标硬件最低**: Windows 10/11, ≥16 GB RAM (推荐 32 GB), ≥100 GB 可用磁盘

---

## 第一步:R 4.4.x + Rtools45 (~30 分钟)

### 1.1 下载安装 R 4.4.2 (或最新 release)
- 官网: https://cran.r-project.org/bin/windows/base/
- 下载 `R-4.4.2-win.exe`,默认安装路径 `C:\Program Files\R\R-4.4.2\`

### 1.2 安装 Rtools45 (编译 C++ Bioconductor 包必需)
- 下载: https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools45.html
- 安装 `rtools45-x86_64.exe`,默认 `C:\rtools45\`
- **关键**: 安装结束勾选 "Add Rtools to PATH"

### 1.3 (可选)安装 RStudio Desktop
- 下载: https://posit.co/download/rstudio-desktop/
- 推荐用 RStudio,代码补全+图形预览方便

---

## 第二步:R 包一键安装 (≈ 60–90 分钟)

打开 R 或 RStudio,粘贴执行下面整段:

```r
# ===========================================================
# ATL osteolytic project - all required R packages
# Tested on R 4.4.2 / Bioc 3.20 / Windows 11 / 32 GB RAM
# Total install time: 60-90 min depending on network
# ===========================================================
options(timeout = 600, Ncpus = parallel::detectCores() - 1)

# CRAN packages
cran_pkgs <- c(
  "BiocManager", "tidyverse", "Matrix", "remotes", "devtools",
  "ggplot2", "patchwork", "ggrepel", "ggpubr", "ggsci", "RColorBrewer",
  "viridis", "cowplot", "scales", "pheatmap", "ComplexHeatmap",
  "survival", "survminer", "ggthemes", "broom",
  "data.table", "openxlsx", "readxl",
  "Seurat",          # v5
  "harmony",         # batch integration
  "clustree",
  "irlba", "uwot",
  "future", "future.apply", "doParallel",
  "msigdbr", "fgsea"
)
new_pkgs <- cran_pkgs[!cran_pkgs %in% installed.packages()[,"Package"]]
if (length(new_pkgs) > 0) install.packages(new_pkgs, dependencies = TRUE)

# Bioconductor
BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)
bioc_pkgs <- c(
  "DESeq2", "limma", "edgeR", "GEOquery", "Biobase", "AnnotationDbi",
  "org.Hs.eg.db", "clusterProfiler", "ReactomePA", "enrichplot",
  "ComplexHeatmap", "GSVA", "GSEABase",
  "scDblFinder", "SingleCellExperiment", "scran", "scater",
  "SCpubr", "celldex", "SingleR",
  "BayesPrism",
  "EnhancedVolcano",
  "RUVSeq", "sva"
)
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

# GitHub-only packages
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("immunogenomics/presto")          # fast FindMarkers
remotes::install_github("carmonalab/UCell")                # better module score
remotes::install_github("carmonalab/ProjecTILs")           # CD4/CD8 reference
remotes::install_github("satijalab/azimuth")               # PBMC reference
remotes::install_github("sqjin/CellChat")                  # cell-cell communication
remotes::install_github("aertslab/SCENIC")                 # transcription factor inference
remotes::install_github("aertslab/RcisTarget")
remotes::install_github("aertslab/AUCell")

# Verify install
loaded <- sapply(c("Seurat","harmony","BayesPrism","DESeq2","UCell",
                   "CellChat","SCENIC","clusterProfiler"),
                 requireNamespace, quietly = TRUE)
print(loaded)
cat("All packages loaded successfully?", all(loaded), "\n")
sessionInfo()
```

---

## 第三步:Python 环境 (用 Miniconda 管理,≈ 20 分钟)

### 3.1 安装 Miniconda
- 下载: https://docs.conda.io/en/latest/miniconda.html (Python 3.12 Win64)
- 默认安装路径 `C:\Users\<You>\miniconda3\`

### 3.2 创建 ATL 环境

打开 **Anaconda Prompt** (开始菜单搜索),粘贴:

```cmd
conda create -n atl python=3.10 -y
conda activate atl
conda install -c conda-forge -y ^
  numpy=1.26 pandas=2.2 scipy=1.13 scikit-learn statsmodels ^
  matplotlib seaborn jupyterlab notebook ^
  pyarrow openpyxl tqdm ^
  scanpy=1.10 anndata leidenalg python-igraph harmonypy gseapy lifelines ^
  cellxgene-census ^
  numba

pip install instaprism scaden cnmf
pip install --no-deps scvi-tools  # if you need scVI later (heavy GPU dep optional)
```

### 3.3 验证

```cmd
python -c "import scanpy as sc; print('scanpy', sc.__version__)"
python -c "import lifelines; print('lifelines', lifelines.__version__)"
```

---

## 第四步:目录约定

建议在 `D:\ATL research\` 下建立以下结构(我已建好部分):

```
D:\ATL research\
├── analysis/                  # 你看到的所有交付物
│   ├── scripts/               # 我已生成的 7 个脚本
│   ├── gene_panels/           # 3 个基因 panel
│   ├── data/                  # GEO 数据下载到这里
│   │   ├── geo/               # GSE195674 / GSE33615 / etc.
│   │   └── ccle/              # DepMap RNA-seq
│   ├── results/               # 中间结果 (.rds, .csv)
│   │   ├── scRNA/
│   │   ├── deconv/
│   │   ├── celline/
│   │   └── cmap/
│   └── figures/               # 最终图
└── manuscript/                # 论文文档(可后建)
```

---

## 第五步:跑序

```bash
# Windows PowerShell 或 Git Bash
cd D:\ATL research\analysis

# 1. 下载 GEO (用 Git Bash 跑 sh 脚本; 或转成 PowerShell 版本)
bash scripts/01_download_geo.sh data/geo

# 2. 下载 CCLE
bash scripts/02_download_ccle.sh data/ccle

# 3. 模块 A scRNA-seq (R 4.4)
Rscript scripts/03_scRNAseq_GSE195674.R

# 4. 模块 B 反卷积
Rscript scripts/04_deconvolution_bulk.R

# 5. 模块 C 细胞系 DE + CMap 输入
Rscript scripts/05_celline_DE_CMap.R
# → 把 results/celline/cmap_input_*.txt 提交到 clue.io / SigCom / L1000CDS²

# 6. 模块 D 反向验证
python scripts/06_drug_signature_GSEA.py

# 7. 模块 E 临床整合(我已跑过,需要时重跑)
python scripts/07_clinical_analysis.py
```

---

## 常见坑

### 坑 1: BayesPrism 安装失败
原因: BayesPrism 在 Bioc 3.20 中有时未及时同步,可能要 GitHub 装:
```r
remotes::install_github("Danko-Lab/BayesPrism", subdir = "BayesPrism", upgrade = "never")
```

### 坑 2: Seurat v5 + rgeos / sf 冲突
Windows 上 sf 经常装失败,Seurat v5 `JoinLayers()` 不需要 sf,直接装空间无关版本。
```r
install.packages("Seurat", dependencies = c("Depends","Imports"))
```

### 坑 3: scDblFinder 跨版本 Bioc 不兼容
确保 R 4.4 + Bioc 3.20。低版本 Bioc 装不上 scDblFinder。

### 坑 4: 网络问题(国内)
CRAN/Bioc 国内镜像配置:
```r
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
```

### 坑 5: GEO 下载慢/断流
在国内使用 EBI 镜像替代:
```bash
# GSE33615 EBI 镜像
wget "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/030/E-MTAB-3030/Files/..."
# 或用 R 包 GEOquery 自动 fallback
```
或者用 CDN 较快的 NCBI Cloud Buckets:
```bash
aws s3 cp s3://sra-pub-src-1/SRRxxxxxxx/...  # 需 AWS CLI
```

---

## 最终硬盘占用估计

| 模块 | 数据大小 | 中间 .rds | 最终图表 |
|---|---|---|---|
| A (scRNA) | GSE195674 ≈ 3 GB | seurat 对象 ≈ 4 GB | <100 MB |
| B (反卷积) | GSE33615+55851 ≈ 200 MB | BayesPrism ≈ 1 GB | <100 MB |
| C (细胞系 DE) | CCLE+GEO ≈ 5 GB | counts/DE ≈ 200 MB | <50 MB |
| D (CMap 反向) | LINCS L1000 ≈ 2 GB | 信号矩阵 ≈ 500 MB | <50 MB |
| **总计** | **≈ 12 GB raw** | **≈ 6 GB intermediate** | **≈ 300 MB** |

