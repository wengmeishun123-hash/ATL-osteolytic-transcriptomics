# ATL-osteolytic-transcriptomics

**Single-cell and bulk transcriptomic dissection of osteolytic programs in adult T-cell leukemia/lymphoma reveals tumour-T-cell-intrinsic and Treg-like contributions to bone destruction, supported by serum bone turnover phenotypes in a Chinese cohort**

Weng M, Zhan X, Chen R, Wang Y, Gao X. *Frontiers in Immunology* 2026 (submitted).

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Zenodo](https://img.shields.io/badge/Zenodo-DOI%20pending-blue)](https://zenodo.org/)

## Overview

This repository contains all code, gene panels, and analysis scripts used in our integrative single-cell, bulk-deconvolution, drug-repurposing, and clinical analysis of ATL osteolytic disease. The published manuscript and Supplementary Materials are linked from the corresponding-author institutional page.

## Repository structure

```
ATL-osteolytic-transcriptomics/
├── R/                              # R analysis scripts
│   ├── module_A_scRNA_GSE195674.R          # Cluster annotation + osteolytic module
│   ├── module_B_BayesPrism_deconv.R        # GSE33615 / GSE55851 deconvolution
│   ├── module_B_MuSiC_crossval.R           # Cross-method validation
│   ├── module_C_DESeq2_celline.R           # ATL+/− DE + CMap input
│   └── clinical_KM_Cox.R                    # 29-patient clinical analysis
├── python/                         # Python analysis scripts
│   ├── clinical_analysis.py                # Clinical descriptive + Wilcoxon
│   ├── cmap_convergence_v4.py              # 3-tool reverser convergence
│   ├── candidate_GSEA_NES_proxy.py         # Section 3.5 NES proxy approximation
│   └── profileid_validation.py             # DepMap OmicsProfiles validation
├── data_panels/                    # Curated gene panels
│   ├── osteolytic_core_panel.txt   # 21 genes
│   ├── osteolytic_extended_panel.txt
│   └── atl_signature_genes.txt
├── results_examples/               # Anonymized example outputs
│   ├── GSE33615_celltype_fractions_refined.csv
│   ├── cross_method_concordance_MuSiC_vs_BayesPrism.csv
│   └── cmap_convergence_results_v4.csv
├── docs/                           # Methods documentation
│   ├── 00_INSTALL_Windows.md       # Environment setup
│   ├── RUN_ORDER.md                # Pipeline execution order
│   └── manuscript_references_FINAL.csv  # Bibliography
├── README.md                       # This file
├── LICENSE                         # MIT
└── CITATION.cff                    # Citation file format
```

## Quick start

### Requirements
- R 4.4.2+ (Bioconductor 3.20+); recommended OS: Ubuntu 22.04 / Windows 11
- Python 3.10+ via Miniconda (atl env)
- Cell-line analysis: ≥ 16 GB RAM (32 GB recommended for GSE195674 single-cell processing)
- Disk: ≥ 50 GB for full GEO + CCLE data

### Installation

Detailed step-by-step instructions are in `docs/00_INSTALL_Windows.md`. Briefly:

```bash
# R packages (CRAN + Bioconductor + GitHub)
Rscript -e "install.packages(c('Seurat','harmony','remotes','BiocManager'))"
Rscript -e "BiocManager::install(c('DESeq2','clusterProfiler','BayesPrism','scDblFinder'))"
Rscript -e "remotes::install_github(c('xuranw/MuSiC','immunogenomics/presto'))"

# Python environment
conda create -n atl python=3.10 -y
conda activate atl
pip install scanpy anndata scipy statsmodels pandas matplotlib seaborn
```

### Pipeline execution order

```bash
# Module A: scRNA-seq cell-of-origin
Rscript R/module_A_scRNA_GSE195674.R

# Module B: BayesPrism + MuSiC deconvolution
Rscript R/module_B_BayesPrism_deconv.R
Rscript R/module_B_MuSiC_crossval.R

# Module C: cell-line DE + CMap input
Rscript R/module_C_DESeq2_celline.R

# Module D: CMap convergence analysis
python python/cmap_convergence_v4.py
python python/candidate_GSEA_NES_proxy.py

# Clinical analysis
python python/clinical_analysis.py
Rscript R/clinical_KM_Cox.R
```

See `docs/RUN_ORDER.md` for detailed prerequisites and expected outputs per module.

## Data

Public datasets analyzed in this study:
- GSE195674 — ATL skin scRNA-seq (1 patient + matched TCR)
- GSE33615 — ATL PBMC microarray (52 ATL + 21 healthy controls)
- GSE55851 — HTLV-1 progression cohort (9 carriers + 12 ATL)
- CCLE / DepMap Public 24Q4 — RNA-seq for 7 T-cell lines

29-patient Chinese ATL clinical cohort (de-identified) is available from the corresponding author upon reasonable request, subject to IRB approval (Ningde Municipal Hospital of Ningde Normal University, approval No. NSYKYLL-2026-145) and a formal data-sharing agreement.

## Citation

If you use code or analyses from this repository, please cite:

> Weng M, Zhan X, Chen R, Wang Y, Gao X. Single-cell and bulk transcriptomic dissection of osteolytic programs in adult T-cell leukemia/lymphoma reveals tumour-T-cell-intrinsic and Treg-like contributions to bone destruction, supported by serum bone turnover phenotypes in a Chinese cohort. *Frontiers in Immunology* (2026, submitted).

A `CITATION.cff` file is provided.

## Contact

**Corresponding author:** Xiaojuan Gao, MD (gyinxin2005@163.com), Ningde Municipal Hospital of Ningde Normal University.

**First-author code maintainer:** Meishun Weng, MS (wengmeishun2026@163.com).

## Acknowledgements

We thank the Joo et al. (Frontiers in Immunology 2022) team for sharing GSE195674, the Choi et al. (Blood 2009) team for GSE33615, the Yamagishi et al. (2012) team for GSE55851, and the CCLE / DepMap Public 24Q4 consortium for the reference T-cell-line transcriptomes.

## License

MIT — see `LICENSE` file.
