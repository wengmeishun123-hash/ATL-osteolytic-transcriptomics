# 🔧 FIX 脚本执行顺序(给 Trae 用)

**执行前**: 把当前 results/scRNA、results/deconv、results/celline 备份重命名为 `*_trae_v1`,以便对照。

```powershell
# Windows PowerShell
cd D:\ATL research\analysis\results
ren scRNA scRNA_trae_v1
ren deconv deconv_trae_v1
ren celline celline_trae_v1
```

## Step 1: 重新下载完整 GSE195674 (~30 min)

```bash
Rscript scripts/FIX/FIX_01_GSE195674_diagnose_redownload.R
```

预期输出: `data/geo/GSE195674/GSM_inventory.csv` 列出所有 GSM,自动区分 GEX vs TCR,GEX 样本组织在 `data/geo/GSE195674/10X_samples/` 标准格式。

**🔴 关键检查点**: 跑完后查看 `GSM_inventory.csv`,**GSE195674 应有 ≥2 个 GEX 样本**(原 Joo 2022 论文报告多个 ATL 患者皮肤样本)。如果还是只 1 个 GEX,说明: (a) GSE195674 本身样本就少,或 (b) 网络下载失败,需手动重试。

## Step 2: 重跑 Module A (~1-2 hr)

```bash
Rscript scripts/FIX/FIX_02_module_A_proper.R
```

脚本会**自动检测样本数 N=1 vs N≥2**:
- N=1: 跑聚类 + module score,但日志会警告"不能 claim cell-of-origin"
- N≥2: 走完整 Harmony 整合 + 多样本互证

预期输出: `results/scRNA/seurat_GSE195674_annotated.rds` (含 `cell_type` 字段) + Fig1/2/2b/2c

## Step 3: 真跑 BayesPrism Module B (~30-60 min)

```bash
Rscript scripts/FIX/FIX_03_module_B_real_BayesPrism.R
```

预期输出:
- `results/deconv/GSE33615_bayesprism.rds` ← **真正的 BayesPrism 结果**
- `results/deconv/GSE33615_celltype_fractions.csv` (73 samples × cell types)
- `results/deconv/GSE33615_celltype_Z.rds` (cell-type-specific expression)
- `results/deconv/Fig3a_*.pdf, Fig3b_*.pdf`

## Step 4: 修复 Module C (~15-30 min)

**前置**: 需要在 DepMap 下载 raw counts 文件:
1. 访问 https://depmap.org/portal/data_page/?tab=allData
2. 找到 "Public 24Q4 → Expression Public 24Q4"
3. 下载 `OmicsExpressionGenesExpectedCountProfile.csv` (注意是 *Counts*, 不是 LogTPM)
4. 保存到 `D:\ATL research\analysis\data\ccle\`

```bash
Rscript scripts/FIX/FIX_04_module_C_proper_CCLE_DE_CMap.R
```

预期输出:
- `results/celline/DE_ATLpos_vs_ATLneg_full.csv`
- ⭐ `results/celline/cmap_input_up.txt` (top 150 upregulated genes)
- ⭐ `results/celline/cmap_input_down.txt` (top 150 downregulated genes)
- `results/celline/Fig4_volcano.pdf`

## Step 5: 手动提交 CMap (用户操作, ~30 min)

把 `cmap_input_up.txt` 和 `cmap_input_down.txt` 分别提交到:

| 工具 | URL | 模式 | 阈值 | 结果保存为 |
|---|---|---|---|---|
| clue.io | https://clue.io/query | Touchstone v2 | tau ≤ -90 | results/cmap/clue_query_results.csv |
| SigCom LINCS | https://maayanlab.cloud/sigcom-lincs/ | Reverser | z < -5 | results/cmap/sigcom_lincs_results.csv |
| L1000CDS² | https://maayanlab.cloud/L1000CDS2/ | Reverse | top 50 | results/cmap/l1000cds2_results.csv |

## Step 6: 跑 Module D 三工具收敛 (~1 min)

```bash
conda activate atl
python scripts/FIX/FIX_05_module_D_CMap_convergence.py
```

预期输出:
- `results/cmap/CMap_convergent_candidates.csv` (≥2/3 工具命中的化合物)
- `results/cmap/Fig5_UpSet_3tools.pdf`

---

## ⚠️ 重要修正 — 标题/Abstract 必须按实际样本数调整

如果跑完 FIX_01 后 GSE195674 仍只有 1 个 GEX 样本:
- ❌ 不能用 "single-cell and bulk deconvolution reveals cell-of-origin"
- ✅ 改用候选 3: "Single-cell and bulk transcriptomic dissection of osteolytic programs in ATL: bulk-deconvolution-based analysis with skin scRNA-seq corroboration"
- ✅ Section 5.1-5.2 (Module A 主图) 降级到 supplementary
- ✅ 主叙事改为 **Bulk 反卷积 (GSE33615 n=52 ATL)** 主图,scRNA 辅证

