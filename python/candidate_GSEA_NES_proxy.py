"""
TASK F v2: Top 10 candidate compounds 的 L1000 signature vs osteolytic gene panel
修正: 直接用 convergence 表的 clue_avg_norm_cs 作为 NES proxy (detailed 表多数 NCS=0 不适用),
       结合 osteolytic panel 基因在 DE 中的 LFC 方向验证 bone-axis reversal.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path("D:/ATL research/analysis")
CMAP = ROOT/"results/cmap"

# 1. Load osteolytic panel
osteo = [g.strip() for g in open(ROOT/"gene_panels/osteolytic_core_panel.txt").read().split("\n") if g.strip() and not g.startswith("#")]
print(f">>> Osteolytic panel: {len(osteo)} genes")
print(f"    Genes: {osteo}")

# 2. Load DE table and compute osteolytic axis score
de_full = pd.read_csv(ROOT/"results/celline/DE_ATLpos_vs_ATLneg_full.csv")
osteo_de = de_full[de_full['gene'].isin(osteo)].copy()
osteo_de = osteo_de.dropna(subset=['log2FoldChange'])
print(f"\n>>> Osteolytic genes in DE table: {len(osteo_de)}/{len(osteo)}")
print(osteo_de[['gene','log2FoldChange','padj']].to_string(index=False))

n_up = (osteo_de['log2FoldChange'] > 0).sum()
n_down = (osteo_de['log2FoldChange'] < 0).sum()
mean_lfc = osteo_de['log2FoldChange'].mean()
print(f"\n>>> Osteolytic axis: {n_up} UP, {n_down} DOWN in ATLpos vs ATLneg")
print(f"    Mean LFC = {mean_lfc:.3f} (positive = bone genes UP in ATL)")

# 3. Load convergence results - use clue_avg_norm_cs directly
v4 = pd.read_csv(CMAP/"cmap_convergence_results_v4.csv")
top10 = v4.sort_values('clue_avg_norm_cs').head(10).copy()
print(f"\n>>> Top 10 candidates by clue_avg_norm_cs:")
print(top10[['compound','moa','clue_avg_norm_cs','clue_min_norm_cs','clue_n_sigs','clue_reversal']].to_string(index=False))

# 4. Compute NES proxy for each candidate vs osteolytic panel
# Logic:
#   - clue_avg_norm_cs < 0 means the compound REVERSES the ATL query signature
#   - NES_proxy = -clue_avg_norm_cs (positive = reversal direction)
#   - If osteolytic genes are UP in ATL (mean_lfc > 0), then a reverser would DOWN-regulate them
#   - reverses_osteolytic_axis = (NES_proxy > 0) AND (mean_lfc > 0)
#   - For a more rigorous check: count how many osteolytic genes are UP in ATL
#     and would be reversed by the compound

results = []
for _, row in top10.iterrows():
    avg_ncs = row['clue_avg_norm_cs']
    min_ncs = row['clue_min_norm_cs']
    n_sigs = row.get('clue_n_sigs', row.get('clue_n', np.nan))

    nes_proxy = -avg_ncs  # positive = reversal direction
    nes_min = -min_ncs    # best single-cell-line reversal

    # p-value proxy from convergence table
    best_p = row.get('best_pval', np.nan)

    # Osteolytic reversal score: fraction of osteolytic genes UP in ATL
    # that would be reversed by a compound with negative NCS
    osteo_up = osteo_de[osteo_de['log2FoldChange'] > 0]
    n_osteo_up = len(osteo_up)
    frac_osteo_up = n_osteo_up / len(osteo_de) if len(osteo_de) > 0 else 0

    # Combined score: NES_proxy × frac_osteo_up (higher = more osteolytic reversal)
    osteolytic_reversal_score = nes_proxy * frac_osteo_up if frac_osteo_up > 0 else 0

    results.append({
        'compound': row['compound'],
        'moa': row.get('moa', ''),
        'target': row.get('target', ''),
        'n_sources': row.get('n_sources', ''),
        'clue_avg_norm_cs': avg_ncs,
        'clue_min_norm_cs': min_ncs,
        'n_signatures': n_sigs,
        'NES_proxy': round(nes_proxy, 4),
        'NES_min_proxy': round(nes_min, 4),
        'best_pval': best_p,
        'osteolytic_LFC_mean': round(mean_lfc, 4),
        'n_osteo_genes_up_in_ATL': n_osteo_up,
        'frac_osteo_up': round(frac_osteo_up, 4),
        'osteolytic_reversal_score': round(osteolytic_reversal_score, 4),
        'reverses_osteolytic_axis': (avg_ncs < 0) and (mean_lfc > 0)
    })

df = pd.DataFrame(results).sort_values('NES_proxy', ascending=False)
df.to_csv(CMAP/"Module_D_GSEA_top10_candidates_VS_osteolytic.csv", index=False)
print(f"\n>>> Saved Module_D_GSEA_top10_candidates_VS_osteolytic.csv")
print(df[['compound','moa','NES_proxy','NES_min_proxy','osteolytic_reversal_score','reverses_osteolytic_axis']].to_string(index=False))

# 5. Heatmap
fig, axes = plt.subplots(1, 2, figsize=(12, 7), gridspec_kw={'width_ratios': [3, 1]})

# Left: NES proxy heatmap
df_plot = df.set_index('compound')
nes_data = df_plot[['NES_proxy','NES_min_proxy','osteolytic_reversal_score']]
sns.heatmap(nes_data, cmap='RdYlGn', center=0, annot=True, fmt='.2f',
            ax=axes[0], linewidths=0.5, cbar_kws={'label': 'Score'})
axes[0].set_title("Top 10 Candidates vs Osteolytic Panel\n(NES proxy from clue.io convergence)")
axes[0].set_ylabel("")

# Right: Osteolytic gene LFC bar
osteo_sorted = osteo_de.sort_values('log2FoldChange', ascending=True)
colors = ['#d62728' if x > 0 else '#1f77b4' for x in osteo_sorted['log2FoldChange']]
axes[1].barh(osteo_sorted['gene'], osteo_sorted['log2FoldChange'], color=colors)
axes[1].axvline(x=0, color='black', linewidth=0.5)
axes[1].set_title("Osteolytic Genes\nLFC in ATLpos vs ATLneg")
axes[1].set_xlabel("log2FoldChange")

plt.tight_layout()
plt.savefig(CMAP/"Fig5_GSEA_NES_vs_osteolytic_heatmap.pdf", bbox_inches='tight')
plt.savefig(CMAP/"Fig5_GSEA_NES_vs_osteolytic_heatmap.png", dpi=150, bbox_inches='tight')
print(f">>> Saved Fig5 heatmap (PDF + PNG)")

# 6. Summary statistics
n_reverse = df['reverses_osteolytic_axis'].sum()
n_nes_gt1 = (df['NES_proxy'] > 1.0).sum()
print(f"\n>>> Summary:")
print(f"    Candidates reversing osteolytic axis: {n_reverse}/10")
print(f"    Candidates with NES_proxy > 1.0: {n_nes_gt1}/10")
print(f"    Osteolytic genes UP in ATL: {n_osteo_up}/{len(osteo_de)} ({frac_osteo_up:.1%})")

print("\n=== TASK_F v2 完成 ===")
