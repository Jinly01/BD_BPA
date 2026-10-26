"""
GSE70403 check. No control arrays here - everything is patient vs a common
reference, so the values are M-values and the test is one-sample against 0,
not a group comparison.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import ttest_1samp
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# --- CONFIG
EXPRESSION_FILE = "results/DEG_WGCNA/candidate_genes_union_matrix.csv"
CORE_GENE_FILE = "results/ml/final_core_genes.csv"
RESULT_DIR = "results/validation"
COHORT = "GSE70403"

SATURATION_PERCENTILE = 95   # colour saturation cutoff for the heatmap

os.makedirs(RESULT_DIR, exist_ok=True)




def significance_symbol(p_value):
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


# --- 1. Load the data
x_val = pd.read_csv(EXPRESSION_FILE, index_col=0)
core_genes = pd.read_csv(CORE_GENE_FILE)["Core_Gene"]
x_val_core = x_val[core_genes]

print(f"Loaded {x_val_core.shape[0]} samples and "
      f"{x_val_core.shape[1]} core genes.")
print(x_val_core.head())

plot_data_long = x_val_core.melt(var_name="Gene", value_name="M_Value")

# red if the mean M-value is positive, blue otherwise
gene_means = plot_data_long.groupby("Gene")["M_Value"].mean()
palette_dict = {gene: ("darkred" if val > 0 else "navy")
                for gene, val in gene_means.items()}


# --- 2. Boxplot with a one-sample t-test against zero
plt.figure(figsize=(6, 6))
plt.axhline(y=0, color="gray", linestyle="--", linewidth=1.5, alpha=0.7,
            zorder=0)
sns.boxplot(data=plot_data_long, x="Gene", y="M_Value", palette=palette_dict,
            width=0.6, fliersize=0)
sns.stripplot(data=plot_data_long, x="Gene", y="M_Value", color="black",
              alpha=0.3, size=3, jitter=True)

y_max = plot_data_long["M_Value"].max()
y_range = y_max - plot_data_long["M_Value"].min()

for i, gene in enumerate(plot_data_long["Gene"].unique()):
    _, p_val = ttest_1samp(x_val_core[gene], 0)
    plt.text(i, y_max + y_range * 0.05,
             f"{significance_symbol(p_val)}\n(p={p_val:.1e})",
             ha="center", va="bottom", fontsize=9, color="black")

plt.title(f"External Validation ({COHORT})", fontsize=14)
plt.ylabel("M-value (log2 Ratio)")
plt.xlabel("")
plt.ylim(top=y_max + y_range * 0.25)

custom_lines = [Line2D([0], [0], color="darkred", lw=4),
                Line2D([0], [0], color="navy", lw=4)]
plt.legend(custom_lines, ["Up (Mean > 0)", "Down (Mean < 0)"],
           loc="upper right")
plt.tight_layout()
plt.savefig(os.path.join(RESULT_DIR, f"{COHORT}_core_gene_boxplot.jpg"), dpi=300)
plt.show()


# ---- 3. PCA of the core genes -------
scaler = StandardScaler()
x_val_core_scaled = scaler.fit_transform(x_val_core)

pca = PCA(n_components=2)
pca_result = pca.fit_transform(x_val_core_scaled)

plt.figure(figsize=(10, 6))
plt.scatter(pca_result[:, 0], pca_result[:, 1], c="b", label="Samples",
            alpha=0.6)
plt.title(f"PCA of Core Genes ({COHORT})")
plt.xlabel(f"Principal Component 1 "
           f"({pca.explained_variance_ratio_[0] * 100:.1f}%)")
plt.ylabel(f"Principal Component 2 "
           f"({pca.explained_variance_ratio_[1] * 100:.1f}%)")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(RESULT_DIR, f"{COHORT}_pca.jpg"), dpi=300)
plt.show()


# ---- 4. Radar plot of the mean M-value profile -------
categories = list(core_genes)
values = x_val_core[categories].mean(axis=0).values
n_axes = len(categories)
angles = [n / float(n_axes) * 2 * np.pi for n in range(n_axes)]
values = np.concatenate((values, [values[0]]))
angles += angles[:1]

fig, ax = plt.subplots(figsize=(6, 6), dpi=100, subplot_kw=dict(polar=True))
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)
ax.plot(angles, values, linewidth=1, linestyle="solid", label="Core Genes")
ax.fill(angles, values, "b", alpha=0.1)
ax.set_yticklabels([])
ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories)
plt.title("Radar Chart of Key Gene Expression")
plt.tight_layout()
plt.savefig(os.path.join(RESULT_DIR, f"{COHORT}_radar.jpg"), dpi=300)
plt.show()


# 5. Clustered heatmap with a robust colour range
plot_data = x_val_core.T        # genes x samples

# saturate at p95 of |M|, otherwise a few outliers wash out the whole map
saturation_val = np.percentile(np.abs(plot_data.values),
                               SATURATION_PERCENTILE)
print(f"Colour saturation set to +/-{saturation_val:.2f}")

fig_width = min(20, max(8, plot_data.shape[1] * 0.3))
fig_height = min(15, max(6, plot_data.shape[0] * 0.5))

g = sns.clustermap(plot_data,
                   cmap="RdBu_r",
                   center=0,
                   vmin=-saturation_val,
                   vmax=saturation_val,
                   figsize=(fig_width, fig_height),
                   row_cluster=True,
                   col_cluster=True,
                   dendrogram_ratio=(0.15, 0.15),
                   cbar_pos=(0.02, 0.8, 0.03, 0.15),
                   linewidths=0.5,
                   linecolor="white",
                   yticklabels=True,
                   xticklabels=True)

plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10,
         fontweight="bold")
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=8)
g.ax_cbar.set_title("M-value", fontsize=10)
g.ax_cbar.tick_params(labelsize=9)
g.fig.suptitle(COHORT, fontweight="bold", fontsize=14, y=0.98)

plt.savefig(os.path.join(RESULT_DIR, f"{COHORT}_heatmap.jpg"), dpi=300, bbox_inches="tight")
plt.show()

print("\nStep 09 complete.")
