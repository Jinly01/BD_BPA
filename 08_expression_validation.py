"""
Expression level checks on the core genes: heatmaps, boxplots with t-tests,
split violin, PCA with 2-sigma ellipses, radar.
"""

import os
from math import pi

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Ellipse
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler, StandardScaler

# ---- CONFIG ----------
TRAIN_SCALED_FILE = "data/train_scaled.csv"
TEST_SCALED_FILE = "data/test_scaled.csv"
CORE_GENE_FILE = "results/ml/final_core_genes.csv"
RESULT_DIR = "results/validation"

ID_COLUMN = "SampleID"
LABEL_COLUMN = "label"

# TODO: pull the cohort names from the metadata instead of hard-coding
TRAIN_COHORT = "GSE17114"
TEST_COHORT = "GSE209567"

GROUP_PALETTE = {"Control": "navy", "Case": "darkred"}
PCA_PALETTE = {"Control": "#1f77b4", "Case": "#d62728"}

os.makedirs(RESULT_DIR, exist_ok=True)


def fig(name):
    return os.path.join(RESULT_DIR, name)


def significance_symbol(p_value):
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


# 1. Load the data
core_genes = pd.read_csv(CORE_GENE_FILE)["Core_Gene"]
core_gene_list = core_genes.tolist()

data_train = pd.read_csv(TRAIN_SCALED_FILE)
sample_ids_train = data_train[ID_COLUMN]
y_train = data_train[LABEL_COLUMN]
x_train_core = data_train.drop(columns=[ID_COLUMN, LABEL_COLUMN])[core_genes]

data_test = pd.read_csv(TEST_SCALED_FILE)
sample_ids_test = data_test[ID_COLUMN]
y_test = data_test[LABEL_COLUMN]
x_test_core = data_test.drop(columns=[ID_COLUMN, LABEL_COLUMN])[core_genes]

print(f"Core genes: {core_gene_list}")


# 2. Clustered heatmaps
def plot_classified_heatmap(expression_df, labels, sample_ids, title,
                            save_name):
    """Clustered heatmap, genes x samples, with a group colour bar on top.
    expression_df is samples x genes."""
    data_for_plot = expression_df.copy()
    data_for_plot.index = sample_ids
    plot_data = data_for_plot.T            # genes x samples

    unique_labels = sorted(labels.unique())
    label_colors_map = {unique_labels[0]: "#80CBC4",
                        unique_labels[1]: "#F06292"}
    sample_colors = pd.Series(labels.values,
                              index=sample_ids).map(label_colors_map)

    # scale the canvas or the sample labels collide
    width = max(10, plot_data.shape[1] * 0.35)
    height = max(8, plot_data.shape[0] * 0.5)

    g = sns.clustermap(plot_data,
                       col_colors=sample_colors,
                       cmap="RdBu_r",
                       center=0,
                       vmin=-2, vmax=2,
                       figsize=(width, height),
                       row_cluster=True,
                       col_cluster=True,
                       dendrogram_ratio=(0.1, 0.15),
                       cbar_pos=(0.02, 0.85, 0.02, 0.12),
                       linewidths=0.5, linecolor="white",
                       xticklabels=True, yticklabels=True)

    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=9)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10,
             fontweight="bold")

    legend_patches = [mpatches.Patch(color=c, label=str(l))
                      for l, c in label_colors_map.items()]
    g.ax_col_dendrogram.legend(handles=legend_patches, title="Group",
                               loc="lower left", bbox_to_anchor=(0, 1.05),
                               ncol=1, frameon=False, fontsize=10)

    g.fig.suptitle(title, fontsize=16, fontweight="bold", y=1.05)
    g.ax_cbar.set_title("Z-Score", fontsize=9)

    plt.savefig(save_name, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"Saved heatmap: {save_name}")


plot_classified_heatmap(x_train_core, y_train, sample_ids_train,
                        TRAIN_COHORT,
                        fig(f"internal_heatmap_{TRAIN_COHORT}.jpg"))
plot_classified_heatmap(x_test_core, y_test, sample_ids_test,
                        TEST_COHORT,
                        fig(f"external_heatmap_{TEST_COHORT}.jpg"))


# --- 3. Grouped boxplots with a two-sample t-test
def plot_core_gene_boxplot(expression_df, labels, title, save_name):
    plot_data = expression_df.copy()
    plot_data["Group"] = np.where(labels == 1, "Case", "Control")
    plot_long = plot_data.melt(id_vars="Group", var_name="Gene",
                               value_name="Z-Score Expression")

    plt.figure(figsize=(10, 6))
    sns.boxplot(data=plot_long, x="Gene", y="Z-Score Expression", hue="Group",
                palette=GROUP_PALETTE, width=0.6, fliersize=0)
    sns.stripplot(data=plot_long, x="Gene", y="Z-Score Expression",
                  hue="Group", dodge=True, color="black", alpha=0.3, size=3,
                  jitter=True, legend=False)

    y_max = plot_long["Z-Score Expression"].max()
    y_range = y_max - plot_long["Z-Score Expression"].min()

    for i, gene in enumerate(core_gene_list):
        v_case = plot_data[plot_data["Group"] == "Case"][gene]
        v_ctrl = plot_data[plot_data["Group"] == "Control"][gene]
        _, p_val = ttest_ind(v_case, v_ctrl)
        plt.text(i, y_max + y_range * 0.05,
                 f"{significance_symbol(p_val)}\n(p={p_val:.1e})",
                 ha="center", va="bottom", fontsize=9, color="black")

    plt.title(title, fontsize=14)
    plt.ylabel("Standardized Expression (Z-score)")
    plt.ylim(top=y_max + y_range * 0.2)      # leave room for the annotations
    plt.xlabel("")
    plt.legend(loc="upper right", title="Group")
    plt.tight_layout()
    plt.savefig(save_name, dpi=300)
    plt.show()


plot_core_gene_boxplot(x_train_core, y_train,
                       "Expression of Core Genes in Internal Validation Cohort",
                       fig("internal_gene_boxplot.jpg"))
plot_core_gene_boxplot(x_test_core, y_test,
                       "Expression of Core Genes in External Validation Cohort",
                       fig("external_gene_boxplot.jpg"))


# ---- 4. Split violin plot on the external cohort ----------
sns.set_context("paper", font_scale=1.5)
sns.set_style("ticks")

X_plot = x_test_core.copy()
y_str = y_test.map({1: "Case", 0: "Control"})
X_plot["Group"] = y_str

df_melt = X_plot.melt(id_vars="Group", var_name="Gene",
                      value_name="Expression")

plt.figure(figsize=(10, 6))
ax = sns.violinplot(x="Gene", y="Expression", hue="Group", data=df_melt,
                    split=True,            # one group per side
                    density_norm="width",
                    inner="quartile",
                    palette=GROUP_PALETTE,
                    saturation=1, linewidth=1.2, alpha=0.7)

# Restyle the quartile lines drawn by seaborn.
for line in ax.lines:
    line.set_linestyle("--")
    line.set_linewidth(1.0)
    line.set_color("gray")
    line.set_alpha(0.8)

sns.stripplot(x="Gene", y="Expression", hue="Group", data=df_melt, dodge=True,
              jitter=0.2, size=4, color="gray", linewidth=0, alpha=0.3, ax=ax)

plt.title("External Validation", fontsize=14, pad=20)
plt.ylabel("Standardized Expression (Z-score)", fontsize=14)
plt.xlabel("")

# drop the duplicated stripplot legend entries
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], title="Group", loc="upper right",
          frameon=False, fontsize=12)
sns.despine(trim=True)
plt.tight_layout()
plt.savefig(fig("external_violin.jpg"), dpi=300, bbox_inches="tight")
plt.show()


# 5. PCA with 2-sigma confidence ellipses
def draw_confidence_ellipse(x, y, ax, n_std=2.0, facecolor="none", **kwargs):
    """n_std confidence ellipse for a 2D cloud."""
    if x.size < 2 or y.size < 2:
        return
    cov = np.cov(x, y)
    eigenvalues, eigenvectors = np.linalg.eig(cov)
    radius_x, radius_y = np.sqrt(eigenvalues)
    ellipse = Ellipse((0, 0), width=radius_x * 2 * n_std,
                      height=radius_y * 2 * n_std, facecolor=facecolor,
                      **kwargs)
    ellipse.set_transform(ax.transData)
    ellipse.set_center((np.mean(x), np.mean(y)))
    ellipse.set_angle(np.degrees(np.arctan2(eigenvectors[1, 0],
                                            eigenvectors[0, 0])))
    ax.add_patch(ellipse)


scaler = StandardScaler()
X_pca_input = scaler.fit_transform(x_test_core)
pca = PCA(n_components=2)
components = pca.fit_transform(X_pca_input)

df_pca = pd.DataFrame(components, columns=["PC1", "PC2"])
df_pca["Group"] = y_str.values

plt.figure(figsize=(8, 7))
sns.scatterplot(x="PC1", y="PC2", hue="Group", style="Group", data=df_pca,
                palette=PCA_PALETTE, s=120, alpha=0.9)

ax = plt.gca()
for group, color in PCA_PALETTE.items():
    subset = df_pca[df_pca["Group"] == group]
    draw_confidence_ellipse(subset["PC1"], subset["PC2"], ax, n_std=2.0,
                            edgecolor=color, linestyle="--", alpha=0.6, lw=2)

plt.title(f"PCA: External Validation Cohort\n"
          f"(Based on {len(core_gene_list)} Core Genes)",
          fontsize=14, fontweight="bold")
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
plt.grid(True, linestyle="--", alpha=0.3)
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(fig("external_pca.jpg"), dpi=300)
plt.show()


# ---- 6. Radar plot of the mean expression profile ----------
# min-max so every spoke sits in 0-1
scaler_mm = MinMaxScaler()
X_radar = pd.DataFrame(scaler_mm.fit_transform(x_test_core),
                       columns=core_gene_list)
X_radar["Group"] = y_str.values
means = X_radar.groupby("Group").mean()

categories = list(means.columns)
n_axes = len(categories)
angles = [n / float(n_axes) * 2 * pi for n in range(n_axes)]
angles += angles[:1]        # close the polygon

plt.figure(figsize=(8, 8))
ax = plt.subplot(111, polar=True)
ax.set_theta_offset(pi / 2)
ax.set_theta_direction(-1)
plt.xticks(angles[:-1], categories, size=14, fontweight="bold")

for group, color in [("Case", "#d62728"), ("Control", "#1f77b4")]:
    if group not in means.index:
        continue
    values = means.loc[group].values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, linewidth=2, linestyle="solid", label=group,
            color=color)
    ax.fill(angles, values, color, alpha=0.15)

plt.legend(loc="upper right", bbox_to_anchor=(0.1, 0.1), fontsize=10)
plt.tight_layout()
plt.savefig(fig("external_radar.jpg"), dpi=300)
plt.show()

print("\nStep 08 complete.")
