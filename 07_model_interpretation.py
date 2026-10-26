"""
SHAP + nomogram.

Refits the model with the same hyper-parameters as step 06 so the SHAP values
and the nomogram describe the same fit. Keep the two in sync if you change C.
"""

import os

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from sklearn.linear_model import LogisticRegression

# CONFIG
TRAIN_SCALED_FILE = "data/train_scaled.csv"
TEST_SCALED_FILE = "data/test_scaled.csv"
CORE_GENE_FILE = "results/ml/final_core_genes.csv"
RESULT_DIR = "results/ml"

ID_COLUMN = "SampleID"
LABEL_COLUMN = "label"

LR_PENALTY = "l2"
LR_C = 0.01
RANDOM_STATE = 2025

FORCE_PLOT_SAMPLE = 0      # change to look at a different test sample

os.makedirs(RESULT_DIR, exist_ok=True)
plt.rcParams["axes.unicode_minus"] = False   # render minus signs correctly


def out(name):
    return os.path.join(RESULT_DIR, name)


# ---- 1. Data and model -------
core_genes = pd.read_csv(CORE_GENE_FILE)["Core_Gene"]
core_gene_list = core_genes.tolist()

data_train = pd.read_csv(TRAIN_SCALED_FILE)
y_train = data_train[LABEL_COLUMN]
x_train_core = data_train.drop(columns=[ID_COLUMN, LABEL_COLUMN])[core_genes]

data_test = pd.read_csv(TEST_SCALED_FILE)
sample_ids_test = data_test[ID_COLUMN]
x_test_core = data_test.drop(columns=[ID_COLUMN, LABEL_COLUMN])[core_genes]

best_model = LogisticRegression(penalty=LR_PENALTY, C=LR_C,
                                class_weight="balanced", solver="liblinear",
                                random_state=RANDOM_STATE)
best_model.fit(x_train_core, y_train)


# ---- 2. SHAP values ----------
explainer = shap.LinearExplainer(best_model, x_train_core)
shap_obj = explainer(x_test_core)


# --- 2a. Feature importance bar plot ---------------------------------
plt.figure(figsize=(8, 6))
shap.summary_plot(shap_obj, x_test_core, plot_type="bar", show=False)

# blue gradient instead of the SHAP default red
ax = plt.gca()
bars = ax.patches
for bar, val in zip(bars, np.linspace(0.4, 0.8, len(bars))):
    bar.set_facecolor(cm.Blues(val))

plt.title("Feature Importance (Mean |SHAP value|)", fontsize=14)
plt.tight_layout()
plt.savefig(out("shap_feature_importance_bar.jpg"), dpi=300)
plt.show()


# --- 2b. Beeswarm summary plot ---------------------------------------
plt.figure(figsize=(8, 6))
shap.summary_plot(shap_obj, x_test_core, show=False)
plt.title("SHAP Summary Plot", fontsize=14)
plt.tight_layout()
plt.savefig(out("shap_beeswarm.jpg"), dpi=300, bbox_inches="tight")
plt.show()


# --- 2c. Dependence plot per core gene -------------------------------
print(">>> Generating dependence plots ...")
for gene in core_gene_list:
    plt.figure(figsize=(6, 5))
    # 'auto' colours by the strongest interacting feature, None turns it off
    shap.dependence_plot(gene, shap_obj.values, x_test_core,
                         interaction_index="auto", show=False,
                         title=f"Dependence Plot: {gene}")
    plt.tight_layout()
    plt.savefig(out(f"shap_dependence_{gene}.jpg"), dpi=300)
    plt.close()


# --- 2d. Force plot for a single sample ------------------------------
current_sample_id = sample_ids_test.iloc[FORCE_PLOT_SAMPLE]
print(f">>> Force plot for sample {current_sample_id}")

x_test_display = x_test_core.round(3)
shap_values_rounded = np.round(shap_obj.values, 3)

plt.figure(figsize=(8, 6))
# matplotlib=True -> static png rather than the JS widget
shap.force_plot(explainer.expected_value,
                shap_values_rounded[FORCE_PLOT_SAMPLE],
                x_test_display.iloc[FORCE_PLOT_SAMPLE],
                matplotlib=True, show=False)
plt.savefig(out(f"shap_force_plot_{current_sample_id}.jpg"), dpi=300,
            bbox_inches="tight")
plt.show()


# ---- 3. Nomogram ----
def create_nomogram(model, feature_names, feature_ranges, save_path,
                    figsize=(12, 8)):
    """Draw a nomogram for a fitted logistic regression.

    Parameters
    ----------
    model : fitted sklearn LogisticRegression
    feature_names : list of feature names, in the order used to fit the model
    feature_ranges : dict mapping feature name -> (min_value, max_value)
    save_path : output image path
    """
    coefs = model.coef_[0]
    intercept = model.intercept_[0]

    print("=" * 50)
    print("Model coefficients (sign check against SHAP):")
    for name, coef in zip(feature_names, coefs):
        direction = ("positive (higher value -> higher risk)" if coef > 0
                     else "negative (higher value -> lower risk)")
        print(f"  {name}: {coef:.4f} --> {direction}")
    print(f"  Intercept: {intercept:.4f}")
    print("=" * 50)

    fig, axes = plt.subplots(len(feature_names) + 3, 1, figsize=figsize)
    fig.suptitle("Nomogram (Based on Logistic Regression)", fontsize=14,
                 fontweight="bold")

    # Points axis, 0-100.
    ax_points = axes[0]
    ax_points.set_xlim(0, 100)
    ax_points.set_yticks([])
    ax_points.set_xlabel("")
    ax_points.set_title("Points", loc="left", fontsize=11)
    ax_points.xaxis.set_ticks(np.arange(0, 101, 10))
    ax_points.xaxis.tick_top()

    # contribution range per feature, normalises the 0-100 point scale
    total_min, total_max = 0.0, 0.0
    feature_contributions = {}
    for name, coef in zip(feature_names, coefs):
        f_min, f_max = feature_ranges[name]
        contrib_min = coef * f_min
        contrib_max = coef * f_max
        feature_contributions[name] = (min(contrib_min, contrib_max),
                                       max(contrib_min, contrib_max))
        total_min += min(contrib_min, contrib_max)
        total_max += max(contrib_min, contrib_max)

    for idx, (name, coef) in enumerate(zip(feature_names, coefs)):
        ax = axes[idx + 1]
        f_min, f_max = feature_ranges[name]

        contrib_range = (feature_contributions[name][1]
                         - feature_contributions[name][0])
        points_range = (contrib_range / (total_max - total_min)) * 100

        # negative coef -> reverse the ticks so points still run left to right
        tick_vals = (np.linspace(f_min, f_max, 6) if coef > 0
                     else np.linspace(f_max, f_min, 6))
        ax.set_xlim(0, points_range)
        ax.set_xticks(np.linspace(0, points_range, 6))
        ax.set_xticklabels([f"{v:.1f}" for v in tick_vals])
        ax.set_yticks([])
        ax.set_title(f"{name} (coef={coef:.3f})", loc="left", fontsize=11)

    ax_total = axes[-2]
    ax_total.set_xlim(0, 100)
    ax_total.set_yticks([])
    ax_total.set_title("Total Points", loc="left", fontsize=11)
    ax_total.xaxis.set_ticks(np.arange(0, 101, 20))

    # total points -> probability
    ax_risk = axes[-1]
    logit_min = total_min + intercept
    logit_max = total_max + intercept
    prob_min = 1 / (1 + np.exp(-logit_min))
    prob_max = 1 / (1 + np.exp(-logit_max))

    ax_risk.set_xlim(0, 100)
    ax_risk.set_yticks([])
    ax_risk.set_title("Risk of Disease", loc="left", fontsize=11)

    valid_ticks, valid_positions = [], []
    for p in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        if prob_min < p < prob_max:
            logit_p = np.log(p / (1 - p))
            valid_positions.append(
                (logit_p - logit_min) / (logit_max - logit_min) * 100)
            valid_ticks.append(p)

    ax_risk.set_xticks(valid_positions)
    ax_risk.set_xticklabels([f"{v:.1f}" for v in valid_ticks])

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()
    return fig


feature_ranges = {
    gene: (x_train_core[gene].min(), x_train_core[gene].max())
    for gene in core_gene_list
}

create_nomogram(model=best_model,
                feature_names=core_gene_list,
                feature_ranges=feature_ranges,
                save_path=out("nomogram.jpg"))

print("\nStep 07 complete.")
