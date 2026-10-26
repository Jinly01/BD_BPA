"""
Internal (LOOCV on GSE17114) and external (GSE209567) validation of the
core gene logistic regression. ROC / PR / DCA plus the Youden threshold.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (auc, average_precision_score,
                             balanced_accuracy_score, classification_report,
                             confusion_matrix, precision_recall_curve,
                             roc_auc_score, roc_curve)
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

# CONFIG
RAW_TRAIN_FILE = "data/intersect_gene_expression_matrix_GSE17114.csv"
TRAIN_SCALED_FILE = "data/train_scaled.csv"
TEST_SCALED_FILE = "data/test_scaled.csv"
CORE_GENE_FILE = "results/ml/final_core_genes.csv"
RESULT_DIR = "results/ml"

RAW_ID_COLUMN = "SampleID_CEL"
ID_COLUMN = "SampleID"
LABEL_COLUMN = "label"

# Logistic regression hyper-parameters (fixed, no grid search).
LR_PENALTY = "l2"
LR_C = 0.01           # 1.0 overfits badly on this sample size
RANDOM_STATE = 2025

COLOR_INTERNAL = "#1f77b4"
COLOR_EXTERNAL = "#d62728"

os.makedirs(RESULT_DIR, exist_ok=True)


def out(name):
    return os.path.join(RESULT_DIR, name)


def oriented_roc(y_true, score):
    """Single-feature ROC. AUC < 0.5 just means the gene runs the other way,
    so flip the sign and report the discrimination."""
    fpr, tpr, _ = roc_curve(y_true, score)
    roc_auc = auc(fpr, tpr)
    if roc_auc < 0.5:
        fpr, tpr, _ = roc_curve(y_true, -score)
        roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc


def net_benefit(y_true, y_proba, thresholds):
    """DCA net benefit across threshold probabilities."""
    y_true = np.asarray(y_true)
    n = len(y_true)
    values = []
    for thresh in thresholds:
        thresh = min(thresh, 0.999)   # thresh == 1 blows up the weight
        y_pred = (y_proba >= thresh).astype(int)
        tp = np.sum((y_true == 1) & (y_pred == 1))
        fp = np.sum((y_true == 0) & (y_pred == 1))
        weight = thresh / (1 - thresh)
        values.append((tp / n) - (fp / n) * weight)
    return np.array(values)


def net_benefit_treat_all(y_true, thresholds):
    """Net benefit of treating every subject."""
    y_true = np.asarray(y_true)
    n = len(y_true)
    tp_all = np.sum(y_true == 1)
    fp_all = np.sum(y_true == 0)
    weights = thresholds / (1 - thresholds)
    return (tp_all / n) - (fp_all / n) * weights


# 1. Internal validation on the unscaled training cohort
print(">>> Internal validation (LOOCV)")
data_train_raw = pd.read_csv(RAW_TRAIN_FILE)
x_train_raw = data_train_raw.drop(columns=[RAW_ID_COLUMN, LABEL_COLUMN])
label_train = data_train_raw[LABEL_COLUMN]

core_genes = pd.read_csv(CORE_GENE_FILE)["Core_Gene"]
X_core_raw = x_train_raw[core_genes]

# scaler inside the pipeline: each split standardises on its own N-1 samples
lr_pipeline = make_pipeline(
    StandardScaler(),
    LogisticRegression(penalty=LR_PENALTY, C=LR_C, class_weight="balanced",
                       random_state=RANDOM_STATE, solver="liblinear"),
)

y_pred_probs_loo = cross_val_predict(
    lr_pipeline, X_core_raw, label_train, cv=LeaveOneOut(),
    method="predict_proba", n_jobs=-1)[:, 1]

fpr_int, tpr_int, _ = roc_curve(label_train, y_pred_probs_loo)
auc_int = auc(fpr_int, tpr_int)
prec_int, rec_int, _ = precision_recall_curve(label_train, y_pred_probs_loo)
ap_int = average_precision_score(label_train, y_pred_probs_loo)
print(f"LOOCV AUC = {auc_int:.3f}, AP = {ap_int:.3f}")

fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(14, 6))

ax_roc.plot(fpr_int, tpr_int, color=COLOR_INTERNAL, lw=2,
            label=f"Internal LOOCV (AUC = {auc_int:.3f})")
ax_roc.plot([0, 1], [0, 1], color="grey", lw=1.5, linestyle="--",
            label="Chance (AUC = 0.500)")
ax_roc.set_xlim([-0.01, 1.0])
ax_roc.set_ylim([0.0, 1.05])
ax_roc.set_xlabel("False Positive Rate", fontsize=12)
ax_roc.set_ylabel("True Positive Rate", fontsize=12)
ax_roc.set_title("Internal Validation ROC Curve (LOOCV)", fontsize=14)
ax_roc.legend(loc="lower right")
ax_roc.grid(alpha=0.3)

baseline_int = label_train.sum() / len(label_train)
ax_pr.plot(rec_int, prec_int, color=COLOR_INTERNAL, lw=2,
           label=f"Internal LOOCV (AP = {ap_int:.3f})")
ax_pr.plot([0, 1], [baseline_int, baseline_int], color="grey", lw=1.5,
           linestyle="--", label=f"No Skill (AP = {baseline_int:.3f})")
ax_pr.set_xlim([0.0, 1.0])
ax_pr.set_ylim([0.4, 1.05])
ax_pr.set_xlabel("Recall", fontsize=12)
ax_pr.set_ylabel("Precision", fontsize=12)
ax_pr.set_title("Internal Validation PR Curve (LOOCV)", fontsize=14)
ax_pr.legend(loc="best")
ax_pr.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(out("internal_validation_roc_pr.jpg"), dpi=300)
plt.show()


# --- 2. Load the z-scored training and test cohorts
print("\n>>> Loading the scaled cohorts")
# The LOOCV probabilities above are indexed by the row order of the raw
# training file, and they are plotted below against labels taken from the
# scaled file, so both files must list the samples in the same order.
data_train_ss = pd.read_csv(TRAIN_SCALED_FILE)
sample_ids_train = data_train_ss[ID_COLUMN]
y_train = data_train_ss[LABEL_COLUMN]
x_train_core = data_train_ss.drop(columns=[ID_COLUMN, LABEL_COLUMN])[core_genes]

data_test_ss = pd.read_csv(TEST_SCALED_FILE)
sample_ids_test = data_test_ss[ID_COLUMN]
y_test = data_test_ss[LABEL_COLUMN]
x_test_core = data_test_ss.drop(columns=[ID_COLUMN, LABEL_COLUMN])[core_genes]

print(f"Core features: {len(core_genes)}")
print(f"Train shape: {x_train_core.shape}  Test shape: {x_test_core.shape}")

palette_name = "tab10" if len(core_genes) <= 10 else "tab20"
gene_colors = sns.color_palette(palette_name, len(core_genes))
gene_color_map = dict(zip(core_genes, gene_colors))


# --- 3. Single-gene ROC curves, internal and external
print("\n>>> Single-gene ROC curves")
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

for ax, (X, y, title) in zip(
    axes,
    [(x_train_core, y_train, "Internal Validation"),
     (x_test_core, y_test, "External Validation")],
):
    for gene in core_genes:
        fpr, tpr, gene_auc = oriented_roc(y, X[gene])
        ax.plot(fpr, tpr, color=gene_color_map[gene], lw=2.5, alpha=0.8,
                label=f"{gene} (AUC = {gene_auc:.3f})")
    ax.plot([0, 1], [0, 1], "k--", lw=1.5, alpha=0.6)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.set_xlabel("False Positive Rate", fontsize=12)
    ax.set_ylabel("True Positive Rate", fontsize=12)
    ax.set_title(title, fontsize=14)
    ax.legend(loc="lower right", fontsize=10, frameon=True, edgecolor="gray")
    ax.grid(True, linestyle="--", alpha=0.4)

plt.tight_layout()
plt.savefig(out("single_gene_roc_panels.jpg"), dpi=300, bbox_inches="tight")
plt.show()

# Same information on one axis: solid = train, dashed = test.
plt.figure(figsize=(10, 9))
contrast_colors = sns.color_palette("Set1", len(core_genes))
for i, gene in enumerate(core_genes):
    color = contrast_colors[i]
    fpr_tr, tpr_tr, auc_tr = oriented_roc(y_train, x_train_core[gene])
    plt.plot(fpr_tr, tpr_tr, color=color, linestyle="-", lw=2.5, alpha=0.9,
             label=f"{gene} (Train AUC={auc_tr:.3f})")
    fpr_te, tpr_te, auc_te = oriented_roc(y_test, x_test_core[gene])
    plt.plot(fpr_te, tpr_te, color=color, linestyle="--", lw=2, alpha=0.9,
             label=f"{gene} (Test AUC={auc_te:.3f})")

plt.plot([0, 1], [0, 1], "k--", lw=1.5, alpha=0.5)
plt.xlim([-0.02, 1.02])
plt.ylim([-0.02, 1.02])
plt.xlabel("False Positive Rate", fontsize=13)
plt.ylabel("True Positive Rate", fontsize=13)
plt.title("ROC Curve: Internal vs External Validation", fontsize=15)
plt.legend(loc="lower right", fontsize=10, frameon=True, edgecolor="gray")
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()
# plt.savefig(out("single_gene_roc_combined.jpg"), dpi=300)
plt.show()


# 4. Train the logistic regression and evaluate it externally
print("\n>>> Training the logistic regression on the core genes")
best_model = LogisticRegression(penalty=LR_PENALTY, C=LR_C,
                                class_weight="balanced", solver="liblinear",
                                random_state=RANDOM_STATE)
best_model.fit(x_train_core, y_train)

y_pred_test = best_model.predict(x_test_core)
y_proba_test = best_model.predict_proba(x_test_core)[:, 1]

bal_acc = balanced_accuracy_score(y_test, y_pred_test)
auc_ext = roc_auc_score(y_test, y_proba_test)
ap_ext = average_precision_score(y_test, y_proba_test)

print(f"Balanced accuracy: {bal_acc:.4f}")
print(f"ROC-AUC:           {auc_ext:.4f}")
print(f"PR-AUC:            {ap_ext:.4f}")
print("\nConfusion matrix:")
print(confusion_matrix(y_test, y_pred_test))
print("\nClassification report:")
print(classification_report(y_test, y_pred_test))

coefficients = pd.DataFrame({
    "Gene": list(core_genes),
    "Coefficient": best_model.coef_[0],
})
coefficients.to_csv(out("logistic_regression_coefficients.csv"), index=False)

fpr_ext, tpr_ext, _ = roc_curve(y_test, y_proba_test)
prec_ext, rec_ext, _ = precision_recall_curve(y_test, y_proba_test)


# --- 5. Internal vs external comparison (model only)
fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(15, 6))

ax_roc.plot(fpr_int, tpr_int, color=COLOR_INTERNAL, lw=2, alpha=0.8,
            label=f"Internal (LOOCV) (AUC = {auc_int:.3f})")
ax_roc.plot(fpr_ext, tpr_ext, color=COLOR_EXTERNAL, lw=2, alpha=0.8,
            label=f"External (Test) (AUC = {auc_ext:.3f})")
ax_roc.plot([0, 1], [0, 1], color="grey", lw=1.5, linestyle="--",
            label="Chance (AUC = 0.500)")
ax_roc.set_xlim([-0.02, 1.0])
ax_roc.set_ylim([0.0, 1.02])
ax_roc.set_xlabel("False Positive Rate", fontsize=12)
ax_roc.set_ylabel("True Positive Rate", fontsize=12)
ax_roc.set_title("ROC Curve Comparison", fontsize=14)
ax_roc.legend(loc="lower right", fontsize=10)
ax_roc.grid(alpha=0.3)

no_skill_ext = y_test.sum() / len(y_test)
ax_pr.plot(rec_int, prec_int, color=COLOR_INTERNAL, lw=2, alpha=0.8,
           label=f"Internal (LOOCV) (AP = {ap_int:.3f})")
ax_pr.plot(rec_ext, prec_ext, color=COLOR_EXTERNAL, lw=2, alpha=0.8,
           label=f"External (Test) (AP = {ap_ext:.3f})")
ax_pr.plot([0, 1], [no_skill_ext, no_skill_ext], color="grey", lw=1.5,
           linestyle="--", label=f"Ext. Baseline (AP = {no_skill_ext:.3f})")
ax_pr.set_xlim([0.0, 1.0])
ax_pr.set_ylim([0.0, 1.02])
ax_pr.set_xlabel("Recall", fontsize=12)
ax_pr.set_ylabel("Precision", fontsize=12)
ax_pr.set_title("Precision-Recall Curve Comparison", fontsize=14)
ax_pr.legend(loc="best", fontsize=10)
ax_pr.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(out("model_internal_vs_external_roc_pr.jpg"), dpi=300)
plt.show()


# 6. Single genes plus the combined model, three panels
print("\n>>> Combined single-gene and model panels")
probs_internal_model = y_pred_probs_loo
probs_external_model = y_proba_test

prec_model_int, rec_model_int, _ = precision_recall_curve(
    y_train, probs_internal_model)
ap_model_int = average_precision_score(y_train, probs_internal_model)
prec_model_ext, rec_model_ext, _ = precision_recall_curve(
    y_test, probs_external_model)
ap_model_ext = average_precision_score(y_test, probs_external_model)

fig, axes = plt.subplots(1, 3, figsize=(20, 6))
lw_gene, lw_model, alpha_gene = 2, 2.5, 0.6

for ax, (X, y, model_fpr, model_tpr, model_auc, title) in zip(
    axes[:2],
    [(x_train_core, y_train, fpr_int, tpr_int, auc_int,
      "ROC Curves (Internal Validation)"),
     (x_test_core, y_test, fpr_ext, tpr_ext, auc_ext,
      "ROC Curves (External Validation)")],
):
    for gene in core_genes:
        fpr, tpr, gene_auc = oriented_roc(y, X[gene])
        ax.plot(fpr, tpr, color=gene_color_map[gene], lw=lw_gene,
                alpha=alpha_gene, label=f"{gene} (AUC={gene_auc:.3f})")
    ax.plot(model_fpr, model_tpr, color="black", lw=lw_model, linestyle="--",
            label=f"Logistic Regression (AUC={model_auc:.3f})")
    ax.plot([0, 1], [0, 1], "k:", lw=1, alpha=0.4)
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.3)

ax = axes[2]
base_int = y_train.sum() / len(y_train)
base_ext = y_test.sum() / len(y_test)
ax.plot(rec_model_int, prec_model_int, color="tab:purple", lw=2,
        label=f"Logistic Regression - Internal (AP={ap_model_int:.3f})")
ax.axhline(y=base_int, color="tab:purple", linestyle=":", alpha=0.6, lw=1,
           label="Internal Baseline")
ax.plot(rec_model_ext, prec_model_ext, color="tab:red", lw=2,
        label=f"Logistic Regression - External (AP={ap_model_ext:.3f})")
ax.axhline(y=base_ext, color="tab:red", linestyle=":", alpha=0.6, lw=1,
           label="External Baseline")
ax.set_title("PR Curves (Model Comparison)", fontsize=14)
ax.set_xlabel("Recall (Sensitivity)")
ax.set_ylabel("Precision (PPV)")
ax.set_xlim([0, 1.02])
ax.set_ylim([0, 1.02])
ax.legend(loc="lower left", fontsize=10, frameon=True)
ax.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
plt.savefig(out("combined_roc_pr_panels.jpg"), dpi=300, bbox_inches="tight")
plt.show()


# --- 7. Decision curve analysis
print("\n>>> Decision curve analysis")
thresholds = np.linspace(0.01, 0.99, 100)
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

for ax, (X_fit, X_eval, y_fit, y_eval, model_probs, title) in zip(
    axes,
    [(x_train_core, x_train_core, y_train, y_train, probs_internal_model,
      "Internal Validation (DCA)"),
     (x_train_core, x_test_core, y_train, y_test, probs_external_model,
      "External Validation (DCA)")],
):
    ax.plot(thresholds, np.zeros_like(thresholds), color="black",
            linestyle=":", lw=1.5, label="Treat None")
    nb_all = net_benefit_treat_all(y_eval, thresholds)
    ax.plot(thresholds, nb_all, color="gray", linestyle=":", lw=1.5,
            label="Treat All")

    nb_model = net_benefit(y_eval, model_probs, thresholds)
    ax.plot(thresholds, nb_model, color="black", linestyle="--", lw=3,
            label="Combined Model")

    # univariate LR to turn each gene into a probability
    for gene in core_genes:
        lr_single = LogisticRegression(solver="liblinear")
        lr_single.fit(X_fit[[gene]], y_fit)
        probs_gene = lr_single.predict_proba(X_eval[[gene]])[:, 1]
        ax.plot(thresholds, net_benefit(y_eval, probs_gene, thresholds),
                color=gene_color_map[gene], lw=2, alpha=0.6, label=gene)

    ax.set_title(title, fontsize=14)
    ax.set_xlabel("Threshold Probability", fontsize=12)
    ax.set_ylabel("Net Benefit", fontsize=12)
    ax.set_ylim(-0.05, max(nb_all.max(), nb_model.max()) + 0.05)
    ax.set_xlim(0, 1.0)
    ax.legend(loc="upper right", fontsize=9, frameon=True, edgecolor="gray")
    ax.grid(True, linestyle="--", alpha=0.4)

plt.tight_layout()
plt.savefig(out("decision_curve_analysis.jpg"), dpi=300, bbox_inches="tight")
plt.show()


# 8. Threshold optimisation (Youden's index)
print("\n>>> Threshold optimisation")
fpr_t, tpr_t, thresh_t = roc_curve(y_test, y_proba_test)
best_thresh = thresh_t[np.argmax(tpr_t - fpr_t)]
print(f"Optimal threshold (Youden's index): {best_thresh:.4f}")

y_pred_opt = (y_proba_test >= best_thresh).astype(int)
cm_opt = confusion_matrix(y_test, y_pred_opt)
print("Confusion matrix at the optimal threshold:")
print(cm_opt)

plt.figure(figsize=(6, 6))
sns.heatmap(cm_opt, annot=True, fmt="d", cmap="Blues",
            annot_kws={"size": 16}, cbar=False)
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.xticks([0.5, 1.5], ["Control", "Case"])
plt.yticks([0.5, 1.5], ["Control", "Case"])
plt.tight_layout()
plt.savefig(out("external_confusion_matrix.jpg"), dpi=300)
plt.show()


# 9. Risk score distribution and waterfall plot
df_risk = pd.DataFrame({"Risk Score": y_proba_test, "Group": y_test})
df_risk["Group"] = df_risk["Group"].map({0: "Control", 1: "Case"})

plt.figure(figsize=(6, 5))
sns.violinplot(x="Group", y="Risk Score", data=df_risk,
               palette={"Control": COLOR_INTERNAL, "Case": COLOR_EXTERNAL},
               inner=None, alpha=0.5, cut=0)
sns.stripplot(x="Group", y="Risk Score", data=df_risk, color="black",
              size=5, jitter=0.2, alpha=0.6)
sns.pointplot(x="Group", y="Risk Score", data=df_risk, estimator=np.mean,
              color="black", linestyles="--", capsize=0.1, markers="D",
              errorbar=None)
plt.title("Predicted Risk Score Distribution", fontsize=14, fontweight="bold")
plt.ylabel("Predicted Probability")
plt.xlabel("")
plt.axhline(0.5, color="gray", linestyle=":", label="Default Threshold")
plt.tight_layout()
plt.savefig(out("risk_score_distribution.jpg"), dpi=300)
plt.show()

df_sort = df_risk.sort_values(by="Risk Score").reset_index(drop=True)
bar_colors = df_sort["Group"].map({"Control": COLOR_INTERNAL,
                                   "Case": COLOR_EXTERNAL})

plt.figure(figsize=(10, 5))
plt.bar(range(len(df_sort)), df_sort["Risk Score"], color=bar_colors,
        alpha=0.8)
plt.axhline(0.5, color="black", linestyle="--", linewidth=1,
            label="Threshold 0.5")
plt.title("Prediction Waterfall Plot", fontsize=14, fontweight="bold")
plt.xlabel("Individual Samples (Sorted by Risk)")
plt.ylabel("Predicted Probability")

from matplotlib.patches import Patch  # noqa: E402  (used only for the legend)

legend_elements = [Patch(facecolor=COLOR_EXTERNAL, label="Case"),
                   Patch(facecolor=COLOR_INTERNAL, label="Control")]
plt.legend(handles=legend_elements, loc="upper left")
plt.tight_layout()
plt.savefig(out("prediction_waterfall.jpg"), dpi=300)
plt.show()

print("\nStep 06 complete.")
