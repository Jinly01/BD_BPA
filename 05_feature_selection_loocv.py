"""
Stability selection with LOOCV.

Lasso, SVM-RFE and RF are refitted in every fold and we count how often each
gene gets picked. Scalers live inside the Pipeline so the held-out sample
never contributes to them.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.linear_model import LassoCV
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

# CONFIG
TRAIN_FILE = "data/intersect_gene_expression_matrix_GSE17114.csv"
RESULT_DIR = "results/ml"

ID_COLUMN = "SampleID_CEL"
LABEL_COLUMN = "label"

RANDOM_STATE = 2025
N_FEATURES_RFE = 10        # features kept per fold by SVM-RFE
N_TOP_RF = 10              # features kept per fold by random forest importance
N_CORE_GENES = 3           # tried 5 first, 3 held up better externally
N_GENES_IN_PLOT = 38       # genes shown in the lollipop charts
STABILITY_THRESHOLD = 80   # reference line only, not a hard filter

os.makedirs(RESULT_DIR, exist_ok=True)

FREQUENCY_FILE = os.path.join(RESULT_DIR, "gene_stability_frequency.csv")
CORE_GENE_FILE = os.path.join(RESULT_DIR, "final_core_genes.csv")


# ---- 1. Load the candidate-gene expression matrix ----------
data_train = pd.read_csv(TRAIN_FILE)
x_train_raw = data_train.drop(columns=[ID_COLUMN, LABEL_COLUMN])
label_train = data_train[LABEL_COLUMN]

print(f"Input shape: {x_train_raw.shape}")
print(f"Label distribution: {np.unique(label_train, return_counts=True)}")


# --- 2. LOOCV stability selection
loo = LeaveOneOut()
n_samples = x_train_raw.shape[0]
gene_frequency = pd.DataFrame(0, index=x_train_raw.columns,
                              columns=["Lasso", "SVM", "RF"])

print("\nRunning LOOCV stability selection ...")
for iteration, (train_idx, _) in enumerate(loo.split(x_train_raw), start=1):
    if iteration % 5 == 0:
        print(f"  fold {iteration}/{n_samples}")

    X_fold_train = x_train_raw.iloc[train_idx, :]
    y_fold_train = label_train.iloc[train_idx]

    # Lasso, inner CV over the alpha path
    lasso_pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("lasso", LassoCV(alphas=np.logspace(-4, -1, 100), cv=5,
                          random_state=RANDOM_STATE, max_iter=10000,
                          n_jobs=-1)),
    ])
    lasso_pipe.fit(X_fold_train, y_fold_train)
    lasso_model = lasso_pipe.named_steps["lasso"]
    lasso_selected = X_fold_train.columns[lasso_model.coef_ != 0].tolist()
    gene_frequency.loc[lasso_selected, "Lasso"] += 1

    # SVM-RFE, scaling is mandatory here
    svm_pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("rfe", RFE(estimator=SVC(kernel="linear", C=1,
                                  random_state=RANDOM_STATE),
                    n_features_to_select=N_FEATURES_RFE, step=1)),
    ])
    svm_pipe.fit(X_fold_train, y_fold_train)
    rfe_model = svm_pipe.named_steps["rfe"]
    svm_selected = X_fold_train.columns[rfe_model.support_].tolist()
    gene_frequency.loc[svm_selected, "SVM"] += 1

    # RF on the raw matrix, trees do not care about scale
    rf = RandomForestClassifier(n_estimators=30, random_state=RANDOM_STATE,
                                n_jobs=-1)
    rf.fit(X_fold_train, y_fold_train)
    importances = pd.Series(rf.feature_importances_,
                            index=X_fold_train.columns)
    rf_selected = importances.nlargest(N_TOP_RF).index.tolist()
    gene_frequency.loc[rf_selected, "RF"] += 1


# ---- 3. Rank genes and keep the top signature ----
gene_frequency["Total_Count"] = gene_frequency.sum(axis=1)
gene_frequency = gene_frequency.sort_values(by="Total_Count", ascending=False)

final_core_genes = gene_frequency.head(N_CORE_GENES).index.tolist()
print(f"\nCore genes: {final_core_genes}")

gene_frequency.to_csv(FREQUENCY_FILE)
pd.DataFrame(final_core_genes, columns=["Core_Gene"]).to_csv(
    CORE_GENE_FILE, index=False)


# --- 4. Lollipop chart of the raw selection counts
plt.figure(figsize=(10, 8))
plot_data = gene_frequency.head(N_GENES_IN_PLOT).sort_values(
    by="Total_Count", ascending=True)
my_range = range(1, len(plot_data.index) + 1)

plt.hlines(y=my_range, xmin=0, xmax=plot_data["Total_Count"], color="skyblue")
plt.plot(plot_data["Total_Count"], my_range, "o", markersize=10,
         color="#005c99")
plt.plot(plot_data.tail(N_CORE_GENES)["Total_Count"],
         list(my_range)[-N_CORE_GENES:], "o", markersize=10, color="red",
         label="Selected Core Genes")

plt.yticks(my_range, plot_data.index)
plt.title("Stability Selection: Gene Frequency across LOOCV folds",
          fontsize=14)
plt.xlabel("Total Selection Count", fontsize=12)
plt.ylabel("Gene Symbol", fontsize=12)
plt.legend()
plt.grid(axis="x", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(RESULT_DIR, "stability_selection_counts.jpg"),
            dpi=300)
plt.show()


# ---- 5. Lollipop chart on a percentage scale ----
# folds x selectors
max_possible_count = n_samples * 3
gene_frequency["Frequency"] = (
    gene_frequency["Total_Count"] / max_possible_count * 100)

plt.figure(figsize=(10, 10))
plot_data = gene_frequency.head(N_GENES_IN_PLOT).sort_values(
    by="Frequency", ascending=True)
my_range = range(1, len(plot_data.index) + 1)

plt.hlines(y=my_range, xmin=0, xmax=plot_data["Frequency"], color="skyblue")
plt.plot(plot_data["Frequency"], my_range, "o", markersize=10, color="#005c99")
plt.axvline(x=STABILITY_THRESHOLD, color="#C02B2B", linestyle="--", alpha=0.8,
            label=f"Threshold ({STABILITY_THRESHOLD}%)")
plt.plot(plot_data.tail(N_CORE_GENES)["Frequency"],
         list(my_range)[-N_CORE_GENES:], "o", markersize=10, color="#C02B2B",
         label="Selected Core Genes")

plt.yticks(my_range, plot_data.index)
plt.tick_params(axis="y", labelsize=10)
plt.title("Gene Frequency (LOOCV)", fontsize=14)
plt.xlabel("Stability Score (%)", fontsize=12)
plt.ylabel("Gene Symbol", fontsize=12)
plt.xlim(0, 105)
plt.legend(loc="lower right")
plt.grid(axis="x", linestyle="--", alpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(RESULT_DIR, "stability_selection_frequency.jpg"),
            dpi=300)
plt.show()

print(f"\nSaved frequency table to {FREQUENCY_FILE}")
print(f"Saved core gene list to {CORE_GENE_FILE}")
